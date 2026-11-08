"""
TDA heuristics against the proven optima, and TDA techniques inside the exact PLNE.

Part `heuristics`: every TDA heuristic on the Petit instances, judged by verify.py.

Part `exact`: the tightened Gurobi model of plne_gurobi.py on the Petit_N60
instances, in five variants, one thread, same deterministic work limit:
  base        tightened model, no starting solution
  warm        + starting solution: best feasible TDA heuristic
  cuts        + rounded capacity cuts on persistence clusters (TDA_algos/plne_cuts.py)
  warm+cuts   both
  pipeline    PLNE on a restricted set of edges, then warm+cuts started from its solution

Runs are compared on Gurobi's work units (m.Work), which do not depend on the load
of the machine, unlike seconds. Each result is checked by verify.py, and every cut is
tested on the solution of the base run: a valid cut cannot be violated by a feasible network.

Written to output/tda_plne/: heuristics.jsonl, exact.jsonl (one line per run, a stopped
run resumes), heuristics.csv, exact.csv.

Usage:
  python run_tda_plne.py heuristics
  python run_tda_plne.py exact --work-limit 120 --workers 2
  python run_tda_plne.py report
"""
import argparse
import contextlib
import csv
import io
import json
import logging
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import gurobipy as gp
from gurobipy import GRB
import networkx as nx

logging.disable(logging.INFO)

import plne_gurobi as P
from run_plne_verified import instance_from_graph, warm_start_values
from src.instances import load_instance
from TDA_algos import bounds, departs, local_search, partition, reductions, spectral
from TDA_algos.decomposition import SubGraph
from TDA_algos.persistence import h1_planar, shortest_distances
from TDA_algos.plne_cuts import add_cuts, h0_sets, tomato_sets
from verify import check_solution, solution_key

HERE = Path(__file__).resolve().parent
INSTANCE_DIR = HERE / "data" / "verified_instances"
OUT = HERE / "output" / "tda_plne"
OPT_FILE = HERE / "output" / "gurobi_lower_bounds_n60" / "results.csv"
CAMPAIGN_FILE = HERE / "output" / "gurobi_campaign" / "summary.csv"
LOAD_LABELS = ('Conso', 'Prod')
STD = ('none', 'capacity', 'sqrt')
CONFIGS = ('base', 'warm', 'cuts', 'warm+cuts', 'pipeline')
N_NEAREST = 3            # restricted model: shortest paths from each load to its 3 nearest terminals
RESTRICTED_WORK = 20.0
LS_TIME = 30.0


def load(path):
    with contextlib.redirect_stdout(io.StringIO()):
        return load_instance(path)


def known_optima() -> dict:
    """{instance: (J*, Gurobi seconds, thread count, solved how)} from experiment 5 and the Gurobi campaign."""
    out = {}
    with open(OPT_FILE, encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if float(r['dp_bound']) >= float(r['gurobi_obj']) - 1e-3 or r['gurobi_status'] == 'OPTIMAL':
                out[r['instance']] = {'J': float(r['gurobi_obj']), 'seconds': float(r['gurobi_time_s']), 'source': 'exp5'}
    if CAMPAIGN_FILE.exists():
        with open(CAMPAIGN_FILE, encoding='utf-8') as fh:
            for r in csv.DictReader(fh):
                if r['status'] in ('OPTIMAL', 'OPTIMAL_BY_DP'):
                    out[r['instance']] = {'J': float(r['obj']), 'seconds': float(r['time_solve']), 'source': 'campaign'}
    return out


# ------------------------------------------------------------------------------ heuristics

def heuristics(g, ls_time: float = LS_TIME) -> dict:
    """{method: {cert, time, trees}} for every TDA heuristic on instance g, plus the portfolio."""
    out = {}

    def add(name, trees, t):
        out[name] = {'cert': check_solution(g, trees), 'time': t, 'trees': trees}

    def solve(name, inits, **kw):
        t = time.perf_counter()
        best = partition.best_of(g, inits, **kw)
        add(name, best['trees'], time.perf_counter() - t)
        return best

    std = {n: partition.standard_inits(g)[n] for n in STD}
    solve('voronoi_variable', std)
    t = time.perf_counter()
    tom = partition.tomato_init_weights(g)['weights']
    t_tom = time.perf_counter() - t
    solve('tomato_init', {**std, 'tomato': tom})
    out['tomato_init']['time'] += t_tom
    solve('grow', std, grow=True)
    t = time.perf_counter()
    attr = spectral.spectral_edge_costs(g.G)
    t_emb = time.perf_counter() - t
    solve('spectral_costs', std, weight=attr)
    out['spectral_costs']['time'] += t_emb
    t = time.perf_counter()
    add('spectral_cells', partition.steiner_builder(g, spectral.anchored_spectral_cells(g.G)), time.perf_counter() - t)
    mt = None
    for m in departs.METHODS:
        run = solve(f'departs_{m}', std, tree_builder=partition.departs_builder(m))
        mt = run if m == 'mergetree' else mt
    t = time.perf_counter()
    h1 = h1_planar(g.G)
    t_h1 = time.perf_counter() - t
    t = time.perf_counter()
    tr, _ = local_search.local_search(g, mt['trees'], local_search.candidate_insert_edges(g, h1),
                                      local_search.candidate_relocations(g, 0.3), time_limit=ls_time)
    add('mergetree+ls_tda', tr, out['departs_mergetree']['time'] + t_h1 + time.perf_counter() - t)
    t = time.perf_counter()
    tr, _ = local_search.local_search(g, mt['trees'], local_search.candidate_insert_edges(g),
                                      local_search.candidate_relocations(g), time_limit=ls_time)
    add('mergetree+ls_all', tr, out['departs_mergetree']['time'] + time.perf_counter() - t)
    best = min(out.values(), key=lambda r: solution_key(r['cert']))
    out['portfolio'] = {'cert': best['cert'], 'time': sum(r['time'] for r in out.values()), 'trees': best['trees']}
    return out


def heuristics_job(path: Path) -> dict:
    g = load(path)
    G = g.G
    res = heuristics(g)
    lb = bounds.lower_bound(g)
    R, st = reductions.reduce_graph(G)
    return {'instance': path.stem, 'n_nodes': len(G), 'n_edges': G.number_of_edges(),
            'n_loads': sum(l in LOAD_LABELS for _, l in G.nodes('label')),
            'n_sources': sum(l == 'Source' for _, l in G.nodes('label')),
            'lb_h0': lb['lb'], 'reduced_nodes': st['nodes_after'], 'reduced_edges': st['edges_after'],
            'methods': {k: {'J': r['cert']['J'], 'feasible': r['cert']['feasible'], 'time': r['time'],
                            'n_violations': sum(v[0] for v in r['cert']['violations'].values()) + len(r['cert']['errors'])}
                        for k, r in res.items()}}


# ------------------------------------------------------------------------------ exact

class Progress:
    """Gurobi callback: root bound, and the work and time of the best solution."""

    def __init__(self):
        self.root = None
        self.work_best = None
        self.first_obj = None

    def __call__(self, model, where):
        if where == GRB.Callback.MIPNODE and model.cbGet(GRB.Callback.MIPNODE_NODCNT) == 0:
            self.root = model.cbGet(GRB.Callback.MIPNODE_OBJBND)
        elif where == GRB.Callback.MIPSOL:
            self.work_best = model.cbGet(GRB.Callback.WORK)
            if self.first_obj is None:
                self.first_obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)


def exact_run(graph, work_limit: float, start=None, cuts=None) -> dict:
    """One Gurobi run on `graph`: dict of metrics, with 'trees' and 'arcs' (arcs used, for the check of the cuts)."""
    inst = instance_from_graph(graph)
    P.limit_feeders(inst)
    t0 = time.perf_counter()
    m, var = P.build_model(inst, lmax_mode="shortest", disjoint=True)
    P.align_with_verify(m, var, inst)
    P.tighten_model(m, var, inst)
    n_cuts = add_cuts(m, var, inst, cuts) if cuts else 0
    m.update()
    t_build = time.perf_counter() - t0
    used_start = False
    if start is not None:
        values = warm_start_values(graph, inst, var, var["f"], start)
        if values is not None:
            for v, val in values.items():
                v.Start = val
            used_start = True
    m.Params.OutputFlag = 0
    m.Params.Threads = 1
    m.Params.WorkLimit = work_limit
    m.Params.MIPGap = 1e-6
    cb = Progress()
    t0 = time.perf_counter()
    m.optimize(cb)
    res = {'status': P.STATUS_NAMES.get(m.Status, str(m.Status)), 'obj': None, 'bound': None, 'gap': None,
           'work': m.Work, 'seconds': time.perf_counter() - t0, 'n_vars': m.NumVars, 'n_cons': m.NumConstrs,
           'nodes': m.NodeCount, 'root_bound': cb.root, 'work_best': cb.work_best, 'first_obj': cb.first_obj,
           'n_cuts': n_cuts, 'used_start': used_start, 'time_build': t_build, 'feasible': None, 'J_verify': None}
    if abs(m.ObjBound) < GRB.INFINITY:
        res['bound'] = m.ObjBound
    if m.SolCount == 0:
        return {**res, 'trees': None, 'arcs': None}
    res['obj'], res['gap'] = m.ObjVal, m.MIPGap
    trees = {r: set() for r in inst["sources"]}
    arcs = set()
    for k in inst["feeders"]:
        trees[inst["src_of"][k]].update(e for e in var["E"] if var["z"][e, k].X > 0.5)
    for (a, k), v in var["x"].items():
        if v.X > 0.5:
            arcs.add(a)
    trees = {r: sorted(es) for r, es in trees.items()}
    return {**res, 'trees': trees, 'arcs': arcs}


def restricted_graph(graph, solutions: list) -> tuple:
    """(restricted graph, stats): edges of the heuristic solutions + shortest paths from each load to its nearest terminals."""
    G = graph.G
    loads = [v for v, l in G.nodes('label') if l in LOAD_LABELS]
    terminals = loads + [v for v, l in G.nodes('label') if l == 'Source']
    keep = {tuple(sorted(e)) for trees in solutions for es in trees.values() for e in es}
    _, index, D = shortest_distances(G, loads)
    for i, l in enumerate(loads):
        near = sorted((t for t in terminals if t != l), key=lambda t: D[i, index[t]])[:N_NEAREST]
        paths = nx.single_source_dijkstra_path(G, l)
        for t in near:
            keep.update(tuple(sorted(e)) for e in zip(paths[t], paths[t][1:]))
    sub = SubGraph(graph, list(G))
    sub.G.remove_edges_from([e for e in G.edges if tuple(sorted(e)) not in keep])
    return sub, {'edges_full': G.number_of_edges(), 'edges_kept': sub.G.number_of_edges()}


def cuts_violated(cuts: dict, arcs: set) -> int:
    """Number of cuts violated by the network made of the arcs `arcs`."""
    return sum(sum(1 for a in arcs if a[1] in S and a[0] not in S) < c['r'] for S, c in cuts.items())


def strip(res: dict) -> dict:
    return {k: v for k, v in res.items() if k not in ('trees', 'arcs')}


def certified(graph, res: dict) -> dict:
    if res['trees'] is not None:
        cert = check_solution(graph, res['trees'])
        res.update(feasible=cert['feasible'], J_verify=cert['J'])
    return res


def exact_job(stem: str, work_limit: float, done: set) -> list:
    """The variants of `CONFIGS` not yet in `done`, on instance `stem`."""
    g = load(INSTANCE_DIR / f"{stem}.json")
    t = time.perf_counter()
    res_h = heuristics(g, ls_time=10.0)
    t_heur = time.perf_counter() - t
    feasible = [r for r in res_h.values() if r['cert']['feasible']]
    start = min(feasible, key=lambda r: r['cert']['J'])['trees'] if feasible else None
    t = time.perf_counter()
    h0 = h0_sets(g)
    cuts = tomato_sets(g)
    cuts.update(h0)
    t_cuts = time.perf_counter() - t
    common = {'instance': stem, 'heuristic_J': min((r['cert']['J'] for r in feasible), default=None),
              'time_heuristics': t_heur, 'time_cuts': t_cuts, 'n_h0_sets': len(h0)}
    rows, base_arcs = [], None
    for cfg in CONFIGS:
        if (stem, cfg) in done:
            continue
        extra = {}
        if cfg == 'base':
            r = certified(g, exact_run(g, work_limit))
            base_arcs = r['arcs']
            if base_arcs is not None:
                extra['cuts_violated_by_base'] = cuts_violated(cuts, base_arcs)
        elif cfg == 'warm':
            r = certified(g, exact_run(g, work_limit, start=start))
        elif cfg == 'cuts':
            r = certified(g, exact_run(g, work_limit, cuts=cuts))
        elif cfg == 'warm+cuts':
            r = certified(g, exact_run(g, work_limit, start=start, cuts=cuts))
        else:
            sub, st = restricted_graph(g, [x['trees'] for x in res_h.values()])
            rr = certified(sub, exact_run(sub, RESTRICTED_WORK, start=start))
            sub_start = rr['trees'] if rr['trees'] is not None and rr['feasible'] else start
            r = certified(g, exact_run(g, max(work_limit - rr['work'], 1.0), start=sub_start, cuts=cuts))
            extra = {**st, 'restricted_status': rr['status'], 'restricted_J': rr['J_verify'],
                     'restricted_feasible': rr['feasible'], 'restricted_work': rr['work'],
                     'restricted_vars': rr['n_vars']}
            r['work'] += rr['work']   # the proof pays for the restricted solve
        rows.append({**common, 'config': cfg, **strip(r), **extra})
    return rows


# ------------------------------------------------------------------------------ main

def read_jsonl(path: Path) -> list:
    if not path.exists():
        return []
    with open(path, encoding='utf-8') as fh:
        return [json.loads(line) for line in fh if line.strip()]


def to_builtin(o):
    return o.item() if hasattr(o, 'item') else str(o)


def write_csv(rows: list, path: Path) -> None:
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter=';')
        w.writeheader()
        w.writerows(rows)


def run_heuristics(a) -> None:
    path = OUT / "heuristics.jsonl"
    done = {r['instance'] for r in read_jsonl(path)}
    files = sorted(p for p in INSTANCE_DIR.glob(f"{a.prefix}*.json") if p.stem not in done)
    print(f"{len(files)} instances to run ({len(done)} already done)", flush=True)
    for i, p in enumerate(files, 1):
        rec = heuristics_job(p)
        with open(path, 'a', encoding='utf-8') as fh:
            fh.write(json.dumps(rec, default=to_builtin) + "\n")
        pf = rec['methods']['portfolio']
        print(f"[{i}/{len(files)}] {p.stem}: portfolio J={pf['J']:.2f} feasible={pf['feasible']} ({pf['time']:.1f} s)", flush=True)


def run_exact(a) -> None:
    path = OUT / "exact.jsonl"
    done = {(r['instance'], r['config']) for r in read_jsonl(path)}
    with open(OPT_FILE, encoding='utf-8') as fh:
        order = sorted(csv.DictReader(fh), key=lambda r: float(r['gurobi_time_s']))
    todo = [r['instance'] for r in order if any((r['instance'], c) not in done for c in CONFIGS)]
    print(f"{len(todo)} instances to run, work limit {a.work_limit:g}", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as pool:
        futs = {pool.submit(exact_job, s, a.work_limit, done): s for s in todo}
        for n, fut in enumerate(as_completed(futs), 1):
            rows = fut.result()
            with open(path, 'a', encoding='utf-8') as fh:
                for r in rows:
                    fh.write(json.dumps(r, default=to_builtin) + "\n")
            print(f"[{n}/{len(todo)}] {futs[fut]}: " + "  ".join(
                f"{r['config']}={r['status'][:5]},w{r['work']:.0f},lb{(r['bound'] or 0):.1f}" for r in rows)
                + f"  [{(time.time() - t0) / 60:.1f} min]", flush=True)


def report(_) -> None:
    write_csv([{**{k: v for k, v in r.items() if k != 'methods'}, 'method': m, **x}
               for r in read_jsonl(OUT / "heuristics.jsonl") for m, x in r['methods'].items()], OUT / "heuristics.csv")
    write_csv(read_jsonl(OUT / "exact.jsonl"), OUT / "exact.csv")
    print("CSV written in", OUT)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("part", choices=("heuristics", "exact", "report"))
    p.add_argument("--prefix", default="Petit")
    p.add_argument("--work-limit", type=float, default=120.0)
    p.add_argument("--workers", type=int, default=2)
    a = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    {'heuristics': run_heuristics, 'exact': run_exact, 'report': report}[a.part](a)


if __name__ == "__main__":
    main()
