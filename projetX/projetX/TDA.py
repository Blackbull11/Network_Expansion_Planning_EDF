"""
Experiments of the Topological Data Analysis tools (TDA_algos/) on the
network planning instances. Every solution is judged by verify.check_solution.

Experiments (python TDA.py <experiment> [options]):
  e1  descriptors  persistence diagrams (H0 loads, H1 planarised graph, H0
                   distance to the sources, ToMATo density) and simple
                   features; do the diagrams separate the instance types?
  e2  conditions   necessary conditions: verify.py vs TDA_algos.bounds
                   (reach filtration, transportation LP, length budget)
  e3  bounds       H0 lower bound on J vs heuristic J and PLNE (bound, J) when
                   output/plne_verified/<instance>.json exists
  e4  reductions   exact graph reductions (size), --plne T: PLNE with and
                   without reduction (T seconds each)
  e5  init         variable distance from the 3 fixed weightings vs from the
                   ToMATo cluster assignment
  e6  partition    Voronoi on lengths vs Voronoi on spectral edge costs vs
                   anchored spectral cells
  e7  departs      Steiner tree per cell vs designed départs (kruskal = H0
                   capacitated single linkage, savings, mergetree)
  e8  local        local search after the heuristic: insertion (all edges /
                   H1-persistent edges), relocation (all loads / small margin)
  e9  grow         variable distance shrinking only vs shrinking + growing
  e10 select       predict the best configuration from simple / TDA features
                   (k-NN, leave one instance type out); --labels <results.csv>
                   of e5, e6, e7, e8, e9 or grid
  e11 scaling      computing time of each TDA component, summer / winter
                   distance between diagrams (verified_instances)
  grid             full ablation: init {std, tda} x distance {graph, spectral}
                   x départs {steiner, mergetree} x local search {none, tda}
  all              e1 to e9, grid, then e10 on the grid results

Instances:
  --set gen        data/instances (g01-g40, default); --petites / --grandes; motifs
  --set verified   verified_instances (sorted by real size); --prefix Petit, --max-nodes, --n

Results in output/tda/<date>_<experiment>[_label]/: results.csv, summary.md
(+ figures, diagrams.npz for e1).

Examples:
  python TDA.py e1 --petites
  python TDA.py e7 g05 g10 --label test
  python TDA.py grid --set verified --prefix Petit --max-nodes 60
  python TDA.py e10 --labels output/tda/<run>_grid/results.csv
  python TDA.py all --grandes --label grandes
"""
import argparse
import contextlib
import io
import json
import logging
import math
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

logging.disable(logging.INFO)   # GraphData and voronoi log every step

from src.instances import load_instance, verified_files, instance_size
from verify import check_solution, necessary_conditions as verify_conditions, solution_key, source_limits
from TDA_algos import bounds, departs, features, local_search, partition, reductions, spectral
from TDA_algos.persistence import h1_planar, terminal_h0_diagram
from TDA_algos.vectorize import wasserstein

ROOT = Path(__file__).resolve().parent
INSTANCES_DIR = ROOT / "data" / "instances"
PLNE_DIR = ROOT / "output" / "plne_verified"
TDA_DIR = ROOT / "output" / "tda"


# ----------------------------------------------------------------------------
# Instances and reports
# ----------------------------------------------------------------------------

def select_instances(args) -> list:
    if args.set == 'verified':
        files = verified_files(args.prefix)
        if args.max_nodes:
            files = [f for f in files if instance_size(f)[0] <= args.max_nodes]
    else:
        files = sorted(INSTANCES_DIR.glob("*.json"))
        if args.petites:
            files = [f for f in files if int(f.stem[1:3]) <= 20]
        if args.grandes:
            files = [f for f in files if int(f.stem[1:3]) > 20]
    if args.motifs:
        files = [f for f in files if any(m in f.stem for m in args.motifs)]
    if args.n:
        files = files[:args.n]
    return files


def instance_group(path: Path) -> str:
    """Instance type: ville / montagne / epars / mixte, or size class + season for verified_instances."""
    parts = path.stem.split('_')
    if path.parent == INSTANCES_DIR:
        return parts[1]
    season = 'Ete' if 'Ete' in path.stem else ('Hiver' if 'Hiver' in path.stem else '')
    return f"{parts[0]}_{season}"


def load(path):
    """GraphData of `path`, without the prints of GraphData."""
    with contextlib.redirect_stdout(io.StringIO()):
        return load_instance(path)


def lmax_of(graph) -> float:
    sources = [n for n, lab in graph.G.nodes('label') if lab == 'Source']
    l = min(source_limits(graph, s)['l_max_depart'] for s in sources)
    return l if math.isfinite(l) and l > 0 else 1.0


def cert_row(cert: dict) -> dict:
    v = cert['violations']
    return {'feasible': bool(cert['feasible']), 'J': cert['J'], 'PD': cert['PD'],
            'n_viol': sum(x[0] for x in v.values()), 'viol_L': v['L'][0], 'viol_P': v['P'][0],
            'viol_source': v['source'][0], 'viol_n': v['n'][0], 'excess': cert['excess'],
            'n_errors': len(cert['errors'])}


def fmt(v, digits=4) -> str:
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        if math.isnan(v):
            return ""
        if math.isinf(v):
            return "inf"
        return f"{v:.{digits}g}"
    return str(v)


def md_table(df: pd.DataFrame, index: bool = True) -> str:
    if index:
        df = df.reset_index()
    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def add_gaps(df: pd.DataFrame, key: str = 'method') -> pd.DataFrame:
    """
    gap = J / J_ref - 1 per instance, J_ref = best J among the feasible rows of
    the instance (among all rows if none is feasible); best = row with the best
    (feasible first, violations, excess, J) of its instance.
    """
    df = df.copy()
    df['gap'] = np.nan
    df['best'] = False
    for inst, sub in df.groupby('instance'):
        feas = sub[sub['feasible']]
        ref = (feas if len(feas) else sub)['J'].min()
        df.loc[sub.index, 'gap'] = sub['J'] / ref - 1 if ref > 0 else 0.0
        order = sub.sort_values(['feasible', 'n_viol', 'excess', 'J'], ascending=[False, True, True, True])
        top = order.iloc[0]
        tie = (sub['feasible'] == top['feasible']) & (sub['n_viol'] == top['n_viol']) & \
              (np.isclose(sub['J'], top['J'], rtol=1e-9))
        df.loc[sub.index[tie.values], 'best'] = True
    return df


def method_summary(df: pd.DataFrame, key: str = 'method', reference: str = None) -> list:
    """Sections of summary.md for an experiment comparing methods (one row per instance x method)."""
    df = add_gaps(df, key)
    agg = df.groupby(key, sort=False).apply(lambda s: pd.Series({
        'instances': len(s), 'feasible': int(s['feasible'].sum()), 'best': int(s['best'].sum()),
        'gap_J_feasible_%': 100 * s.loc[s['feasible'], 'gap'].mean() if s['feasible'].any() else math.nan,
        'gap_J_all_%': 100 * s['gap'].mean(), 'violations': s['n_viol'].mean(),
        'viol_L': s['viol_L'].sum(), 'viol_P': s['viol_P'].sum(), 'time_s': s['time'].mean()}),
        include_groups=False)
    sections = [("Par méthode", md_table(agg))]
    piv = df.pivot_table(index='group', columns=key, values='feasible', aggfunc='sum', sort=False)
    sections.append(("Solutions faisables par type d'instance", md_table(piv)))
    reference = reference or df[key].iloc[0]
    sections.append((f"Comparaison appariée à '{reference}' (instances faisables pour les deux)",
                     md_table(paired_tests(df, key, reference))))
    return sections


def paired_tests(df: pd.DataFrame, key: str, reference: str) -> pd.DataFrame:
    """Wilcoxon signed-rank test of J against the reference method, on the instances feasible for both."""
    from scipy.stats import wilcoxon
    ref = df[df[key] == reference].set_index('instance')
    rows = []
    for m, sub in df.groupby(key, sort=False):
        if m == reference:
            continue
        sub = sub.set_index('instance')
        both = [i for i in sub.index if i in ref.index and sub.loc[i, 'feasible'] and ref.loc[i, 'feasible']]
        rel = np.array([sub.loc[i, 'J'] / ref.loc[i, 'J'] - 1 for i in both]) * 100
        p = math.nan
        if len(rel) >= 5 and np.any(rel != 0):
            p = float(wilcoxon(rel).pvalue)
        rows.append({key: m, 'pairs': len(both), 'median_dJ_%': float(np.median(rel)) if len(rel) else math.nan,
                     'mean_dJ_%': float(rel.mean()) if len(rel) else math.nan,
                     'better': int((rel < -1e-9).sum()), 'worse': int((rel > 1e-9).sum()), 'p_wilcoxon': p})
    return pd.DataFrame(rows).set_index(key) if rows else pd.DataFrame()


def write_report(name: str, out: Path, df: pd.DataFrame, sections: list, args, files: list) -> None:
    df.to_csv(out / "results.csv", index=False)
    lines = [f"# TDA — {name}", "",
             f"Run `{out.name}` : {len(files)} instances ({args.set}), "
             f"commande `python TDA.py {' '.join(args.argv)}`.", ""]
    for title, body in sections:
        lines += [f"## {title}", "", body, ""]
    (out / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"\nRésultats : {out / 'results.csv'}\nSynthèse  : {out / 'summary.md'}")


# ----------------------------------------------------------------------------
# Solving helpers
# ----------------------------------------------------------------------------

def init_weights(graph, names) -> dict:
    """{name: weights} for names among none, capacity, sqrt (voronoi.py) and tomato."""
    std = partition.standard_inits(graph)
    out = {n: std[n] for n in names if n in std}
    if 'tomato' in names:
        out['tomato'] = partition.tomato_init_weights(graph)['weights']
    return out


def solve(graph, inits: dict, **kwargs) -> tuple:
    """(best run of partition.best_of, time)."""
    t = time.time()
    best = partition.best_of(graph, inits, **kwargs)
    return best, time.time() - t


STD = ('none', 'capacity', 'sqrt')


def baseline(graph, args) -> dict:
    """Variable distance from the 3 fixed weightings, Steiner trees (the heuristic of main.py)."""
    best, _ = solve(graph, init_weights(graph, STD), n_iter=args.n_iter)
    return best


# ----------------------------------------------------------------------------
# Experiments
# ----------------------------------------------------------------------------

def exp_descriptors(files, args, out):
    rows, dgms, groups = [], {}, {}
    for f in files:
        g = load(f)
        timings = {}
        feats, d = features.tda_features(g, spectral=True, resistance=args.resistance, timings=timings)
        scale = lmax_of(g)
        dgms[f.stem] = {k: v / (scale if k != 'tomato' else max(v[:, 0].max(), 1e-12)) for k, v in d.items()}
        groups[f.stem] = instance_group(f)
        rows.append({'instance': f.stem, 'group': groups[f.stem], **features.simple_features(g), **feats,
                     **{f"t_{k}": v for k, v in timings.items()}})
        print(f"{f.stem}: {len(feats)} descripteurs ({sum(timings.values()):.2f} s)")
    df = pd.DataFrame(rows)
    np.savez(out / "diagrams.npz", **{f"{i}__{k}": v for i, dd in dgms.items() for k, v in dd.items()})
    keys = ['h1_total_pers', 'h1_max_pers', 'h0dist_n_bars', 'h0dist_max_pers', 'h0term_total_pers',
            'tomato_n_clusters', 'lambda2', 'lb_over_lmax', 'net_over_cap', 'mean_dist_over_lmax']
    keys = [k for k in keys if k in df.columns]
    sections = [("Moyenne des descripteurs par type", md_table(df.groupby('group')[keys].mean()))]
    sep = diagram_separation(dgms, groups)
    if sep is not None:
        sections.append(("Séparation des types par les diagrammes (distance W1 moyenne, intra / inter type)",
                         md_table(sep)))
    plot_diagrams(dgms, groups, out)
    sections.append(("Figures", "`diagrams_<type de diagramme>.png` : diagrammes superposés par type d'instance "
                                "(longueurs divisées par l_max, densités par leur maximum)."))
    write_report("e1 descripteurs", out, df, sections, args, files)


def diagram_separation(dgms: dict, groups: dict, max_instances: int = 80):
    names = list(dgms)[:max_instances]
    if len(set(groups[n] for n in names)) < 2:
        return None
    rows = []
    for kind in ('h1', 'h0_dist', 'h0_term', 'tomato'):
        intra, inter = [], []
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                d = wasserstein(dgms[a][kind], dgms[b][kind], p=1, max_points=200)
                if not math.isfinite(d):
                    fa, fb = (x[np.isfinite(x[:, 1])] for x in (dgms[a][kind], dgms[b][kind]))
                    d = wasserstein(fa, fb, p=1, max_points=200)
                (intra if groups[a] == groups[b] else inter).append(d)
        rows.append({'diagramme': kind, 'W1_intra': np.mean(intra) if intra else math.nan,
                     'W1_inter': np.mean(inter) if inter else math.nan,
                     'ratio_inter_intra': (np.mean(inter) / np.mean(intra)) if intra and inter else math.nan})
    return pd.DataFrame(rows).set_index('diagramme')


def plot_diagrams(dgms: dict, groups: dict, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    types = sorted(set(groups.values()))
    for kind in ('h1', 'h0_dist', 'tomato'):
        fig, axes = plt.subplots(1, len(types), figsize=(4 * len(types), 4), squeeze=False)
        for ax, t in zip(axes[0], types):
            for name, dd in dgms.items():
                if groups[name] != t:
                    continue
                d = dd[kind]
                d = d[np.isfinite(d[:, 1])]
                ax.scatter(d[:, 0], d[:, 1], s=6, alpha=0.4)
            lim = ax.get_xlim() + ax.get_ylim()
            lo, hi = min(lim), max(lim)
            ax.plot([lo, hi], [lo, hi], color='grey', lw=0.8)
            ax.set_title(t)
            ax.set_xlabel('naissance')
            ax.set_ylabel('mort')
        fig.suptitle(f"Diagrammes {kind} par type d'instance")
        fig.tight_layout()
        fig.savefig(out / f"diagrams_{kind}.png", dpi=110)
        plt.close(fig)


def exp_conditions(files, args, out):
    rows = []
    for f in files:
        g = load(f)
        t = time.time()
        old = verify_conditions(g)
        t_old = time.time() - t
        t = time.time()
        new = bounds.necessary_conditions(g)
        t_new = time.time() - t
        ratio = max((lb / cap for lb, cap in new['budget'].values() if cap > 0), default=0.0)
        rows.append({'instance': f.stem, 'group': instance_group(f), 'proof_verify': bool(old),
                     'proof_tda': bool(new['reasons']), 'lp_status': new['lp_status'],
                     'n_reach_components': new['n_reach_components'],
                     'n_forced': sum(new['n_forced'].values()), 'max_budget_ratio': ratio,
                     'reasons_verify': " / ".join(old), 'reasons_tda': " / ".join(new['reasons']),
                     't_verify': t_old, 't_tda': t_new})
        print(f"{f.stem}: verify {'preuve' if old else '-'}, TDA {'preuve' if new['reasons'] else '-'}")
    df = pd.DataFrame(rows)
    agg = df.groupby('group')[['proof_verify', 'proof_tda']].sum()
    agg.loc['total'] = agg.sum()
    proofs = df[df['proof_tda'] | df['proof_verify']][['instance', 'reasons_verify', 'reasons_tda']]
    sections = [("Instances prouvées infaisables", md_table(agg)),
                ("Preuves", md_table(proofs, index=False) if len(proofs) else "Aucune."),
                ("Budget de longueur", "max_budget_ratio = max sur les sources de (borne des charges forcées) / "
                 "(n_departs_max x l_max) : > 1 prouve l'infaisabilité, proche de 1 signale une source tendue.\n\n"
                 + md_table(df[['instance', 'max_budget_ratio', 'n_forced', 'n_reach_components']]
                            .sort_values('max_budget_ratio', ascending=False).head(15), index=False))]
    write_report("e2 conditions nécessaires", out, df, sections, args, files)


def plne_result(stem: str) -> dict:
    p = PLNE_DIR / f"{stem}.json"
    if not p.exists():
        return {}
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def exp_bounds(files, args, out):
    rows = []
    for f in files:
        g = load(f)
        t = time.time()
        lb = bounds.lower_bound(g)
        t_lb = time.time() - t
        row = {'instance': f.stem, 'group': instance_group(f),
               'lb': lb['lb'], 'lb_mst': lb['lb_mst'], 'lb_nn': lb['lb_nn'], 't_lb': t_lb}
        best = baseline(g, args)
        row.update({'J_heur': best['cert']['J'], 'heur_feasible': best['cert']['feasible'],
                    'lb_over_heur': lb['lb'] / best['cert']['J'] if best['cert']['J'] else math.nan})
        pl = plne_result(f.stem)
        if pl:
            row.update({'plne_J': pl.get('obj'), 'plne_bound': pl.get('bound'),
                        'plne_status': pl.get('highs_status', pl.get('status')),
                        'lb_over_plne_bound': lb['lb'] / pl['bound'] if pl.get('bound') else math.nan,
                        'lb_over_plne_J': lb['lb'] / pl['obj'] if pl.get('obj') else math.nan})
        rows.append(row)
        print(f"{row['instance']}: LB = {lb['lb']:.1f}")
    df = pd.DataFrame(rows)
    cols = [c for c in ('lb_over_heur', 'lb_over_plne_bound', 'lb_over_plne_J') if c in df]
    agg = df.groupby('group')[cols].mean()
    sections = [("Rapport borne / référence (moyenne par type)", md_table(agg)),
                ("Lecture", "lb_over_heur < 1 toujours ; 1 - lb_over_heur majore l'écart relatif de l'heuristique "
                 "à l'optimum *non capacité*. lb_over_plne_bound > 1 : la borne H0 est meilleure que celle de "
                 "HiGHS à la fin de la limite de temps.")]
    write_report("e3 bornes inférieures", out, df, sections, args, files)


def exp_reductions(files, args, out):
    from run_plne_verified import solve_graph
    import copy
    rows = []
    for f in files:
        g = load(f)
        t = time.time()
        R, st = reductions.reduce_graph(g.G)
        row = {'instance': f.stem, 'group': instance_group(f), **st, 't_reduce': time.time() - t,
               'node_reduction_%': 100 * (1 - st['nodes_after'] / st['nodes_before']),
               'edge_reduction_%': 100 * (1 - st['edges_after'] / st['edges_before'])}
        if args.plne:
            res0, trees0 = solve_graph(g, args.plne, warm=False)
            gr = copy.copy(g)
            gr.G = R
            res1, trees1 = solve_graph(gr, args.plne, warm=False)
            for tag, res, trees, exp in (('orig', res0, trees0, False), ('red', res1, trees1, True)):
                cert = None
                if trees is not None:
                    trees = {s: reductions.expand_edges(R, es) for s, es in trees.items()} if exp else trees
                    cert = check_solution(g, trees)
                row.update({f'{tag}_vars': res['n_vars'], f'{tag}_status': res.get('highs_status', res['status']),
                            f'{tag}_J': cert['J'] if cert else None,
                            f'{tag}_feasible': cert['feasible'] if cert else None,
                            f'{tag}_bound': res['bound'], f'{tag}_time': res['time_solve']})
        rows.append(row)
        print(f"{f.stem}: -{row['node_reduction_%']:.1f} % nœuds, -{row['edge_reduction_%']:.1f} % arêtes")
    df = pd.DataFrame(rows)
    cols = ['node_reduction_%', 'edge_reduction_%', 'removed_deg1', 'contracted_deg2', 'parallel_dropped']
    cols += [c for c in ('orig_time', 'red_time', 'orig_bound', 'red_bound', 'orig_J', 'red_J') if c in df]
    sections = [("Réduction par type", md_table(df.groupby('group')[cols].mean()))]
    write_report("e4 réductions exactes", out, df, sections, args, files)


def method_rows(f, g, results: dict) -> list:
    """Rows {instance, group, method, time, ...cert, extras} of {method: (cert, time, extras)}."""
    rows = []
    for m, (cert, t, extra) in results.items():
        rows.append({'instance': f.stem, 'group': instance_group(f), 'method': m, 'time': t,
                     **cert_row(cert), **extra})
    return rows


def run_methods(files, args, out, name, methods_of, reference=None):
    """Generic experiment: methods_of(graph) -> {method: (cert, time, extras)}."""
    rows = []
    for f in files:
        g = load(f)
        res = methods_of(g)
        rows += method_rows(f, g, res)
        print(f"{f.stem}: " + ", ".join(f"{m} {'F' if c['feasible'] else 'inf'} J={c['J']:.1f}"
                                        for m, (c, _, _) in res.items()))
    df = pd.DataFrame(rows)
    write_report(name, out, add_gaps(df), method_summary(df, reference=reference), args, files)
    return df


def exp_init(files, args, out):
    def methods(g):
        res = {}
        for name, w in init_weights(g, STD + ('tomato',)).items():
            t = time.time()
            r = partition.variable_distance(g, w, n_iter=args.n_iter)
            res[name] = (r['cert'], time.time() - t, {'n_iter': r['n_iter'], 'iteration': r['iteration']})
        std = min((res[n] for n in STD), key=lambda x: solution_key(x[0]))
        res['best3'] = (std[0], sum(res[n][1] for n in STD), {})
        allb = min((res[n] for n in STD + ('tomato',)), key=lambda x: solution_key(x[0]))
        res['best3+tomato'] = (allb[0], res['best3'][1] + res['tomato'][1], {})
        return res
    run_methods(files, args, out, "e5 initialisation des poids", methods, reference='best3')


def exp_partition(files, args, out):
    def methods(g):
        res = {}
        inits = init_weights(g, STD)
        best, t = solve(g, inits, n_iter=args.n_iter)
        res['voronoi_longueurs'] = (best['cert'], t, {})
        t = time.time()
        attr = spectral.spectral_edge_costs(g.G, beta=args.beta)
        t_emb = time.time() - t
        best, t = solve(g, inits, weight=attr, n_iter=args.n_iter)
        res['voronoi_spectral'] = (best['cert'], t + t_emb, {})
        t = time.time()
        cells = spectral.anchored_spectral_cells(g.G)
        trees = partition.steiner_builder(g, cells)
        res['cellules_spectrales'] = (check_solution(g, trees), time.time() - t, {})
        return res
    run_methods(files, args, out, "e6 partition spectrale", methods, reference='voronoi_longueurs')


def exp_departs(files, args, out):
    def methods(g):
        res = {}
        inits = init_weights(g, STD)
        best, t = solve(g, inits, n_iter=args.n_iter)
        res['steiner'] = (best['cert'], t, {})
        for m in departs.METHODS:
            best, t = solve(g, inits, tree_builder=partition.departs_builder(m, args.tau), n_iter=args.n_iter)
            res[f'departs_{m}'] = (best['cert'], t, {})
        return res
    run_methods(files, args, out, "e7 conception des départs", methods, reference='steiner')


def ls_variants(g, trees, args, h1=None) -> dict:
    """{variant: (cert, time, extras)} of the local search variants from `trees`."""
    t = time.time()
    h1 = h1 or h1_planar(g.G)
    t_h1 = time.time() - t
    t = time.time()
    margin_loads = local_search.candidate_relocations(g, args.margin_quantile)
    t_margin = time.time() - t
    cands = {'insert_all': (local_search.candidate_insert_edges(g), None, 0.0),
             'insert_h1': (local_search.candidate_insert_edges(g, h1, args.h1_quantile), None, t_h1),
             'relocate_all': (None, local_search.candidate_relocations(g), 0.0),
             'relocate_margin': (None, margin_loads, t_margin),
             'both_all': (local_search.candidate_insert_edges(g), local_search.candidate_relocations(g), 0.0),
             'both_tda': (local_search.candidate_insert_edges(g, h1, args.h1_quantile), margin_loads,
                          t_h1 + t_margin)}
    res = {}
    for name, (ins, rel, t_cand) in cands.items():
        tr, st = local_search.local_search(g, trees, ins, rel, time_limit=args.ls_time)
        res[name] = (check_solution(g, tr), st['time'] + t_cand,
                     {'candidates': len(ins or []) + len(rel or []), 'tried': st['tried'],
                      'accepted': st['accepted_insert'] + st['accepted_relocate'], 't_candidates': t_cand})
    return res


def exp_local(files, args, out):
    def methods(g):
        best = baseline(g, args)
        res = {'sans': (best['cert'], 0.0, {'candidates': 0, 'tried': 0, 'accepted': 0, 't_candidates': 0.0})}
        res.update(ls_variants(g, best['trees'], args))
        return res
    run_methods(files, args, out, "e8 recherche locale", methods, reference='sans')


def exp_grow(files, args, out):
    def methods(g):
        inits = init_weights(g, STD)
        res = {}
        for name, grow in (('retrecit', False), ('retrecit+agrandit', True)):
            best, t = solve(g, inits, grow=grow, n_iter=args.n_iter)
            res[name] = (best['cert'], t, {'n_iter': best['n_iter']})
        return res
    run_methods(files, args, out, "e9 distance variable bidirectionnelle", methods, reference='retrecit')


def exp_grid(files, args, out):
    def methods(g):
        res = {}
        t = time.time()
        attr = spectral.spectral_edge_costs(g.G, beta=args.beta)
        t_spec = time.time() - t
        t = time.time()
        h1 = h1_planar(g.G)
        t_h1 = time.time() - t
        std = init_weights(g, STD)
        t = time.time()
        tda = {**std, 'tomato': partition.tomato_init_weights(g)['weights']}
        t_tom = time.time() - t
        for init_name, inits, t_init in (('std', std, 0.0), ('tda', tda, t_tom)):
            for part, weight, t_part in (('graph', 'weight', 0.0), ('spectral', attr, t_spec)):
                for dep, builder in (('steiner', partition.steiner_builder),
                                     ('mergetree', partition.departs_builder('mergetree', args.tau))):
                    best, t = solve(g, inits, weight=weight, tree_builder=builder, n_iter=args.n_iter)
                    key = f"init={init_name}|dist={part}|dep={dep}"
                    t_base = t + t_init + t_part
                    res[key + "|ls=none"] = (best['cert'], t_base, {})
                    ls = ls_variants(g, best['trees'], args, h1)['both_tda']
                    res[key + "|ls=tda"] = (ls[0], t_base + ls[1] + t_h1, {})
        return res
    return run_methods(files, args, out, "grid ablation", methods,
                       reference="init=std|dist=graph|dep=steiner|ls=none")


def exp_select(files, args, out, labels: Path = None):
    labels = Path(labels or args.labels)
    lab = add_gaps(pd.read_csv(labels))
    best = lab[lab['best']].groupby('instance').first()['method']
    best_set = lab[lab['best']].groupby('instance')['method'].apply(set)   # ties: several best methods
    by_name = {f.stem: f for f in files}
    names = [n for n in best.index if n in by_name]
    if not names:
        raise SystemExit(f"Aucune instance de {labels} dans la sélection courante")
    feats_simple, feats_tda, groups = [], [], []
    for n in names:
        g = load(by_name[n])
        feats_simple.append(features.simple_features(g))
        feats_tda.append(features.tda_features(g, spectral=True)[0])
        groups.append(instance_group(by_name[n]))
        print(f"{n}: descripteurs calculés")
    S, T = pd.DataFrame(feats_simple), pd.DataFrame(feats_tda)
    y = best.loc[names].values
    rows, preds = [], {}
    for set_name, X in (('simples', S), ('TDA', T), ('simples+TDA', pd.concat([S, T], axis=1))):
        r = features.knn_select(X.values, list(y), groups, k=args.k)
        preds[set_name] = r['predictions']
        gaps, feas = [], []
        for n, p in zip(names, r['predictions']):
            row = lab[(lab['instance'] == n) & (lab['method'] == p)]
            if len(row):
                gaps.append(float(row['gap'].iloc[0]))
                feas.append(bool(row['feasible'].iloc[0]))
        hits = [p in best_set[n] for n, p in zip(names, r['predictions']) if p is not None]
        rows.append({'descripteurs': set_name, 'n_features': X.shape[1],
                     'accuracy': float(np.mean(hits)) if hits else math.nan,
                     'majority_accuracy': r['majority_accuracy'],
                     'regret_J_%': 100 * np.mean(gaps) if gaps else math.nan,
                     'feasible_predicted': int(np.sum(feas))})
    res = pd.DataFrame(rows)
    df = pd.DataFrame({'instance': names, 'group': groups, 'best': y,
                       **{f"pred_{k}": v for k, v in preds.items()}})
    oracle = int(lab[lab['best']].groupby('instance')['feasible'].first().loc[names].sum())
    sections = [("Sélection de configuration (k-NN, une famille d'instances laissée de côté)", md_table(res, index=False)),
                ("Lecture", f"Étiquettes : {labels}. Meilleure configuration faisable pour {oracle}/{len(names)} "
                 "instances (oracle). accuracy = part des instances où la configuration prédite est l'une des "
                 "meilleures (ex aequo compris) ; majority_accuracy = toujours prédire l'étiquette la plus "
                 "fréquente (référence, ex aequo non comptés). regret_J = écart moyen de J de la configuration "
                 "prédite à la meilleure.")]
    write_report("e10 sélection de configuration", out, df, sections, args, files)


def exp_scaling(files, args, out):
    rows, h0 = [], {}
    for f in files:
        g = load(f)
        timings = {}
        feats, d = features.tda_features(g, spectral=True, resistance=args.resistance, timings=timings)
        t = time.time()
        bounds.necessary_conditions(g)
        timings['conditions'] = time.time() - t
        t = time.time()
        reductions.reduce_graph(g.G)
        timings['reductions'] = time.time() - t
        t = time.time()
        departs.departs_by_cells(g, partition.variable_distance(g, partition.standard_inits(g)['none'],
                                                                n_iter=1)['cells'])
        timings['departs_mergetree'] = time.time() - t
        n_nodes, n_edges, n_src = instance_size(f) if f.parent != INSTANCES_DIR else \
            (g.G.number_of_nodes(), g.G.number_of_edges(), None)
        h0[f.stem] = (d['h0_term'] / lmax_of(g), d['tomato'] / max(d['tomato'][:, 0].max(), 1e-12))
        rows.append({'instance': f.stem, 'group': instance_group(f), 'n_nodes': n_nodes, 'n_edges': n_edges,
                     **{f"t_{k}": v for k, v in timings.items()}})
        print(f"{f.stem}: {g.G.number_of_nodes()} nœuds, {sum(timings.values()):.2f} s")
    df = pd.DataFrame(rows)
    tcols = [c for c in df.columns if c.startswith('t_')]
    df['size_class'] = df['group'].str.split('_').str[0]
    sections = [("Temps moyen (s) par classe de taille", md_table(df.groupby('size_class')[['n_nodes'] + tcols].mean()))]
    pairs = []
    for name, (dh0, dtom) in h0.items():
        if 'Ete' in name:
            other = name.replace('Ete', 'Hiver')
            if other in h0:
                pairs.append({'paire': name.replace('_Ete', ''),
                              'W1_h0_charges': wasserstein(dh0, h0[other][0], p=1),
                              'W1_tomato': wasserstein(dtom[np.isfinite(dtom[:, 1])],
                                                       h0[other][1][np.isfinite(h0[other][1][:, 1])], p=1)})
    if pairs:
        sections.append(("Distance été / hiver entre diagrammes", md_table(pd.DataFrame(pairs), index=False)))
    write_report("e11 passage à l'échelle", out, df, sections, args, files)


EXPERIMENTS = {'e1': exp_descriptors, 'e2': exp_conditions, 'e3': exp_bounds, 'e4': exp_reductions,
               'e5': exp_init, 'e6': exp_partition, 'e7': exp_departs, 'e8': exp_local, 'e9': exp_grow,
               'e10': exp_select, 'e11': exp_scaling, 'grid': exp_grid}
NAMES = {'e1': 'descriptors', 'e2': 'conditions', 'e3': 'bounds', 'e4': 'reductions', 'e5': 'init',
         'e6': 'partition', 'e7': 'departs', 'e8': 'local', 'e9': 'grow', 'e10': 'select',
         'e11': 'scaling', 'grid': 'grid'}


def out_dir(exp: str, label: str, stamp: str) -> Path:
    out = TDA_DIR / f"{stamp}_{exp}_{NAMES[exp]}{'_' + label if label else ''}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def main(argv=None):
    import sys
    p = argparse.ArgumentParser(description="Expériences TDA (voir le docstring de TDA.py)")
    p.add_argument("experiment", choices=list(EXPERIMENTS) + ['all'])
    p.add_argument("motifs", nargs="*", help="ne garder que les instances dont le nom contient un de ces motifs")
    p.add_argument("--set", choices=('gen', 'verified'), default='gen')
    p.add_argument("--petites", action="store_true", help="g01 à g20")
    p.add_argument("--grandes", action="store_true", help="g21 à g40")
    p.add_argument("--prefix", default="", help="verified : Petit, Moyen, Grand ou TrèsGrand")
    p.add_argument("--max-nodes", type=int, help="verified : au plus ce nombre de nœuds")
    p.add_argument("--n", type=int, help="seulement les n premières instances")
    p.add_argument("--label", help="suffixe du dossier de résultats")
    p.add_argument("--n-iter", type=int, default=60, help="itérations de la distance variable")
    p.add_argument("--tau", type=float, default=0.0, help="simplification de l'arbre de fusion (km)")
    p.add_argument("--beta", type=float, default=1.0, help="poids du saut spectral dans le coût des arêtes")
    p.add_argument("--h1-quantile", type=float, default=None,
                   help="fraction des arêtes H1 les plus persistantes (défaut : plus grand écart)")
    p.add_argument("--margin-quantile", type=float, default=0.3, help="fraction des charges de plus petite marge")
    p.add_argument("--ls-time", type=float, default=60.0, help="limite de temps de la recherche locale (s)")
    p.add_argument("--resistance", action="store_true", help="e1/e11 : résistance effective (plus lent)")
    p.add_argument("--plne", type=float, help="e4 : temps limite de la PLNE avec / sans réduction")
    p.add_argument("--labels", help="e10 : results.csv d'une expérience (colonnes instance, method, feasible, J)")
    p.add_argument("--k", type=int, default=3, help="e10 : nombre de voisins")
    args = p.parse_args(argv)
    args.argv = sys.argv[1:] if argv is None else list(argv)
    files = select_instances(args)
    if not files:
        raise SystemExit("Aucune instance sélectionnée")
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    exps = [e for e in EXPERIMENTS if e not in ('e10', 'e11', 'grid')] + ['grid'] \
        if args.experiment == 'all' else [args.experiment]
    if args.experiment == 'e10' and not args.labels:
        raise SystemExit("e10 demande --labels <results.csv>")
    print(f"{len(files)} instances, expériences {exps}")
    grid_df = None
    for e in exps:
        print(f"\n=== {e} ({NAMES[e]}) ===")
        res = EXPERIMENTS[e](files, args, out_dir(e, args.label, stamp))
        if e == 'grid':
            grid_df = res
    if args.experiment == 'all' and grid_df is not None:
        print("\n=== e10 (select) sur les résultats de grid ===")
        exp_select(files, args, out_dir('e10', args.label, stamp),
                   labels=out_dir('grid', args.label, stamp) / "results.csv")


if __name__ == "__main__":
    main()
