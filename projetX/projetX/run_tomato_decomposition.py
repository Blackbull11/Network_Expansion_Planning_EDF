"""
ToMATo decomposition (TDA_algos/decomposition.py) of one instance into smaller
subproblems, solved separately and put back together. Compared with the same
solver on the whole instance, and with the same merging rules started from the
plain Voronoi cells of the sources (control: what does ToMATo itself bring?).

Written to output/tomato_decomposition/<instance>/ :
  report.md             results and figures
  results.csv           one line per configuration
  solution.json         trees of the best decomposition {source: edges}
  persistence.png       persistence diagram with the bootstrap noise band
  regions.png           regions and trees of the best decomposition
  subproblems/*.json    the regions of the best decomposition as standalone
                        instances (readable by src.instances.load_instance)

Usage:
  python run_tomato_decomposition.py [motif of the instance name]
  python run_tomato_decomposition.py TrèsGrand_N2000_S15_0.25_Hiver --workers 12 --ls-time 20
"""
import argparse
import csv
import json
import logging
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from itertools import product
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.collections import LineCollection

from output import INK, INK_2, GRID, SURFACE, _new_axes, _save
from run_voronoi_verified import PALETTE, graph_area, region_polygons
from src.instances import load_instance
from TDA_algos.decomposition import (bootstrap_tau, decompose, export_subinstances, model_size, prominences,
                                     solve_regions, solve_subproblem, voronoi_control)
from verify import check_solution, necessary_conditions, solution_key

HERE = Path(__file__).resolve().parent
INSTANCE_DIR = HERE / "data" / "verified_instances"
OUT_ROOT = HERE / "output" / "tomato_decomposition"
BANDWIDTHS = (3.0, 5.0, 7.0, 10.0)
MASSES = ('power', 'count')
MAX_SOURCES = (1, 2, 3)

_graph = None


def _init(path):
    global _graph
    logging.disable(logging.INFO)
    _graph = load_instance(path)


def summary(cert: dict) -> dict:
    return {'J': cert['J'], 'PD': cert['PD'], 'feasible': cert['feasible'], 'n_errors': len(cert['errors']),
            'n_violations': sum(v[0] for v in cert['violations'].values()), 'excess': cert['excess'],
            **{f"viol_{k}": v[0] for k, v in cert['violations'].items()}}


def run_config(cfg: dict, n_iter: int, ls_time: float, polish_time: float) -> dict:
    """One configuration, in a worker: {cfg, ...metrics, trees, regions}."""
    graph = _graph
    t0 = time.perf_counter()
    if cfg['kind'] == 'full':
        res = solve_subproblem(graph, n_iter, ls_time)
        wall = time.perf_counter() - t0
        return {'cfg': cfg, 'first': summary(res['cert']), 'best': summary(res['cert']), 'trees': res['trees'],
                'regions': None, 'n_regions': 1, 'n_atoms': None, 'time_cpu': wall, 'time_parallel': wall, 'stats': None}
    if cfg['kind'] == 'tomato':
        dec = decompose(graph, cfg['h'], cfg['mass'], max_sources=cfg['max_sources'])
    else:
        dec = voronoi_control(graph, cfg['max_sources'])
    t_dec = time.perf_counter() - t0
    sol = solve_regions(graph, dec['regions'], n_iter, ls_time, polish_time)
    return {'cfg': cfg, 'first': summary(sol['cert']), 'best': summary(sol['polished_cert']),
            'trees': sol['polished_trees'], 'regions': dec['regions'], 'n_regions': len(dec['regions']),
            'n_atoms': dec['n_atoms'], 'unresolved': dec['unresolved'], 'stats': dec['stats'],
            'time_cpu': t_dec + sol['time_sum'] + sol['time_polish'],
            'time_parallel': t_dec + sol['time_max'] + sol['time_polish'], 'per_region': sol['regions']}


def rank_key(r: dict) -> tuple:
    b = r['best']
    return (b['n_errors'], b['n_violations'], b['excess'], b['J'])


def label_of(cfg: dict) -> str:
    if cfg['kind'] == 'full':
        return "whole instance"
    if cfg['kind'] == 'voronoi':
        return f"Voronoi control, ≤{cfg['max_sources']} src"
    return f"ToMATo h={cfg['h']:g} {cfg['mass']}, ≤{cfg['max_sources']} src"


# ---------------------------------------------------------------- figures

def plot_diagram(boot: dict, h: float, mass: str, path: Path):
    dgm = boot['dgm']
    fin = dgm[np.isfinite(dgm[:, 1])]
    f, ax = _new_axes(width=6, height=5.4)
    top = dgm[:, 0].max() * 1.05
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.plot([0, top], [0, top], color=INK_2, lw=1, zorder=2)
    ax.fill_between([0, top], [0, top], [-boot['tau'], top - boot['tau']], color=GRID, alpha=0.8, zorder=1, linewidth=0)
    ax.scatter(fin[:, 0], fin[:, 1], s=34, color="#2a78d6", edgecolors=SURFACE, linewidths=1, zorder=3)
    ax.annotate(f"shaded band = noise: prominence < τ = {boot['tau']:.1f}", (top * 0.03, top * 0.95),
                ha='left', va='center', fontsize=9, color=INK_2)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_xlabel("density at the peak (birth)", color=INK_2, fontsize=9)
    ax.set_ylabel("density at the merge (death)", color=INK_2, fontsize=9)
    _save(f, ax, f"Persistence diagram, h = {h:g} km, {mass} density: "
                 f"{boot['n_significant']} significant peak(s) of {len(fin) + 1}", path)


def plot_regions(graph, regions: dict, trees: dict, title: str, path: Path):
    G = graph.G
    owner = {v: i for i, ns in regions.items() for v in ns}
    touch = nx.Graph()
    touch.add_nodes_from(regions)
    touch.add_edges_from((owner[u], owner[v]) for u, v in G.edges if owner[u] != owner[v])
    colors = {i: PALETTE[c % len(PALETTE)] for i, c in nx.greedy_color(touch, strategy='largest_first').items()}
    pos = dict(G.nodes('pos'))
    f, ax = plt.subplots(figsize=(11, 11))
    f.patch.set_facecolor(SURFACE)
    for i, poly in region_polygons(G, regions, graph_area(G)).items():
        if not poly.is_empty:
            gpd.GeoSeries([poly]).plot(ax=ax, color=colors[i], alpha=0.22, linewidth=0, zorder=1)
    ax.add_collection(LineCollection([(pos[u], pos[v]) for u, v in G.edges], colors=[GRID], linewidths=0.3, zorder=2))
    ax.add_collection(LineCollection([(pos[u], pos[v]) for es in trees.values() for u, v in es], colors=[INK],
                                     linewidths=1.3, zorder=3))
    loads = [v for v, l in G.nodes('label') if l in ('Conso', 'Prod')]
    sources = [v for v, l in G.nodes('label') if l == 'Source']
    ax.scatter(*np.array([pos[v] for v in loads]).T, s=12, c=[colors[owner[v]] for v in loads], edgecolors=INK,
               linewidths=0.3, zorder=4)
    ax.scatter(*np.array([pos[v] for v in sources]).T, s=90, marker='s', c=[colors[owner[v]] for v in sources],
               edgecolors=INK, linewidths=1.3, zorder=5)
    for i, ns in regions.items():
        x, y = np.mean([pos[v] for v in ns], axis=0)
        ax.annotate(str(i), (x, y), ha='center', va='center', fontsize=12, fontweight='bold', color=INK, zorder=6,
                    bbox=dict(boxstyle='round,pad=0.2', fc=SURFACE, ec='none', alpha=0.85))
    ax.set_aspect('equal')
    ax.autoscale_view()
    ax.axis('off')
    ax.set_title(title, color=INK, fontsize=11, loc='left')
    f.savefig(path, dpi=130, bbox_inches='tight', facecolor=SURFACE)
    plt.close(f)


# ---------------------------------------------------------------- report

def write_csv(results: list, path: Path):
    cols = ['configuration', 'n_regions', 'n_atoms', 'J', 'feasible', 'n_violations', 'viol_L', 'viol_P', 'viol_source',
            'viol_n', 'excess', 'J_before_polish', 'violations_before_polish', 'max_region_loads', 'max_region_vars',
            'time_cpu_s', 'time_if_parallel_s']
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter=';')
        w.writerow(cols)
        for r in results:
            st = r['stats'] or {}
            b = r['best']
            w.writerow([label_of(r['cfg']), r['n_regions'], r['n_atoms'], f"{b['J']:.2f}", b['feasible'], b['n_violations'],
                        b['viol_L'], b['viol_P'], b['viol_source'], b['viol_n'], f"{b['excess']:.4f}",
                        f"{r['first']['J']:.2f}", r['first']['n_violations'],
                        max((s['n_loads'] for s in st.values()), default=''), max((s['n_vars'] for s in st.values()), default=''),
                        f"{r['time_cpu']:.1f}", f"{r['time_parallel']:.1f}"])


def row(r: dict) -> str:
    b, st = r['best'], r['stats'] or {}
    max_loads = max((s['n_loads'] for s in st.values()), default='—')
    return (f"| {label_of(r['cfg'])} | {r['n_regions']} | {max_loads} | {b['J']:.1f} | {b['n_violations']} "
            f"| {b['viol_L']} / {b['viol_P']} / {b['viol_source']} / {b['viol_n']} | {r['time_cpu']:.0f} "
            f"| {r['time_parallel']:.0f} |")


def write_report(out: Path, stem: str, graph, results: list, boots: list, best: dict, full: dict, proofs: list, a):
    G = graph.G
    n_load = sum(l in ('Conso', 'Prod') for _, l in G.nodes('label'))
    n_feeders = int(sum(graph.nodes.loc[s, 'n_departs_max'] for s, l in G.nodes('label') if l == 'Source'))
    full_vars = model_size(n_load, G.number_of_edges(), len(G), n_feeders)
    head = ("| Configuration | Regions | Max loads in a region | J (km) | Violations | L / P / source / n | CPU (s) "
            "| Time if regions run in parallel (s) |\n|---|---|---|---|---|---|---|---|")
    st = best['stats']
    md = [f"# ToMATo decomposition of {stem}", "",
          f"{datetime.now():%d/%m/%Y %H:%M} · n_iter = {a.n_iter}, local search {a.ls_time:g} s per region "
          f"({a.ls_time * a.full_ls_factor:g} s for the whole instance), polish {a.polish_time:g} s.", "",
          "## Instance", "",
          f"{len(G)} nodes, {G.number_of_edges()} edges, {n_load} loads, {n_feeders} feeders. "
          f"Necessary conditions of verify.py: {'; '.join(proofs) if proofs else 'none violated'}.", "",
          "## Is there a natural clustering?", "",
          "Bottleneck bootstrap of AuToMATo on the persistence diagrams (α = 0.35, "
          f"B = {a.n_boot}, loads resampled with replacement):", "",
          "| h (km) | Density | Peaks | Largest finite prominence | τ bootstrap | Significant peaks |", "|---|---|---|---|---|---|"]
    for b in boots:
        p = prominences(b['boot']['dgm'])
        md.append(f"| {b['h']:g} | {b['mass']} | {len(p) + 1} | {p[0]:.2f} | {b['boot']['tau']:.2f} "
                  f"| {b['boot']['n_significant']} |")
    md += ["", "![Persistence diagram](persistence.png)", "",
           "## Results", "", "Ranked by verify.py on the whole graph: errors, violations, excess, then J.", "", head]
    md += [row(r) for r in results]
    md += ["", "## Best decomposition", "", f"**{label_of(best['cfg'])}**, against the whole-instance run: "
           f"J = {best['best']['J']:.1f} km with {best['best']['n_violations']} violation(s), "
           f"against J = {full['best']['J']:.1f} km with {full['best']['n_violations']} violation(s).", "",
           f"Full PLNE model (plne_reseau.build_model): {full_vars:,} variables. Largest subproblem: "
           f"{max(s['n_vars'] for s in st.values()):,}; sum over the subproblems: {sum(s['n_vars'] for s in st.values()):,}.", "",
           "![Regions](regions.png)", "",
           "| Region | Nodes | Edges | Sources | Loads | Net power (MW) | Usage | Variables of the PLNE | J (km) | Violations |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for i, s in st.items():
        pr = best['per_region'][i]
        md.append(f"| {i} | {s['n_nodes']} | {s['n_edges']} | {s['n_sources']} | {s['n_loads']} | {s['power']:g} "
                  f"| {s['usage']:.2f} | {s['n_vars']:,} | {pr['J']:.1f} | {sum(v[0] for v in pr['violations'].values())} |")
    (out / "report.md").write_text("\n".join(md) + "\n", encoding='utf-8')


def to_builtin(o):
    return o.item() if hasattr(o, 'item') else str(o)


def main():
    p = argparse.ArgumentParser(description="ToMATo decomposition of an instance")
    p.add_argument("motif", nargs="?", default="TrèsGrand_N2000_S15_0.25_Hiver")
    p.add_argument("--workers", type=int, default=12)
    p.add_argument("--n-iter", type=int, default=60)
    p.add_argument("--ls-time", type=float, default=20.0, help="local search per region (s)")
    p.add_argument("--full-ls-factor", type=int, default=15, help="whole instance: local search of ls-time x this")
    p.add_argument("--polish-time", type=float, default=60.0, help="final local search on the whole graph (s)")
    p.add_argument("--n-boot", type=int, default=200)
    a = p.parse_args()
    logging.disable(logging.INFO)
    matches = sorted(f for f in INSTANCE_DIR.glob(f"*{a.motif}*.json") if ' (' not in f.name)
    if len(matches) != 1:
        p.error(f"{len(matches)} instances match {a.motif!r}: {[m.name for m in matches]}")
    path = matches[0]
    out = OUT_ROOT / path.stem
    out.mkdir(parents=True, exist_ok=True)
    graph = load_instance(path)
    proofs = necessary_conditions(graph)

    cfgs = [{'kind': 'tomato', 'h': h, 'mass': m, 'max_sources': k} for h, m, k in product(BANDWIDTHS, MASSES, MAX_SOURCES)]
    cfgs += [{'kind': 'voronoi', 'max_sources': k} for k in MAX_SOURCES] + [{'kind': 'full'}]
    results = []
    t0 = time.time()
    # load_instance rewrites a file that is not normalized: the workers read the copy made above
    with ProcessPoolExecutor(a.workers, initializer=_init, initargs=(graph.filepath,)) as pool:
        futures = {pool.submit(run_config, c, a.n_iter,
                               a.ls_time * a.full_ls_factor if c['kind'] == 'full' else a.ls_time, a.polish_time): c
                   for c in cfgs}
        for fut in as_completed(futures):
            r = fut.result()
            results.append(r)
            b = r['best']
            print(f"[{len(results)}/{len(cfgs)}] {label_of(r['cfg'])}: regions={r['n_regions']} J={b['J']:.1f} "
                  f"violations={b['n_violations']} (L {b['viol_L']}, P {b['viol_P']}, src {b['viol_source']}, n {b['viol_n']}) "
                  f"cpu={r['time_cpu']:.0f}s  [{(time.time() - t0) / 60:.1f} min]", flush=True)
    results.sort(key=rank_key)
    full = next(r for r in results if r['cfg']['kind'] == 'full')
    best = next(r for r in results if r['cfg']['kind'] == 'tomato')
    check = check_solution(graph, best['trees'])   # recomputed in this process, on the whole graph
    assert abs(check['J'] - best['best']['J']) < 1e-6, "the certificate of the worker is not reproduced"

    boots = [{'h': h, 'mass': m, 'boot': bootstrap_tau(graph.G, h, m, n_boot=a.n_boot)} for h, m in product(BANDWIDTHS, MASSES)]
    h, mass = best['cfg']['h'], best['cfg']['mass']
    plot_diagram(next(b['boot'] for b in boots if (b['h'], b['mass']) == (h, mass)), h, mass, out / "persistence.png")
    plot_regions(graph, best['regions'], best['trees'],
                 f"{path.stem}: {label_of(best['cfg'])}. {best['n_regions']} regions, J = {best['best']['J']:.1f} km, "
                 f"{best['best']['n_violations']} violation(s)", out / "regions.png")
    export_subinstances(path, best['regions'], out / "subproblems")
    with open(out / "solution.json", 'w', encoding='utf-8') as fh:
        json.dump({'configuration': label_of(best['cfg']), 'certificate': best['best'],
                   'trees': {str(s): es for s, es in best['trees'].items()}}, fh, default=to_builtin)
    write_csv(results, out / "results.csv")
    write_report(out, path.stem, graph, results, boots, best, full, proofs, a)
    print(f"\nBest: {label_of(best['cfg'])}: {best['best']}\nWhole instance: {full['best']}\nReport: {out / 'report.md'}", flush=True)


if __name__ == "__main__":
    main()
