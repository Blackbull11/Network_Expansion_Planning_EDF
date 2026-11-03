"""
Heuristique Voronoi + Steiner (les 4 distances de main.py : non pondérée,
capacité, racine capacité, variable) sur les instances de verified_instances/,
avec un rapport compact. Tout est écrit dans output/voronoi_verified/ :

  results.jsonl         une ligne par instance, écrite après chaque instance :
                        un run interrompu reprend là où il s'est arrêté
  results.csv           une ligne par instance et par distance
  figures/<inst>.png    une image par instance, un panneau par distance :
                        cellules de Voronoi (zones colorées), arbres de Steiner
                        (en gras), capacité des sources (MW), aire des charges
                        proportionnelle à |demande| ; J, temps et faisabilité
                        dans le titre du panneau
  charts/*.png          graphiques de synthèse
  report.md             tableaux et graphiques (aperçu Markdown de VS Code)

Temps mesuré : Voronoi + Steiner (+ ajustement des poids pour la distance
variable), sans la vérification ni le dessin.

Usage :
  python run_voronoi_verified.py                  toutes les instances (reprise si results.jsonl existe)
  python run_voronoi_verified.py Petit_ Moyen_    seulement les instances dont le nom contient un de ces motifs
  python run_voronoi_verified.py --max-nodes 250  seulement les instances d'au plus 250 nœuds
  python run_voronoi_verified.py --restart        efface les résultats précédents et recommence
  python run_voronoi_verified.py --report         refait seulement le rapport à partir de results.jsonl
  python run_voronoi_verified.py --no-figures     sans les images par instance
  python run_voronoi_verified.py --dir <dossier>  instances lues dans un autre dossier
  python run_voronoi_verified.py --out <dossier>  résultats écrits dans un autre dossier
"""
import argparse
import csv
import json
import logging
import re
import time
import traceback
from datetime import datetime
from math import sqrt
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from shapely import MultiPoint, voronoi_polygons
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

from main import solve, evaluate
from nx_graph import cell_colors
from output import (METHODS, WEIGHTING_NAMES, WEIGHTING_COLORS, GREEN, RED, INK, INK_2, GRID, SURFACE,
                    SEQ_BLUE, colored, md_weighting, _new_axes, _save)
from src.instances import SIZE_CLASSES, VERIFIED_DIR, instance_size, is_instance, load_instance
from verify import necessary_conditions

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "output" / "voronoi_verified"
# le dossier des instances a été à ces deux endroits : le premier qui contient des .json
INSTANCE_DIRS = (VERIFIED_DIR, HERE / "verified_instances", HERE / "data" / "verified_instances")
DEFAULT_DIR = next((d for d in INSTANCE_DIRS if any(d.glob("*.json"))), VERIFIED_DIR)
NAME_RE = re.compile(r"^(?P<classe>[^_]+)_N(?P<N>\d+)_S(?P<S>\d+)_(?P<ratio>[\d.]+)_(?P<saison>[^_ ]+)")
CONSTRAINTS = {'L': "L_max départ", 'P': "P_max départ", 'source': "capacité source",
               'n': "nb départs", 'errors': "structure"}
LOAD_LABELS = ('Conso', 'Prod')


# ---------------------------------------------------------------- instances

def instance_files(folder: Path, motifs: list, max_nodes: int = None) -> list:
    """Instances JSON de `folder` (nom commençant par une classe de SIZE_CLASSES), triées par taille."""
    files = [f for f in folder.glob("*.json") if f.stem.split('_')[0] in SIZE_CLASSES and is_instance(f)]
    files = [f for f in files if not motifs or any(m in f.stem for m in motifs)]
    sizes = {f: instance_size(f) for f in files}
    if max_nodes:
        files = [f for f in files if sizes[f][0] <= max_nodes]
    return sorted(files, key=lambda f: (sizes[f], f.name))


def name_params(stem: str) -> dict:
    """Paramètres lus dans le nom : classe, S nominal, ratio, saison (None si le nom n'a pas ce format)."""
    m = NAME_RE.match(stem)
    if not m:
        return {'classe': stem.split('_')[0], 'S_nom': None, 'ratio': None, 'saison': None}
    return {'classe': m['classe'], 'S_nom': int(m['S']), 'ratio': float(m['ratio']), 'saison': m['saison']}


def run_instance(path: Path, figures_dir: Path = None) -> dict:
    """Les METHODS sur l'instance `path` : enregistrement de results.jsonl, image dans figures_dir."""
    graph = load_instance(path)
    G = graph.G
    loads = [v for v, lab in G.nodes('label') if lab in LOAD_LABELS]
    rec = {'instance': path.stem, **name_params(path.stem), 'n_nodes': len(G), 'n_edges': G.number_of_edges(),
           'n_sources': sum(lab == 'Source' for _, lab in G.nodes('label')), 'n_loads': len(loads),
           'n_prod': sum(G.nodes[v]['conso_ete'] < 0 for v in loads),
           'proofs': necessary_conditions(graph), 'methods': {}, 'figure': None}
    snaps = {}
    for method in METHODS:
        t0 = time.perf_counter()
        solve(graph, method)
        seconds = time.perf_counter() - t0
        ev = evaluate(graph)
        rec['methods'][method] = {'J': ev['J'], 'PD': ev['PD'], 'feasible': ev['feasible'], 'time': seconds,
                                  'violations': {k: list(v) for k, v in ev['violations'].items()},
                                  'errors': [str(e) for e in ev['errors']]}
        snaps[method] = snapshot(graph, ev, seconds)
    if figures_dir:
        filepath = figures_dir / f"{path.stem}.png"
        draw_instance(graph, snaps, path.stem, filepath)
        rec['figure'] = f"{figures_dir.name}/{filepath.name}"
    return rec


def n_violations(m: dict) -> int:
    return sum(v[0] for v in m['violations'].values()) + len(m['errors'])


# ---------------------------------------------------------------- figure par instance

PALETTE = [c for i, c in enumerate(plt.get_cmap('tab10').colors) if i != 7]   # tab10 sans le gris
ZERO_AREA = 14   # aire du marqueur des charges de demande nulle


def snapshot(graph, ev: dict, seconds: float) -> dict:
    """Ce dont draw_instance a besoin, une fois une distance résolue et vérifiée."""
    return {'cells': {k: set(v) for k, v in graph.cells.items()},
            'edge_cells': {k: list(v) for k, v in graph.edge_cells.items() if k != 'split'},
            'pos': dict(graph.G_split.nodes('pos')), 'trees': graph.steiner_edges,
            'J': ev['J'], 'feasible': ev['feasible'], 'time': seconds,
            'n_violations': sum(v[0] for v in ev['violations'].values()) + len(ev['errors'])}


def neighbour_colors(G, snaps: dict, sources: list) -> dict:
    """{source: couleur} : coloration gloutonne des cellules voisines dans au moins un panneau (cell_colors sinon)."""
    touch = nx.Graph()
    touch.add_nodes_from(sources)
    for snap in snaps.values():
        cell = {n: k for k, ns in snap['cells'].items() for n in ns}
        touch.add_edges_from((cell[u], cell[v]) for u, v in G.edges
                             if cell.get(u) in touch and cell.get(v) in touch and cell[u] != cell[v])
    color = nx.greedy_color(touch, strategy='largest_first')
    if max(color.values(), default=0) >= len(PALETTE):
        return cell_colors(sources)
    return {s: PALETTE[c] for s, c in color.items()}


def graph_area(G):
    """Bande autour des arêtes de G (demi-longueur médiane d'arête), sans ses trous."""
    lines = [LineString([G.nodes[u]['pos'], G.nodes[v]['pos']]) for u, v in G.edges]
    width = 0.5 * float(np.median([d for _, _, d in G.edges.data('weight')]) or 1)
    band = unary_union(lines).buffer(width, quad_segs=4)
    return unary_union([Polygon(p.exterior) for p in getattr(band, 'geoms', [band])])


def region_polygons(G, cells: dict, area) -> dict:
    """{cellule: polygone} union des polygones de Voronoi plans de ses nœuds, coupée à `area`."""
    nodes = list(G)
    polys = voronoi_polygons(MultiPoint([G.nodes[n]['pos'] for n in nodes]),
                             extend_to=area.envelope, ordered=True).geoms
    cell = {n: k for k, ns in cells.items() for n in ns}
    by_cell = {}
    for n, poly in zip(nodes, polys):
        by_cell.setdefault(cell.get(n), []).append(poly)
    return {k: unary_union(ps).intersection(area) for k, ps in by_cell.items() if k is not None}


def draw_panel(ax, G, snap: dict, colors: dict, area, sizes: dict, src_font: float):
    """Une distance : zones, arêtes par cellule, arbre de Steiner, sources et charges."""
    pos = snap['pos']
    for k, poly in region_polygons(G, snap['cells'], area).items():
        if not poly.is_empty:
            gpd.GeoSeries([poly]).plot(ax=ax, color=colors.get(k, 'r'), alpha=0.16, linewidth=0, zorder=1)
    tree = {e for es in snap['trees'].values() for e in es}
    for k, edges in snap['edge_cells'].items():
        thin = [(pos[u], pos[v]) for u, v in edges if (min(u, v), max(u, v)) not in tree]
        bold = [(pos[u], pos[v]) for u, v in edges if (min(u, v), max(u, v)) in tree]
        ax.add_collection(LineCollection(thin, colors=[colors[k]], linewidths=0.4, alpha=0.45, zorder=2))
        ax.add_collection(LineCollection(bold, colors=[colors[k]], linewidths=1.6, zorder=3))
    cell = {n: k for k, ns in snap['cells'].items() for n in ns}
    color_of = lambda n: colors.get(cell.get(n), 'r')
    by_label = {}
    for n, lab in G.nodes('label'):
        by_label.setdefault(lab, []).append(n)
    if by_label.get('Transit'):
        ax.scatter(*np.array([pos[n] for n in by_label['Transit']]).T, s=2, color=INK_2, alpha=0.5,
                   linewidths=0, zorder=4)
    for lab, marker in (('Conso', 'o'), ('Prod', '^')):
        ns = [n for n in by_label.get(lab, []) if sizes[n] > 0]
        if ns:
            ax.scatter(*np.array([pos[n] for n in ns]).T, s=[sizes[n] for n in ns], marker=marker,
                       c=[color_of(n) for n in ns], edgecolors=INK, linewidths=0.4, alpha=0.9, zorder=5)
    zero = [n for n in sizes if sizes[n] == 0]
    if zero:
        ax.scatter(*np.array([pos[n] for n in zero]).T, s=ZERO_AREA, marker='o', facecolors=SURFACE,
                   edgecolors=[color_of(n) for n in zero], linewidths=1, zorder=5)
    for s in by_label.get('Source', []):
        x, y = pos[s]
        ax.scatter([x], [y], s=90, marker='s', color=color_of(s), edgecolors=INK, linewidths=1.2, zorder=6)
        ax.annotate(f"{G.nodes[s]['conso_ete']:g} MW", (x, y), xytext=(5, 5), textcoords='offset points',
                    fontsize=src_font, color=INK, fontweight='bold', zorder=7,
                    bbox=dict(boxstyle='round,pad=0.15', fc=SURFACE, ec='none', alpha=0.8))
    ax.set_aspect('equal')
    ax.autoscale_view()
    ax.axis('off')


def draw_instance(graph, snaps: dict, name: str, filepath: Path):
    """Image de l'instance : un panneau par distance de `snaps` = {distance: snapshot}."""
    G = graph.G
    n = len(G)
    sources = [s for s, lab in G.nodes('label') if lab == 'Source']
    loads = [v for v, lab in G.nodes('label') if lab in LOAD_LABELS]
    demand = {v: abs(G.nodes[v]['conso_ete']) for v in loads}
    d_max = max(demand.values(), default=1) or 1
    max_area = min(300, max(25, 300 * sqrt(60 / n)))
    sizes = {v: max_area * d / d_max for v, d in demand.items()}
    colors = neighbour_colors(G, snaps, sources)
    area = graph_area(G)

    f, axes = plt.subplots(2, 2, figsize=(14, 14))
    f.patch.set_facecolor(SURFACE)
    for ax, (method, snap) in zip(axes.ravel(), snaps.items()):
        draw_panel(ax, G, snap, colors, area, sizes, 8 if len(sources) <= 6 else 6)
        verdict = "faisable" if snap['feasible'] else f"infaisable ({snap['n_violations']} violations)"
        ax.set_title(f"{WEIGHTING_NAMES[method]}\nJ = {snap['J']:.1f} km  ·  {snap['time']:.2f} s  ·  {verdict}",
                     fontsize=11, color=GREEN if snap['feasible'] else RED, loc='left')
    legend = [Line2D([], [], marker='s', ls='', mfc='w', mec=INK, ms=9, label="source (capacité en MW)"),
              Line2D([], [], color=INK_2, lw=1.6, label="arbre de Steiner"),
              Line2D([], [], color=INK_2, lw=0.4, label="autres arêtes")]
    legend += [Line2D([], [], marker='o', ls='', mfc=SURFACE, mec=INK, ms=sqrt(max_area * fr),
                      label=f"conso {d_max * fr:.3g} MW") for fr in (1, 0.25)]
    legend += [Line2D([], [], marker='^', ls='', mfc=SURFACE, mec=INK, ms=sqrt(max_area * 0.25),
                      label=f"prod {d_max * 0.25:.3g} MW"),
               Line2D([], [], marker='o', ls='', mfc=SURFACE, mec=INK, ms=sqrt(ZERO_AREA), label="charge 0 MW")]
    f.legend(handles=legend, loc='lower center', ncol=len(legend), frameon=False, fontsize=9, labelcolor=INK)
    n_prod = sum(G.nodes[v]['conso_ete'] < 0 for v in loads)
    f.suptitle(f"{name}\n{n} nœuds · sources : {len(sources)} · charges : {len(loads)} dont {n_prod} productions"
               "  —  zones colorées = cellules de Voronoi, aire des charges ∝ |demande|",
               fontsize=12, color=INK, x=0.02, ha='left')
    f.tight_layout(rect=(0, 0.04, 1, 0.96))
    f.savefig(filepath, dpi=110, facecolor=SURFACE)
    plt.close(f)


# ---------------------------------------------------------------- graphiques de synthèse

def grouped_bars(groups: list, values: dict, title: str, ylabel: str, path: Path, fmt=None, ymax=None):
    """Barres groupées : un groupe par élément de `groups`, une barre par distance ({distance: [valeur par groupe]})."""
    f, ax = _new_axes(width=max(8, 1.3 * len(groups) + 3))
    width = 0.8 / len(METHODS)
    for j, m in enumerate(METHODS):
        xs = np.arange(len(groups)) + (j - (len(METHODS) - 1) / 2) * width
        ax.bar(xs, values[m], width, color=WEIGHTING_COLORS[m], edgecolor=SURFACE, linewidth=2,
               label=WEIGHTING_NAMES[m], zorder=2)
        if fmt:
            for x, v in zip(xs, values[m]):
                ax.annotate(fmt(v), (x, v), xytext=(0, 2), textcoords='offset points', ha='center',
                            va='bottom', fontsize=7, color=INK_2)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, color=INK)
    ax.set_ylabel(ylabel, color=INK_2, fontsize=9)
    ax.set_ylim(0, ymax)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, title="Distance de Voronoi", title_fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1))
    _save(f, ax, title, path)


def chart_time(recs: list, path: Path):
    """Temps de calcul en fonction du nombre de nœuds (log-log), médiane par taille."""
    f, ax = _new_axes(width=9)
    markers = dict(zip(METHODS, ('o', 's', 'D', '^')))
    for m in METHODS:
        x = np.array([r['n_nodes'] for r in recs])
        t = np.maximum([r['methods'][m]['time'] for r in recs], 1e-3)
        ax.scatter(x, t, s=16, marker=markers[m], color=WEIGHTING_COLORS[m], alpha=0.35, linewidths=0, zorder=2)
        sizes = sorted(set(x))
        med = [np.median(t[x == n]) for n in sizes]
        ax.plot(sizes, med, color=WEIGHTING_COLORS[m], lw=2, marker=markers[m], ms=6,
                label=WEIGHTING_NAMES[m], zorder=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("nombre de nœuds", color=INK_2, fontsize=9)
    ax.set_ylabel("temps Voronoi + Steiner (s)", color=INK_2, fontsize=9)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, title="Distance (ligne = médiane)", title_fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1))
    _save(f, ax, "Temps de calcul selon la taille de l'instance", path)


def chart_feasibility_matrix(recs: list, path: Path):
    """Carte nombre de nœuds x distance de la part d'instances faisables."""
    sizes = sorted({r['n_nodes'] for r in recs})
    M = np.array([[np.mean([r['methods'][m]['feasible'] for r in recs if r['n_nodes'] == n]) for m in METHODS]
                  for n in sizes])
    counts = {n: sum(r['n_nodes'] == n for r in recs) for n in sizes}
    f, ax = plt.subplots(figsize=(7, 0.42 * len(sizes) + 1.4))
    f.patch.set_facecolor(SURFACE)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq", SEQ_BLUE)
    ax.pcolormesh(M, cmap=cmap, vmin=0, vmax=1, edgecolors=SURFACE, linewidth=2)
    for (i, j), v in np.ndenumerate(M):
        k = round(v * counts[sizes[i]])
        ax.text(j + 0.5, i + 0.5, f"{k}/{counts[sizes[i]]}", ha='center', va='center', fontsize=8,
                color="white" if v >= 0.55 else INK)
    ax.set_xticks(np.arange(len(METHODS)) + 0.5)
    ax.set_xticklabels([WEIGHTING_NAMES[m] for m in METHODS], color=INK, fontsize=9)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(sizes)) + 0.5)
    ax.set_yticklabels([f"{n} nœuds" for n in sizes], color=INK, fontsize=8)
    ax.invert_yaxis()
    for side in ax.spines.values():
        side.set_visible(False)
    ax.tick_params(length=0)
    _save(f, ax, "Instances faisables par taille et par distance", path)


def write_charts(recs: list, charts_dir: Path) -> list:
    """Graphiques de synthèse dans charts_dir ; renvoie les lignes Markdown qui les affichent."""
    charts_dir.mkdir(parents=True, exist_ok=True)
    classes = [c for c in SIZE_CLASSES if any(r['classe'] == c for r in recs)]
    pct = lambda v: f"{v:.0f} %"
    md = ["## Graphiques", ""]

    def add(key, title, caption, draw):
        path = charts_dir / f"{key}.png"
        draw(path)
        md.extend([f"### {title}", "", caption, "", f"![{title}]({charts_dir.name}/{path.name})", ""])

    feas = {m: [100 * np.mean([r['methods'][m]['feasible'] for r in recs if r['classe'] == c]) for c in classes]
            for m in METHODS}
    add("faisables_classe", "Instances faisables par classe de taille",
        "Part des instances où la solution respecte toutes les contraintes.",
        lambda p: grouped_bars([f"{c}\n({sum(r['classe'] == c for r in recs)} inst.)" for c in classes], feas,
                               "Instances faisables par classe de taille", "instances faisables (%)", p, pct, 105))
    add("faisables_taille", "Instances faisables par nombre de nœuds",
        "Nombre d'instances faisables / nombre d'instances, pour chaque taille et chaque distance.",
        lambda p: chart_feasibility_matrix(recs, p))
    combos = sorted({(r['ratio'], r['saison']) for r in recs if r['ratio'] is not None})
    if combos:
        feas_rs = {m: [100 * np.mean([r['methods'][m]['feasible'] for r in recs
                                      if (r['ratio'], r['saison']) == c]) for c in combos] for m in METHODS}
        add("faisables_ratio_saison", "Instances faisables par ratio et saison",
            "Ratio et saison lus dans le nom de l'instance (…_<ratio>_<saison>).",
            lambda p: grouped_bars([f"{ra:.2f}\n{s}" for ra, s in combos], feas_rs,
                                   "Instances faisables par ratio et saison", "instances faisables (%)", p, pct, 105))
    gaps = {m: [np.mean([gap_J(r, m) for r in recs if r['classe'] == c]) for c in classes] for m in METHODS}
    add("ecart_J", "Écart de J à la meilleure distance",
        "Pour chaque instance : 100 × (J − J_min) / J_min, J_min = plus petit J des 4 distances (faisable ou non) ;"
        " moyenne par classe.",
        lambda p: grouped_bars(classes, gaps, "Écart moyen de J à la meilleure distance", "écart (%)", p,
                               lambda v: f"{v:.1f}"))
    viol = {m: [100 * np.mean([(r['methods'][m]['violations'].get(k, [0])[0] if k != 'errors'
                                else len(r['methods'][m]['errors'])) > 0 for r in recs]) for k in CONSTRAINTS]
            for m in METHODS}
    add("violations", "Contraintes violées",
        "Part des instances où au moins une contrainte du type est violée (structure = arbre invalide).",
        lambda p: grouped_bars(list(CONSTRAINTS.values()), viol, "Instances avec au moins une violation, par type",
                               "instances (%)", p, pct, 105))
    add("temps", "Temps de calcul",
        "Un point par instance, ligne = médiane par taille. Échelles logarithmiques.",
        lambda p: chart_time(recs, p))
    return md


# ---------------------------------------------------------------- rapport

def gap_J(rec: dict, m: str) -> float:
    J_min = min(x['J'] for x in rec['methods'].values())
    return 100 * (rec['methods'][m]['J'] - J_min) / J_min if J_min > 0 else 0.0


def best_feasible(rec: dict):
    """Distance de plus petit J parmi les solutions faisables (None si aucune)."""
    ok = [m for m in METHODS if rec['methods'][m]['feasible']]
    return min(ok, key=lambda m: rec['methods'][m]['J']) if ok else None


def write_csv(recs: list, path: Path):
    cols = ['instance', 'classe', 'n_nodes', 'n_edges', 'n_sources', 'S_nom', 'ratio', 'saison', 'n_loads',
            'n_prod', 'distance', 'J', 'PD', 'time_s', 'feasible', 'n_violations'] + [f"viol_{k}" for k in CONSTRAINTS]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter=";")
        w.writerow(cols)
        for r in recs:
            for m in METHODS:
                x = r['methods'][m]
                w.writerow([r[c] for c in cols[:10]] + [m, f"{x['J']:.3f}", f"{x['PD']:.3f}", f"{x['time']:.4f}",
                            x['feasible'], n_violations(x)]
                           + [x['violations'].get(k, [0])[0] if k != 'errors' else len(x['errors'])
                              for k in CONSTRAINTS])


def write_report(recs: list, errors: list, out_dir: Path) -> Path:
    """report.md, results.csv et charts/ à partir des enregistrements de results.jsonl."""
    write_csv(recs, out_dir / "results.csv")
    n = len(recs)
    md = [f"# Voronoi + Steiner sur les instances vérifiées", "",
          f"{datetime.now():%d/%m/%Y %H:%M} · {n} instances · 4 distances : "
          + ", ".join(md_weighting(m) for m in METHODS) + ".", "",
          "J = L_tot, longueur totale des arbres (km), à minimiser. Une solution est faisable si les arbres forment"
          " un réseau radial valide couvrant toutes les charges et respectent L_max et P_max par départ, la capacité"
          " des sources et le nombre de départs (verify.py). Temps = Voronoi + Steiner (+ ajustement des poids pour"
          " la distance variable).", "",
          "## Synthèse", "",
          "| Distance | Faisables | Meilleure solution faisable | Écart moyen de J (%) | Temps moyen (s)"
          " | Temps médian (s) | Temps max (s) |", "|---|---|---|---|---|---|---|"]
    best = {r['instance']: best_feasible(r) for r in recs}
    for m in METHODS:
        t = [r['methods'][m]['time'] for r in recs]
        k = sum(r['methods'][m]['feasible'] for r in recs)
        md.append(f"| {md_weighting(m)} | {k}/{n} ({100 * k / n:.0f} %) | {sum(b == m for b in best.values())}"
                  f" | {np.mean([gap_J(r, m) for r in recs]):.1f} | {np.mean(t):.2f} | {np.median(t):.2f}"
                  f" | {max(t):.2f} |")
    any_ok = sum(b is not None for b in best.values())
    proven = sum(bool(r['proofs']) for r in recs)
    md += ["", f"Au moins une distance faisable : **{any_ok}/{n}** instances. "
           f"Instances prouvées infaisables (conditions nécessaires de verify.py) : {proven}.", "",
           "### Par classe de taille", "",
           "| Classe | Instances | Nœuds | " + " | ".join(f"Faisables {md_weighting(m)}" for m in METHODS)
           + " | Temps moyen variable (s) |", "|---|---|---|" + "---|" * len(METHODS) + "---|"]
    for c in SIZE_CLASSES:
        rc = [r for r in recs if r['classe'] == c]
        if not rc:
            continue
        sizes = sorted({r['n_nodes'] for r in rc})
        md.append(f"| {c} | {len(rc)} | {', '.join(map(str, sizes))} | "
                  + " | ".join(f"{sum(r['methods'][m]['feasible'] for r in rc)}/{len(rc)}" for m in METHODS)
                  + f" | {np.mean([r['methods']['variable']['time'] for r in rc]):.2f} |")
    md += [""] + write_charts(recs, out_dir / "charts")

    md += ["## Détail par instance", "",
           f"J (km) et temps (s) de chaque distance : {colored('vert', GREEN)} = faisable, {colored('rouge', RED)}"
           " = infaisable, gras = plus petit J des 4. Le nom renvoie à l'image de l'instance.", "",
           "| Instance | Nœuds | Sources | Charges (prod) | " + " | ".join(md_weighting(m) for m in METHODS)
           + " | Meilleure faisable |", "|---|---|---|---|" + "---|" * len(METHODS) + "---|"]
    for r in recs:
        J_min = min(x['J'] for x in r['methods'].values())
        cells = []
        for m in METHODS:
            x = r['methods'][m]
            cells.append(colored(f"{x['J']:.1f}", GREEN if x['feasible'] else RED, bold=x['J'] == J_min)
                         + f" ({x['time']:.2f} s)")
        name = f"[{r['instance']}]({r['figure'].replace(' ', '%20')})" if r.get('figure') else r['instance']
        b = best[r['instance']]
        md.append(f"| {name} | {r['n_nodes']} | {r['n_sources']} | {r['n_loads']} ({r['n_prod']}) | "
                  + " | ".join(cells) + f" | {md_weighting(b) if b else '—'} |")
    if errors:
        md += ["", "## Instances en erreur", ""] + [f"- {e['instance']} : `{e['error']}`" for e in errors]
    path = out_dir / "report.md"
    path.write_text("\n".join(md) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------- main

def to_builtin(o):
    """Pour json.dumps : scalaires numpy -> nombres Python, le reste en texte."""
    return o.item() if hasattr(o, 'item') else str(o)


def read_results(path: Path) -> dict:
    """{instance: enregistrement} de results.jsonl (la dernière ligne d'une instance l'emporte)."""
    recs = {}
    if path.exists():
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    r = json.loads(line)
                    recs[r['instance']] = r
    return recs


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("motifs", nargs="*", help="seulement les instances dont le nom contient un de ces motifs")
    p.add_argument("--dir", type=Path, default=DEFAULT_DIR, help=f"dossier des instances (défaut : {DEFAULT_DIR})")
    p.add_argument("--out", type=Path, default=OUT_DIR, help=f"dossier des résultats (défaut : {OUT_DIR})")
    p.add_argument("--max-nodes", type=int, help="seulement les instances d'au plus ce nombre de nœuds")
    p.add_argument("--no-figures", action="store_true", help="pas d'image par instance")
    p.add_argument("--restart", action="store_true", help="effacer les résultats précédents")
    p.add_argument("--report", action="store_true", help="refaire seulement le rapport")
    a = p.parse_args()
    logging.disable(logging.INFO)

    out = a.out
    out.mkdir(parents=True, exist_ok=True)
    results_path = out / "results.jsonl"
    if a.restart:
        results_path.unlink(missing_ok=True)
    recs = read_results(results_path)
    if not a.report:
        if not a.dir.is_dir():
            p.error(f"dossier d'instances introuvable : {a.dir}")
        figures_dir = None if a.no_figures else out / "figures"
        if figures_dir:
            figures_dir.mkdir(exist_ok=True)
        files = instance_files(a.dir, a.motifs, a.max_nodes)
        todo = [f for f in files if 'error' in recs.get(f.stem, {'error': 1})]
        print(f"{len(files)} instances dans {a.dir}, {len(files) - len(todo)} déjà faites, {len(todo)} à faire",
              flush=True)
        t_start = time.time()
        for i, f in enumerate(todo, 1):
            t0 = time.time()
            try:
                rec = run_instance(f, figures_dir)
                ms = rec['methods']
                status = " ".join(f"{m[:3]}={'OK' if ms[m]['feasible'] else 'KO'}" for m in METHODS)
                msg = f"J_var = {ms['variable']['J']:.1f} km  {status}"
            except Exception as exc:   # une instance qui plante n'arrête pas le run
                rec = {'instance': f.stem, 'error': repr(exc), 'traceback': traceback.format_exc()}
                msg = f"ERREUR {exc!r}"
            with open(results_path, "a", encoding="utf-8") as fh:
                fh.write(json.dumps(rec, default=to_builtin) + "\n")
            recs[f.stem] = rec
            print(f"[{i}/{len(todo)}] {f.stem} ({instance_size(f)[0]} nœuds) : {msg}  "
                  f"[{time.time() - t0:.1f} s, total {(time.time() - t_start) / 60:.1f} min]", flush=True)
    done = [r for r in recs.values() if 'error' not in r]
    done.sort(key=lambda r: (r['n_nodes'], r['instance']))
    errors = [r for r in recs.values() if 'error' in r]
    if not done:
        print("Aucun résultat à mettre dans le rapport.")
        return
    report = write_report(done, errors, out)
    print(f"\nRapport : {report}\nRésultats : {out / 'results.csv'}\nImages : {out / 'figures'}")


if __name__ == "__main__":
    main()
