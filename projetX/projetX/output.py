"""
Reports of a run of main.py, written in output/:
  recap_<stamp>.md / .csv   summary tables, one row per graph and distance
  charts/<stamp>_*.png      charts embedded in recap_<stamp>.md
  details_<stamp>.txt       per-source and per-départ details of each evaluation
"""
import csv
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from voronoi import WEIGHTINGS

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
CHARTS_DIR = OUTPUT_DIR / "charts"

# the 3 fixed distances of voronoi.py, then the variable distance of voronoi_variable_distance.py
METHODS = WEIGHTINGS + ('variable',)
WEIGHTING_NAMES = {'none': "non pondérée", 'capacity': "capacité", 'sqrt': "racine capacité",
                   'variable': "variable"}


def format_evaluation(name: str, weighting: str, ev: dict, weights: dict = None) -> list:
    """Lines describing each source and départ of an evaluation (`weights`: {source: w} of the distance)."""
    ok = lambda b: "OK" if b else "KO"
    lines = [f"=== {name} | distance {WEIGHTING_NAMES[weighting]} ===",
             f"J = L_tot = {ev['J']:.1f} km"
             f"  PD = {ev['PD']:.1f} MW.km"
             f"  -> {'toutes contraintes satisfaites' if ev['feasible'] else 'contraintes violées'}"]
    lines += [f"  ERREUR : {e}" for e in ev['errors']]
    for s, src in ev['sources'].items():
        lines.append(f"  Source {s}: P = {src['power']:7.2f} / {src['capacity']:6.2f} MW {ok(src['source_ok'])}"
                     f" | départs {src['n_departs']} / {src['n_departs_max']:g} {ok(src['n_ok'])}"
                     f" | PD = {src['pd']:7.1f} MW.km | {src['n_cell']} noeuds dans l'arbre"
                     + (f" | w = {weights[s]:.3g}" if weights else ""))
        for i, d in enumerate(src['departs']):
            lines.append(f"    départ {i} (tête {d['head']:3d}, {d['n_nodes']:3d} noeuds):"
                         f" L = {d['length']:6.1f} / {src['l_max_depart']:g} km {ok(d['L_ok'])}"
                         f" | P = {d['power']:6.2f} / {src['p_max_depart']:g} MW {ok(d['P_ok'])}"
                         f" | PD = {d['pd']:7.1f} MW.km")
    return lines


def constraint_cell(v: list) -> str:
    """'OK' if no violation, else 'KO (violated/total)'."""
    return "OK" if v[0] == 0 else f"KO ({v[0]}/{v[1]})"


# Colors of the markdown summary (HTML spans, shown e.g. by the VS Code preview)
GREEN, RED = "#1a7f37", "#cf222e"
WEIGHTING_COLORS = {'none': "#2a78d6", 'capacity': "#eb6834", 'sqrt': "#4a3aa7", 'variable': "#1baf7a"}


def colored(text, color: str, bold: bool = False) -> str:
    text = f"<b>{text}</b>" if bold else text
    return f'<span style="color:{color}">{text}</span>'


def md_constraint(v: list) -> str:
    return colored("OK", GREEN) if v[0] == 0 else colored(f"KO ({v[0]}/{v[1]})", RED)


def md_weighting(w: str) -> str:
    return colored(WEIGHTING_NAMES[w], WEIGHTING_COLORS[w], bold=True)


# Charts of the summary (static PNG, light theme): ink, grid and surface colors
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SEQ_BLUE = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
TYPE_ORDER = ['ville', 'montagne', 'epars', 'mixte']


def relative_gaps(results: list, key: str) -> dict:
    """{graph: {weighting: 100 * (X_w - X_best) / X_best}}, X = ev[key], best = min of the METHODS."""
    by_graph = {}
    for r in results:
        by_graph.setdefault(r['name'], {})[r['weighting']] = r['ev'][key]
    return {name: {w: 100 * (x - min(xs.values())) / min(xs.values()) for w, x in xs.items()}
            for name, xs in by_graph.items()}


def _new_axes(width=8, height=4.2):
    f, ax = plt.subplots(figsize=(width, height))
    f.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    return f, ax


def _save(f, ax, title, path):
    ax.set_title(title, color=INK, fontsize=11, loc='left', pad=12)
    f.savefig(path, dpi=150, bbox_inches='tight', facecolor=SURFACE)
    plt.close(f)


def chart_gaps(gaps: dict, group_of: dict, groups: list, group_label, title, ylabel, path,
               feasible: dict):
    """
    Grouped bars: mean relative gap of each distance in each group of graphs,
    with one mark per graph: a dot if the solution is feasible, a cross otherwise.
    `group_of` = {graph: group}, `feasible` = {(graph, weighting): bool}.
    """
    f, ax = _new_axes(width=max(8, 1.2 * len(groups) + 3))
    width = 0.8 / len(METHODS) - 0.02
    rng = np.random.default_rng(0)
    for i, g in enumerate(groups):
        names = [n for n in gaps if group_of[n] == g]
        for j, w in enumerate(METHODS):
            x = i + (j - (len(METHODS) - 1) / 2) * (width + 0.02)
            values = [gaps[n][w] for n in names]
            mean = sum(values) / len(values)
            # light bar (mean) so that the dots (one per graph) stay visible on it
            ax.bar(x, mean, width, color=WEIGHTING_COLORS[w], alpha=0.35, zorder=2,
                   edgecolor=WEIGHTING_COLORS[w], linewidth=1,
                   label=WEIGHTING_NAMES[w] if i == 0 else None)
            xs = x + rng.uniform(-0.07, 0.07, len(values))
            ok = np.array([feasible[n, w] for n in names], dtype=bool)
            vs = np.array(values)
            ax.scatter(xs[ok], vs[ok], s=24, color=WEIGHTING_COLORS[w], edgecolor=SURFACE,
                       linewidth=1.2, zorder=3)
            ax.scatter(xs[~ok], vs[~ok], s=30, marker='x', color=WEIGHTING_COLORS[w],
                       linewidth=1.5, zorder=3)
            ax.annotate(f"{mean:.0f} %", (x, mean), xytext=(0, 3), textcoords='offset points',
                        ha='center', va='bottom', fontsize=8, color=INK)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([group_label(g, sum(group_of[n] == g for n in gaps)) for g in groups], color=INK)
    ax.set_ylabel(ylabel, color=INK_2, fontsize=9)
    ax.set_ylim(bottom=0)
    ax.scatter([], [], s=24, color=INK_2, label="graphe, solution faisable")
    ax.scatter([], [], s=30, marker='x', color=INK_2, linewidth=1.5, label="graphe, solution infaisable")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, title="Distance de Voronoi", title_fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1))
    _save(f, ax, title, path)


def chart_violations(results: list, path):
    """Heatmap graph x distance of the number of violated constraints."""
    names = list(dict.fromkeys(r['name'] for r in results))
    count = {(r['name'], r['weighting']): sum(v[0] for v in r['ev']['violations'].values())
             + len(r['ev']['errors'])
             for r in results}
    M = np.array([[count[n, w] for w in METHODS] for n in names])
    f, ax = plt.subplots(figsize=(7, 0.32 * len(names) + 1.2))
    f.patch.set_facecolor(SURFACE)
    vmax = max(int(M.max()), 1)
    # one ramp step per count while it fits, then a continuous ramp
    cmap = matplotlib.colors.ListedColormap(SEQ_BLUE[:vmax + 1]) if vmax < len(SEQ_BLUE) \
        else matplotlib.colors.LinearSegmentedColormap.from_list("seq", SEQ_BLUE)
    ax.pcolormesh(M, cmap=cmap, vmin=-0.5, vmax=vmax + 0.5, edgecolors=SURFACE, linewidth=2)
    for (i, j), v in np.ndenumerate(M):
        ax.text(j + 0.5, i + 0.5, str(v), ha='center', va='center', fontsize=8,
                color="white" if v > 0 and v >= 0.55 * vmax else INK)
    ax.set_xticks(np.arange(len(METHODS)) + 0.5)
    ax.set_xticklabels([WEIGHTING_NAMES[w] for w in METHODS], color=INK, fontsize=9)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(names)) + 0.5)
    ax.set_yticklabels(names, color=INK, fontsize=8)
    ax.invert_yaxis()
    for side in ax.spines.values():
        side.set_visible(False)
    ax.tick_params(length=0)
    _save(f, ax, "Nombre de contraintes violées (0 = solution faisable)", path)


def write_charts(results: list, stamp: str) -> list:
    """Save the charts of a run in CHARTS_DIR, return the markdown lines embedding them."""
    CHARTS_DIR.mkdir(parents=True, exist_ok=True)
    gaps_J, gaps_PD = relative_gaps(results, 'J'), relative_gaps(results, 'PD')
    feasible = {(r['name'], r['weighting']): r['ev']['feasible'] for r in results}
    n_loads = {r['name']: r['n_loads'] for r in results}
    types = {r['name']: r['type'] for r in results}
    load_groups = sorted(set(n_loads.values()))
    type_groups = [t for t in TYPE_ORDER if t in types.values()]
    by_loads = lambda g, k: f"{g} charges\n({k} graphes)"
    by_type = lambda g, k: f"{g} ({k} graphes)"
    charts = [
        ("gap_J_charges", "Écart relatif à l'optimum de J, par nombre de charges",
         lambda p: chart_gaps(gaps_J, n_loads, load_groups, by_loads,
                              "Écart de J à la meilleure distance, par nombre de charges",
                              "écart à l'optimum (%)", p, feasible),
         "J = L_tot. Pour chaque graphe, l'optimum est la plus petite valeur des 4 distances (0 %),"
         " faisable ou non. Barres : écart moyen ; un marqueur par graphe : rond si la solution"
         " respecte les contraintes, croix sinon."),
        ("gap_PD_charges", "Écart relatif à l'optimum de PD, par nombre de charges",
         lambda p: chart_gaps(gaps_PD, n_loads, load_groups, by_loads,
                              "Écart de PD à la meilleure distance, par nombre de charges",
                              "écart à l'optimum (%)", p, feasible),
         "Même lecture, pour le produit PD (puissance × distance)."),
        ("gap_J_type", "Écart relatif à l'optimum de J, par type de répartition",
         lambda p: chart_gaps(gaps_J, types, type_groups, by_type,
                              "Écart de J à la meilleure distance, par type de répartition",
                              "écart à l'optimum (%)", p, feasible),
         "Quelle distance convient à quelle répartition des charges."),
        ("violations", "Contraintes violées",
         lambda p: chart_violations(results, p),
         "Nombre total de contraintes violées (L_max, P_max départ, P_max source, nb départs, erreurs de structure)"
         " pour chaque graphe et chaque méthode."),
    ]
    md = ["", "## Graphiques", ""]
    for key, title, draw, caption in charts:
        path = CHARTS_DIR / f"{stamp}_{key}.png"
        draw(path)
        md += [f"### {title}", "", caption, "",
               f"![{title}]({path.relative_to(OUTPUT_DIR).as_posix()})", ""]
    return md


def md_synthesis_by_loads(results: list, best_J: dict, gaps_J: dict) -> list:
    """Markdown table: for each number of loads, feasible runs, best J and mean gap of J per distance."""
    n_loads = {r['name']: r['n_loads'] for r in results}
    feasible = {(r['name'], r['weighting']): r['ev']['feasible'] for r in results}
    md = ["", "## Synthèse par nombre de charges", "",
          "Par groupe de graphes de même nombre de charges : runs faisables, nombre de graphes où la"
          " distance donne le plus petit J, écart moyen de J à l'optimum.", "",
          "| Charges | Graphes | " + " | ".join(f"Faisables {md_weighting(w)}" for w in METHODS)
          + " | " + " | ".join(f"Meilleur J {md_weighting(w)}" for w in METHODS)
          + " | " + " | ".join(f"Écart J (%) {md_weighting(w)}" for w in METHODS) + " |",
          "|---|---|" + "---|" * 3 * len(METHODS)]
    for k in sorted(set(n_loads.values())):
        names = [n for n in n_loads if n_loads[n] == k]
        cells = [f"{sum(feasible[n, w] for n in names)}/{len(names)}" for w in METHODS]
        cells += [str(sum(best_J[n] == w for n in names)) for w in METHODS]
        cells += [f"{sum(gaps_J[n][w] for n in names) / len(names):.1f}" for w in METHODS]
        md.append(f"| {k} | {len(names)} | " + " | ".join(cells) + " |")
    return md


def write_summary(results: list, stamp: str) -> Path:
    """Write the summary table (markdown and csv) of a run, return the markdown path."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    by_graph = {}
    for r in results:
        by_graph.setdefault(r['name'], {})[r['weighting']] = r['ev']
    # all the solutions are compared, feasible or not (feasibility is only shown by the colors)
    best_J = {name: min(evs, key=lambda w: evs[w]['J']) for name, evs in by_graph.items()}
    best_PD = {name: min(evs, key=lambda w: evs[w]['PD']) for name, evs in by_graph.items()}

    columns = ["Graphe", "Type", "Charges (dont <0)", "Sources", "Distance", "J = L_tot (km)",
               "PD (MW·km)", "L_max départ", "P_max départ", "P_max source",
               "Nb départs", "Faisable"]
    rows, md_rows = [], []
    for r in results:
        ev, viol, w, name = r['ev'], r['ev']['violations'], r['weighting'], r['name']
        rows.append([name, r['type'], f"{r['n_loads']} ({r['n_neg']})", r['n_sources'],
                     WEIGHTING_NAMES[w], f"{ev['J']:.1f}", f"{ev['PD']:.1f}", constraint_cell(viol['L']),
                     constraint_cell(viol['P']), constraint_cell(viol['source']),
                     constraint_cell(viol['n']), "oui" if ev['feasible'] else "non"])
        first = w == METHODS[0]   # graph columns only on the first row of each graph
        md_rows.append([f"**{name}**" if first else "",
                        r['type'] if first else "",
                        f"{r['n_loads']} ({r['n_neg']})" if first else "",
                        r['n_sources'] if first else "",
                        md_weighting(w),
                        colored(f"{ev['J']:.1f}", GREEN if ev['feasible'] else RED, bold=w == best_J[name]),
                        colored(f"{ev['PD']:.1f}", GREEN, bold=True) if w == best_PD[name] else f"{ev['PD']:.1f}",
                        md_constraint(viol['L']), md_constraint(viol['P']),
                        md_constraint(viol['source']), md_constraint(viol['n']),
                        colored("oui", GREEN, bold=True) if ev['feasible'] else colored("non", RED, bold=True)])
    with open(OUTPUT_DIR / f"recap_{stamp}.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(columns)
        writer.writerows(rows)

    md = [f"# Récapitulatif du run {stamp}", "",
          "Objectif minimisé : J = L_tot, longueur totale des arbres (km), comme dans"
          " data/exemple_ete_4_PS/model.lp. Les contraintes (L_max et P_max par départ, capacité"
          " des sources, nombre de départs) ne sont pas dans J : elles sont seulement vérifiées."
          " Toutes les solutions sont comparées, même celles qui violent une contrainte"
          " (infaisables : en rouge dans les tableaux, croix sur les graphiques).", "",
          "PD = Σ sur les arêtes des arbres de |puissance transitée| × longueur (MW·km), "
          "non inclus dans J.", "",
          f"Lecture : {colored('OK', GREEN)} contrainte satisfaite partout, "
          f"{colored('KO (violées/testées)', RED)} sinon ; "
          f"J en {colored('vert', GREEN)} si la solution est faisable, en {colored('rouge', RED)} sinon,"
          " en gras pour la plus courte solution du graphe (faisable ou non) ;"
          f" PD : {colored('valeur', GREEN, bold=True)} = plus petit PD des 4 distances ; "
          "distances : " + ", ".join(md_weighting(w) for w in METHODS) + ".", "",
          "## Détail par graphe et par distance", "",
          "| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    md += ["| " + " | ".join(str(x) for x in row) + " |" for row in md_rows]

    # J and PD per graph, the METHODS side by side for each quantity
    md += ["", "## Fonction objectif J = L_tot (km) et produit PD (MW·km) par distance", "",
           f"J : {colored('vert', GREEN)} si toutes les contraintes sont satisfaites, "
           f"{colored('rouge', RED)} sinon. En gras : plus petite valeur des 4 distances,"
           " solutions faisables ou non.", "",
           "| Graphe | " + " | ".join(f"J {md_weighting(w)}" for w in METHODS) + " | "
           + " | ".join(f"PD {md_weighting(w)}" for w in METHODS) + " | Meilleur J | Meilleur PD |",
           "|---|" + "---|" * (2 * len(METHODS) + 2)]
    wins = {'J': {w: 0 for w in METHODS}, 'PD': {w: 0 for w in METHODS}}
    for name, evs in by_graph.items():
        wins['J'][best_J[name]] += 1
        wins['PD'][best_PD[name]] += 1
        J_cells = [colored(f"{evs[w]['J']:.1f}", GREEN if evs[w]['feasible'] else RED,
                           bold=w == best_J[name]) for w in METHODS]
        PD_cells = [f"**{evs[w]['PD']:.1f}**" if w == best_PD[name] else f"{evs[w]['PD']:.1f}"
                    for w in METHODS]
        md.append(f"| {name} | " + " | ".join(J_cells + PD_cells)
                  + f" | {md_weighting(best_J[name])} | {md_weighting(best_PD[name])} |")
    gaps = {'J': relative_gaps(results, 'J'), 'PD': relative_gaps(results, 'PD')}
    md += ["",
           "## Synthèse", "",
           "Écart à l'optimum : 100 × (X − X_min) / X_min, X_min = plus petite valeur des 4 distances"
           " sur le graphe ; moyenne sur les graphes (comparable entre graphes de tailles différentes,"
           " contrairement aux moyennes de J et PD).", "",
           "| Distance | Meilleur J (nb graphes) | Meilleur PD (nb graphes) | Runs faisables"
           " | Écart moyen J (%) | Écart moyen PD (%) | J moyen (km) | PD moyen (MW·km) |",
           "|---|---|---|---|---|---|---|---|"]
    for w in METHODS:
        rw = [r['ev'] for r in results if r['weighting'] == w]
        mean = lambda k: sum(e[k] for e in rw) / len(rw)
        mean_gap = lambda k: sum(g[w] for g in gaps[k].values()) / len(gaps[k])
        md.append(f"| {md_weighting(w)} | {wins['J'][w]} | {wins['PD'][w]}"
                  f" | {sum(e['feasible'] for e in rw)}/{len(rw)}"
                  f" | {mean_gap('J'):.1f} | {mean_gap('PD'):.1f}"
                  f" | {mean('J'):.1f} | {mean('PD'):.1f} |")
    md += md_synthesis_by_loads(results, best_J, gaps['J'])
    md += write_charts(results, stamp)
    path = OUTPUT_DIR / f"recap_{stamp}.md"
    path.write_text("\n".join(md) + "\n", encoding="utf-8")
    return path


def write_details(details: list, stamp: str) -> Path:
    """Write the lines of format_evaluation of a run, return the path."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    path = OUTPUT_DIR / f"details_{stamp}.txt"
    path.write_text("\n".join(details), encoding="utf-8")
    return path
