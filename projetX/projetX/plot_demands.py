"""
Carte des demandes de chaque instance de verified_instances/ : une image par
instance dans data/verified_instances/demand_maps/<instance>.png.

L'aire d'une charge est proportionnelle à |demande|, avec la même échelle sur
toutes les images, et sa demande (MW, négative pour une production) est écrite
à côté du nœud. Les sources n'ont pas de demande (leur champ `demande` est leur
capacité) : carrés de taille fixe. Le panneau de droite liste les contraintes
de l'instance : sources, départs, nœuds (lues dans le JSON) et structure
(celle imposée par verify.py). Au-delà de REF_NODES nœuds, l'image grandit avec
sqrt(N) : même densité de marqueurs et d'étiquettes qu'à REF_NODES.

Usage :
  python plot_demands.py                  instances Petit_ et Moyen_
  python plot_demands.py Grand_           Grand et TrèsGrand (seulement les noms contenant un de ces motifs)
"""
import json
import sys
from math import cos, pi, sin, sqrt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
IN_DIR = HERE / "data" / "verified_instances"
OUT_DIR = IN_DIR / "demand_maps"

INK, INK_2, EDGE, GRID, SURFACE = "#0b0b0b", "#52514e", "#c9c8c3", "#e4e3df", "#fcfcfb"
CONSO, PROD = "#2a78d6", "#eb6834"

AREA_PER_MW = 120     # aire (pt²) d'un disque de 1 MW
TRIANGLE_FACTOR = 1.571   # à `s` égal, le triangle de matplotlib couvre 0.5 et le disque pi/4 : même aire réelle
ZERO_AREA = 22        # charge de demande nulle : petit cercle vide
SOURCE_AREA = 130
LABEL_SIZE = 9        # police des demandes (pt)
FIG_W, FIG_H = 13.0, 9.4   # pour au plus REF_NODES nœuds
REF_NODES = 250
GRAPH_RECT = (0.015, 0.02, 0.715, 0.875)
PANEL_RECT = (0.745, 0.02, 0.24, 0.875)


def area(mw: float, marker: str = 'o') -> float:
    return AREA_PER_MW * mw * (TRIANGLE_FACTOR if marker == '^' else 1)


def marker_of(n: dict) -> tuple:
    """(marqueur, couleur, aire) d'une charge non nulle."""
    marker, color = ('^', PROD) if n['type'] == 'Prod' else ('o', CONSO)
    return marker, color, area(abs(n['demande']), marker)


def legend_handles(nodes: list) -> list:
    """Légende limitée à ce que l'instance contient."""
    def values(kind):
        return sorted({abs(n['demande']) for n in nodes if n['type'] == kind and n['demande'] != 0})

    handles = [Line2D([], [], marker='s', ls='', mfc=INK, mec=INK, ms=sqrt(SOURCE_AREA), label="source (taille fixe)")]
    for kind, marker, color, name in (('Conso', 'o', CONSO, "conso"), ('Prod', '^', PROD, "prod")):
        handles += [Line2D([], [], marker=marker, ls='', mfc=color, mec=SURFACE, ms=sqrt(area(v, marker)),
                           label=f"{name} {v:g} MW") for v in values(kind)]
    if any(n['type'] == 'Conso' and n['demande'] == 0 for n in nodes):
        handles.append(Line2D([], [], marker='o', ls='', mfc=SURFACE, mec=CONSO, ms=sqrt(ZERO_AREA), label="conso 0 MW"))
    handles.append(Line2D([], [], marker='o', ls='', mfc=INK_2, mec='none', ms=2.5, label="transit"))
    return handles


def overlap(box, boxes: np.ndarray) -> np.ndarray:
    """Aire commune de `box` = (x0, y0, x1, y1) avec chaque ligne de `boxes`."""
    w = np.minimum(box[2], boxes[:, 2]) - np.maximum(box[0], boxes[:, 0])
    h = np.minimum(box[3], boxes[:, 3]) - np.maximum(box[1], boxes[:, 1])
    return np.clip(w, 0, None) * np.clip(h, 0, None)


def place_labels(f, ax, loads: list, sources: list, transit: list) -> int:
    """
    Écrit la demande de chaque charge à côté de son nœud, à la place libre parmi
    8 directions. Renvoie le nombre d'étiquettes qui touchent un marqueur ou une
    autre étiquette malgré tout.
    """
    ax.apply_aspect()    # l'aspect 'equal' appliqué, transData est définitif
    to_pt = lambda xys: ax.transData.transform(np.array(xys, dtype=float).reshape(-1, 2)) * 72 / f.dpi
    frame = (ax.get_window_extent().extents * 72 / f.dpi)[None]

    def boxes(xys, r):
        c = to_pt(xys)
        return np.column_stack((c[:, 0] - r, c[:, 1] - r, c[:, 0] + r, c[:, 1] + r))

    radius = [sqrt(marker_of(n)[2] if n['demande'] else ZERO_AREA) / 2 for n in loads]
    centers = to_pt([n['geometry'] for n in loads])
    marks = np.vstack((boxes([n['geometry'] for n in loads], np.array(radius)), boxes(sources, sqrt(SOURCE_AREA) / 2)))
    dots = boxes(transit, 1.5)
    placed = np.zeros((len(loads), 4))
    collisions = 0
    halo = [path_effects.withStroke(linewidth=2.5, foreground=SURFACE)]
    angles = [i * pi / 4 for i in (1, 0, 2, 3, 4, 7, 6, 5)]    # haut-droite d'abord

    for i, n in enumerate(loads):
        text = f"{n['demande']:g}".replace('-', '−')
        x, y = centers[i]
        r, w, h = radius[i], 0.62 * LABEL_SIZE * len(text), 0.95 * LABEL_SIZE
        best = None
        for rank, angle in enumerate(angles):
            c, s = cos(angle), sin(angle)
            dx, dy = (r + 2) * c, (r + 2) * s
            ha = 'left' if c > 0.3 else 'right' if c < -0.3 else 'center'
            va = 'bottom' if s > 0.3 else 'top' if s < -0.3 else 'center'
            x0 = x + dx - {'left': 0, 'center': w / 2, 'right': w}[ha]
            y0 = y + dy - {'bottom': 0, 'center': h / 2, 'top': h}[va]
            lab = (x0, y0, x0 + w, y0 + h)
            on_marks = overlap(lab, marks)
            on_marks[i] = 0    # son propre marqueur ne compte pas
            hard = on_marks.sum() + overlap(lab, placed[:i]).sum()
            outside = w * h - overlap(lab, frame)[0]
            cost = 20 * hard + 1000 * outside + overlap(lab, dots).sum() + 0.01 * rank
            if best is None or cost < best[0]:
                best = (cost, hard, lab, dx, dy, ha, va)
        _, hard, lab, dx, dy, ha, va = best
        collisions += int(hard > 0)
        placed[i] = lab
        ax.annotate(text, n['geometry'], xytext=(dx, dy), textcoords='offset points', ha=ha, va=va,
                    fontsize=LABEL_SIZE, color=INK, path_effects=halo, annotation_clip=False, zorder=5)
    return collisions


def constraints_panel(f, ax, nodes: list, elec: dict, k: float):
    """Colonne de droite : contraintes de l'instance, puis légende (`k` = échelle du texte et des interlignes)."""
    src = [n for n in nodes if n['type'] == 'Source']
    values = lambda key: ' / '.join(f'{v:g}' for v in sorted({n[key] for n in src}))
    sections = [
        ("Sources", [("Capacité  |P_source|", f"≤ {values('demande')} MW"),
                     ("Nombre de départs", f"≤ {values('n_departs_max')}")]),
        ("Départs", [("Puissance  |P_départ|", f"≤ {elec['puissance_max_ligne']:g} MW"),
                     ("Longueur  L_départ", f"≤ {elec['longueur_max_depart']:g} km")]),
        ("Nœuds", [("Produit P×L", f"≤ {elec['PL_max_noeud']:g}")]),
        ("Structure", ["Réseau radial : un arbre par source", "Arbres disjoints en nœuds",
                       "Pas de passage par une autre source", "Chaque Conso / Prod est raccordée"]),
    ]
    height = ax.get_position().height * FIG_H * k * 72    # l'axe est gradué en points depuis le bas
    ax.set_xlim(0, 1)
    ax.set_ylim(0, height)
    ax.axis('off')
    y = height - 8 * k
    ax.text(0, y, "Contraintes", fontsize=12 * k, color=INK, fontweight='bold', va='center')
    y -= 24 * k
    for title, rows in sections:
        ax.text(0, y, title, fontsize=10 * k, color=INK_2, fontweight='bold', va='center')
        y -= 8 * k
        ax.plot([0, 1], [y, y], color=GRID, lw=k)
        y -= 12 * k
        for row in rows:
            if isinstance(row, tuple):
                ax.text(0, y, row[0], fontsize=10 * k, color=INK, va='center')
                ax.text(1, y, row[1], fontsize=10 * k, color=INK, va='center', ha='right', fontweight='bold')
            else:
                ax.text(0, y, row, fontsize=10 * k, color=INK, va='center')
            y -= 18 * k
        y -= 8 * k
    handles = legend_handles(nodes)
    f.legend(handles=handles, loc='lower left', bbox_to_anchor=(PANEL_RECT[0] - 0.008, PANEL_RECT[1]), frameon=False,
             fontsize=10 * k, labelcolor=INK, labelspacing=1.75, handletextpad=1.2, handlelength=2.4, borderaxespad=0)


def draw(path: Path, out: Path) -> int:
    """Dessine la carte de `path` ; renvoie le nombre d'étiquettes en collision."""
    data = json.loads(path.read_text(encoding='utf8'))
    nodes = data['nodes']

    k = max(1.0, sqrt(len(nodes) / REF_NODES))
    f = plt.figure(figsize=(FIG_W * k, FIG_H * k))
    f.patch.set_facecolor(SURFACE)
    ax = f.add_axes(GRAPH_RECT)
    ax.add_collection(LineCollection([e['geometry'] for e in data['edges']], colors=EDGE, linewidths=0.5, zorder=1))

    transit = [n['geometry'] for n in nodes if n['type'] == 'Transit']
    if transit:
        ax.scatter(*zip(*transit), s=6, color=INK_2, alpha=0.55, linewidths=0, zorder=2)

    # les plus grosses d'abord : une petite charge n'est jamais cachée par une grosse
    loads = sorted((n for n in nodes if n['type'] in ('Conso', 'Prod')), key=lambda n: -abs(n['demande']))
    for n in loads:
        x, y = n['geometry']
        if n['demande'] == 0:
            ax.scatter([x], [y], s=ZERO_AREA, facecolors=SURFACE, edgecolors=CONSO, linewidths=1, zorder=3)
        else:
            marker, color, s = marker_of(n)
            ax.scatter([x], [y], s=s, marker=marker, color=color, edgecolors=SURFACE, linewidths=1, zorder=3)

    sources = [n['geometry'] for n in nodes if n['type'] == 'Source']
    ax.scatter(*zip(*sources), s=SOURCE_AREA, marker='s', color=INK, edgecolors=SURFACE, linewidths=1, zorder=4)

    ax.set_aspect('equal')
    ax.margins(0.05)
    ax.autoscale_view()
    ax.axis('off')
    collisions = place_labels(f, ax, loads, sources, transit)

    conso = sum(n['demande'] for n in loads if n['type'] == 'Conso')
    prod = -sum(n['demande'] for n in loads if n['type'] == 'Prod')
    f.text(0.02, 0.965, path.stem, fontsize=14 * k, color=INK, ha='left', va='center', fontweight='bold')
    f.text(0.02, 0.925, f"{len(nodes)} nœuds · {len(sources)} source{'s' if len(sources) > 1 else ''} · "
                        f"{len(loads)} charges · conso {conso:g} MW · prod {prod:g} MW  —  aire ∝ |demande| · "
                        "valeur écrite = demande en MW (production < 0)",
           fontsize=10 * k, color=INK_2, ha='left', va='center')
    constraints_panel(f, f.add_axes(PANEL_RECT), nodes, data['electrotechnique'], k)
    f.savefig(out, dpi=max(100, 150 / k), facecolor=SURFACE)
    plt.close(f)
    return collisions


def main():
    motifs = sys.argv[1:] or ['Petit_', 'Moyen_']
    files = sorted(p for p in IN_DIR.glob("*.json") if any(m in p.stem for m in motifs))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    collisions = {path.stem: draw(path, OUT_DIR / f"{path.stem}.png") for path in files}
    print(f"{len(files)} images dans {OUT_DIR}")
    print(f"étiquettes en collision : {sum(collisions.values())}", {k: v for k, v in collisions.items() if v})


if __name__ == "__main__":
    main()
