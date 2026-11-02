# Network Expansion Planning — EDF

Planification d'un réseau de distribution radial. Il faut relier toutes les charges (consommateurs et producteurs) aux postes sources par une forêt d'arbres de longueur totale **J = L_tot** minimale. Contraintes :

- longueur maximale d'un départ : `l_max_depart` ;
- puissance maximale d'un départ : `p_max_depart` ;
- capacité de chaque source ;
- nombre maximal de départs par source : `n_departs_max`.

Le code est dans `projetX/projetX/` et les scripts se lancent depuis ce dossier. Les résultats de l'heuristique sont dans **[experiment_1.md](projetX/projetX/experiment_1.md)**.

---

## 1. Algorithmes, fichier par fichier

### Heuristique Voronoi + Steiner

Chaîne : `main.py` → `voronoi*.py` → `steiner.py` → `verify.py` → `output.py`.

| Fichier | Algorithme |
|---|---|
| [main.py](projetX/projetX/main.py) | Boucle sur les instances × 4 distances. Pour chaque couple : (1) test des conditions nécessaires d'infaisabilité, (2) Voronoi + Steiner, (3) vérification, (4) rapport et figure. |
| [voronoi.py](projetX/projetX/voronoi.py) | **Cellules de Voronoi sur le graphe** par un Dijkstra multi-sources, adapté de networkx. Le coût d'une arête est divisé par le poids `w_s` de la source d'où part le chemin. Poids possibles : `w_s = 1`, `w_s = capacité` ou `w_s = √capacité`. Une arête partagée entre deux cellules est coupée au point équidistant des deux sources (nœud fictif, utilisé seulement pour l'affichage). |
| [voronoi_variable_distance.py](projetX/projetX/voronoi_variable_distance.py) | **Distance variable** : ajustement itératif des poids `w_s` (60 itérations au plus). À chaque itération : Voronoi, Steiner, vérification, puis calcul de `r_s` (pire ratio valeur/limite parmi les contraintes de la source). Le poids d'une source en surcharge est divisé par `r_s^(pas·0.95^k)`, au plus par 2. Sa cellule rétrécit et ses nœuds de bord partent vers les voisines. Arrêt dès que la solution est faisable. L'algorithme est lancé depuis les 3 pondérations fixes et garde la meilleure solution : faisable d'abord, puis plus petit J. |
| [steiner.py](projetX/projetX/steiner.py) | **Arbre de Steiner dans chaque cellule** (approximation de Mehlhorn, facteur 2, via networkx). Les terminaux sont la source et les charges ; les nœuds de transit ne sont utilisés que s'ils raccourcissent l'arbre. |
| [verify.py](projetX/projetX/verify.py) | **Vérificateur indépendant**, seul juge de la faisabilité. Il contrôle la structure : arbres disjoints, sans cycle, ne passant par aucune autre source, couvrant toutes les charges. Il découpe chaque arbre en **départs** (un départ = le sous-arbre accroché à un voisin de la source), puis calcule par parcours postfixe le flux sur chaque arête, L et P par départ, P par source, J et PD. Il fournit aussi `necessary_conditions`, deux preuves d'infaisabilité : bilan de puissance global, et charge hors de portée `l_max_depart` de toutes les sources. |
| [output.py](projetX/projetX/output.py) | Rapports : `recap_*.md/.csv` (tableau par graphe et par distance, écarts relatifs à la meilleure distance, synthèse par taille), `details_*.txt` (détail par départ), graphiques matplotlib. |
| [nx_graph.py](projetX/projetX/nx_graph.py) | Classe `GraphData` : lit le JSON (nœuds, arêtes), construit des GeoDataFrames et un graphe networkx, vérifie que tout nœud est atteignable et dessine les cellules et les arbres. |

### Modèle exact (PLNE)

| Fichier | Algorithme |
|---|---|
| [plne_reseau.py](projetX/projetX/plne_reseau.py) | **PLNE « forêt de Steiner enracinée avec capacités »**, en flot mono-commodité par départ. Variables : `y[v,k]` (charge v affectée au départ k), `z[e,k]` (arête construite), `x[a,k]` (arc utilisé) et `f[a,k]` (flot). Objectif : Σ L_e·z_e. Contraintes de base : affectation unique, `P ≤ Pmax` par départ, couplage arc/arête et flot/arc (big-M), conservation du flot, au plus un arc entrant par nœud (arborescence). Contraintes optionnelles : `--ltree` (longueur du départ), `--pd` (puissance × distance), `--lmax` (fixation par plus court chemin ; `--lmax-exact` pour la version exacte par potentiels), et des inégalités valides sur les transits. Cassage de symétrie entre les départs d'une même source. Une solution initiale gloutonne est donnée au solveur : insertion au plus proche avec capacités, puis arbre couvrant minimal et élagage des feuilles. La solution est revérifiée de façon indépendante. Solveurs : HiGHS (par défaut), CBC, Gurobi, CPLEX, SCIP. |

### Données

| Fichier / dossier | Contenu |
|---|---|
| [data/generate_instances.py](projetX/projetX/data/generate_instances.py) | Générateur des 40 instances g01–g40. Positions selon 4 scénarios : **ville** (nuages gaussiens), **montagne** (routes en marche aléatoire + vallées), **épars** (uniforme), **mixte**. Graphe de Gabriel rendu connexe, chaque source reliée à ses k plus proches voisins, longueurs multipliées par un facteur de détour. Consommations log-normales, ~10 % de productions. Pour les grandes instances : sources placées par k-moyennes et capacités dimensionnées sur la charge de la zone. |
| `data/instances/` | Les 40 instances JSON, avec un aperçu PNG dans `graphics/`. |
| `data/exemple_ete_4_PS/` | Exemple de référence EDF à 4 postes sources : modèle `.lp`, solution optimale (J = 320, HiGHS). |
| `data/JCgraph/json/` | 50 instances au format `noeuds`. |
| `verified_instances/` | ~430 instances Petit (60–90 nœuds) → TrèsGrand (2000–4000 nœuds), scénarios été et hiver. |
| [src/outils.py](projetX/projetX/src/outils.py) | Utilitaires : logger, extraction de solution linopy. |

---

## 2. État d'avancement

### Fait ✅

- Générateur de 40 instances couvrant 4 typologies et 20 à 500 charges.
- Heuristique Voronoi + Steiner avec 4 distances, dont la distance variable qui corrige les surcharges.
- Vérificateur indépendant (structure, contraintes, J, PD) et conditions nécessaires d'infaisabilité.
- Campagne complète sur les 40 instances, avec rapports et graphiques : la distance variable est faisable sur 31/40 instances (voir [experiment_1.md](projetX/projetX/experiment_1.md)).
- Modèle PLNE complet, avec solution initiale gloutonne, inégalités valides et vérification de la solution.

### Reste à faire 🔲

**PLNE**

- [ ] **La lancer sur les instances g01–g20** et comparer à l'heuristique (écart à l'optimum, temps de calcul). Aucun résultat PLNE n'est encore enregistré dans `output/`.
- [ ] **Aligner le modèle sur `verify.py`**, sinon on ne compare pas les mêmes contraintes :
  - `l_max_depart` est lu dans l'instance mais pas utilisé : il faut passer `--ltree` à la main ;
  - `|P| ≤ Pmax` (`--pmax-sym`), la capacité des sources (`--cap-source`) et les départs disjoints (`--disjoint`) sont optionnels et désactivés par défaut ;
  - un départ PLNE peut quitter la source par plusieurs arcs. L'option `un_arc_source` existe mais n'est pas exposée en ligne de commande.
- [ ] Ajouter le module `jcgraph` (importé par `plne_reseau.py` pour lire `data/JCgraph/json`), qui manque dans le dépôt.

**Heuristique**

- [ ] Statuer sur les **9 instances sans solution faisable** (g02–g05, g14, g16, g19, g20, g30) : infaisabilité réelle, ou limite de l'heuristique ? La PLNE peut trancher sur les petites.
- [ ] Laisser la distance variable **agrandir** une cellule, et pas seulement la rétrécir.
- [ ] Réparer les départs trop longs ou trop chargés en déplaçant des sous-arbres entre départs, au lieu de passer uniquement par les poids de Voronoi.
- [ ] Ajouter une recherche locale après Steiner (échanges de nœuds de bord entre cellules).

**Données et validation**

- [ ] Lancer l'heuristique sur `verified_instances/` (passage à l'échelle jusqu'à 4000 nœuds, scénario hiver) après avoir vérifié la compatibilité du format avec `GraphData`.
- [ ] Comparer l'heuristique à l'optimum connu de `data/exemple_ete_4_PS` (J = 320).
- [ ] Les scripts de `exemple_monosource/` (glouton et PLNE) n'existent plus qu'en `.pyc` dans `__pycache__`. Il faut les restaurer ou supprimer l'exemple.

---

## 3. Utilisation

```bash
pip install networkx geopandas shapely pandas numpy matplotlib   # heuristique
pip install "pulp>=2.8,<4" highspy                               # PLNE

python main.py --petites --label petites     # heuristique sur g01–g20 (--grandes : g21–g40)
python main.py g05 --show                    # une instance, figures à l'écran
python plne_reseau.py data/instances/g06_ville_c20_s3.json --ltree 15 --pmax-sym --cap-source --disjoint --plot sol.png
python data/generate_instances.py            # régénère les 40 instances
```

Sorties de `main.py` dans `output/` : `recap_<date>_<label>.md/.csv`, `details_*.txt`, `charts/`, `figures/`.
