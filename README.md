# Network Expansion Planning — EDF

Planification d'un réseau de distribution radial. Il faut relier toutes les charges (consommateurs et producteurs) aux postes sources par une forêt d'arbres de longueur totale **J = L_tot** minimale. Contraintes :

- longueur maximale d'un départ : `l_max_depart` ;
- puissance maximale d'un départ : `p_max_depart` ;
- capacité de chaque source ;
- nombre maximal de départs par source : `n_departs_max`.

Le code est dans `projetX/projetX/` et les scripts se lancent depuis ce dossier. Résultats de l'heuristique :

- sur les 40 instances g01–g40 : **[experiment_1.md](projetX/projetX/experiment_1.md)** ;
- test sur 8 instances de `verified_instances/` : **[output/test_verified.md](projetX/projetX/output/test_verified.md)**.

---

## 1. Algorithmes, fichier par fichier

### Heuristique Voronoi + Steiner

Chaîne : `main.py` (ou `run_verified.py`) → `voronoi*.py` → `steiner.py` → `verify.py` → `output.py`.

| Fichier | Algorithme |
|---|---|
| [main.py](projetX/projetX/main.py) | Boucle sur les instances de `data/instances/` × 4 distances. Pour chaque couple : (1) test des conditions nécessaires d'infaisabilité, (2) Voronoi + Steiner, (3) vérification, (4) rapport et figure. |
| [run_verified.py](projetX/projetX/run_verified.py) | Même boucle sur les instances de `verified_instances/`, avec les métadonnées de `ratios_charge_nette.csv` (taille, densité, saison, charge nette). Ajoute, par instance et par distance : le temps de calcul, le flux maximal sur une ligne, le PD maximal d'un départ. Sorties : les rapports de `main.py`, plus `verified_<date>.csv` (écrit au fil du run) et `solutions_<date>.json` (arbres). Pas de figure par instance par défaut (`--figures variable` ou `all` pour en avoir). |
| [voronoi.py](projetX/projetX/voronoi.py) | **Cellules de Voronoi sur le graphe** par un Dijkstra multi-sources, adapté de networkx. Le coût d'une arête est divisé par le poids `w_s` de la source d'où part le chemin. Poids possibles : `w_s = 1`, `w_s = capacité` ou `w_s = √capacité`. Une arête partagée entre deux cellules est coupée au point équidistant des deux sources (nœud fictif, utilisé seulement pour l'affichage). |
| [voronoi_variable_distance.py](projetX/projetX/voronoi_variable_distance.py) | **Distance variable** : ajustement itératif des poids `w_s` (60 itérations au plus). À chaque itération : Voronoi, Steiner, vérification, puis calcul de `r_s` (pire ratio valeur/limite parmi les contraintes de la source). Le poids d'une source en surcharge est divisé par `r_s^(pas·0.95^k)`, au plus par 2. Sa cellule rétrécit et ses nœuds de bord partent vers les voisines. Arrêt dès que la solution est faisable. L'algorithme est lancé depuis les 3 pondérations fixes et garde la meilleure solution : faisable d'abord, puis plus petit J. |
| [steiner.py](projetX/projetX/steiner.py) | **Arbre de Steiner dans chaque cellule** (approximation de Mehlhorn, facteur 2, via networkx). Les terminaux sont la source et les charges ; les nœuds de transit ne sont utilisés que s'ils raccourcissent l'arbre. Seule la longueur est minimisée : la source est un terminal comme les autres, et les départs ne sont pas construits mais constatés ensuite par `verify.py` (voir § 3). |
| [verify.py](projetX/projetX/verify.py) | **Vérificateur indépendant**, seul juge de la faisabilité. Il contrôle la structure : arbres disjoints, sans cycle, ne passant par aucune autre source, couvrant toutes les charges. Il découpe chaque arbre en **départs** (un départ = le sous-arbre accroché à un voisin de la source), puis calcule par parcours postfixe le flux sur chaque arête, L et P par départ, P par source, J et PD. Il fournit aussi `necessary_conditions`, deux preuves d'infaisabilité : bilan de puissance global, et charge hors de portée `l_max_depart` de toutes les sources. |
| [output.py](projetX/projetX/output.py) | Rapports : `recap_*.md/.csv` (tableau par graphe et par distance, écarts relatifs à la meilleure distance, synthèse par taille), `details_*.txt` (détail par départ), graphiques matplotlib. Pour `verified_instances`, le « type » des graphiques est la classe de taille (Petit → TrèsGrand). |
| [nx_graph.py](projetX/projetX/nx_graph.py) | Classe `GraphData` : lit le JSON (nœuds, arêtes), construit des GeoDataFrames et un graphe networkx, vérifie que tout nœud est atteignable et dessine les cellules et les arbres. Lit aussi le format `verified_instances` (bloc `electrotechnique`, voir ci-dessous). |

### Modèle exact (PLNE)

| Fichier | Algorithme |
|---|---|
| [plne_reseau.py](projetX/projetX/plne_reseau.py) | **PLNE « forêt de Steiner enracinée avec capacités »**, en flot mono-commodité par départ. Variables : `y[v,k]` (charge v affectée au départ k), `z[e,k]` (arête construite), `x[a,k]` (arc utilisé) et `f[a,k]` (flot). Objectif : Σ L_e·z_e. Contraintes de base : affectation unique, `P ≤ Pmax` par départ, couplage arc/arête et flot/arc (big-M), conservation du flot, au plus un arc entrant par nœud (arborescence). Contraintes optionnelles : `--ltree` (longueur du départ), `--pd` (puissance × distance), `--lmax` (fixation par plus court chemin ; `--lmax-exact` pour la version exacte par potentiels), et des inégalités valides sur les transits. Cassage de symétrie entre les départs d'une même source. Une solution initiale gloutonne est donnée au solveur : insertion au plus proche avec capacités, puis arbre couvrant minimal et élagage des feuilles. La solution est revérifiée de façon indépendante. Solveurs : HiGHS (par défaut), CBC, Gurobi, CPLEX, SCIP. |

### Données

| Fichier / dossier | Contenu |
|---|---|
| [data/generate_instances.py](projetX/projetX/data/generate_instances.py) | Générateur des 40 instances g01–g40. Positions selon 4 scénarios : **ville** (nuages gaussiens), **montagne** (routes en marche aléatoire + vallées), **épars** (uniforme), **mixte**. Graphe de Gabriel rendu connexe, chaque source reliée à ses k plus proches voisins, longueurs multipliées par un facteur de détour. Consommations log-normales, ~10 % de productions. Pour les grandes instances : sources placées par k-moyennes et capacités dimensionnées sur la charge de la zone. |
| `data/instances/` | Les 40 instances JSON, avec un aperçu PNG dans `graphics/`. Limites stockées sur chaque source (`p_max_depart`, `l_max_depart`). |
| `data/exemple_ete_4_PS/` | Exemple de référence EDF à 4 postes sources : modèle `.lp`, solution optimale (J = 320, HiGHS). |
| `data/JCgraph/json/` | 50 instances au format `noeuds`. |
| `verified_instances/` | **214 instances** `<taille>_N<nœuds>_S<sources>_<densité>_<saison>.json` : Petit (60–90 nœuds), Moyen (150–250), Grand (500–700), TrèsGrand (2000–4000) ; S de 1 à 18 sources ; densité 0,15 / 0,20 / 0,25 ; été (consommateurs et producteurs) et hiver (consommateurs seuls). Chaque JSON a son aperçu `.pdf` (graphe d'entrée). `ratios_charge_nette.csv` : consommation, production, charge nette, capacités et ratios de chaque instance. `TrèsGrand_N4000_S18_0.25` (été et hiver) n'existe qu'en PDF. |
| [src/outils.py](projetX/projetX/src/outils.py) | Utilitaires : logger, extraction de solution linopy. |

### Format `verified_instances` : le bloc `electrotechnique`

Ces instances, comme l'exemple EDF, ne stockent pas les limites sur chaque source mais dans un bloc global du JSON :

```json
"electrotechnique": {"puissance_max_ligne": 8, "longueur_max_depart": 50, "PL_max_noeud": 200}
```

`GraphData` recopie ces valeurs sur chaque source au chargement, et prend `demande` comme `demande_hivers` (une seule saison par fichier). `verify.py` et le reste du code sont inchangés.

| Paramètre | Valeur | Utilisé comme | Contrôlé par `verify.py` |
|---|---|---|---|
| `demande` d'une Source | 70 MW | capacité de la source | oui |
| `n_departs_max` | 8 | nombre de départs max par source | oui |
| `puissance_max_ligne` | 8 MW | `p_max_depart` : \|P_départ\| ≤ 8 MW (P_max du sujet, § 2.4) | oui, sur la ligne de tête du départ |
| `longueur_max_depart` | 50 km | `l_max_depart` : longueur totale des lignes d'un départ ≤ 50 km | oui |
| `PL_max_noeud` | 200 | budget puissance × longueur | non, indicateur seulement (PD de chaque départ) |

Toutes les sources ont la même capacité (70 MW) : les trois distances fixes (non pondérée, capacité, racine capacité) donnent les mêmes cellules et les mêmes arbres sur ces instances.

---

## 2. État d'avancement

### Fait ✅

- Générateur de 40 instances couvrant 4 typologies et 20 à 500 charges.
- Heuristique Voronoi + Steiner avec 4 distances, dont la distance variable qui corrige les surcharges.
- Vérificateur indépendant (structure, contraintes, J, PD) et conditions nécessaires d'infaisabilité.
- Campagne complète sur les 40 instances, avec rapports et graphiques : la distance variable est faisable sur 31/40 instances (voir [experiment_1.md](projetX/projetX/experiment_1.md)).
- Modèle PLNE complet, avec solution initiale gloutonne, inégalités valides et vérification de la solution.
- Environnement Python reproductible : [requirements.txt](requirements.txt) et `.gitignore` (voir § 4).
- **`verified_instances/` remis en ordre.** Les fichiers téléchargés ne contenaient pas l'instance de leur nom : chaque fichier contenait l'instance nommée deux fichiers avant lui (décalage d'une paire été/hiver), le tableau des ratios était stocké sous un nom `.json` et deux `.json` étaient des PDF. Les fichiers ont été renommés d'après leur contenu (N, S, consommation, production comparés à `ratios_charge_nette.csv` ; titre « N nœuds et E arêtes » pour les PDF) et le suffixe « (1) » a été retiré.
- Lecture du format `verified_instances` par `GraphData` et script [run_verified.py](projetX/projetX/run_verified.py).
- **Test sur 8 instances** (une par taille, densité 0,20, été et hiver) : 1 seule faisable (`Petit_N60_S2_0.20_Ete`). Les sources ne sont jamais saturées, mais les départs dépassent 50 km et 8 MW (voir [output/test_verified.md](projetX/projetX/output/test_verified.md) et § 3).

### Reste à faire 🔲

**PLNE**

- [ ] **La lancer sur les instances g01–g20** et comparer à l'heuristique (écart à l'optimum, temps de calcul). Aucun résultat PLNE n'est encore enregistré dans `output/`.
- [ ] **Aligner le modèle sur `verify.py`**, sinon on ne compare pas les mêmes contraintes :
  - `l_max_depart` est lu dans l'instance mais pas utilisé : il faut passer `--ltree` à la main ;
  - `|P| ≤ Pmax` (`--pmax-sym`), la capacité des sources (`--cap-source`) et les départs disjoints (`--disjoint`) sont optionnels et désactivés par défaut ;
  - un départ PLNE peut quitter la source par plusieurs arcs. L'option `un_arc_source` existe mais n'est pas exposée en ligne de commande.
- [ ] Ajouter le module `jcgraph` (importé par `plne_reseau.py` pour lire `data/JCgraph/json`), qui manque dans le dépôt.

**Heuristique** (détails au § 3)

- [ ] Statuer sur les **9 instances sans solution faisable** (g02–g05, g14, g16, g19, g20, g30) : infaisabilité réelle, ou limite de l'heuristique ? La PLNE peut trancher sur les petites.
- [ ] **Scinder les départs trop gros en sous-départs** (§ 3.1).
- [ ] Réparer les départs trop longs ou trop chargés en déplaçant des sous-arbres entre départs, au lieu de passer uniquement par les poids de Voronoi (§ 3.2).
- [ ] Construire des arbres enracinés qui tiennent compte des départs dès la construction (§ 3.3).
- [ ] Laisser la distance variable **agrandir** une cellule, et pas seulement la rétrécir ; ne pas relancer des points de départ identiques (§ 3.4).
- [ ] Ajouter une recherche locale après Steiner (échanges de nœuds de bord entre cellules).

**Données et validation**

- [ ] Lancer l'heuristique sur **toutes** les instances de `verified_instances/` (212 restantes, jusqu'à 4000 nœuds) et en faire la synthèse façon `experiment_1.md`.
- [ ] Comparer l'heuristique à l'optimum connu de `data/exemple_ete_4_PS` (J = 320).
- [ ] Les scripts de `exemple_monosource/` (glouton et PLNE) n'existent plus qu'en `.pyc` dans `__pycache__`. Il faut les restaurer ou supprimer l'exemple.

---

## 3. Pistes d'amélioration de l'heuristique

**Diagnostic.** Sur `verified_instances`, les violations portent sur les **départs**, pas sur les sources : chaque source reste sous sa capacité, mais n'utilise que 1 à 3 départs sur 8. Exemple (`Grand_N600_S6_0.20_Hiver`) : une source a un seul départ de 110 km et 22 MW, pour des limites de 50 km et 8 MW. La cause est dans Steiner (Mehlhorn) : il minimise la longueur seule, la source y est un terminal comme les autres et se relie à son terminal le plus proche ; le nombre de départs n'est que le nombre de voisins de la source dans l'arbre. La distance variable ne peut pas corriger cela : elle ne fait que rétrécir des cellules.

Point d'entrée commun aux pistes 3.1 à 3.3 : `steiner_trees_by_cells(G, cells)` dans [steiner.py](projetX/projetX/steiner.py), qui renvoie `{source: [(u, v), ...] avec u < v}` (un arbre par cellule). Tout ce qui respecte ce contrat fonctionne avec les distances fixes, la distance variable, `verify.py` et les rapports.

### 3.1 Scinder les départs en sous-départs (arêtes partagées)

Comme dans la formulation du sujet (variable `z_eκ` par départ), une même arête peut porter les lignes de plusieurs départs. Un départ trop gros peut donc être scindé en **sous-départs** qui partagent le début de son tracé, chacun respectant `P_max` et `L_max`. Méthode, en post-traitement des arbres de Steiner (sans modifier le code existant) :

1. Enraciner l'arbre de chaque source en s. Pour chaque départ, lister ses charges dans l'ordre d'un parcours en profondeur depuis sa tête : des charges proches dans l'arbre sont consécutives.
2. Un sous-départ = un intervalle [i..j] de cette liste. Ses lignes sont l'union des chemins s → v de l'arbre pour v dans l'intervalle. Sa puissance est la somme des demandes de **ses** charges (pas de celles qu'il traverse sans les desservir) ; sa longueur se calcule de façon incrémentale (en ajoutant v, remonter vers s jusqu'au premier nœud déjà couvert). La longueur croît avec j, la puissance non (producteurs et consommateurs se compensent).
3. Découpage optimal par programmation dynamique (procédure « Split » de Prins pour les tournées) : `C[j][k]` = longueur minimale pour desservir les j premières charges avec k sous-départs admissibles (|P| ≤ 8 MW, L ≤ 50 km). Un sous-départ ne mélange pas deux départs d'origine, et le total par source reste ≤ `n_departs_max`.
4. Ce que la scission ne peut pas corriger : une charge dont le chemin dans l'arbre dépasse déjà 50 km, une source au-delà de sa capacité, un besoin de plus de 8 sous-départs. Ces cas sont détectés par la DP et relèvent d'une autre piste.

Coût : chaque sous-départ supplémentaire paie à nouveau le tronc qu'il partage (objectif Σ_κ Σ_e L_e z_eκ du sujet). Variante : pour chaque intervalle, comparer ce coût à celui d'un rebranchement direct sur la source par un plus court chemin dans la cellule (piste 3.2), et garder le moins cher.

Vérification : les arbres de chaque run sont déjà enregistrés dans `output/solutions_<date>.json`. La longueur de chaque sous-départ peut être contrôlée avec `verify.departs_of_tree` ; la puissance doit être recalculée sur les charges affectées, et le modèle « nœuds disjoints » de `check_structure` ne s'applique plus tel quel.

### 3.2 Rebrancher des sous-arbres sur la source

Pour un départ trop long ou trop chargé, couper un sous-arbre et le relier directement à la source par un plus court chemin qui ne touche pas les autres départs. Cela ouvre un nouveau départ (5 à 7 restent libres par source) en gardant des départs disjoints, au prix d'un peu de longueur.

### 3.3 Construire des arbres enracinés avec capacités

Le problème est proche de l'**arbre couvrant minimal avec capacités** (CMST) : chaque sous-arbre accroché à la racine a une capacité limitée. L'heuristique classique d'**Esau-Williams** construit directement de tels arbres. Il faut l'adapter à deux limites (P ≤ 8 MW, L ≤ 50 km), au plafond de 8 départs, et aux nœuds de transit (Steiner plutôt qu'arbre couvrant).

### 3.4 Distance variable

- Lui permettre d'**agrandir** une cellule (une source sous-chargée attire les nœuds de bord de ses voisines), et pas seulement de rétrécir les cellules en surcharge.
- Quand les sources ont toutes la même capacité, les trois points de départ (poids 1, capacité, √capacité) sont proportionnels : les cellules et toute la suite des itérations sont identiques. Ne pas relancer un point de départ proportionnel à un autre diviserait le temps de calcul par 3 à résultats identiques (non implémenté).

### 3.5 Contraintes non contrôlées

- `puissance_max_ligne` n'est contrôlée que sur la ligne de tête de chaque départ. En été, une ligne interne peut porter plus (producteurs et consommateurs se compensent en tête) : voir les colonnes `max_flow` et `lines_over` de `verified_*.csv`.
- `PL_max_noeud` n'est pas une contrainte. Deux lectures possibles : par départ (Σ |flux| × longueur, le budget puissance-distance du sujet, § 2.7, mesuré aujourd'hui) ou par nœud (Σ |flux| × longueur sur le chemin source → nœud, approximation de la chute de tension, ce que suggère le nom).

---

## 4. Utilisation

Environnement (Python 3.13, depuis la racine du dépôt) :

```bash
uv venv --python 3.13 .venv                    # ou : python -m venv .venv
uv pip install --python .venv/Scripts/python.exe -r requirements.txt   # ou : .venv/Scripts/pip install -r requirements.txt
.venv\Scripts\Activate.ps1                     # Windows PowerShell (Linux/macOS : source .venv/bin/activate)
```

Puis depuis `projetX/projetX/` :

```bash
python main.py --petites --label petites     # heuristique sur g01–g20 (--grandes : g21–g40)
python main.py g05 --show                    # une instance, figures à l'écran
python run_verified.py --label verified      # heuristique sur verified_instances/ (sans figures)
python run_verified.py Petit_N60_S2_0.20 --figures all --label test   # quelques instances, avec figures
python plne_reseau.py data/instances/g06_ville_c20_s3.json --ltree 15 --pmax-sym --cap-source --disjoint --plot sol.png
python data/generate_instances.py            # régénère les 40 instances
```

Sorties dans `output/` : `recap_<date>_<label>.md/.csv`, `details_*.txt`, `charts/`, `figures/` ; plus `verified_*.csv` et `solutions_*.json` pour `run_verified.py`.
