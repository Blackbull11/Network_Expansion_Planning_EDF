# Expérimentation 1 — Heuristique Voronoi + Steiner sur 40 instances

> Runs du 4 octobre 2026 (`main.py`). Ce fichier en résume l'essentiel. Les 40 instances `data/instances/`, les résultats bruts de `output/` et `model.lp` ont été retirés du dépôt au commit 1ec1776 (retrouvables dans l'historique git, `data/generate_instances.py` régénère les instances).

## 1. Objectif

Construire rapidement un réseau de distribution radial (une forêt d'arbres, un arbre par poste source) qui relie toutes les charges, puis mesurer :

- **J = L_tot** : longueur totale des arbres (km). C'est l'objectif, comme dans `data/exemple_ete_4_PS/model.lp`.
- **PD** = Σ |puissance transitée| × longueur (MW·km). Indicateur secondaire, hors de J.
- **Faisabilité** : contraintes dures, vérifiées a posteriori par [verify.py](../verify.py) :

| Contrainte | Portée |
|---|---|
| `L_départ ≤ l_max_depart` | chaque départ (sous-arbre accroché à un voisin de la source) |
| `\|P_départ\| ≤ p_max_depart` | chaque départ |
| `\|P_source\| ≤ capacité` | chaque source |
| `nb départs ≤ n_departs_max` | chaque source |

## 2. Méthode

1. **Cellules de Voronoi** sur le graphe : chaque nœud va à la source la plus « proche », avec un coût d'arête divisé par un poids `w_s` de la source ([voronoi.py](../voronoi.py)).
2. **Arbre de Steiner** (Mehlhorn, 2-approximation) dans chaque cellule, qui relie la source et ses charges ([steiner.py](../steiner.py)).
3. **Vérification** de la solution (structure + contraintes) par [verify.py](../verify.py).

Quatre choix de distance (`w_s`) sont comparés :

| Distance | Poids `w_s` |
|---|---|
| non pondérée | 1 |
| capacité | capacité de la source |
| racine capacité | √capacité |
| **variable** | poids ajustés itérativement : on part des 3 précédents, et à chaque itération on rétrécit la cellule des sources en surcharge (`w_s ← w_s / r_s^pas`, `r_s` = pire ratio valeur/limite). Arrêt dès que la solution est faisable ([voronoi_variable_distance.py](../voronoi_variable_distance.py)). |

## 3. Instances

40 graphes générés par [data/generate_instances.py](../data/generate_instances.py) (graphe de Gabriel, ~10 % de producteurs) :

| Lot | Instances | Charges | Sources |
|---|---|---|---|
| Petites | g01 – g20 | 20 à 100 | 2 à 8 |
| Grandes | g21 – g40 | 200 à 500 | 13 à 33 |

Quatre répartitions des charges : **ville** (concentrées), **montagne** (le long de vallées, graphe quasi linéaire), **épars** (rural dispersé), **mixte** (ville + campagne).

## 4. Historique des runs

| Run (`output/recap_<run>.md`) | Instances | Distances | Remarque |
|---|---|---|---|
| `20261004_120244`, `_120650` | g01–g20 | 3 fixes | premiers essais, seulement des graphiques dans `charts/` |
| `20261004_122324` → `_124751` (×4) | g01–g20 | 3 fixes | mêmes résultats (CSV identiques), seule la mise en forme du rapport change |
| `20261004_131243_grandes` | g21–g40 | 3 fixes | |
| `20261004_131451_agrege` | g01–g40 | 3 fixes | synthèse des 40 instances |
| `20261004_152826_test` | g01 | 4 | test de la distance variable |
| **`20261004_152958_petites_variable`** | g01–g20 | 4 | **run de référence, petites** |
| **`20261004_215518_grandes_variable`** | g21–g40 | 4 | **run de référence, grandes** |

Les chiffres ci-dessous viennent des deux runs de référence.

## 5. Résultats

### 5.1 Faisabilité (nombre d'instances où la solution respecte toutes les contraintes)

| Distance | Petites (/20) | Grandes (/20) | **Total (/40)** |
|---|---|---|---|
| non pondérée | 7 | 13 | 20 |
| capacité | 5 | 4 | 9 |
| racine capacité | 7 | 7 | 14 |
| **variable** | **12** | **19** | **31** |

- La distance variable trouve une solution faisable **dès qu'une des 4 méthodes en trouve une** (31/31), et elle est la seule à en trouver une sur 7 instances : g10, g28, g31, g33, g35, g36, g39.
- **9 instances restent sans solution faisable** : g02, g03, g04, g05, g14, g16, g19, g20, g30. Les conditions nécessaires de `verify.necessary_conditions` ne prouvent l'infaisabilité d'aucune d'elles, donc on ne sait pas si c'est l'heuristique ou l'instance qui bloque.

Par type d'instance (distance variable) : ville 7/11, montagne 9/10, épars 8/10, mixte 7/9.

### 5.2 Qualité de J (écart moyen à la plus petite valeur des 4 distances, faisable ou non)

| Distance | Écart J petites | Écart J grandes | Écart PD petites | Écart PD grandes |
|---|---|---|---|---|
| non pondérée | 3.5 % | 1.0 % | 4.7 % | **1.4 %** |
| capacité | **1.0 %** | 1.8 % | 10.9 % | 8.1 % |
| racine capacité | 2.9 % | 0.9 % | 5.3 % | 3.6 % |
| variable | 6.7 % | **0.8 %** | 10.3 % | 2.8 % |

Lecture :

- **Petites instances** : la pondération par la capacité donne souvent le J le plus court (12/20), mais au prix de violations (5/20 faisables seulement). La distance variable paie sa faisabilité par des arbres ~7 % plus longs.
- **Grandes instances** : les écarts sont tous faibles (< 2 %). La distance variable est à la fois la plus souvent faisable et la plus proche de l'optimum.
- **Parmi les solutions faisables**, la distance variable donne toujours le plus petit J (31/31 instances), et le gagne strictement sur 7 d'entre elles (jusqu'à −2.6 % sur g15). C'est attendu : elle part des 3 distances fixes et garde la meilleure solution faisable.
- PD est minimisé le plus souvent par la distance non pondérée. La pondération par la capacité le dégrade nettement (+8 à +11 %).

### 5.3 Contraintes violées (nombre d'instances, sur 40)

| Distance | L_max départ | P_max départ | P_max source | Nb départs |
|---|---|---|---|---|
| non pondérée | 14 | 14 | 10 | 1 |
| capacité | 26 | 21 | 14 | 1 |
| racine capacité | 20 | 16 | 12 | 0 |
| variable | 6 | 6 | 2 | 0 |

Les contraintes **par départ** (longueur et puissance) sont les plus souvent violées. Agrandir la cellule des grosses sources (pondération par la capacité) allonge les départs, d'où le plus grand nombre de violations de `L_max`.

## 6. Conclusions

1. **La distance variable est la meilleure méthode** : 31/40 instances faisables contre 20/40 pour la meilleure distance fixe, sans perte sur J parmi les solutions faisables.
2. Les distances fixes ne se valent pas selon le critère : la capacité raccourcit J mais casse les contraintes de départ, la distance non pondérée est la plus sûre et donne le meilleur PD.
3. Les petites instances sont les plus difficiles (12/20 faisables contre 19/20 pour les grandes) : peu de sources, donc peu de marge pour réaffecter les charges.

## 7. Limites et suites

- Voronoi + Steiner ne choisit pas les départs : ils découlent de la forme de l'arbre. Une violation de `L_max` ou `P_max` sur un départ ne peut se corriger qu'en déplaçant la frontière des cellules.
- Pas de borne inférieure : les écarts sont mesurés entre heuristiques, pas à l'optimum. Prochaine étape : comparer à la PLNE ([plne_reseau.py](../plne_reseau.py)) sur les petites instances, et statuer sur les 9 instances sans solution.
- Le mode variable rétrécit seulement les cellules en surcharge. Il ne peut pas agrandir une cellule ni découper un départ trop long.

## 8. Reproduire

```bash
python main.py --petites --label petites_variable   # g01–g20
python main.py --grandes --label grandes_variable   # g21–g40
python main.py g05 --show                           # une instance, figures à l'écran
```

| Sortie | Contenu |
|---|---|
| `output/recap_<run>.md` / `.csv` | tableau par graphe et distance + synthèses |
| `output/details_<run>.txt` | détail par source et par départ (L, P, PD, poids `w`) |
| `output/charts/<run>_*.png` | écarts à l'optimum (par taille, par type), violations |
| `output/figures/<graphe>_<distance>.png` | cellules de Voronoi et arbres de Steiner |
