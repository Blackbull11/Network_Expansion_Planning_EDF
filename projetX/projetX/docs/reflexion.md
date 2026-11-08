# Réflexion : transcription des tableaux

Transcription fidèle des tableaux, dans leur langue d'origine (français et anglais mélangés). Les lectures incertaines sont marquées `[?]` et listées à la fin.

## 1. Notations et formulation

**Schéma du réseau d'exemple.** Graphe à 12 sommets, avec le sommet 1 (carré) comme source. Annotations : « longueur » → $L_{48}$ sur l'arête 4–8, « demande » → $D_8$ sur le sommet 8, et $D_{12} =$ sur le sommet 12.

- Sources : $V^s = \{1\}$
- Charges : $V^c = \{2, 5, 6, 7, 8, 10, 12\}$ `[?]`
- Transit : $V^t = \{3, 9, 11\}$

**Ensembles.**

- $E$ : set of edges
- $V = V^s \cup V^c \cup V^t$ : set of vertices
- $C = \bigcup_{s \in V^s} C_s$ : set of available feeders

**Programme (encadré).**

$$\min \sum_{k \in C,\, e \in E} L_e\, z_e^k$$

s.c.

1. $\forall k, \quad \sum_{v \in V^c} D_v\, y_v^k \le P_{\max}^k$
2. $\forall k, \quad \sum_{e \in E} L_e\, z_e^k \le L_{\max}^k$
3. $\sum_{k \in C} y_v^k = 1 \quad \forall v \in V^c$
4. $\sum_{k \in C} \sum_{v \in V^c} D_v\, \delta_v^k\, y_v^k \le PD$ `[?]`

## 2. Calculer des minorants

- Pas de dual pour le PLNE exact
  - ↳ Dual pour la relaxation
    - ↳ Branch and Bound sur la relax° ?
- Dreyfus-Wagner algo : sous-pb dépendant uniquement de $(S \subseteq V^c, s)$ `[?]`
- Set-partitioning PL → 1 var par feeder (autorisé) — annotation sous « autorisé » : *? Wdym ?*
  - ↳ ~11k allowed feeders, ~33k groupes
  - ⊕ DUAL : $\max \sum_i \pi_{v_i}$ s.c. $\sum_{i \in I} \pi_{v_i} \le \text{cost}(S) \quad \forall S$ allowed feeder

**Encadré à droite.**

∀ valid network :

$$\text{length} \ge \sum_{S \text{ allowed feeder}} \text{cost}(S) \ge \sum_{i \in V^c} \pi_i$$

- $\text{cost}(S)$ est entouré : *computed by D-W*
- Sous $\sum_{i \in V^c} \pi_i$ : *sur toutes les charges*
- $\pi$ = coût

## 3. Gurobi

- Hard step : proving optimality
  - ↳ Pb : Dyn. Prog in $O(3^n)$ (↳ partitioning) → Useless over 20 loads
    - ↗ Challenge : Use DP on subproblems to get a fast lower bound
- Accelerates with feasible network as input
  - ↳ Challenge : find a feasible solution (the more optimal the better but any admissible solution helps)
- ↗ Efficient clustering in subproblems

## 4. Clustering into subproblems

$L_{\max}^{\text{tree}} = 50$

→ TDA methods or hierarchic reduction ?

- ↳ TDA not so efficient apparently

**Schéma.** Réseau maillé avec une source (carré) en bas, entouré d'un contour qui délimite un cluster autour de la source.

- Frontier at $\delta_v^k = L_{\max}^{\text{tree}}$
- ⇒ Natural exact subproblem
- → How to deal with nodes in 2 clusters ?
  - → Posterior treatment

## 5. Need to identify unfeasible problems

*(« identify unfeasible problems » encadré, précédé d'un symbole entouré `[?]`)*

⚠ unfeasible ⇒ time limit

↳ if $S \subset V^c$ intersects with one unique cluster AND $\sum_{v \in S} D_v \ge P_{\max}^{\text{tree}}$, then UNFEASIBLE

↗ Exact check taking more intersecting clusters is hard → Different level of generalisation possible

i.e. clusters = {feeder_id : [list of reachable charges]}

\+ feeders = *(laissé vide au tableau)*

## 6. Unfeasability checking

*(orthographe du tableau : « UNFEASABILITY CHECKING »)*

**Schéma.** Graphe biparti entre 10 charges et 3 boîtes (feeders) de capacité 8 chacune.

| Charge | c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8 | c9 | c10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Demande | 3 | 2 | 3 | 4 | 2 | 3 | 4 | 2 | 3 | 3 |

Liens dessinés `[?]` :

- Boîte 1 (capacité 8) : c2, c3, c5
- Boîte 2 (capacité 8) : c1, c4, c7
- Boîte 3 (capacité 8) : c5, c6, c8
- c9 et c10 : aucun lien dessiné

↳ Pb d'attribution des charges dans des boîtes de capacité limité ∀ ensemble de feeders dont `[?]` qui partagent des charges atteignables

## Lectures incertaines

| Passage | Ce qui est incertain |
|---|---|
| $V^c$ | Le tableau semble indiquer « 9 » et « K » en 5ᵉ et 7ᵉ positions. Lu 8 et 12, car 9 est déjà dans $V^t$ et le schéma porte $D_8$ et $D_{12}$. |
| $V^t$ | Le sommet 4 du schéma n'apparaît dans aucun des trois ensembles tels qu'ils sont écrits. |
| Contrainte 4 | La fin de la ligne chevauche le bord de l'encadré. Le membre de droite se lit « PD », peut-être $PD_{\max}$. |
| Dreyfus-Wagner | Le sous-problème peut se lire $(S_I \subseteq V^c, s)$. |
| Section 5 | Le symbole entouré devant « Need to identify » ressemble à ⊕ ou à un 4. Une petite paire d'ovales avec un « ? » est dessinée à côté et n'a pas été interprétée. |
| Section 6, liens | Reconstitués à partir de lignes qui se croisent sur la photo. |
| Section 6, dernière phrase | « dont » semble remplacé par « qui », la zone est partiellement effacée. |
