# Récapitulatif du run 20261102_193405_test_verified

Objectif minimisé : J = L_tot, longueur totale des arbres (km), comme dans data/exemple_ete_4_PS/model.lp. Les contraintes (L_max et P_max par départ, capacité des sources, nombre de départs) ne sont pas dans J : elles sont seulement vérifiées. Toutes les solutions sont comparées, même celles qui violent une contrainte (infaisables : en rouge dans les tableaux, croix sur les graphiques).

PD = Σ sur les arêtes des arbres de |puissance transitée| × longueur (MW·km), non inclus dans J.

Lecture : <span style="color:#1a7f37">OK</span> contrainte satisfaite partout, <span style="color:#cf222e">KO (violées/testées)</span> sinon ; J en <span style="color:#1a7f37">vert</span> si la solution est faisable, en <span style="color:#cf222e">rouge</span> sinon, en gras pour la plus courte solution du graphe (faisable ou non) ; PD : <span style="color:#1a7f37"><b>valeur</b></span> = plus petit PD des 4 distances ; distances : <span style="color:#2a78d6"><b>non pondérée</b></span>, <span style="color:#eb6834"><b>capacité</b></span>, <span style="color:#4a3aa7"><b>racine capacité</b></span>, <span style="color:#1baf7a"><b>variable</b></span>.

## Détail par graphe et par distance

| Graphe | Type | Charges (dont <0) | Sources | Distance | J = L_tot (km) | PD (MW·km) | L_max départ | P_max départ | P_max source | Nb départs | Faisable |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Petit_N60_S2_0.20_Ete** | Petit | 12 (6) | 2 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#1a7f37"><b>104.1</b></span> | <span style="color:#1a7f37"><b>173.5</b></span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37"><b>oui</b></span> |
|  |  |  |  | <span style="color:#eb6834"><b>capacité</b></span> | <span style="color:#1a7f37">104.1</span> | 173.5 | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37"><b>oui</b></span> |
|  |  |  |  | <span style="color:#4a3aa7"><b>racine capacité</b></span> | <span style="color:#1a7f37">104.1</span> | 173.5 | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37"><b>oui</b></span> |
|  |  |  |  | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#1a7f37">104.1</span> | 173.5 | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37"><b>oui</b></span> |
| **Petit_N60_S2_0.20_Hiver** | Petit | 12 (0) | 2 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#cf222e">103.6</span> | <span style="color:#1a7f37"><b>136.8</b></span> | <span style="color:#cf222e">KO (1/2)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#eb6834"><b>capacité</b></span> | <span style="color:#cf222e">103.6</span> | 136.8 | <span style="color:#cf222e">KO (1/2)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#4a3aa7"><b>racine capacité</b></span> | <span style="color:#cf222e">103.6</span> | 136.8 | <span style="color:#cf222e">KO (1/2)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#cf222e"><b>103.2</b></span> | 136.8 | <span style="color:#cf222e">KO (1/2)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
| **Moyen_N200_S4_0.20_Ete** | Moyen | 40 (20) | 4 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#cf222e">322.2</span> | 828.2 | <span style="color:#cf222e">KO (3/10)</span> | <span style="color:#cf222e">KO (1/10)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#eb6834"><b>capacité</b></span> | <span style="color:#cf222e">322.2</span> | 828.2 | <span style="color:#cf222e">KO (3/10)</span> | <span style="color:#cf222e">KO (1/10)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#4a3aa7"><b>racine capacité</b></span> | <span style="color:#cf222e">322.2</span> | 828.2 | <span style="color:#cf222e">KO (3/10)</span> | <span style="color:#cf222e">KO (1/10)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#cf222e"><b>315.3</b></span> | <span style="color:#1a7f37"><b>703.4</b></span> | <span style="color:#cf222e">KO (2/10)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
| **Moyen_N200_S4_0.20_Hiver** | Moyen | 40 (0) | 4 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#cf222e">325.0</span> | 1119.0 | <span style="color:#cf222e">KO (1/10)</span> | <span style="color:#cf222e">KO (2/10)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#eb6834"><b>capacité</b></span> | <span style="color:#cf222e">325.0</span> | 1119.0 | <span style="color:#cf222e">KO (1/10)</span> | <span style="color:#cf222e">KO (2/10)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#4a3aa7"><b>racine capacité</b></span> | <span style="color:#cf222e">325.0</span> | 1119.0 | <span style="color:#cf222e">KO (1/10)</span> | <span style="color:#cf222e">KO (2/10)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#cf222e"><b>321.8</b></span> | <span style="color:#1a7f37"><b>1117.0</b></span> | <span style="color:#cf222e">KO (2/10)</span> | <span style="color:#cf222e">KO (1/10)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
| **Grand_N600_S6_0.20_Ete** | Grand | 120 (60) | 6 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#cf222e"><b>620.3</b></span> | 2010.9 | <span style="color:#cf222e">KO (5/11)</span> | <span style="color:#cf222e">KO (3/11)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#eb6834"><b>capacité</b></span> | <span style="color:#cf222e">620.3</span> | 2010.9 | <span style="color:#cf222e">KO (5/11)</span> | <span style="color:#cf222e">KO (3/11)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#4a3aa7"><b>racine capacité</b></span> | <span style="color:#cf222e">620.3</span> | 2010.9 | <span style="color:#cf222e">KO (5/11)</span> | <span style="color:#cf222e">KO (3/11)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#cf222e">662.2</span> | <span style="color:#1a7f37"><b>1666.6</b></span> | <span style="color:#cf222e">KO (6/12)</span> | <span style="color:#cf222e">KO (2/12)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
| **Grand_N600_S6_0.20_Hiver** | Grand | 120 (0) | 6 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#cf222e"><b>641.3</b></span> | <span style="color:#1a7f37"><b>3419.5</b></span> | <span style="color:#cf222e">KO (7/15)</span> | <span style="color:#cf222e">KO (7/15)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#eb6834"><b>capacité</b></span> | <span style="color:#cf222e">641.3</span> | 3419.5 | <span style="color:#cf222e">KO (7/15)</span> | <span style="color:#cf222e">KO (7/15)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#4a3aa7"><b>racine capacité</b></span> | <span style="color:#cf222e">641.3</span> | 3419.5 | <span style="color:#cf222e">KO (7/15)</span> | <span style="color:#cf222e">KO (7/15)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#cf222e">665.0</span> | 4756.7 | <span style="color:#cf222e">KO (6/13)</span> | <span style="color:#cf222e">KO (7/13)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
| **TrèsGrand_N2000_S12_0.20_Ete** | TrèsGrand | 400 (200) | 12 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#cf222e"><b>1610.3</b></span> | <span style="color:#1a7f37"><b>6480.7</b></span> | <span style="color:#cf222e">KO (13/25)</span> | <span style="color:#cf222e">KO (10/25)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#eb6834"><b>capacité</b></span> | <span style="color:#cf222e">1610.3</span> | 6480.7 | <span style="color:#cf222e">KO (13/25)</span> | <span style="color:#cf222e">KO (10/25)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#4a3aa7"><b>racine capacité</b></span> | <span style="color:#cf222e">1610.3</span> | 6480.7 | <span style="color:#cf222e">KO (13/25)</span> | <span style="color:#cf222e">KO (10/25)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#cf222e">1661.3</span> | 11522.3 | <span style="color:#cf222e">KO (12/26)</span> | <span style="color:#cf222e">KO (9/26)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
| **TrèsGrand_N2000_S12_0.20_Hiver** | TrèsGrand | 400 (0) | 12 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#cf222e"><b>1681.5</b></span> | <span style="color:#1a7f37"><b>14078.3</b></span> | <span style="color:#cf222e">KO (16/25)</span> | <span style="color:#cf222e">KO (20/25)</span> | <span style="color:#cf222e">KO (2/12)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#eb6834"><b>capacité</b></span> | <span style="color:#cf222e">1681.5</span> | 14078.3 | <span style="color:#cf222e">KO (16/25)</span> | <span style="color:#cf222e">KO (20/25)</span> | <span style="color:#cf222e">KO (2/12)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#4a3aa7"><b>racine capacité</b></span> | <span style="color:#cf222e">1681.5</span> | 14078.3 | <span style="color:#cf222e">KO (16/25)</span> | <span style="color:#cf222e">KO (20/25)</span> | <span style="color:#cf222e">KO (2/12)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |
|  |  |  |  | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#cf222e">1764.9</span> | 14964.1 | <span style="color:#cf222e">KO (14/25)</span> | <span style="color:#cf222e">KO (18/25)</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#1a7f37">OK</span> | <span style="color:#cf222e"><b>non</b></span> |

## Fonction objectif J = L_tot (km) et produit PD (MW·km) par distance

J : <span style="color:#1a7f37">vert</span> si toutes les contraintes sont satisfaites, <span style="color:#cf222e">rouge</span> sinon. En gras : plus petite valeur des 4 distances, solutions faisables ou non.

| Graphe | J <span style="color:#2a78d6"><b>non pondérée</b></span> | J <span style="color:#eb6834"><b>capacité</b></span> | J <span style="color:#4a3aa7"><b>racine capacité</b></span> | J <span style="color:#1baf7a"><b>variable</b></span> | PD <span style="color:#2a78d6"><b>non pondérée</b></span> | PD <span style="color:#eb6834"><b>capacité</b></span> | PD <span style="color:#4a3aa7"><b>racine capacité</b></span> | PD <span style="color:#1baf7a"><b>variable</b></span> | Meilleur J | Meilleur PD |
|---|---|---|---|---|---|---|---|---|---|---|
| Petit_N60_S2_0.20_Ete | <span style="color:#1a7f37"><b>104.1</b></span> | <span style="color:#1a7f37">104.1</span> | <span style="color:#1a7f37">104.1</span> | <span style="color:#1a7f37">104.1</span> | **173.5** | 173.5 | 173.5 | 173.5 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#2a78d6"><b>non pondérée</b></span> |
| Petit_N60_S2_0.20_Hiver | <span style="color:#cf222e">103.6</span> | <span style="color:#cf222e">103.6</span> | <span style="color:#cf222e">103.6</span> | <span style="color:#cf222e"><b>103.2</b></span> | **136.8** | 136.8 | 136.8 | 136.8 | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#2a78d6"><b>non pondérée</b></span> |
| Moyen_N200_S4_0.20_Ete | <span style="color:#cf222e">322.2</span> | <span style="color:#cf222e">322.2</span> | <span style="color:#cf222e">322.2</span> | <span style="color:#cf222e"><b>315.3</b></span> | 828.2 | 828.2 | 828.2 | **703.4** | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#1baf7a"><b>variable</b></span> |
| Moyen_N200_S4_0.20_Hiver | <span style="color:#cf222e">325.0</span> | <span style="color:#cf222e">325.0</span> | <span style="color:#cf222e">325.0</span> | <span style="color:#cf222e"><b>321.8</b></span> | 1119.0 | 1119.0 | 1119.0 | **1117.0** | <span style="color:#1baf7a"><b>variable</b></span> | <span style="color:#1baf7a"><b>variable</b></span> |
| Grand_N600_S6_0.20_Ete | <span style="color:#cf222e"><b>620.3</b></span> | <span style="color:#cf222e">620.3</span> | <span style="color:#cf222e">620.3</span> | <span style="color:#cf222e">662.2</span> | 2010.9 | 2010.9 | 2010.9 | **1666.6** | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#1baf7a"><b>variable</b></span> |
| Grand_N600_S6_0.20_Hiver | <span style="color:#cf222e"><b>641.3</b></span> | <span style="color:#cf222e">641.3</span> | <span style="color:#cf222e">641.3</span> | <span style="color:#cf222e">665.0</span> | **3419.5** | 3419.5 | 3419.5 | 4756.7 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#2a78d6"><b>non pondérée</b></span> |
| TrèsGrand_N2000_S12_0.20_Ete | <span style="color:#cf222e"><b>1610.3</b></span> | <span style="color:#cf222e">1610.3</span> | <span style="color:#cf222e">1610.3</span> | <span style="color:#cf222e">1661.3</span> | **6480.7** | 6480.7 | 6480.7 | 11522.3 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#2a78d6"><b>non pondérée</b></span> |
| TrèsGrand_N2000_S12_0.20_Hiver | <span style="color:#cf222e"><b>1681.5</b></span> | <span style="color:#cf222e">1681.5</span> | <span style="color:#cf222e">1681.5</span> | <span style="color:#cf222e">1764.9</span> | **14078.3** | 14078.3 | 14078.3 | 14964.1 | <span style="color:#2a78d6"><b>non pondérée</b></span> | <span style="color:#2a78d6"><b>non pondérée</b></span> |

## Synthèse

Écart à l'optimum : 100 × (X − X_min) / X_min, X_min = plus petite valeur des 4 distances sur le graphe ; moyenne sur les graphes (comparable entre graphes de tailles différentes, contrairement aux moyennes de J et PD).

| Distance | Meilleur J (nb graphes) | Meilleur PD (nb graphes) | Runs faisables | Écart moyen J (%) | Écart moyen PD (%) | J moyen (km) | PD moyen (MW·km) |
|---|---|---|---|---|---|---|---|
| <span style="color:#2a78d6"><b>non pondérée</b></span> | 5 | 5 | 1/8 | 0.4 | 4.8 | 676.0 | 3530.9 |
| <span style="color:#eb6834"><b>capacité</b></span> | 0 | 0 | 1/8 | 0.4 | 4.8 | 676.0 | 3530.9 |
| <span style="color:#4a3aa7"><b>racine capacité</b></span> | 0 | 0 | 1/8 | 0.4 | 4.8 | 676.0 | 3530.9 |
| <span style="color:#1baf7a"><b>variable</b></span> | 3 | 3 | 1/8 | 2.3 | 15.4 | 699.7 | 4380.0 |

## Synthèse par nombre de charges

Par groupe de graphes de même nombre de charges : runs faisables, nombre de graphes où la distance donne le plus petit J, écart moyen de J à l'optimum.

| Charges | Graphes | Faisables <span style="color:#2a78d6"><b>non pondérée</b></span> | Faisables <span style="color:#eb6834"><b>capacité</b></span> | Faisables <span style="color:#4a3aa7"><b>racine capacité</b></span> | Faisables <span style="color:#1baf7a"><b>variable</b></span> | Meilleur J <span style="color:#2a78d6"><b>non pondérée</b></span> | Meilleur J <span style="color:#eb6834"><b>capacité</b></span> | Meilleur J <span style="color:#4a3aa7"><b>racine capacité</b></span> | Meilleur J <span style="color:#1baf7a"><b>variable</b></span> | Écart J (%) <span style="color:#2a78d6"><b>non pondérée</b></span> | Écart J (%) <span style="color:#eb6834"><b>capacité</b></span> | Écart J (%) <span style="color:#4a3aa7"><b>racine capacité</b></span> | Écart J (%) <span style="color:#1baf7a"><b>variable</b></span> |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 12 | 2 | 1/2 | 1/2 | 1/2 | 1/2 | 1 | 0 | 0 | 1 | 0.2 | 0.2 | 0.2 | 0.0 |
| 40 | 2 | 0/2 | 0/2 | 0/2 | 0/2 | 0 | 0 | 0 | 2 | 1.6 | 1.6 | 1.6 | 0.0 |
| 120 | 2 | 0/2 | 0/2 | 0/2 | 0/2 | 2 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | 5.2 |
| 400 | 2 | 0/2 | 0/2 | 0/2 | 0/2 | 2 | 0 | 0 | 0 | 0.0 | 0.0 | 0.0 | 4.1 |

## Graphiques

### Écart relatif à l'optimum de J, par nombre de charges

J = L_tot. Pour chaque graphe, l'optimum est la plus petite valeur des 4 distances (0 %), faisable ou non. Barres : écart moyen ; un marqueur par graphe : rond si la solution respecte les contraintes, croix sinon.

![Écart relatif à l'optimum de J, par nombre de charges](charts/20261102_193405_test_verified_gap_J_charges.png)

### Écart relatif à l'optimum de PD, par nombre de charges

Même lecture, pour le produit PD (puissance × distance).

![Écart relatif à l'optimum de PD, par nombre de charges](charts/20261102_193405_test_verified_gap_PD_charges.png)

### Écart relatif à l'optimum de J, par type de répartition

Quelle distance convient à quelle répartition des charges.

![Écart relatif à l'optimum de J, par type de répartition](charts/20261102_193405_test_verified_gap_J_type.png)

### Contraintes violées

Nombre total de contraintes violées (L_max, P_max départ, P_max source, nb départs, erreurs de structure) pour chaque graphe et chaque méthode.

![Contraintes violées](charts/20261102_193405_test_verified_violations.png)

