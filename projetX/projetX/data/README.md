




### 1. Génération de graphe maillé de départ : *generation_graphe_maille.py*

Le graphe maillé généré se décompose en quatre parties :  **l'ensemble de noeuds**, **l'ensemble d'arêtes**, **l'ensemble d'arcs** et **les données electrotechniques**.



#### a) Ensemble de noeuds

Les noeuds sont définis dans un *DataFrame*, dont les colonnes contiennent les différents paramètres du problème associé aux noeuds. Les colonnes sont :

- *id_noeud* (dtype=int) : un indice unique permettant d'identifier le noeud.
- *type* (dtype=str) : le type du noeud parmis la liste suivante 
    - "Source" -> poste source,
    - "Conso" -> poste charge de type consommateur,
    - "Prod" -> poste charge de type producteur,
    - "Transit" -> noeud de transit.
- *id_selon_type* (dtype = int) : un indice permettant d'itérer sur de 0 à n_type -1 en parcourant tous les noeuds de ce type. Notamment utile pour la déclaration de l'association des postes charges à des départs. 
- *n_departs_max* : le nombre de départs maximal si le noeud est un poste source.
- *demande* (dtype=float) : la demande du poste source (positive si c'est un consommateur, négative si c'est un producteur) & égal à la capacité du poste source pour les neouds de ce type. S'il s'agit d'un noeud de transit, la demande vaut 0.
- *geometry* (dtype=shapely.Point) : coordonnées du noeud dans le plan
- *neighbors* (dtype=list[int]) : la liste d'indice des noeuds voisins au noeud considéré.
- *connected_edges* (dtype=list[int]) : la liste d'indice des arêtes connectées au noeud considéré.

#### b) Ensemble des arêtes

Les arêtes sont également définies dans un *DataFrame* dont les colonnes sont :
- *id_arete* (dtype=int) : l'indice unique permettant d'identifier l'arête. 
- *id_noeud_depart* (dtype=int) : l'indice du noeud de départ. 
- *id_noeud_arrivee* (dtype=int) : l'indice du noeud d'arrivée. 
- *longueur* (dtype=float) : la longueur de l'arête
- *geometry* (dtype=shapely.LineString) : la position de l'arête. 


**Remarque :** Les noeuds et arêtes sont en fait rangés dans des **Geo**DataFrame. Cela facilite l'affichage, car *geopandas* met à disposition une fonction plot permettant d'afficher les données du GeoDataFrame en respectant leurs coordoonnées, spécifiées dans la colonne "geometry" (GeoSeries) du GeoDataFrame. 

#### c) Données électrotechniques
- **"electrotechnique"** :  spécifie les limites électrotechniques du problème. Il contient :
    - *puissance_ligne* : la puissance maximale pouvant transiter sur une ligne électrique,
    - *longueur_max_depart* : la longueur maximale de lignes électriques pour un départ,
    - *PL_max_noeud* : la borne que le PL ne doit pas dépasser.

### 2. Instances vérifiées : *verified_instances/*

Chaque instance existe en `<Classe>_N<nœuds>_S<sources>_<ratio>_<Ete|Hiver>.json` (même format que ci-dessus) et en `.pdf` (dessin du graphe).

- `unfeasible/` : les 26 instances prouvées infaisables, retirées du jeu de test (voir [Expérience 3](../experiments/experiment_3_unfeasible_instances.md)). Il en reste 190 : 54 Petit, 50 Moyen, 49 Grand, 37 TrèsGrand.
- `demand_maps/` : une carte des demandes par instance (216 images, instances infaisables comprises).
- `ratios_charge_nette.csv` : bilans de puissance de chaque instance.
