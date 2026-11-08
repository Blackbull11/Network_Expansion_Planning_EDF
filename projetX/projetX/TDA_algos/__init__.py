"""
Topological Data Analysis tools for the network planning heuristic.

Modules (each one is independent of the others except through persistence.py):
  union_find   : union-find with path halving, used by every filtration
  persistence  : H0 of vertex / edge filtrations of a graph, H0 of the terminals
                 under the shortest-path metric (Mehlhorn), H1 of the planarised
                 graph (primal / dual union-find)
  vectorize    : diagram summaries, persistence images, Betti curves,
                 Wasserstein and bottleneck distances
  clustering   : ToMATo on the graph (persistence-based clustering of the power density)
  merge_tree   : split tree of the distance to a source inside a cell
  bounds       : lower bound on J and necessary conditions of feasibility
  spectral     : Laplacian embedding, spectral edge costs and cells, effective resistance
  reductions   : exact graph reductions (degree 1 / degree 2 transit nodes)
  departs      : design of the départs of a cell (merge tree, capacitated Kruskal)
  partition    : variable-distance Voronoi with pluggable tree builder, init and growth
  local_search : edge insertion (H1-filtered) and boundary relocation moves
  features     : instance descriptors and k-NN configuration selection
  decomposition: subproblems from the ToMATo clusters merged until self-sufficient

TDA.py (next to main.py) runs the experiments.
"""
