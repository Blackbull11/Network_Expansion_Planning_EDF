"""
Exact reductions of the graph (they keep an optimal solution of the problem
with all the constraints of verify.py):
  - a transit node of degree 1 is never needed (dropping it from a tree keeps
    it feasible and shorter): removed;
  - a transit node v of degree 2, between a and b, is used by a tree only
    through both its edges: the path a - v - b is replaced by one edge of
    length l(a, v) + l(v, b). If a - b already exists, the shorter of the two
    is kept (a tree never uses both, they would close a cycle).
Repeated until no node can be reduced. Sources and loads are never reduced.

Each edge of the reduced graph keeps in 'path' the list of the original
nodes it stands for, so that solutions can be mapped back (expand_edges).
"""
import networkx as nx

TERMINAL_LABELS = ('Source', 'Conso', 'Prod')


def reduce_graph(G, weight: str = "weight") -> tuple:
    """(reduced graph, stats). Node attributes are copied, edges get 'path'."""
    R = nx.Graph()
    R.add_nodes_from(G.nodes(data=True))
    for u, v, w in G.edges(data=weight):
        R.add_edge(u, v, **{weight: w, 'path': (u, v)})
    terminal = {v for v, lab in G.nodes('label') if lab in TERMINAL_LABELS}
    stack = [v for v in R if v not in terminal]
    n_deg1 = n_deg2 = n_parallel = 0
    while stack:
        v = stack.pop()
        if v not in R or v in terminal:
            continue
        nbrs = list(R[v])
        if len(nbrs) <= 1:
            R.remove_node(v)
            n_deg1 += 1
            stack += nbrs
        elif len(nbrs) == 2:
            a, b = nbrs
            pa, pb = R[a][v]['path'], R[v][b]['path']
            pa = pa if pa[0] == a else pa[::-1]
            pb = pb if pb[0] == v else pb[::-1]
            w = R[a][v][weight] + R[v][b][weight]
            R.remove_node(v)
            n_deg2 += 1
            if R.has_edge(a, b):
                n_parallel += 1
                if w < R[a][b][weight]:
                    R[a][b].update({weight: w, 'path': tuple(pa) + tuple(pb[1:])})
            else:
                R.add_edge(a, b, **{weight: w, 'path': tuple(pa) + tuple(pb[1:])})
            stack += [a, b]
    stats = {'nodes_before': G.number_of_nodes(), 'edges_before': G.number_of_edges(),
             'nodes_after': R.number_of_nodes(), 'edges_after': R.number_of_edges(),
             'removed_deg1': n_deg1, 'contracted_deg2': n_deg2, 'parallel_dropped': n_parallel}
    return R, stats


def expand_edges(R, edges) -> list:
    """Original edges (u, v), u < v, of a list of edges of the reduced graph R."""
    out = []
    for a, b in edges:
        p = R[a][b]['path']
        out += [(min(x, y), max(x, y)) for x, y in zip(p, p[1:])]
    return sorted(set(out))
