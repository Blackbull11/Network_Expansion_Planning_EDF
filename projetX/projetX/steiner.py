"""
Steiner trees inside the Voronoi cells: in each cell, connect the source and
all the nodes that have a demand or a production (Conso / Prod), using
transit nodes only where they help shorten the tree.
"""
import networkx as nx
from networkx.algorithms.approximation import steiner_tree

from nx_graph import GraphData
from src.outils import init_stream_logger

logger = init_stream_logger()

TERMINAL_LABELS = ('Source', 'Conso', 'Prod')


def steiner_trees_by_cells(G: nx.Graph, cells: dict, weight: str = "weight", log: bool = True) -> dict:
    """
    Approximate Steiner tree (Mehlhorn, 2-approximation) of each cell.

    Parameters
    ----------
    G : nx.Graph
        graph with a 'label' attribute on nodes and `weight` on edges.
    cells : dict
        {source: set of nodes}. Each cell must be connected (true for Voronoi
        cells, since shortest paths stay in their cell). 'unreachable' is ignored.
    log : bool
        log the length of each tree.

    Returns
    -------
    trees : dict
        {source: list of edges (u, v) with u < v}
    """
    trees = {}
    for k, nodes in cells.items():
        if k == 'unreachable':
            continue
        G_k = G.subgraph(nodes)
        terminals = [n for n in nodes if G.nodes[n]['label'] in TERMINAL_LABELS]
        tree = steiner_tree(G_k, terminals, weight=weight, method='mehlhorn')
        trees[k] = sorted((min(u, v), max(u, v)) for u, v in tree.edges)
        if not log:
            continue
        length = sum(G[u][v][weight] for u, v in trees[k])
        logger.info("Cellule %s : arbre de Steiner de %d arêtes, longueur %.1f (%d terminaux)",
                    k, len(trees[k]), length, len(terminals))
    return trees


def compute_steiner_trees(graph: GraphData) -> None:
    """
    Build the Steiner tree of each Voronoi cell of `graph` (compute_voronoi_cells
    must have been called) and store them in graph.steiner_edges, in the
    format {source: list of edges (u, v), u < v}.
    """
    graph.steiner_edges = steiner_trees_by_cells(graph.G, graph.cells)


def show_steiner_trees(graph: GraphData, title: str = "Steiner trees by Voronoi cell", **kwargs):
    """
    Plot the cells, with the Steiner tree edges solid and other edges dotted.
    `kwargs` (e.g. filepath, show) are passed to draw_graph_from_dataframes.
    """
    # the drawing modifies the DataFrames in place, so give it copies
    return graph.draw_graph_from_dataframes(
        nodes=graph.nodes_split.copy(), edges=graph.edges_split.copy(),
        used_edges_by_cells=graph.steiner_edges, title=title, **kwargs)
