import networkx as nx
import pandas as pd
import json
from pathlib import Path
import os
import geopandas as gpd

import matplotlib
import matplotlib.pyplot as plt

from shapely.geometry import Point, LineString
from src.outils import init_stream_logger

from collections import namedtuple
logger = init_stream_logger()

BASE_COLORS = ['b', 'hotpink', 'orange', 'mediumseagreen', 'purple', 'saddlebrown',
               'c', 'olive', 'gold', 'slategrey']


def cell_colors(cell_ids) -> dict:
    """
    {cell_id: color}, one distinct color per cell ('unreachable' / 'split' excluded).
    Cells are sorted, then take BASE_COLORS in turn; if there are more cells than
    BASE_COLORS, colors are taken evenly in the 'hsv' colormap instead.
    """
    ids = sorted({k for k in cell_ids if k not in ('unreachable', 'split')}, key=str)
    if len(ids) <= len(BASE_COLORS):
        return dict(zip(ids, BASE_COLORS))
    cmap = plt.get_cmap('hsv', len(ids) + 1)
    return {k: cmap(i) for i, k in enumerate(ids)}


class GraphData:

    def __init__(self, folderpath:str="data", filename:str="graphe_maille.json"):
        self.filepath = os.path.join(folderpath, filename)
        self.nodes = None                 # GeoDataFrame of nodes
        self.edges = None                 # GeoDataFrame of edges
        self.G = None                     # networkx undirected graph
        # Optional results, filled by other modules and used only for plotting
        self.cells: dict = None                  # {cell_id: nodes}
        self.edge_cells: dict = None             # {cell_id: edges (u, v)}
        self.nodes_close_to_source: dict = None  # {cell_id: node ids}
        self.used_edges_by_cells: dict = None    # {cell_id: edges on shortest paths}
        self.principal_edges: list = None        # main edges (u, v), u < v
        # Filled by voronoi.compute_voronoi_cells
        self.dist: dict = None                   # {node: scaled distance to its source}
        self.paths: dict = None                  # {node: shortest path from its source}
        self.G_split = None                      # copy of G with 'dummy' split nodes
        self.nodes_split = None                  # nodes GeoDataFrame incl. dummy nodes
        self.edges_split = None                  # edges GeoDataFrame incl. half edges
        # Filled by steiner.compute_steiner_trees
        self.steiner_edges: dict = None          # {cell_id: Steiner tree edges (u, v), u < v}
        # Filled by voronoi_variable_distance.compute_voronoi_variable
        self.voronoi_weights: dict = None        # {source: weight w_s of the distance}

        self.lire_json_graphe(folderpath, filename)
        self.from_geodf_to_undirected_networkx(self.nodes, self.edges)

    def __str__(self):
        # Summary: node, edge and source counts
        return f"GraphData with {len(self.nodes)} nodes, {len(self.edges)} edges and {len(self.nodes[self.nodes['type'] == 'Source'])} sources from {self.filepath}"

    def show(self):
        # Plot with the instance's current cells / edges data
        self.draw_graph_from_dataframes()

    def lire_json_graphe(self, folderpath:str, filename:str="graphe_maille.json"):
        """
        Read the graph stored in `folderpath`/`filename`.

        Parameters
        ----------
        folderpath : str
            path of the folder containing the graph file
        filename : str, optional
            name of the graph file (same format as graphe_maille.json).
            The default is "graphe_maille.json".

        Raises
        ------
        FileNotFoundError
            If `filename` is not found in `folderpath`.

        Returns
        -------
        nodes : geopandas.GeoDataFrame
            the nodes of the graph, one row per node, with a shapely Point
            "geometry" column built from the node coordinates
        edges : geopandas.GeoDataFrame
            the edges of the graph (undirected), one row per edge, with a
            shapely LineString "geometry" column

        """
        # Lecture du graphe depuis un fichier JSON au format graphe_maille.json
        logger.info("Lecture du graphe : %s/%s", folderpath, filename)
        graphe_filename = os.path.join(folderpath, filename)
        if Path(graphe_filename).exists():
            with open(graphe_filename, "r") as f:
                graphe_data = json.load(f)
            # Conversion en DataFrame
            nodes = pd.DataFrame(graphe_data["nodes"])
            edges = pd.DataFrame(graphe_data["edges"])
            #arcs = pd.DataFrame(graphe_data["arcs"]) # arcs are directed, whereas edges are not
            # Conversion explicite en GeoDataFrame pour l'affichage
            nodes["geometry"] = nodes["geometry"].apply(lambda coords: Point(coords))
            nodes = gpd.GeoDataFrame(nodes, geometry="geometry")
            edges["geometry"] = edges["geometry"].apply(lambda coords: LineString(coords))
            edges = gpd.GeoDataFrame(edges, geometry="geometry")
            logger.info("Graphe maille chargé et converti en GeoDataFrame.")
        else:
            raise FileNotFoundError(f"Fichier <{filename}> introuvable dans le dossier {folderpath}. \n veuillez run 'generate_instances.py' ")
        self.nodes = nodes
        self.edges = edges

    def from_geodf_to_undirected_networkx(self, nodes: gpd.GeoDataFrame, edges: gpd.GeoDataFrame) -> nx.Graph :
        """
        Convert nodes and edges GeoDataFrames to a networkx graph.

        Exits the program (through check_reachable_nodes) if a node has no
        incident edge.

        Parameters
        ----------
        nodes : geopandas.GeoDataFrame
            the nodes of the graph, with columns "geometry", "type", "demande"
            and "demande_hivers"
        edges : geopandas.GeoDataFrame
            the edges of the graph (undirected), with columns "start_node_idx",
            "end_node_idx" and "length_km"

        Returns
        -------
        G : nx.Graph
            networkx undirected Graph. Nodes are the index of `nodes`, with
            attributes:
                - pos : (x, y) coordinates
                - label : node type ('Source', 'Transit', 'Conso', 'Prod')
                - conso_hiver : - demande_hivers
                - conso_ete : demande (used as the source capacity in
                the Voronoi computations)
            Edges have a single attribute "weight" = length_km.

        """
        G = nx.Graph()
        # One graph node per row, keyed by the DataFrame index
        for idx, row in nodes.iterrows():
            point = row['geometry']
            G.add_node(idx, pos=(point.x, point.y),
                    label=row['type'],
                    conso_hiver = - row['demande_hivers'], 
                    conso_ete = row['demande']
                    )
        # Edge weight is the line length in km
        for line in edges.itertuples():
            start_idx = line.start_node_idx
            end_idx = line.end_node_idx
            G.add_edge(start_idx, end_idx, 
                        weight=line.length_km)
        self.check_reachable_nodes(G)
        self.G = G
        logger.info("Objet GraphData créé.")

    def check_reachable_nodes(self, G: nx.Graph):
        """
        Check that every node of `G` has at least one incident edge.

        Prints the list of isolated nodes (empty if there are none), then logs
        an error and exits the program with code 1 if the list is not empty.
        Note that this does not check that the graph is connected.

        Parameters
        ----------
        G : nx.Graph
            graph to check

        Returns
        -------
        None.

        """
        unreachable = []
        # A node with an empty adjacency dict has no incident edge
        for n, nbrsdict in G.adjacency():
            if len(nbrsdict)==0:
                unreachable += [n]
        print(unreachable)
        if len(unreachable) > 0:
            logger.error(f'Some nodes have no connected arcs: {unreachable}')
            exit(1)
        # print(len(G))
        # G.remove_nodes_from(unreachable)
        # print(len(G))
        # for n in unreachable:
        #     nodes.drop(n, inplace=True)
        # print(len(nodes))

    def draw_graph_from_dataframes(self, nodes: gpd.GeoDataFrame=None,
                                edges: gpd.GeoDataFrame=None,
                                cells: dict=None,
                                edge_cells: dict=None,
                                nodes_close_to_source: dict=None,
                                used_edges_by_cells: dict=None,
                                principal_edges: list=None,
                                title: str=None,
                                filepath: str=None,
                                show: bool=True,
                                ) -> matplotlib.axes._axes.Axes:
        """
        Plot graph from GeoDataFrames nodes and edges, and draw nodes/edges
        in the color associated with their cell, then save and/or show the figure.

        Each cell gets its own color (see cell_colors), whatever the number of
        cells. The special keys 'unreachable' / 'split' are drawn in red / not drawn.

        Parameters
        ----------
        nodes : gpd.GeoDataFrame
            nodes of the graph (and data associated with nodes)
        edges : gpd.GeoDataFrame
            edges of the graph (and data associated with edges)
        cells : dict, optional
            Nodes cells in the form {cell_id: list of nodes in the cell}.
            If None, nothing is drawn for the nodes. The default is None.
        edge_cells : dict, optional
            Edges cells in the form {cell_id: list of edges (u, v) in the cell},
            as returned by voronoi_cells_edges.
            Note that initial edges may be split in two if their starting and ending points
            are in different cells. In this case, a point is inserted in the initial segment,
            and two edges are created and associated with their respective cells.
            The initial edges are listed under the key 'split' and are not drawn.
            `nodes` and `edges` must then contain the added 'dummy' nodes and edges.
            If None, only the edges with no existing line are drawn (in grey),
            and only if `cells` is given. The default is None.
        nodes_close_to_source : dict, optional
            {cell_id: collection of node ids ("id_noeud" column)}. These nodes are
            marked with a black cross. Only used if `cells` is given.
            The default is None.
        used_edges_by_cells : dict, optional
            {cell_id: list of edges (u, v) with u < v}, as returned by
            edges_on_shortest_paths. Used edges are drawn as solid lines, the other
            edges of the cell as dotted lines. Only used if `edge_cells` is given.
            The default is None.
        principal_edges : list, optional
            list of edges (u, v) with u < v, drawn as black dotted lines.
            The default is None.
        title : str, optional
            title of the figure. The default is None.
        filepath : str, optional
            if given, the figure is saved to this file. The default is None.
        show : bool, optional
            call plt.show() if True, else close the figure. The default is True.

        Returns
        -------
        ax : matplotlib.axes._axes.Axes
            the axes the graph was drawn on

        Notes
        -----
        The input DataFrames are modified in place: a "cell" column is added to
        `nodes` when `cells` is given, and `edges` is re-indexed by
        (node1, node2) = (min, max) of its end nodes when `edge_cells` or
        `principal_edges` is given (plus a "cell" column when `edge_cells` is given).

        """
        # Defaults come from the instance (can't be written in the signature)
        nodes = self.nodes if nodes is None else nodes
        edges = self.edges if edges is None else edges
        cells = self.cells if cells is None else cells
        edge_cells = self.edge_cells if edge_cells is None else edge_cells
        nodes_close_to_source = self.nodes_close_to_source if nodes_close_to_source is None else nodes_close_to_source
        used_edges_by_cells = self.used_edges_by_cells if used_edges_by_cells is None else used_edges_by_cells
        principal_edges = self.principal_edges if principal_edges is None else principal_edges

        f, ax = plt.subplots()

        if not (cells or edge_cells or principal_edges):
            # nothing to colour: draw the raw graph
            edges.plot(ax=ax, color="grey", linewidth=0.5)
            nodes.plot(ax=ax, color="b", markersize=5)

        colors = cell_colors(list(cells or {}) + list(edge_cells or {}))
        
        # --- Nodes: coloured by cell, marker depends on node type ---
        if cells:
            inv_cells = { n : k for k in cells for n in cells[k]}  # node -> cell
            nodes['cell'] = pd.Series(inv_cells)
            cell_color = { k : colors[k] if k != 'unreachable' else 'r' for k in cells }
            #print(cell_color)
            if not edge_cells:
                edges[edges['nb_lignes_existantes'] == 0].plot(ax=ax, color="grey", linewidth=0.2)
                #edges[edges['nb_lignes_existantes'] != 0].plot(ax=ax, color="red", linewidth=0.5)
            Plot_data = namedtuple('Plot_data', ['marker', 'markersize'])
            node_types = {'Source': Plot_data('s', 30),
                        'Transit': Plot_data('.', 5),
                        'Conso': Plot_data('o', 12),
                        'Prod': Plot_data('s', 12)}
            for k in cells:
                nodes_k = nodes[nodes['cell'] == k]
                color_k = cell_color[k]            
                for nt, nt_plot in node_types.items():
                    nodes_nt = nodes_k[nodes_k['type'] == nt]
                    if len(nodes_nt):
                        nodes_nt.plot(ax=ax, color=color_k, marker=nt_plot.marker, markersize=nt_plot.markersize)
                if nodes_close_to_source:
                    criterion = nodes['id_noeud'].map(lambda i: i in nodes_close_to_source[k])
                    nodes[criterion].plot(ax=ax, color='k', marker='x', markersize=10)

        # --- Edges: coloured by cell; used edges solid, others dotted ---
        if edge_cells:
            inv_cells = { edge : k for k in edge_cells for edge in edge_cells[k]}  # edge -> cell
            # Index edges by (min, max) end nodes to match the (u, v) keys
            edges['node1'] = edges[['start_node_idx','end_node_idx']].min(axis=1)
            edges['node2'] = edges[['start_node_idx','end_node_idx']].max(axis=1)
            edges.set_index(['node1', 'node2'], inplace=True)
            edges['cell']=pd.Series(inv_cells) #.loc[edges['nodes']]
            if used_edges_by_cells:
                inv_used_edges = {edge : 1 for k in used_edges_by_cells for edge in used_edges_by_cells[k]}
                s = pd.Series(inv_used_edges)
                s.name='used'
                edges = pd.concat([edges, s], axis=1)
            cell_color = { k : colors[k] if k != 'split' else 'r' for k in edge_cells }
            for k in edge_cells:
                if k != 'split':
                    edges_k = edges[edges['cell'] == k]
                    if used_edges_by_cells:
                        edges_k[edges_k['used'] == 1].plot(ax=ax, color=cell_color[k], linewidth=1)
                        edges_k[edges_k['used'] != 1].plot(ax=ax, color=cell_color[k], linewidth=0.5, linestyle = "dotted")
                    else :
                        edges_k.plot(ax=ax, color=cell_color[k], linewidth=0.2)
            nodes[nodes['type'] == 'dummy'].plot(ax=ax, color='k', marker=".", markersize=1)
        
        # --- Principal edges: black dotted overlay ---
        if principal_edges:
            edges['node1'] = edges[['start_node_idx','end_node_idx']].min(axis=1)
            edges['node2'] = edges[['start_node_idx','end_node_idx']].max(axis=1)
            edges.set_index(['node1', 'node2'], inplace=True)
            edges.loc[principal_edges].plot(ax=ax, color='black', linewidth=1, linestyle = "dotted")
        if not (cells or edge_cells or principal_edges):
            edges.plot(ax=ax, color="grey", linewidth=0.5)
            nodes.plot(ax=ax, color="b", markersize=5)
        # Set margins for the axes so that nodes aren't clipped
        ax.plot()
        ax = plt.gca()
        ax.margins(0.20)
        plt.axis("off")
        if title:
            plt.title(title, fontsize=6)
        if filepath:
            plt.savefig(filepath, dpi=200, bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.close(f)
        
        return ax

    def plot_with_plotly(G):
        """
        Draw the graph with plotly, without the 'dummy' nodes and edges, and
        show it. Hovering over an edge shows its end nodes "u-v".

        Parameters
        ----------
        G : nx.Graph
            graph with a "pos" attribute on nodes

        Returns
        -------
        None.

        """
        import plotly.graph_objects as go
        # Real (non-dummy) edges and their endpoint coordinates
        graph_edges = [(u,v) for u,v,d in G.edges.data('label') if d != 'dummy']
        coord =[(G.nodes[u]['pos'], G.nodes[v]['pos']) for u,v in graph_edges]
        # None separates segments so plotly draws each edge on its own
        edge_x = []
        edge_y = []
        # nodes in the middle for hover
        mnode_x = []
        mnode_y = []
        edge_label = []
        for p,q in coord:
            x_0, y_0 = p
            x_1, y_1 = q
            edge_x.extend([x_0, x_1, None])
            edge_y.extend([y_0, y_1, None])
            # nodes for hover
            mnode_x.extend([(x_0 + x_1)/2])
            mnode_y.extend([(y_0 + y_1)/2])
            mnode_x.extend([(3*x_0 + x_1)/4])
            mnode_y.extend([(3*y_0 + y_1)/4])
            mnode_x.extend([(x_0 + 3*x_1)/4])
            mnode_y.extend([(y_0 + 3*y_1)/4])
        # 3 "middle" points for each edge
        for u,v in graph_edges:
            edge_label += [f'{u}-{v}']*3
        
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            mode='lines')
        
        mnode_trace = go.Scatter(x = mnode_x, y = mnode_y, mode = "markers", showlegend = False,
                            #hovertemplate = "Edge %{hovertext}<extra></extra>",
                            hoverinfo='text',
                            hovertext = edge_label, marker = go.Marker(opacity = 0))
        
        graph_nodes = [n for (n,d) in G.nodes.data('label') if d != 'dummy']
        node_x = [G.nodes[n]['pos'][0] for n in graph_nodes]
        node_y = [G.nodes[n]['pos'][1] for n in graph_nodes]
        
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers', #"markers+text
            hoverinfo='text',
            #hovertext=node_labels,
            marker=dict(
                # showscale=True,
                # colorscale options
                #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
                #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
                #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
                # colorscale='YlGnBu',
                # reversescale=True,
                color='orange',
                size=10, #node_sizes
                # colorbar=dict(
                #     thickness=15,
                #     title=dict(
                #       text='Node Connections',
                #       side='right'
                #     ),
                #     xanchor='left',
                # ),
                line_width=2), #line=dict(width=1, color="black")
            text=list(G.nodes()),
            textposition="bottom center",
            textfont=dict(size=12),
            )
        
        fig = go.Figure(layout=go.Layout(showlegend=False,
                                        hovermode='closest',
                                        xaxis=dict(showgrid=True, zeroline=False, showticklabels=True),
                                        yaxis=dict(showgrid=True, zeroline=False, showticklabels=True)))

        fig.add_trace(edge_trace)
        fig.add_trace(node_trace)
        fig.add_trace(mnode_trace)
        fig.show()




