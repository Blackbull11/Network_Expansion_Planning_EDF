"""
Split tree (merge tree of the superlevel sets) of f = distance to the source
inside a cell.

Sweeping the nodes of the cell by decreasing f, a component is born at each
local maximum of f (the end of a corridor of the network) and two components
merge at a saddle (where two corridors join, on their way to the source). The
source (f = 0) is the root. Only the components that contain terminals create
tree nodes: the tree is the corridor structure of the loads, as seen from the
source.

Node of the tree: {'value': f at birth (leaf) or at the merge (internal node),
'children': [...], 'own': terminals attached to this branch below its birth}.
Persistence of a branch = value(node) - value(parent): the length of the
corridor before it joins another one. simplify(tau) merges the branches of
persistence < tau into their parent.
"""
import math

from .union_find import UnionFind


class MergeTree:

    def __init__(self):
        self.nodes = {}
        self.root = None

    def new(self, value, children=(), own=()):
        k = len(self.nodes)
        self.nodes[k] = {'value': value, 'children': list(children), 'own': list(own)}
        return k

    def parent_map(self) -> dict:
        return {c: k for k, n in self.nodes.items() for c in n['children']}

    def terminals(self, k) -> list:
        """Terminals of the subtree of node k."""
        out, stack = [], [k]
        while stack:
            n = self.nodes[stack.pop()]
            out += n['own']
            stack += n['children']
        return out

    def postorder(self) -> list:
        out, stack = [], [(self.root, False)]
        while stack:
            k, done = stack.pop()
            if k is None:
                continue
            if done:
                out.append(k)
            else:
                stack.append((k, True))
                stack += [(c, False) for c in self.nodes[k]['children']]
        return out

    def branch_persistence(self) -> dict:
        """{node: value(node) - value(parent)} (inf for the root)."""
        par = self.parent_map()
        return {k: (n['value'] - self.nodes[par[k]]['value']) if k in par else math.inf
                for k, n in self.nodes.items()}

    def simplify(self, tau: float) -> "MergeTree":
        """Merge, leaves first, the branches of persistence < tau into their parent (in place)."""
        if tau <= 0:
            return self
        changed = True
        while changed:
            changed = False
            par = self.parent_map()
            pers = self.branch_persistence()
            for k in self.postorder():
                if k in par and pers[k] < tau and not self.nodes[k]['children']:
                    p = par[k]
                    self.nodes[p]['own'] += self.nodes[k]['own']
                    self.nodes[p]['children'].remove(k)
                    del self.nodes[k]
                    changed = True
                    break
        # contract internal nodes left with a single child
        for k in list(self.postorder()):
            n = self.nodes.get(k)
            if n is not None and len(n['children']) == 1 and k != self.root:
                c = n['children'][0]
                self.nodes[c]['own'] += n['own']
                p = self.parent_map()[k]
                self.nodes[p]['children'] = [c if x == k else x for x in self.nodes[p]['children']]
                del self.nodes[k]
        return self


def split_tree(H, root, f: dict, terminals) -> MergeTree:
    """
    Split tree of f on the graph H (nodes with finite f only), with `root`
    (f = 0) at the root. `terminals`: nodes that create / feed branches.
    """
    terminals = set(terminals)
    nodes = [v for v in H if math.isfinite(f.get(v, math.inf))]
    order = sorted(nodes, key=lambda v: (-f[v], str(v)))
    mt = MergeTree()
    uf = UnionFind()
    tnode = {}    # {component root: merge tree node or None}
    for v in order:
        uf.add(v)
        roots = {uf.find(u) for u in H[v] if u in uf} - {v}
        tree_nodes = [tnode[r] for r in roots if tnode.get(r) is not None]
        for r in roots:
            uf.union(v, r)
        if len(tree_nodes) >= 2:
            current = mt.new(f[v], children=tree_nodes)
        else:
            current = tree_nodes[0] if tree_nodes else None
        if v in terminals:
            if current is None:
                current = mt.new(f[v], own=[v])
            else:
                mt.nodes[current]['own'].append(v)
        tnode[uf.find(v)] = current
    comp = tnode.get(uf.find(root)) if root in uf else None
    if comp is None:
        comp = mt.new(0.0)
    if mt.nodes[comp]['value'] != f.get(root, 0.0) or mt.nodes[comp]['own']:
        comp = mt.new(f.get(root, 0.0), children=[comp])
    mt.root = comp
    return mt
