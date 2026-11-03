"""Union-find (disjoint sets) on hashable items, with path halving and union by size."""


class UnionFind:

    def __init__(self, items=()):
        self.parent = {}
        self.size = {}
        for x in items:
            self.add(x)

    def add(self, x) -> None:
        if x not in self.parent:
            self.parent[x] = x
            self.size[x] = 1

    def __contains__(self, x) -> bool:
        return x in self.parent

    def find(self, x):
        parent = self.parent
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(self, a, b):
        """Merge the sets of a and b, return the new root (None if already merged)."""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return None
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        return ra

    def groups(self) -> dict:
        """{root: list of items}."""
        out = {}
        for x in self.parent:
            out.setdefault(self.find(x), []).append(x)
        return out
