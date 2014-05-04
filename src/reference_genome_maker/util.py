"""
Contains utility functions for algorithmic problems.
"""

###############################################################################
# Rolling Hash
###############################################################################

PRIME, MOD = 97654321, 2 ** 32  # constants for rolling hash

# Computes the forward rolling hash of a string.
# Implemented as a generator that outputs (index, hash) tuples
#
def forward_hash(s):
    roll = 0
    for i, c in enumerate(s):
        roll = (roll * PRIME + ord(c)) % MOD
        yield i + 1, roll

# Computes the backward rolling hash of a string.
# Implemented as a generator that outputs (index, hash) tuples
#
def backward_hash(s):
    roll, exp = 0, 1
    for i, c in enumerate(s[::-1]):
        roll = (roll + exp * ord(c)) % MOD
        yield i + 1, roll
        exp = (exp * PRIME) % MOD


###############################################################################
# Union Find
###############################################################################

class UnionFind:
    
    # This Union Find data structure assumes that any key not in the
    #   self.parents dictionary is in its own set, so dictionary accesses
    #   to self.parents are always of the form self.parents.get(x, x).
    #
    def __init__(self):
        self.parents = {}

    # Find the representative element of node's set.
    # We use the path compression heuristic.
    def find(self, node):
        parent = node
        while self.parents.get(parent, parent) != parent:
            parent = self.parents[parent]
        while self.parents.get(node, node) != node:
            self.parents[node] = parent
            node = self.parents[node]
        return parent

    # Union the set of node1 and the set of node2.
    def union(self, node1, node2):
        parent1 = self.find(node1)
        parent2 = self.find(node2)
        self.parents[parent1] = parent2


###############################################################################
# Intervals
###############################################################################

class Interval:

    # An interval is from start inclusive, to end exclusive.
    def __init__(self, start, end):
        self.start = start
        self.end = end

    @staticmethod
    def create(start, length):
        return Interval(start, start + length)

    def __str__(self):
        return '(%d, %d)' % (self.start, self.end)

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.end - self.start

    def endpoints(self):
        return self.start, self.end

    def contains(self, pos):
        return pos >= self.start and pos < self.end

    def sliceInto(self, s):
        return s[self.start : self.end]


###############################################################################
# Interval Mapper
###############################################################################

class IntervalMapper:

    # This data structure stores a mapping from intervals to intervals
    #   (we say from "ref intervals" to "alt intervals")
    #   that is similar to Runtime Liftover, but provides optimizations
    #   for in-order accesses.
    #
    def __init__(self):
        self.mapper = []
        self.generator = None
        self.flags = (False, 0)  # internal state for asserting invariants

    # Appends a new mapping from the ref interval to the alt interval.
    # IMPORTANT: ref intervals must be added in order of interval start pos.
    #
    def add_interval(self, ref, alt):
        assert (False, ref.start) >= self.flags, \
                'add_interval() called out of order'
        self.flags = False, ref.start

        self.mapper.append((ref, alt))

    # Returns where the given position will map to.
    # IMPORTANT: calls to get_mapping must be in order of position,
    #   and after all calls to add_interval have finished.
    #
    def get_mapping(self, pos):
        assert (True, pos) >= self.flags, 'get_mapping() called out of order'
        self.flags = True, pos

        if not self.generator:
            self.generator = enumerate(self.mapper)

        while True:
            index, mapping = self.generator.next()
            if mapping[0].contains(pos):
                # If position is the nth character of the interval,
                #   then return the nth character of the mapped interval.
                return min(mapping[1].end,
                        pos - mapping[0].start + mapping[1].start)

