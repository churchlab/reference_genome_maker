"""
Contains algorithmic helper functions.
"""

from collections import defaultdict
from util import UnionFind

# Computes the length of the largest overlap between the end of s1 and
#   the beginning of s2
#
def overlapLen(s1, s2):
    for i in range(len(s1) + 1):
        if s2.startswith(s1[i:]):
            return len(s1) - i


# Computes the approximate longest Hamiltonian path of the given directed
#   graph, which is a |V| x |V| matrix with entry ij equal to the length
#   of edge ij. Uses the greedy algorithm to join the longest edges first.
#
# Returns two dictionaries, ins and outs
#   ins is a map from vertex index to its unique in-vertex, and
#   outs is a map from vertex index to its unique out-vertex.
#
def longest_hamiltonian_path(graph):
    # First, sort the edges from largest to smallest.
    # Bias the ties by preserving order (using i > j).
    sorted_edges = [(i, j, edge)
            for i, row in enumerate(graph)
            for j, edge in enumerate(row)]
    sorted_edges.sort(key=lambda (i, j, edge): (-edge, i > j))

    # Go through the edges, combining any pair of vertices as long as
    #   1) neither vertex has been directly connected to, and
    #   2) the two vertices are not connected by any path.
    # For (1) we maintain two flag tables for "in vertices" and "out vertices",
    #   and for (2) we maintain a Union Find data structure.
    ins = {}
    outs = {}
    uf = UnionFind()
    for i, j, edge in sorted_edges:
        if i == j or i in outs or j in ins or uf.find(i) == uf.find(j):
            continue
        outs[i] = j
        ins[j] = i
        uf.union(i, j)

    return ins, outs


# Finds the longest contiguous substring between the two strings.
#
# Returns a tuple (I, J, L) where:
#   I is the index of the substring in s1
#   J is the index of the substring in s2
#   L is the length of the substring
def longest_common_substring(s1, s2):
    table = [[0 for i in range(len(s2) + 1)] for j in range(len(s1) + 1)]
    best = (0, 0, 0)
    for i, c1 in enumerate(s1):
        for j, c2 in enumerate(s2):
            if c1 == c2:
                newLen = 1 + table[i][j]
                table[i + 1][j + 1] = newLen
                if newLen > best[2]:
                    best = (i - newLen + 1, j - newLen + 1, newLen)
    return best


# Finds a short superstring that contains all specified strings.
#   strings: a list of strings
#   returns the shortest common superstring (SCS)
#
def shortest_common_superstring(strings):
    # Use a greedy algorithm based on the following paper by Turner:
    # http://www.math.wsu.edu/math/faculty/bkrishna/FilesMath574/Papers/
    #   ApproxAlgosSCS.pdf
    #
    # The running time of the following algorithm is O(max(m, n^2 lg n)), where
    #   n is the number of strings and m is the total number of characters.
    # The implementation is simpler than that given in the paper, but is
    #   slower than the O(m lg n) approach detailed in the paper. This
    #   algorithm is guaranteed to produce a 3-approximation.

    # Trivial case.
    if len(strings) == 0:
        return [], ''

    # First remove all strings which are substrings of other strings.
    # This can be theoretically done in O(m) with a suffix tree, but
    #   is done naively now in order to not have dependencies.
    strings = filter(lambda s: all(
        [s == string or s not in string for string in strings]), strings)

    # Find longest overlap between every pair of strings
    overlaps = [[overlapLen(s1, s2) for s2 in strings] for s1 in strings]

    # Interpret the overlaps matrix as a graph, and find the longest
    #   Hamiltonian Path.
    ins, outs = longest_hamiltonian_path(overlaps)

    # Find the starting string, the one with no in-vertex.
    # Note: index and superstring gets stored with this initial string.
    for index, superstring in enumerate(strings):
        if index not in ins:
            break

    # Now compute the longest superstring and the locations of its substrings.
    while index in outs:
        next_index = outs[index]
        overlap = overlaps[index][next_index]
        superstring += strings[next_index][overlap:]
        index = next_index

    return superstring


# Find all overlaps between intervals in A with intervals in B.
# Returns a dictionary, with keys equal to indices in A and
#   values equal to lists consisting of indices of intervals in B that
#   overlap with the interval of A.
# An optional limit parameter gives an upper bound on the list sizes.
def find_overlapping_segments(A, B, limit=None):
    # Sort all interval endpoints.
    # Store tuples of the form (position, isB, isEnd, index)
    # Sort biased first on position, then isEnd
    endpoints = []
    for isB, group in enumerate((A, B)):
        for index, interval in enumerate(group):
            for isEnd, pos in enumerate(interval.endpoints()):
                endpoints.append((pos, isB, isEnd, index))
    endpoints.sort(key=lambda (pos, isB, isEnd, index): (pos, isEnd))

    # Line sweep.
    ans = defaultdict(list)
    currents = (set(), set())  # stores current intervals of A and B
    for pos, isB, isEnd, index in endpoints:
        if isEnd:
            if index in currents[isB]:
                currents[isB].remove(index)
        else:
            currents[isB].add(index)
            if isB:
                # Add this B interval to all A intervals that are not full.
                for aIndex in list(currents[0]):
                    ans[aIndex].append(index)
                    if limit and len(ans[aIndex]) >= limit:
                        currents[0].remove(aIndex)
            else:
                # Add at most limit intervals of B to this A interval.
                for counter, bIndex in enumerate(currents[1]):
                    if limit and counter >= limit:
                        currents[0].remove(index)
                        break
                    ans[index].append(bIndex)

    return ans

