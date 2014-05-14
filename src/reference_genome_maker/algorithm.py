"""
Contains algorithmic helper functions.
"""

from collections import defaultdict
from util import forward_hash
from util import backward_hash
from util import UnionFind

# Computes all strings that are substrings of another string
# Returns a list of tuples (I, J, K) where:
#   I is the index of the substring
#   J is the index of the superstring
#   K is the index of the substring in the superstring
def compute_substrings(strings):
    substrings = []
    for i, s in enumerate(strings):
        # Bias looking for nearby strings first, so increment from di = 1.
        # For strings to the right, use index, but
        #   for strings to the left, use rindex to find closest match.
        for di in range(1, len(strings)):
            j = i + di
            if j >= 0 and j < len(strings):
                other = strings[j]
                if s != other and s in other:
                    substrings.append((i, j, other.index(s)))
                    break
            j = i - di
            if j >= 0 and j < len(strings):
                other = strings[j]
                if s != other and s in other:
                    substrings.append((i, j, other.rindex(s)))
                    break
    return substrings


# Computes the longest overlaps between all pairs of strings
#
def compute_overlaps(strings):
    overlaps = [[0] * len(strings) for s in enumerate(strings)]
    hashes = defaultdict(list)
    for i, s in enumerate(strings):
        for index, hash in backward_hash(s):
            hashes[hash].append(i)
    for i, s in enumerate(strings):
        for index, hash in forward_hash(s):
            for j in hashes[hash]:
                if strings[j][-index:] == strings[i][:index]:
                    overlaps[j][i] = index
    return overlaps


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
#
# Returns a tuple (L, S) where:
#   S is the shortest common superstring (SCS)
#   L is a list of indices [A_1, A_2, ... A_n] of the strings in the SCS
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
    filtered_strings = compute_substrings(strings)
    for i, j, index in reversed(filtered_strings):
        strings.pop(i)

    # Find longest overlap between every pair of strings
    overlaps = compute_overlaps(strings)

    # Interpret the overlaps matrix as a graph, and find the longest
    #   Hamiltonian Path.
    ins, outs = longest_hamiltonian_path(overlaps)

    # Find the starting string, the one with no in-vertex.
    # Note: index and superstring gets stored with this initial string.
    for index, superstring in enumerate(strings):
        if index not in ins:
            break

    # Now compute the longest superstring and the locations of its substrings.
    locations = [0] * len(strings)
    while index in outs:
        next_index = outs[index]
        overlap = overlaps[index][next_index]
        locations[next_index] = len(superstring) - overlap
        superstring += strings[next_index][overlap:]
        index = next_index

    # Finally, compute the locations of substrings removed in the first step.
    # First we add them back with placeholders, then we set the correct value.
    for i, j, index in filtered_strings:
        locations.insert(i, -1)  # placeholder
    for i, j, index in filtered_strings:
        locations[i] = locations[j] + index

    return locations, superstring


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

