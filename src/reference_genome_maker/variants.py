"""
The Variant class encapsulates a ref sequence and alt sequence.
"""

from util import Interval
from algorithm import longest_common_substring
from Bio.Seq import Seq

# Constants for calculating primitive variants.
MAX_CALC_SIZE = 100 * 100
MIN_INVERSION_SIZE = 10

# Helper function to find the reverse complement for a string.
def inversion(s):
    return str(Seq(s).reverse_complement())


###############################################################################
# Variants
###############################################################################

class Variant:

    def __init__(self, position, ref, alt, note=None):
        self.position = position
        self.ref = ref
        self.alt = alt
        self.note = note
        self.type = None

    def __eq__(self, other):
        return self.position == other.position and \
                self.ref == other.ref and \
                self.alt == other.alt

    def __cmp__(self, other):
        return self.position.__cmp__(other.position) or \
                len(self.ref).__cmp__(len(other.ref))

    def __str__(self):
        return '{POS: %d, REF: %s, ALT: %s}' % \
                (self.position, self.ref, self.alt)

    def __repr__(self):
        return self.__str__()

    def start(self):
        return self.position

    def end(self):
        return self.position + len(self.ref)

    def interval(self):
        return Interval(self.start(), self.end())

    def primitive_variants(self):
        """Returns a list of nonoverlapping, "primitive" variants that
        are equivalent to this one. Primitive variants include basic
        insertions, deletions, duplications, inversions."""

        # LCS algorithm is O(NM), so abort if input is too large.
        if len(self.ref) * len(self.alt) > MAX_CALC_SIZE:
            return [self]

        variants = self.primitive_variants_helper()

        # Hack: for each insertion, check if it is actually a duplication
        for variant in variants:
            if variant.type == 'INS':
                pos = variant.position - self.position
                if self.ref[pos-len(variant.alt):pos] == variant.alt or \
                        self.ref[pos:pos+len(variant.alt)] == variant.alt:
                            variant.type = 'DUP:TANDEM'

        return variants

    def primitive_variants_helper(self):
        # Find LCS that matches between ref and alt.
        I, J, L = longest_common_substring(self.ref, self.alt)

        # Also try LCS of reverse complement, to check for large inversion.
        # If a larger (and significant) inversion is found, use it instead.
        II, IJ, IL = longest_common_substring(self.ref, inversion(self.alt))
        use_inversion = IL > MIN_INVERSION_SIZE and IL > L
        if use_inversion:
            I, J, L = II, len(self.alt) - (IJ + IL), IL

        # No common strings; add a type if possible, and return this variant.
        if L == 0:
            if self.ref or self.alt:
                if not self.alt:
                    self.type = 'DEL'
                elif not self.ref:
                    self.type = 'INS'
                return [self]
            return []

        variants = []
        variants.extend(Variant(self.position, self.ref[:I], self.alt[:J])
                .primitive_variants())
        if use_inversion:
            inv_ref = self.ref[I:I+L]
            inv = Variant(self.position + I, inv_ref, inversion(inv_ref))
            inv.type = 'INV'
            variants.append(inv)
        variants.extend(Variant(self.position + I + L, self.ref[I+L:],
                self.alt[J+L:]).primitive_variants())

        return variants

    def get_type(self):
        variants = self.primitive_variants()
        if len(variants) == 0:
            return 'NONE'
        elif len(variants) > 1:
            return 'COMPLEX'
        else:
            return variants[0].type

