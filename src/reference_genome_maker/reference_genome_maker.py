"""
This module processes a reference genome (SeqRecord)
and a variant set (an object in this module),
and outputs a modified genome (SeqRecord).

Methods also exist for parsing either a csv or vcf
into a variant set.
"""

import csv
import vcf
from util import Interval
from util import IntervalMapper
from util import feature_interval
from algorithm import shortest_common_superstring
from algorithm import longest_common_substring

from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition

###############################################################################
# Variants
###############################################################################

class Variant:

    def __init__(self, position, ref, alt, note=None):
        self.position = position
        self.ref = ref
        self.alt = alt
        self.note = note

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
        I, J, L = longest_common_substring(self.ref, self.alt)
        if L <= 1:
            return [self] if self.ref or self.alt else []
        variants = Variant(self.position, self.ref[:I], self.alt[:J]) \
                .primitive_variants() + \
                Variant(self.position + I + L, self.ref[I+L:], self.alt[J+L:]) \
                .primitive_variants()
        return variants


###############################################################################
# Reference Genome Maker
###############################################################################

class RefGenomeMaker:

    def parse_from_csv(variant_set_csv):
        """Updates the Ref Genome Maker with a list of variants in a csv.

        Args:
            variant_data_csv: Path to .csv file containing the following cols:
                Required:
                    * position - Position in the starting genome.
                    * ref - Reference sequence at that position.
                    * alt - The alternative to replace with.
                Optional:
                    * note - A note to add to the data.
        """
        self.variants = []

        with open(variant_set_csv) as csv_fh:
            csv_reader = csv.DictReader(csv_fh)
            for row in csv_reader:
                self.variants.append(Variant(
                    int(row['position']) - 1,
                    row['ref'],
                    row['alt'],
                    row.get('note', None)
                ))

    def apply_variants(self):
        """Applies stored variants to the given ref genome (string).

        Returns a tuple (A, M) where:
            A is the modified genome (string)
            M is a mapper object (for determining mappings from ref to mod)
        """

        # First, sort the variants by position
        genome_seq = str(self.genome.seq)
        self.variants.sort()

        # Find groups of variants whose ref sequences overlap.
        # Find the shortest common superstring of the alts of these variants.
        # Replace the refs with the shortest common superstring.
        self.primitive_variants = []
        generator = iter(self.variants)
        variant = next(generator, None)
        while variant:
            start, end = variant.interval().endpoints()
            variant_alts = []
            assert variant.ref == variant.interval().sliceInto(genome_seq)

            while True:
                variant_alts.append(variant.alt)
                variant = next(generator, None)
                if not variant or variant.start() >= end:
                    locations, superstring = shortest_common_superstring(variant_alts)
                    self.primitive_variants.extend(Variant(
                        start, genome_seq[start:end], superstring
                        ).primitive_variants())
                    break
                end = max(end, variant.end())

        # Add sentinel.
        self.primitive_variants.append(Variant(len(self.genome), '', ''))

        # Now apply the primitive variants in order
        self.mapper = IntervalMapper()
        alt_genome = ''
        ref_index, alt_index = 0, 0  # index of where we are in ref and alt
        for variant in self.primitive_variants:
            # First append unchanged part of ref genome
            ref_interval = Interval(ref_index, variant.start())
            alt_interval = Interval.create(alt_index, len(ref_interval))
            alt_genome += ref_interval.sliceInto(genome_seq)
            self.mapper.add_interval(ref_interval, alt_interval)
            alt_index += len(ref_interval)

            alt_interval = Interval.create(alt_index, len(variant.alt))
            alt_genome += variant.alt
            self.mapper.add_interval(variant.interval(), alt_interval)
            ref_index = variant.end()
            alt_index += len(variant.alt)

        # Update genome
        self.genome.seq = Seq(alt_genome, self.genome.seq.alphabet)
        self._update_genome_features()


    # Updates the genome features with the mapper.
    #
    def _update_genome_features(self):
        # Get all feature endpoints to be mapped
        feature_endpoints = sorted(sum([feature_interval(feature).endpoints()
            for feature in self.genome.features], ()))
        mapped_positions = {}
        for endpoint in feature_endpoints:
            mapped_positions[endpoint] = self.mapper.get_mapping(endpoint)

        # Apply the mapping to all the features
        def map_feature_endpoints(feature):
            featureStart, featureEnd = feature_interval(feature).endpoints()
            feature.location = FeatureLocation(
                    ExactPosition(mapped_positions[featureStart]),
                    ExactPosition(mapped_positions[featureEnd]))
        [map_feature_endpoints(feature) for feature in self.genome.features]



###############################################################################
# Main procedure entrypoint
###############################################################################

def run(seq_record, output_root, vcf_path, **kwargs):
    rgm = RefGenomeMaker()

    # Store genome
    rgm.genome = seq_record

    # Parse variants
    rgm.variants = []
    with open(vcf_path) as fh:
        vcf_reader = vcf.Reader(fh)
        for index, record in enumerate(vcf_reader):
            position = record.POS - 1
            ref = record.REF

            # Validate sample
            sample = record.samples[0]
            if not sample.gt_bases:
                continue

            # HACK: For now, if we see a record that we should handle, we
            # assume that we should take the first alt that is different than
            # the ref.
            phase_char = sample.gt_phase_char()
            alts = sample.gt_bases.split(phase_char)
            for alt in alts:
                if alt != ref:
                    break

            if alt != ref:
                rgm.variants.append(Variant(position, ref, alt))

    rgm.apply_variants()

    # Output modified genome
    output_path = output_root + '.genbank'
    rgm.genome.seq.alphabet = DNAAlphabet()  # needed?
    SeqIO.write(rgm.genome, output_path, 'genbank')

    return rgm.genome

