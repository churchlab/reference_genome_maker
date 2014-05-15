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
from algorithm import find_overlapping_segments
from variants import Variant

from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation

# Maximum number of qualifiers for variants that will be attached to
#   each feature.
MAX_VARIANT_QUALIFIERS = 5

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

        # Get a list of all primitive variants to be applied.
        primitive_variants = self._find_all_primitive_variants(
                genome_seq, self.variants)

        # Add a qualifier to each feature for every variant intersecting
        #   the feature.
        self._qualify_genome_features(primitive_variants)

        # Now apply the primitive variants in order
        alt_genome, mapper = self._apply_primitive_variants(
                genome_seq, primitive_variants)

        # Update genome
        self.genome.seq = Seq(alt_genome, self.genome.seq.alphabet)
        self._update_genome_features(mapper)
        self._add_variant_features(primitive_variants)


    # Get all primitive elements, and store into self.primitive_variants
    def _find_all_primitive_variants(self, genome_seq, variants):
        # Find groups of variants whose ref sequences overlap.
        # Find the shortest common superstring of the alts of these variants.
        # Create a variant out of (overlapping refs, SCS)
        primitive_variants = []
        generator = iter(variants)
        variant = next(generator, None)
        while variant:
            start, end = variant.interval().endpoints()
            variant_alts = []
            assert variant.ref == variant.interval().sliceInto(genome_seq)

            while True:
                variant_alts.append(variant.alt)
                variant = next(generator, None)
                if not variant or variant.start() >= end:
                    superstring = shortest_common_superstring(variant_alts)
                    primitive_variants.extend(Variant(
                        start, genome_seq[start:end], superstring
                        ).primitive_variants())
                    break
                end = max(end, variant.end())
        return primitive_variants

    def _qualify_genome_features(self, primitive_variants):
        # Perform a line sweep of the variant and feature endpoints,
        #   in order of position.
        overlaps = find_overlapping_segments(
                [feature_interval(f) for f in self.genome.features],
                [variant.interval() for variant in primitive_variants],
                MAX_VARIANT_QUALIFIERS)

        # Add the qualifiers.
        for index, variantIndices in overlaps.items():
            feature = self.genome.features[index]
            if 'variants' not in feature.qualifiers:
                feature.qualifiers['variants'] = []
            for variantIndex in variantIndices:
                feature.qualifiers['variants'].append(
                        str(primitive_variants[variantIndex]))

    def _apply_primitive_variants(self, genome_seq, primitive_variants):
        # Add a sentinel variant at the end.
        primitive_variants = primitive_variants[:]
        primitive_variants.append(Variant(len(genome_seq), '', ''))

        # Repeatedly add the identically mapped parts of the genome,
        #   then the primitive variant alt, until we reach the end.
        alt_genome = ''
        mapper = IntervalMapper()
        ref_index, alt_index = 0, 0  # index of where we are in ref and alt
        for variant in primitive_variants:
            # First append unchanged part of ref genome
            # REF: ******.............
            # ALT: ******...............
            ref_interval = Interval(ref_index, variant.start())
            alt_interval = Interval.create(alt_index, len(ref_interval))
            alt_genome += ref_interval.sliceInto(genome_seq)
            mapper.add_interval(ref_interval, alt_interval)
            alt_index += len(ref_interval)

            # Now apply the variant
            # REF: ......****.........
            # ALT: ......******.........
            alt_interval = Interval.create(alt_index, len(variant.alt))
            alt_genome += variant.alt
            mapper.add_interval(variant.interval(), alt_interval)
            ref_index = variant.end()
            alt_index += len(variant.alt)

        return alt_genome, mapper

    # Updates the genome features with the mapper.
    #
    def _update_genome_features(self, mapper):
        # Get all feature endpoints to be mapped
        feature_endpoints = sorted(sum([feature_interval(feature).endpoints()
            for feature in self.genome.features], ()))
        mapped_positions = {}
        for endpoint in feature_endpoints:
            mapped_positions[endpoint] = mapper.get_mapping(endpoint)

        # Apply the mapping to all the features
        for feature in self.genome.features:
            feature.location = FeatureLocation(
                    *map(lambda endpoint: mapped_positions[endpoint],
                        feature_interval(feature).endpoints()))

    # Add the variants as annotations to the reference genome
    #
    def _add_variant_features(self, variants):
        for variant in variants:
            self.genome.features.append(SeqFeature(
                type='variation',
                location=FeatureLocation(*variant.interval().endpoints()),
                strand=1,
                ))


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

