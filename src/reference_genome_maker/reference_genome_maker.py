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
from variants import Variant

from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition

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

        # Now apply the primitive variants in order
        alt_genome, mapper = self._apply_primitive_variants(
                genome_seq, primitive_variants)

        # Update genome
        self.genome.seq = Seq(alt_genome, self.genome.seq.alphabet)
        self._update_genome_features(mapper)


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
                    locations, superstring = shortest_common_superstring(variant_alts)
                    primitive_variants.extend(Variant(
                        start, genome_seq[start:end], superstring
                        ).primitive_variants())
                    break
                end = max(end, variant.end())

        # Add sentinel.
        primitive_variants.append(Variant(len(genome_seq), '', ''))
        return primitive_variants

    def _apply_primitive_variants(self, genome_seq, primitive_variants):
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

        # Add sentinel.
        mapper.add_interval(Interval.create(ref_index, 1),
                Interval.create(alt_index, 1))
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

