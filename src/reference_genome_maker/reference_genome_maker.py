"""
This module processes a reference genome (SeqRecord)
and a variant set (an object in this module),
and outputs a modified genome (SeqRecord).

Methods also exist for parsing either a csv or vcf
into a variant set.
"""

import re
import vcf

from util import Interval
from util import IntervalMapper
from util import feature_interval
from algorithm import shortest_common_superstring

from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
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


###############################################################################
# Reference Genome Maker
###############################################################################

class RefGenomeMaker:

    def apply_variants(self):
        """Applies stored variants to the given ref genome (string).

        Returns a tuple (A, M) where:
            A is the modified genome (string)
            M is a mapper object (for determining mappings from ref to mod)
        """

        # First, sort the variants by position
        self.variants.sort()

        # Now we use the following algorithm to modify the genome.
        # Find groups of variants whose ref sequences overlap.
        # Find the shortest common superstring of the alts of these variants.
        # Replace the refs with the shortest common superstring.

        # Iterate through each group of variants, gradually appending
        #   to alt_genome, which contains the bases of the modified genome.
        alt_genome = ''
        ref_index = 0  # index of where we are in the ref genome
        alt_index = 0  # index of where we are in the alt genome
        self.mapper = IntervalMapper()
        for variant_group in self._get_variant_groups(self.variants):
            # First append the unchanged part of ref
            variant_begin = variant_group[0].start()
            ref_interval = Interval(ref_index, variant_begin)
            alt_genome += ref_interval.sliceInto(str(self.genome.seq))
            self.mapper.add_interval(ref_interval,
                    Interval.create(alt_index, len(ref_interval)))

            # Find shortest common superstring of all variant alts
            locations, superstring = shortest_common_superstring(
                    [variant.alt for variant in variant_group])
            alt_genome += superstring
            for index, variant in enumerate(variant_group):
                alt_begin = variant_begin + locations[index]
                self.mapper.add_interval(variant.interval(),
                        Interval.create(alt_begin, len(variant.alt)))

            # Update ref and alt indices
            ref_index = max([variant.end() for variant in variant_group])
            alt_index += len(ref_interval) + len(superstring)

        # Update genome
        self.genome.seq = Seq(alt_genome, self.genome.seq.alphabet)
        self._update_genome_features()


    # Generator function to produce a group of variants with
    #   overlapping ref sequences.
    # variants must be in order of start position.
    # This function also yields a sentinel variant group at the end.
    #
    def _get_variant_groups(self, variants):
        variant_group = []
        end = -1  # max of variant.ends in variant_group
        for variant in variants:
            if variant_group and variant.start() >= end and \
                    variant_group[-1] != variant:  # check for duplicates
                yield variant_group
                variant_group = []
            end = max(end, variant.end())
            variant_group.append(variant)
        yield variant_group
        yield [Variant(len(self.genome), '$', '')]  # sentinel

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

            # Special handling for long ref.
            maybe_long_ref_match = re.match('LONG:(\d+)', ref)
            if maybe_long_ref_match:
                ref_size = int(maybe_long_ref_match.group(1))
                ref = str(seq_record.seq[position:position + ref_size])

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
                alt = alt if not alt == '<DEL>' else ''
                rgm.variants.append(Variant(position, ref, alt))

    rgm.apply_variants()

    # Output modified genome
    output_path = output_root + '.genbank'
    rgm.genome.seq.alphabet = DNAAlphabet()  # needed?
    SeqIO.write(rgm.genome, output_path, 'genbank')

    return rgm.genome

