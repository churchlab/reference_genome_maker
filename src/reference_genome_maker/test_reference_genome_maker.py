from random import randint
from util import feature_interval
from reference_genome_maker import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

r = RefGenomeMaker()

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(3, 'CCC', ''),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'AAAGGG'

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(0, 'AAA', 'TA'),
    Variant(6, 'GGG', 'AT'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'TACCCAT'

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(0, 'AAA', 'TTT'),
    Variant(3, 'CCC', 'TTT'),
    Variant(6, 'GGG', 'TTT'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'TTTTTTTTT'

r.genome = SeqRecord(Seq('AAAAA'))
r.variants = [
    Variant(1, 'A', 'C'),
    Variant(2, 'A', 'C'),
    Variant(3, 'A', 'C'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'ACCCA'

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(3, 'CCC', ''),
    Variant(4, 'CC', 'TT'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'AAATTGGG'

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(3, 'CCC', 'TT'),
    Variant(6, '', 'CCC'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'AAATTCCCGGG'

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(2, 'AC', 'TT'),
    Variant(3, 'CCC', ''),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'AATTGGG'

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(2, 'AC', 'AT'),
    Variant(3, 'CCC', ''),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'AAATGGG'

r.genome = SeqRecord(Seq('AAAGGG'))
r.variants = [
    Variant(2, 'AG', 'TT'),
    Variant(2, 'AG', 'CCC'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'AATTCCCGG'

r.genome = SeqRecord(Seq('AAAGGG'))
r.variants = [
    Variant(1, 'AAG', 'CCC'),
    Variant(2, 'AGG', 'TTT'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'ACCCTTTG'

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(0, 'AAACCC', 'GGGTTT'),
    Variant(3, 'CCCGGG', 'TTTAAA'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'GGGTTTAAA'

r.genome = SeqRecord(Seq('TTTT'))
r.variants = [
    Variant(2, '', 'AA'),
    Variant(2, '', 'AA'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'TTAAAATT'

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(3, '', 'CCC'),
    Variant(4, 'CC', 'TT'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'AAACCCCTTGGG'

r.genome = SeqRecord(Seq('AAAGGG'))
r.variants = [
    Variant(2, 'AG', 'TT'),
    Variant(3, '', 'CCC'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'AATTCCCGG'

r.genome = SeqRecord(Seq('AAACCCGGG'))
r.variants = [
    Variant(0, 'AAACCC', 'GGGTTT'),
    Variant(3, 'CCCGGG', 'AAATTT'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'GGGTTTAAATTT'

r.genome = SeqRecord(Seq('GGGAAACCCGGG'))
r.variants = [
    Variant(3, 'AAACCC', 'GGGTTT'),
    Variant(6, 'CCCGGG', 'TTTAAA'),
    ]
r.apply_variants()
assert str(r.genome.seq) == 'GGGGGGTTTAAA'

# Multiple variants, some overlapping
#   and test features
r.genome = SeqRecord(Seq('AAACCCGGGTTT'))
r.genome.features = [
    SeqFeature(FeatureLocation(ExactPosition(1), ExactPosition(5))),
    SeqFeature(FeatureLocation(ExactPosition(8), ExactPosition(10))),
    SeqFeature(FeatureLocation(ExactPosition(9), ExactPosition(11))),
    ]
r.variants = [
    Variant(1, 'AA', 'TT'),
    Variant(2, 'AC', 'TG'),
    Variant(3, 'CC', 'GG'),
    Variant(6, 'GG', 'AC'),
    Variant(7, 'GGT', 'CCCC'),
    Variant(9, 'TT', ''),
    ]
r.apply_variants()
new_features = map(feature_interval, r.genome.features)
assert str(r.genome.seq) == 'ATTGGCACCCCT'
assert new_features[0].endpoints() == (1, 5)
assert new_features[1].endpoints() == (8, 10)
assert new_features[2].endpoints() == (9, 11)

# Large test case; test efficiency
REFLEN, VARLEN, NUMVARS = 100000, 100, 1000
r.genome = SeqRecord(Seq('A' * REFLEN))
r.variants = []
for i in range(NUMVARS):
    varlen = randint(0, VARLEN)
    pos = randint(0, REFLEN - varlen)
    r.variants.append(
            Variant(pos, 'A' * varlen, 'T' * varlen))
r.apply_variants()

