from reference_genome_maker import *

# Test different types of variants:
# insertions, deletions, tandem duplications, inversions

insertions = [
    Variant(0, 'A', 'AC'),
    Variant(0, 'A', 'CA'),
    Variant(0, 'A', 'ACC'),
    Variant(0, 'A', 'ACGT'),
    Variant(0, '', 'G'),
    Variant(0, '', 'CCAGAT'),
    Variant(0, 'AC', 'ACG'),
    Variant(0, 'AC', 'AGC'),
    Variant(5, 'CCAGAT', 'CCAGTTACGAT'),
    ]
for variant in insertions:
    assert variant.get_type() == 'INS'

deletions = [
    Variant(0, 'A', ''),
    Variant(0, 'AC', 'A'),
    Variant(0, 'AC', 'C'),
    Variant(0, 'ACGC', 'GC'),
    Variant(0, 'ACGC', 'A'),
    Variant(0, 'AAAAA', 'AAA'),
    Variant(0, 'CCAGAT', 'CCAT'),
    Variant(5, 'CCAGAT', ''),
    Variant(7, 'CCAGTTACGAT', 'CCAGAT'),
    ]
for variant in deletions:
    assert variant.get_type() == 'DEL'

duplications = [
    Variant(0, 'A', 'AA'),
    Variant(0, 'AC', 'ACAC'),
    Variant(0, 'AC', 'ACC'),
    Variant(0, 'CGG', 'CCGG'),
    Variant(0, 'ACGC', 'ACGCACGC'),
    Variant(5, 'CCAGAT', 'CCAGATAGAT'),
    Variant(7, 'AAAA', 'AAAAAA'),
    ]
for variant in duplications:
    assert variant.get_type() == 'DUP:TANDEM'

inversions = [
    Variant(0, 'CCCAGTGTATCCATTTG', 'CCCCAAATGGATACACT'),
    Variant(5, 'AAAGGGGGGGGGGGGGAA', 'AAACCCCCCCCCCCCCAA'),
    ]
for variant in inversions:
    assert variant.get_type() == 'INV'

complexs = [
    Variant(0, 'A', 'CAT'),
    Variant(0, 'ACCA', 'CC'),
    Variant(0, 'ACCAG', 'ACACCAGG'),
    Variant(0, 'TGTGACTGATC', 'A')
    ]
for variant in complexs:
    assert variant.get_type() == 'COMPLEX'

