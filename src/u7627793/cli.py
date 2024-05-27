from cogent3 import make_unaligned_seqs
from cogent3.core.alignment import SequenceCollection


__author__ = "YUTONG SHAO"
__copyright__ = "Copyright 2016-2021, YUTONG SHAO"
__credits__ = ["YUTONG SHAO"]
__license__ = "BSD"
__version__ = "2020.6.5"  # A DATE BASED VERSION
__maintainer__ = "YUTONG SHAO"
__email__ = "u7627793@anu.edu.au"
__status__ = "alpha"


def unique_kmers(seqs, k):
    """Returns a {seqname: set(str, ...)}"""
    result = {}
    kmer_to_seqs = {}
    seqs = make_unaligned_seqs(seqs)

    for s in seqs.seqs:
        name = s.name
        s = str(s)
        result[name] = set()
        for i in range(len(s) - k + 1):
            kmer = s[i : i + k]
            if kmer in kmer_to_seqs:
                kmer_to_seqs[kmer].add(name)
            else:
                kmer_to_seqs[kmer] = {name}

    for kmer, names in kmer_to_seqs.items():
        if len(names) == 1:
            unique_name = list(names)[0]
            result[unique_name].add(kmer)

    return result