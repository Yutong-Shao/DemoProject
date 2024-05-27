import pytest

from cogent3 import make_unaligned_seqs, load_unaligned_seqs
from u7627793.cli import unique_kmers


def test_find_unique_kmers_size_2():
    sequences = {"s1": "ATAATCC", "s2": "ATGATCC"}
    kmer_size = 2
    expected_result = {"s1": {"AA", "TA"}, "s2": {"GA", "TG"}}

    result = unique_kmers(sequences, kmer_size)
    assert result == expected_result


def test_find_unique_kmers_size_3():
    sequences = {"s1": "ATAATCC", "s2": "ATGATCC"}
    kmer_size = 3
    expected_result = {"s1": {"AAT", "ATA", "TAA"}, "s2": {"ATG", "GAT", "TGA"}}

    result = unique_kmers(sequences, kmer_size)
    assert result == expected_result


def test_find_unique_kmers_three_sequences():
    sequences = {"s1": "ATAATCC", "s2": "ATGATCC", "s3": "CTGATCC"}
    kmer_size = 2
    expected_result = {"s1": {"AA", "TA"}, "s2": set(), "s3": {"CT"}}

    result = unique_kmers(sequences, kmer_size)
    assert result == expected_result
    
def test_find_unique_kmers_fasta():
    path = "tests/data/sample.fasta"
    kmer_size = 2
    expected_result = {"I": set(), "L": {"CA"}}

    sequences = load_unaligned_seqs(path, moltype="dna")
    result = unique_kmers(sequences, kmer_size)
    assert result == expected_result
