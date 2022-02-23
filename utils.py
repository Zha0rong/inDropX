import os
import csv
import sys
from collections import defaultdict, OrderedDict
from itertools import product, combinations
import gzip
import bz2



def seq_neighborhood(seq, n_subs=1):
    """
    Given a sequence, yield all sequences within n_subs substitutions of
    that sequence by looping through each combination of base pairs within
    each combination of positions.
    """
    for positions in combinations(range(len(seq)), n_subs):
        # yields all unique combinations of indices for n_subs mutations
        for subs in product(*("ATGCN",) * n_subs):
            # yields all combinations of possible nucleotides for strings of length
            # n_subs
            seq_copy = list(seq)
            for p, s in zip(positions, subs):
                seq_copy[p] = s
            yield ''.join(seq_copy)


def build_barcode_neighborhoods(barcodelist):
    """
    Given a set of barcodes, produce sequences which can unambiguously be
    mapped to these barcodes, within 2 substitutions. If a sequence maps to
    multiple barcodes, get rid of it. However, if a sequences maps to a bc1 with
    1change and another with 2changes, keep the 1change mapping.
    """

    # contains all mutants that map uniquely to a barcode
    clean_mapping = dict()
    # contain single or double mutants
    mapping1 = defaultdict(set)
    mapping2 = defaultdict(set)
    # Build the full neighborhood and iterate through barcodes
    for line in barcodelist:
        barcode = line
        # each barcode obviously maps to itself uniquely
        clean_mapping[barcode] = barcode
        for n in seq_neighborhood(barcode, 1):
            mapping1[n].add(barcode)

        for n in seq_neighborhood(barcode, 2):
            mapping2[n].add(barcode)

            # take all single-mutants and find those that could only have come from one
    # specific barcode
    for k, v in mapping1.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]

    for k, v in mapping2.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]
    del mapping1
    del mapping2
    return clean_mapping


def reverse_compliment(sequence):
    reference = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    reversecompliment = ''.join(reference[x] for x in sequence[::-1])
    return reversecompliment


def average_phred_score_calculation(string):
    qualityscore = 0
    for i in range(len(string)):
        qualityscore += (ord(string[i]) - 33)
    return qualityscore / len(string)


def hammingdistance(string, reference):
    answer = 0
    if (len(string) == len(reference)):
        for i in range(len(string)):
            if string[i] != reference[i]:
                answer += 1
    else:
        answer = -1
    return answer


def ParseFastq(pathstofastqs):
    if pathstofastqs[0].endswith('.gz'):
        processes = [gzip.open(fastq) for fastq in pathstofastqs]
    elif pathstofastqs[0].endswith('.bz2'):
        processes = [bz2.open(fastq) for fastq in pathstofastqs]
    elif pathstofastqs[0].endswith('.fastq'):
        processes = [open(fastq) for fastq in pathstofastqs]
    else:
        sys.exit('The format of the file %s is not recognized.' % (str(pathstofastqs)))
    while True:
        names = [next(read).decode().split(' ')[0] for read in processes]
        Sequence = [next(read).decode() for read in processes]
        Blank = [next(read).decode() for read in processes]
        qualityscore = [next(read).decode() for read in processes]
        assert all(name == names[0] for name in names)
        if names:
            yield [names[0], Sequence, qualityscore]
        else:
            break
    for read in processes:
        read.close()

def write_fastq(file,ID, seq, quality_score):
    file.write('%s\n' % ID)
    file.write('%s\n' % seq)
    file.write('+\n')
    file.write('%s\n' % quality_score)
