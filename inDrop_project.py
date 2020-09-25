import os
import csv
import subprocess
import sys
from collections import defaultdict, OrderedDict
from itertools import product, combinations
import gzip
import bz2


def seq_neighborhood(seq, n_subs=1):
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


def write_fastq(ID, file_index, file, seq, quality_score):
    file.write('%s\n' % ID)
    file.write('%s\n' % seq)
    file.write('+\n')
    file.write('%s\n' % quality_score)


def ParseFastq( pathstofastqs):
    if pathstofastqs[0].endswith('.gz'):
        totalreads = [gzip.open(files) for files in pathstofastqs]
    elif pathstofastqs[0].endswith('.bz2'):
        totalreads = [bz2.open(files) for files in pathstofastqs]
    elif pathstofastqs[0].endswith('.fastq'):
        totalreads = [open(files) for files in pathstofastqs]
    else:
        sys.exit('The format of the file %s is not recognized.' % (str(pathstofastqs)))
    while True:
        try:
            names = [next(read).decode().split(' ')[0] for read in totalreads]
            Sequence = [next(read).decode() for read in totalreads]
            Blank = [next(read).decode() for read in totalreads]
            qualityscore = [next(read).decode() for read in totalreads]
            assert all(name == names[0] for name in names)
            if names:
                try:
                    yield [names[0], Sequence, qualityscore]
                except:
                    return
            else:
                break
        except StopIteration:
            break
    for read in totalreads:
        read.close()



class inDrop_Data_processing:
    pathtolibraryindex = ''
    pathtocellbarcode1 = ''
    pathtocellbarcode2umi = ''
    pathtorna = ''
    libraryindex = ''
    # Library index will be a dictionary with sample name as keys
    outputdir = ''
    cellbarcodewhitelist = []
    file_location = {}  # Record the location of every file in every sample: CB1CB2UMI, Read

    def __init__(self, pathtolibraryindex, pathtocellbarcode1, pathtocellbarcode2umi, pathtorna, libraryindex,
                 outputdir):
        if type(libraryindex) is dict:
            self.pathtolibraryindex = pathtolibraryindex
            self.pathtocellbarcode1 = pathtocellbarcode1
            self.pathtocellbarcode2umi = pathtocellbarcode2umi
            self.pathtorna = pathtorna
            self.libraryindex = libraryindex
            self.outputdir = outputdir
            with open('whitelist/cellbarcode.txt') as f:
                for line in f:
                    self.cellbarcodewhitelist.append(line.rstrip())
            f.close()
            # Check for the correct direction of library index.
            for sample in list(self.libraryindex.keys()):
                self.file_location[sample] = [
                        '%s/%s_CBUMI.fastq.gz' % (self.outputdir, sample + '_' + str(self.libraryindex[sample])),
                        '%s/%s_read.fastq.gz' % (self.outputdir, sample + '_' + str(self.libraryindex[sample]))]
        else:
            sys.exit('The library index needs to be a dictionary.')
        if os.path.isfile(self.pathtolibraryindex) is False:
            sys.exit('InDropX pipeline exiting, the library index fastq file does not exist.')
        if os.path.isfile(self.pathtocellbarcode1) is False:
            sys.exit('InDropX pipeline exiting, the Cellbarcode1 fastq file does not exist.')
        if os.path.isfile(self.pathtocellbarcode2umi) is False:
            sys.exit('InDropX pipeline exiting, the Cellbarcode2 and UMI fastq file does not exist.')
        if os.path.isfile(self.pathtorna) is False:
            sys.exit('InDropX pipeline exiting, the RNA read fastq file does not exist.')
        if os.path.isfile(self.pathtorna) is False:
            sys.exit('InDropX pipeline exiting, the RNA read fastq file does not exist.')
        if os.path.isdir(self.outputdir) is False:
            sys.exit('InDropX pipeline exiting, the RNA read fastq file does not exist.')