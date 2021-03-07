import os
import csv
import subprocess
import sys
from collections import defaultdict, OrderedDict
from itertools import product, combinations
import gzip
import bz2
import utils_functions.py


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
    version=''

    def __init__(self, pathtolibraryindex, pathtocellbarcode1, pathtocellbarcode2umi, pathtorna, libraryindex,outputdir,version):
        if type(libraryindex) is dict:
            self.pathtolibraryindex = pathtolibraryindex
            self.pathtocellbarcode1 = pathtocellbarcode1
            self.pathtocellbarcode2umi = pathtocellbarcode2umi
            self.pathtorna = pathtorna
            self.libraryindex = libraryindex
            self.outputdir = outputdir
            self.version = version
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
            sys.exit('InDrop toolkit pipeline exiting, the library index fastq file does not exist.')
        if os.path.isfile(self.pathtocellbarcode1) is False:
            sys.exit('InDrop toolkit pipeline exiting, the Cellbarcode1 fastq file does not exist.')
        if os.path.isfile(self.pathtocellbarcode2umi) is False:
            sys.exit('InDrop toolkit pipeline exiting, the Cellbarcode2 and UMI fastq file does not exist.')
        if os.path.isfile(self.pathtorna) is False:
            sys.exit('InDrop toolkit pipeline exiting, the RNA read fastq file does not exist.')
        if os.path.isfile(self.pathtorna) is False:
            sys.exit('InDrop toolkit pipeline exiting, the RNA read fastq file does not exist.')
        if os.path.isdir(self.outputdir) is False:
            sys.exit('InDrop toolkit pipeline exiting, the RNA read fastq file does not exist.')
        if version not in ['V1','V2','V3']:
            sys.exit('InDrop toolkit only support V1, V2 and V3 for now.')
