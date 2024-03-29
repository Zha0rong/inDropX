import csv
import os
from utils import *
import yaml
import multiprocessing
import time

class inDrop_Data_processing:
    pathtolibraryindex = ''
    pathtocellbarcode1 = ''
    pathtocellbarcode2umi = ''
    pathtorna = ''
    libraryindex = ''
    # Library index will be a dictionary with sample name as keys
    outputdir = ''
    cellbarcodewhitelist = []
    cellbarcodewhitelist_rev = []

    file_location = {}  # Record the location of every file in every sample: CB1, CB2, Read
    unfiltered_file_location = {}
    sd = False
    output_central={}
    Barcode_correction_dict = {}
    Barcode_correction_dict_rev = {}
    Demultiplexing_statistics=''
    post_trimming_length=20
    Total_Read = 0
    Invalid_Library_Index = 0
    Read_statistics = {}
    dictionary_for_fast_index_sample_search = {}
    strict=False
    process_to_use=2
    version=''
    V2toV3conversion=''
    W1='AAGGCGTCACAAGCAATCACT'
    rev_W1=reverse_compliment(W1)
    segregation_point=1000000
    def __init__(self, pathtolibraryindex=None, pathtocellbarcode1=None, pathtocellbarcode2umi=None, pathtorna=None, libraryindex=None,
                 outputdir=None,version = 'V3',post_trimming_length=20):
        try:
            if type(libraryindex) is dict:
                self.pathtolibraryindex = pathtolibraryindex
                self.pathtocellbarcode1 = pathtocellbarcode1
                self.pathtocellbarcode2umi = pathtocellbarcode2umi
                self.pathtorna = pathtorna
                self.libraryindex = libraryindex
                self.outputdir = outputdir
                self.Demultiplexing_statistics = '%s/Read_statistics_for_demultiplexing.tsv' % self.outputdir
                self.post_trimming_length=post_trimming_length
                if len(self.pathtolibraryindex) != len(self.pathtocellbarcode1) != len(self.pathtocellbarcode2umi) != len(self.pathtorna):
                    sys.exit('The fastq files number is wrong.')
                if version == 'V3':
                    with open('whitelist/cellbarcode.txt') as f:
                        for line in f:
                            self.cellbarcodewhitelist.append(line.rstrip())
                            self.cellbarcodewhitelist_rev.append(reverse_compliment(line.rstrip()))
                    f.close()
                    self.Barcode_correction_dict = build_barcode_neighborhoods(barcodelist=self.cellbarcodewhitelist)
                    self.Barcode_correction_dict_rev = build_barcode_neighborhoods(barcodelist=self.cellbarcodewhitelist_rev)

                    for sample in list(self.libraryindex.keys()):
                        CHECK_FOLDER = os.path.isdir('%s/%s' % (self.outputdir, sample))
                        if not CHECK_FOLDER:
                            os.mkdir('%s/%s' % (self.outputdir, sample))
                        self.output_central[sample] = {
                            'Unfiltered_CB1': gzip.open('%s/%s/%s.barcode1.fastq.gz' % (self.outputdir, sample, sample),
                                                        'wt'),
                            'Unfiltered_CB2': gzip.open('%s/%s/%s.barcode2.fastq.gz' % (self.outputdir, sample, sample),
                                                        'wt'),
                            'Unfiltered_RNA': gzip.open('%s/%s/%s.read.fastq.gz' % (self.outputdir, sample, sample),
                                                        'wt'),
                            'Filtered_CB': gzip.open(
                                '%s/%s/%s.filtered.barcodes.umi.fastq.gz' % (self.outputdir, sample, sample), 'wt'),
                            'Filtered_RNA': gzip.open(
                                '%s/%s/%s.filtered.read.fastq.gz' % (self.outputdir, sample, sample), 'wt'),
                            'Filtering.Statistics': {
                                'Total_read': 0,
                                'Valid_read': 0,
                                'Invalid_CB1': 0,
                                'Invalid_CB2': 0,
                                'Invalid_Both_CB': 0,
                                'Read_Too_Short_after_Trimming': 0},
                            'Filtering.Statistics.file': '%s/%s/%s.filtering_statistics.tsv' % (self.outputdir, sample, sample)
                            #'Cell.statistics': {},
                            #'Cell.statistics.file': '%s/%s/%s.cell_statistics.tsv' % (self.outputdir, sample, sample)

                        }

                elif version == 'V2':
                    with open('whitelist/version_2.cell.barcode.txt') as f:
                        for line in f:
                            self.cellbarcodewhitelist.append(line.rstrip())
                    f.close()
                    self.Barcode_correction_dict = build_barcode_neighborhoods(barcodelist=self.cellbarcodewhitelist)
                    self.V2toV3conversion={}
                    with open('whitelist/V2.to.V3.Conversion.Table.tsv') as f:
                        for line in f:
                            split = line.split('\t')
                            self.V2toV3conversion[split[0].strip('\n')] = split[1].strip('\n')

                    for sample in list(self.libraryindex.keys()):
                        CHECK_FOLDER = os.path.isdir('%s/%s' % (self.outputdir, sample))
                        if not CHECK_FOLDER:
                            os.mkdir('%s/%s' % (self.outputdir, sample))
                        self.output_central[sample] = {
                            'Filtered_CB': gzip.open(
                                '%s/%s/%s.filtered.barcodes.umi.fastq.gz' % (self.outputdir, sample, sample), 'wt'),
                            'Filtered_RNA': gzip.open(
                                '%s/%s/%s.filtered.read.fastq.gz' % (self.outputdir, sample, sample), 'wt'),
                            'Filtering.Statistics': {
                                'Total_read': 0,
                                'Valid_read': 0,
                                'Invalid_CB1': 0,
                                'Invalid_CB2': 0,
                                'Invalid_Both_CB': 0,
                                'Read_Too_Short_after_Trimming': 0},
                            'Filtering.Statistics.file': '%s/%s/%s.filtering_statistics.tsv' % (self.outputdir, sample, sample)
                            #'Cell.statistics': {},
                            #'Cell.statistics.file': '%s/%s/%s.cell_statistics.tsv' % (self.outputdir, sample, sample)
                        }
                else:
                    sys.exit('The version information is not recognized, the input should either be \'V2\' or \'V3\'. For data generated by other protocols please contact author Zhaorong Li for support.')
                # Check for the correct direction of library index.
            else:
                sys.exit('The library index needs to be a dictionary.')
        except Exception as e:
            print(str(e))

    def Demultiplexing_and_Correcting(self, strict=False):
        self.strict=strict
        Total_Read = 0
        self.Total_Read=Total_Read
        Invalid_Library_Index = 0
        self.Invalid_Library_Index=Invalid_Library_Index
        Read_statistics = {}
        for sample in self.libraryindex.keys():
            Read_statistics[sample] = 0
        self.Read_statistics=Read_statistics
        dictionary_for_fast_index_sample_search = {}
        for sample in self.libraryindex.keys():
            dictionary_for_fast_index_sample_search[self.libraryindex[sample]] = sample
        self.dictionary_for_fast_index_sample_search=dictionary_for_fast_index_sample_search
        for i in range(len(self.pathtolibraryindex)):
            read_processed=0
            time_period_start=time.time()
            for read in ParseFastq([[self.pathtolibraryindex[i], self.pathtocellbarcode1[i], self.pathtocellbarcode2umi[i],self.pathtorna[i]]]):
                name = read[0].strip('\n')
                librarybarcode = read[1][0].strip('\n')
                CB1read = read[1][1].strip('\n')
                CB2read = read[1][2].strip('\n')
                rnaread = read[1][3].strip('\n')
                CB1_qual = read[2][1].strip('\n')
                CB2_qual = read[2][2].strip('\n')
                rnaread_qual = read[2][3].strip('\n')
                informative_name = '%s %s:%s:%s' % (name, CB1read, CB2read[0:8], CB2read[8:])
                if str(librarybarcode) in dictionary_for_fast_index_sample_search:
                    Total_Read += 1
                    sample = dictionary_for_fast_index_sample_search[str(librarybarcode)]
                    Read_statistics[sample] += 1
                    write_fastq(self.output_central[sample]['Unfiltered_CB1'],informative_name,  CB1read,CB1_qual)
                    write_fastq(self.output_central[sample]['Unfiltered_CB2'],informative_name,  CB2read,CB2_qual)
                    write_fastq(self.output_central[sample]['Unfiltered_RNA'],informative_name,  rnaread,rnaread_qual)
                    '''Insert Correcting and Filtering here'''
                    self.output_central[sample]['Filtering.Statistics']['Total_read'] += 1
                    UMI = CB2read[8:]
                    trimmed_RNA = Trimmer(rnaread, rnaread_qual, min_length=self.post_trimming_length)
                    if trimmed_RNA[2] is True:
                        if CB1read in self.Barcode_correction_dict and (CB2read[0:8]) in self.Barcode_correction_dict_rev:
                            writeCB1 = self.Barcode_correction_dict[CB1read]
                            writeCB2 = (self.Barcode_correction_dict_rev[(CB2read[0:8])])
                            writename = ''
                            if name.startswith('@'):
                                writename = '%s %s:%s:%s' % (name.split(' ')[0], writeCB1, writeCB2, UMI)
                            else:
                                writename = '@%s %s:%s:%s' % (name.split(' ')[0], writeCB1, writeCB2, UMI)
                            truecellname = writeCB1 + writeCB2
                            writeCB = writeCB1 + writeCB2 + UMI
                            writeQual = CB1_qual + CB2_qual
                            #if truecellname in self.output_central[sample]['Cell.statistics']:
                            #    if UMI not in self.output_central[sample]['Cell.statistics'][truecellname][0]:
                            #        self.output_central[sample]['Cell.statistics'][truecellname][0].append(UMI)
                            #        self.output_central[sample]['Cell.statistics'][truecellname][1] += 1
                            #        self.output_central[sample]['Cell.statistics'][truecellname][2] += 1
                            #    else:
                            #        self.output_central[sample]['Cell.statistics'][truecellname][2] += 1
                            #else:
                            #    self.output_central[sample]['Cell.statistics'][truecellname] = [[], 0, 0]
                            #    self.output_central[sample]['Cell.statistics'][truecellname][0].append(UMI)
                            #    self.output_central[sample]['Cell.statistics'][truecellname][1] += 1
                            #    self.output_central[sample]['Cell.statistics'][truecellname][2] += 1
                            self.output_central[sample]['Filtering.Statistics']['Valid_read'] += 1
                            write_fastq(self.output_central[sample]['Filtered_CB'], writename, writeCB, writeQual)
                            write_fastq(self.output_central[sample]['Filtered_RNA'], writename, trimmed_RNA[0], trimmed_RNA[1])
                        else:
                            if CB1read not in self.Barcode_correction_dict and reverse_compliment(
                                    CB2read[0:8]) not in self.Barcode_correction_dict:
                                self.output_central[sample]['Filtering.Statistics']['Invalid_Both_CB'] += 1
                            else:
                                if CB1read not in self.Barcode_correction_dict:
                                    self.output_central[sample]['Filtering.Statistics']['Invalid_CB1'] += 1
                                else:
                                    self.output_central[sample]['Filtering.Statistics']['Invalid_CB2'] += 1
                    else:
                        self.output_central[sample]['Filtering.Statistics']['Read_Too_Short_after_Trimming'] += 1
                else:
                    if strict is False:
                        if min([hammingdistance(librarybarcode, libraryindex) for libraryindex in dictionary_for_fast_index_sample_search]) == 1:
                            correct_barcode = [libraryindex for libraryindex in dictionary_for_fast_index_sample_search if hammingdistance(librarybarcode, libraryindex) == 1][0]

                            sample = dictionary_for_fast_index_sample_search[str(correct_barcode)]
                            Total_Read += 1
                            Read_statistics[sample] += 1
                            write_fastq(self.output_central[sample]['Unfiltered_CB1'], informative_name, CB1read, CB1_qual)
                            write_fastq(self.output_central[sample]['Unfiltered_CB2'], informative_name, CB2read, CB2_qual)
                            write_fastq(self.output_central[sample]['Unfiltered_RNA'], informative_name, rnaread, rnaread_qual)
                            '''Insert Correcting and Filtering here'''
                            self.output_central[sample]['Filtering.Statistics']['Total_read'] += 1
                            UMI = CB2read[8:]
                            trimmed_RNA = Trimmer(rnaread, rnaread_qual, min_length=self.post_trimming_length)
                            if trimmed_RNA[2] is True:
                                if CB1read in self.Barcode_correction_dict and (CB2read[0:8]) in self.Barcode_correction_dict_rev:
                                    writeCB1 = self.Barcode_correction_dict[CB1read]
                                    writeCB2 = (self.Barcode_correction_dict_rev[(CB2read[0:8])])
                                    writename = ''
                                    if name.startswith('@'):
                                        writename = '%s %s:%s:%s' % (name.split(' ')[0], writeCB1, writeCB2, UMI)
                                    else:
                                        writename = '@%s %s:%s:%s' % (name.split(' ')[0], writeCB1, writeCB2, UMI)
                                    truecellname = writeCB1 + writeCB2
                                    writeCB = writeCB1 + writeCB2 + UMI
                                    writeQual = CB1_qual + CB2_qual
                                    #if truecellname in self.output_central[sample]['Cell.statistics']:
                                    #    if UMI not in self.output_central[sample]['Cell.statistics'][truecellname][0]:
                                    #        self.output_central[sample]['Cell.statistics'][truecellname][0].append(UMI)
                                    #        self.output_central[sample]['Cell.statistics'][truecellname][1] += 1
                                    #        self.output_central[sample]['Cell.statistics'][truecellname][2] += 1
                                    #    else:
                                    #        self.output_central[sample]['Cell.statistics'][truecellname][2] += 1
                                    #else:
                                    #    self.output_central[sample]['Cell.statistics'][truecellname] = [[], 0, 0]
                                    #    self.output_central[sample]['Cell.statistics'][truecellname][0].append(UMI)
                                    #    self.output_central[sample]['Cell.statistics'][truecellname][1] += 1
                                    #    self.output_central[sample]['Cell.statistics'][truecellname][2] += 1
                                    self.output_central[sample]['Filtering.Statistics']['Valid_read'] += 1
                                    write_fastq(self.output_central[sample]['Filtered_CB'], writename, writeCB,writeQual)
                                    write_fastq(self.output_central[sample]['Filtered_RNA'], writename, trimmed_RNA[0],trimmed_RNA[1])
                                else:
                                    if CB1read not in self.Barcode_correction_dict and (CB2read[0:8]) not in self.Barcode_correction_dict_rev:
                                        self.output_central[sample]['Filtering.Statistics']['Invalid_Both_CB'] += 1
                                    else:
                                        if CB1read not in self.Barcode_correction_dict:
                                            self.output_central[sample]['Filtering.Statistics']['Invalid_CB1'] += 1
                                        else:
                                             self.output_central[sample]['Filtering.Statistics']['Invalid_CB2'] += 1
                            else:
                                self.output_central[sample]['Filtering.Statistics']['Read_Too_Short_after_Trimming'] += 1
                        else:
                            Total_Read += 1
                            Invalid_Library_Index += 1
                    else:
                        Total_Read += 1
                        Invalid_Library_Index += 1
                read_processed +=1
                if read_processed != 0 and read_processed%self.segregation_point==0:
                    print('Time take to finish %s reads: %s seconds'%(self.segregation_point,time.time()-time_period_start))
                    time_period_start=time.time()


        for sample in self.libraryindex.keys():
            self.output_central[sample]['Unfiltered_CB1'].close()
            self.output_central[sample]['Unfiltered_CB2'].close()
            self.output_central[sample]['Unfiltered_RNA'].close()
            self.output_central[sample]['Filtered_CB'].close()
            self.output_central[sample]['Filtered_RNA'].close()
            with open(self.output_central[sample]['Filtering.Statistics.file'], "w",
                      newline="") as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow(['Category', 'Read Number', 'Percentage'])
                writer.writerow(['Total_read', self.output_central[sample]['Filtering.Statistics']['Total_read'],
                                 self.output_central[sample]['Filtering.Statistics']['Total_read'] / self.output_central[sample]['Filtering.Statistics']['Total_read']])
                writer.writerow(['Valid_read', self.output_central[sample]['Filtering.Statistics']['Valid_read'],
                                 self.output_central[sample]['Filtering.Statistics']['Valid_read'] / self.output_central[sample]['Filtering.Statistics']['Total_read']])
                writer.writerow(['Invalid_CB1', self.output_central[sample]['Filtering.Statistics']['Invalid_CB1'],
                                 self.output_central[sample]['Filtering.Statistics']['Invalid_CB1'] / self.output_central[sample]['Filtering.Statistics']['Total_read']])
                writer.writerow(['Invalid_CB2', self.output_central[sample]['Filtering.Statistics']['Invalid_CB2'],
                                 self.output_central[sample]['Filtering.Statistics']['Invalid_CB2'] / self.output_central[sample]['Filtering.Statistics']['Total_read']])
                writer.writerow(['Invalid_Both_CB', self.output_central[sample]['Filtering.Statistics']['Invalid_Both_CB'],
                                 self.output_central[sample]['Filtering.Statistics']['Invalid_Both_CB'] / self.output_central[sample]['Filtering.Statistics']['Total_read']])
                writer.writerow(['Read_Too_Short_after_Trimming', self.output_central[sample]['Filtering.Statistics']['Read_Too_Short_after_Trimming'],
                                 self.output_central[sample]['Filtering.Statistics']['Read_Too_Short_after_Trimming'] / self.output_central[sample]['Filtering.Statistics']['Total_read']])
            csvfile.close()
            #with open(self.output_central[sample]['Cell.statistics.file'], "w", newline="") as csvfile:
            #    writer = csv.writer(csvfile, delimiter='\t')
            #    writer.writerow(['Cellname', 'Number of unique UMI', 'Number of reads'])
            #    for cell in self.output_central[sample]['Cell.statistics'].keys():
            #        writer.writerow([cell, self.output_central[sample]['Cell.statistics'][cell][1], self.output_central[sample]['Cell.statistics'][cell][2]])
            #csvfile.close()
        with open(self.Demultiplexing_statistics, "w", newline="") as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(['Category', 'Read Number', 'Percentage'])
            writer.writerow(['Total number of Reads', Total_Read, Total_Read / Total_Read])
            for sample in Read_statistics.keys():
                writer.writerow([sample, Read_statistics[sample], Read_statistics[sample] / Total_Read])
            writer.writerow(
                ['Reads not associated with any sample', Invalid_Library_Index, Invalid_Library_Index / Total_Read])
        csvfile.close()
    def Correcting_V2(self):
        return 0
