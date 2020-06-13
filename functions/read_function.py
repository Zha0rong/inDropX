import os
import csv
import subprocess
import sys
from collections import defaultdict, OrderedDict
from itertools import product, combinations
import gzip



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
    reference={'A':'T','C':'G','T':'A','G':'C','N':'N'}
    reversecompliment=''.join(reference[x] for x in sequence[::-1])
    return reversecompliment

def average_phred_score_calculation(string):
    qualityscore=0
    for i in range(len(string)):
        qualityscore+=(ord(string[i])-33)
    return qualityscore/len(string)

def hammingdistance(string,reference):
    answer=0
    if (len(string)==len(reference)):
        for i in range(len(string)):
            if string[i]!=reference[i]:
                answer+=1
    else:
        answer=-1
    return answer

class inDrop_Data_processing:
    pathtolibraryindex=''
    pathtocellbarcode1=''
    pathtocellbarcode2umi=''
    pathtorna=''
    libraryindex=''
    #Library index will be a dictionary with sample name as keys
    outputdir=''
    cellbarcodewhitelist = []
    file_location={} #Record the location of every file in every sample: CB1, CB2, Read
    unfiltered_file_location={}
    mutli_Lane=False
    def __init__(self,pathtolibraryindex,pathtocellbarcode1,pathtocellbarcode2umi,pathtorna,libraryindex,outputdir):
        try:
            if type(libraryindex) is dict:
                self.pathtolibraryindex=pathtolibraryindex
                self.pathtocellbarcode1=pathtocellbarcode1
                self.pathtocellbarcode2umi=pathtocellbarcode2umi
                self.pathtorna=pathtorna
                self.libraryindex=libraryindex
                self.outputdir=outputdir
                with open('whitelist/cellbarcode.txt') as f:
                    for line in f:
                        self.cellbarcodewhitelist.append(line.rstrip())
                f.close()
                #Check for the correct direction of library index.
                for sample in list(self.libraryindex.keys()):
                    self.file_location[sample]=['%s/%s_read.fastq'%(self.outputdir,sample+'_'+str(self.libraryindex[sample])),
                                                '%s/%s_barcode1.fastq'%(self.outputdir,sample+'_'+str(self.libraryindex[sample])),
                                                '%s/%s_barcode2.fastq'%(self.outputdir,sample+'_'+str(self.libraryindex[sample]))]
                    self.unfiltered_file_location[sample]={'RNA':'%s/%s_read.fastq.gz'%(self.outputdir,sample+'_'+str(self.libraryindex[sample])),
                                                'CB1':'%s/%s_barcode1.fastq.gz'%(self.outputdir,sample+'_'+str(self.libraryindex[sample])),
                                                'CB2':'%s/%s_barcode2.fastq.gz'%(self.outputdir,sample+'_'+str(self.libraryindex[sample]))}
                if type(pathtolibraryindex) is list:
                    self.mutli_Lane=True
            else:
                sys.exit('The sample index is wrong.')
        except Exception as e:
            print(str(e))
    def _write_fastq(self,ID,file_index,file,seq,quality_score):
        if file_index=='R1':
            file.write('%s\n'%ID)
            file.write('%s\n' % seq)
            file.write('+\n')
            file.write('%s\n' % quality_score)
        else:
            file.write('%s\n'%ID)
            file.write('%s\n' % seq)
            file.write('+\n')
            file.write('%s\n' % quality_score)
    def _ParseFastq(self,pathstofastqs):
        if pathstofastqs[0].endswith('.gz'):
            processes=[subprocess.Popen(['zcat',(fastq)],stdout=subprocess.PIPE) for fastq in pathstofastqs]
            totalreads = [r.stdout for r in processes]
        elif pathstofastqs[0].endswith('.bz2'):
            processes=[subprocess.Popen(['bzcat',(fastq)],stdout=subprocess.PIPE) for fastq in pathstofastqs]
            totalreads = [r.stdout for r in processes]
        elif pathstofastqs[0].endswith('.fastq'):
            processes=[subprocess.Popen(['cat',(fastq)],stdout=subprocess.PIPE) for fastq in pathstofastqs]
            totalreads = [r.stdout for r in processes]
        else:
            sys.exit('The format of the file %s is not recognized.'%(str(pathtofastq)))
        while True:
            names=[next(read).decode().split(' ')[0] for read in totalreads]
            Sequence=[next(read).decode() for read in totalreads]
            Blank=[next(read).decode() for read in totalreads]
            qualityscore= [next(read).decode() for read in totalreads]
            assert all(name==names[0] for name in names)
            if names:
                yield [names[0], Sequence, qualityscore]
            else:
                break
        for read in totalreads:
            read.close()
    def Demultiplex(self,strict=True):
        Total_Read=0
        Invalid_Library_Index=0
        Read_statistics={}
        #Parse fastq files into a dictionary
        readfile={}
        for sample in self.libraryindex.keys():
            readfile[sample]=[gzip.open('%s/%s_read.fastq.gz'%(self.outputdir,sample+'_'+str(self.libraryindex[sample])),'wt'),
                              gzip.open('%s/%s_barcode1.fastq.gz'%(self.outputdir,sample+'_'+str(self.libraryindex[sample])),'wt'),
                              gzip.open('%s/%s_barcode2.fastq.gz'%(self.outputdir,sample+'_'+str(self.libraryindex[sample])),'wt')]
            #Read, CB1, CB2
            Read_statistics[sample]=0
        dictionary_for_fast_index_sample_search={}
        for sample in self.libraryindex.keys():
            dictionary_for_fast_index_sample_search[self.libraryindex[sample]]=sample
        if self.mutli_Lane is True:
            for i in range(len(self.pathtolibraryindex)):
                for read in self._ParseFastq([self.pathtolibraryindex[i],self.pathtocellbarcode1[i],self.pathtocellbarcode2umi[i],self.pathtorna[i]]):
                    name=read[0].strip('\n')
                    librarybarcode=read[1][0].strip('\n')
                    CB1read=read[1][1].strip('\n')
                    CB2read=read[1][2].strip('\n')
                    rnaread=read[1][3].strip('\n')
                    CB1_qual=read[2][1].strip('\n')
                    CB2_qual=read[2][2].strip('\n')
                    rnaread_qual=read[2][3].strip('\n')
                    informative_name='%s %s:%s:%s'%(name,CB1read,CB2read[0:8],CB2read[8:])
                    if str(librarybarcode) in list(dictionary_for_fast_index_sample_search.keys()):
                        Total_Read+=1
                        sample=dictionary_for_fast_index_sample_search[str(librarybarcode)]
                        Read_statistics[sample]+=1
                        self._write_fastq(informative_name,'R1',readfile[sample][1],CB1read,CB1_qual)#ID,file_index,file,seq,quality_score
                        self._write_fastq(informative_name,'R1',readfile[sample][2],CB2read,CB2_qual)#ID,file_index,file,seq,quality_score
                        self._write_fastq(informative_name,'R1',readfile[sample][0],rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                    else:
                        if strict is False:
                            if min([hammingdistance(librarybarcode,libraryindex) for libraryindex in list(dictionary_for_fast_index_sample_search.keys())])==1:
                                correct_barcode=[libraryindex for libraryindex in list(dictionary_for_fast_index_sample_search.keys()) if hammingdistance(librarybarcode,libraryindex)==1][0]
                                sample=dictionary_for_fast_index_sample_search[str(correct_barcode)]
                                Total_Read+=1
                                Read_statistics[sample]+=1
                                self._write_fastq(informative_name,'R1',readfile[sample][1],CB1read,CB1_qual)#ID,file_index,file,seq,quality_score
                                self._write_fastq(informative_name,'R1',readfile[sample][2],CB2read,CB2_qual)#ID,file_index,file,seq,quality_score
                                self._write_fastq(informative_name,'R1',readfile[sample][0],rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                            else:
                                Total_Read+=1
                                Invalid_Library_Index+=1
                        else:
                            Total_Read+=1
                            Invalid_Library_Index+=1
        else:
            for read in self._ParseFastq([self.pathtolibraryindex,self.pathtocellbarcode1,self.pathtocellbarcode2umi,self.pathtorna]):
                name=read[0].strip('\n')
                librarybarcode=read[1][0].strip('\n')
                CB1read=read[1][1].strip('\n')
                CB2read=read[1][2].strip('\n')
                rnaread=read[1][3].strip('\n')
                CB1_qual=read[2][1].strip('\n')
                CB2_qual=read[2][2].strip('\n')
                rnaread_qual=read[2][3].strip('\n')
                informative_name='%s %s:%s:%s'%(name,CB1read,CB2read[0:8],CB2read[8:])
                if str(librarybarcode) in list(dictionary_for_fast_index_sample_search.keys()):
                    Total_Read+=1
                    sample=dictionary_for_fast_index_sample_search[str(librarybarcode)]
                    Read_statistics[sample]+=1
                    self._write_fastq(informative_name,'R1',readfile[sample][1],CB1read,CB1_qual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(informative_name,'R1',readfile[sample][2],CB2read,CB2_qual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(informative_name,'R1',readfile[sample][0],rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                else:
                    if strict is False:
                        if min([hammingdistance(librarybarcode,libraryindex) for libraryindex in list(dictionary_for_fast_index_sample_search.keys())])==1:
                            correct_barcode=[libraryindex for libraryindex in list(dictionary_for_fast_index_sample_search.keys()) if hammingdistance(librarybarcode,libraryindex)==1][0]
                            sample=dictionary_for_fast_index_sample_search[str(correct_barcode)]
                            Total_Read+=1
                            Read_statistics[sample]+=1
                            self._write_fastq(informative_name,'R1',readfile[sample][1],CB1read,CB1_qual)#ID,file_index,file,seq,quality_score
                            self._write_fastq(informative_name,'R1',readfile[sample][2],CB2read,CB2_qual)#ID,file_index,file,seq,quality_score
                            self._write_fastq(informative_name,'R1',readfile[sample][0],rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                        else:
                            Total_Read+=1
                            Invalid_Library_Index+=1
                    else:
                        Total_Read+=1
                        Invalid_Library_Index+=1
        for sample in readfile:
            for file in readfile[sample]:
                file.close()
        with open('%s/Read_statistics_for_demuliplexing.tsv'%(self.outputdir), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Category','Read Number','Percentage'])
            writer.writerow(['Total number of Reads',Total_Read,Total_Read/Total_Read])
            for sample in Read_statistics.keys():
                writer.writerow([sample,Read_statistics[sample],Read_statistics[sample]/Total_Read])
            writer.writerow(['Reads not associated with any sample',Invalid_Library_Index,Invalid_Library_Index/Total_Read])
        csvfile.close()
    def correct_and_filter_kallisto(self):
        Barcode_correction_dict=build_barcode_neighborhoods(barcodelist=self.cellbarcodewhitelist)
        for sample in self.libraryindex.keys():
            filtering_output_directory='%s/%s_kallisto'%(self.outputdir,sample)
            try:
                os.mkdir(filtering_output_directory)
            except FileExistsError:
                # directory already exists
                pass
            filtering_statistics={
            'Total_read':0,
            'Valid_read':0,
            'Invalid_CB1':0,
            'Invalid_CB2':0,
            'Invalid_Both_CB':0
            }
            Cell_statistics={}#Structure: Cell_barcode:[[UMI_list],number of unique umis per cell, number of reads per cell]
            CB1=gzip.open('%s/%s_filtered_CB1.fastq.gz'%(filtering_output_directory,sample),'wt')
            CB2=gzip.open('%s/%s_filtered_CB2.fastq.gz'%(filtering_output_directory,sample),'wt')
            RNA_read=gzip.open('%s/%s_filtered_Reads.fastq.gz'%(filtering_output_directory,sample),'wt')
            for read in self._ParseFastq([self.unfiltered_file_location[sample]['CB1'],self.unfiltered_file_location[sample]['CB2'],self.unfiltered_file_location[sample]['RNA']]):
                name=read[0].strip('\n')
                CB1read=read[1][0].strip('\n')
                CB2read=read[1][1].strip('\n')
                rnaread=read[1][2].strip('\n')
                CB1_qual=read[2][0].strip('\n')
                CB2_qual=read[2][1].strip('\n')
                rnaread_qual=read[2][2].strip('\n')
                filtering_statistics['Total_read']+=1
                umi=CB2read[8:]
                if CB1read in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) in Barcode_correction_dict:
                    writeCB1=Barcode_correction_dict[CB1read]
                    writeCB2=reverse_compliment(Barcode_correction_dict[reverse_compliment(CB2read[0:8])])
                    writename=''
                    if name.startswith('@'):
                        writename='%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    else:
                        writename='@%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    truecellname=writeCB1+writeCB2
                    if truecellname in Cell_statistics:
                        if umi not in Cell_statistics[truecellname][0]:
                            Cell_statistics[truecellname][0].append(umi)
                            Cell_statistics[truecellname][1]+=1
                            Cell_statistics[truecellname][2]+=1
                        else:
                            Cell_statistics[truecellname][2]+=1
                    else:
                        Cell_statistics[truecellname]=[[],0,0]
                        Cell_statistics[truecellname][0].append(umi)
                        Cell_statistics[truecellname][1]+=1
                        Cell_statistics[truecellname][2]+=1
                    filtering_statistics['Valid_read']+=1
                    self._write_fastq(writename,'R1',CB1,writeCB1,CB1_qual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(writename,'R1',CB2,writeCB2+umi,CB2_qual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(writename,'R1',RNA_read,rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                else:
                    if CB1read not in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) not in Barcode_correction_dict:
                        filtering_statistics['Invalid_Both_CB']+=1
                    else:
                        if CB1read not in Barcode_correction_dict:
                            filtering_statistics['Invalid_CB1']+=1
                        else:
                            filtering_statistics['Invalid_CB2']+=1
            CB1.close()
            CB2.close()
            RNA_read.close()
            #os.system('gzip %s'%('%s/%s_filtered_CB1.fastq'%(filtering_output_directory,sample)))
            #os.system('gzip %s'%('%s/%s_filtered_CB2.fastq'%(filtering_output_directory,sample)))
            #os.system('gzip %s'%('%s/%s_filtered_Reads.fastq'%(filtering_output_directory,sample)))
            for cell in Cell_statistics:
                Cell_statistics[cell][0]=0
            with open('%s/%s_filtering_statistics.tsv'%(filtering_output_directory,sample), "w", newline="") as csvfile:
                writer=csv.writer(csvfile,delimiter='\t')
                writer.writerow(['Category','Read Number','Percentage'])
                writer.writerow(['Total_read',filtering_statistics['Total_read'],filtering_statistics['Total_read']/filtering_statistics['Total_read']])
                writer.writerow(['Valid_read',filtering_statistics['Valid_read'],filtering_statistics['Valid_read']/filtering_statistics['Total_read']])
                writer.writerow(['Invalid_CB1',filtering_statistics['Invalid_CB1'],filtering_statistics['Invalid_CB1']/filtering_statistics['Total_read']])
                writer.writerow(['Invalid_CB2',filtering_statistics['Invalid_CB2'],filtering_statistics['Invalid_CB2']/filtering_statistics['Total_read']])
                writer.writerow(['Invalid_Both_CB',filtering_statistics['Invalid_Both_CB'],filtering_statistics['Invalid_Both_CB']/filtering_statistics['Total_read']])
            csvfile.close()
            with open('%s/%s_Cell_statistics.tsv'%(filtering_output_directory,sample), "w", newline="") as csvfile:
                writer=csv.writer(csvfile,delimiter='\t')
                writer.writerow(['Cellname','Number of unique UMI','Number of reads'])
                for cell in Cell_statistics.keys():
                    writer.writerow([cell,Cell_statistics[cell][1],Cell_statistics[cell][2]])
            csvfile.close()
    def Demultiplex_and_correct_filter(self,strict=False):
        Total_Read=0
        Invalid_Library_Index=0
        Read_statistics={}
        #Parse fastq files into a dictionary
        readfile={}
        for sample in self.libraryindex.keys():
            readfile[sample]=[gzip.open('%s/%s_read.fastq.gz'%(self.outputdir,sample+'_'+str(self.libraryindex[sample])),'wt'),
                              gzip.open('%s/%s_barcode.fastq.gz'%(self.outputdir,sample+'_'+str(self.libraryindex[sample])),'wt')]
            #Read, CB1, CB2
            Read_statistics[sample]=0
        dictionary_for_fast_index_sample_search={}
        for sample in self.libraryindex.keys():
            dictionary_for_fast_index_sample_search[self.libraryindex[sample]]=sample		
	def correct_and_filter_solo(self):
        Barcode_correction_dict=build_barcode_neighborhoods(barcodelist=self.cellbarcodewhitelist)
        for sample in self.libraryindex.keys():
            filtering_output_directory='%s/%s_solo'%(self.outputdir,sample)
            try:
                os.mkdir(filtering_output_directory)
            except FileExistsError:
                # directory already exists
                pass
            filtering_statistics={
            'Total_read':0,
            'Valid_read':0,
            'Invalid_CB1':0,
            'Invalid_CB2':0,
            'Invalid_Both_CB':0
            }
            Cell_statistics={}#Structure: Cell_barcode:[[UMI_list],number of unique umis per cell, number of reads per cell]
            CB=gzip.open('%s/%s_STAR_CB.fastq.gz'%(filtering_output_directory,sample),'wt')
            RNA_read=gzip.open('%s/%s_STAR_Reads.fastq.gz'%(filtering_output_directory,sample),'wt')
            for read in self._ParseFastq([self.unfiltered_file_location[sample]['CB1'],self.unfiltered_file_location[sample]['CB2'],self.unfiltered_file_location[sample]['RNA']]):
                name=read[0].strip('\n')
                CB1read=read[1][0].strip('\n')
                CB2read=read[1][1].strip('\n')
                rnaread=read[1][2].strip('\n')
                CB1_qual=read[2][0].strip('\n')
                CB2_qual=read[2][1].strip('\n')
                rnaread_qual=read[2][2].strip('\n')
                filtering_statistics['Total_read']+=1
                umi=CB2read[8:]
                if CB1read in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) in Barcode_correction_dict:
                    writeCB1=Barcode_correction_dict[CB1read]
                    writeCB2=reverse_compliment(Barcode_correction_dict[reverse_compliment(CB2read[0:8])])
                    writename=''
                    if name.startswith('@'):
                        writename='%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    else:
                        writename='@%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    truecellname=writeCB1+writeCB2
                    writeCB=writeCB1+writeCB2+umi
                    writeQual=CB1_qual+CB2_qual
                    if truecellname in Cell_statistics:
                        if umi not in Cell_statistics[truecellname][0]:
                            Cell_statistics[truecellname][0].append(umi)
                            Cell_statistics[truecellname][1]+=1
                            Cell_statistics[truecellname][2]+=1
                        else:
                            Cell_statistics[truecellname][2]+=1
                    else:
                        Cell_statistics[truecellname]=[[],0,0]
                        Cell_statistics[truecellname][0].append(umi)
                        Cell_statistics[truecellname][1]+=1
                        Cell_statistics[truecellname][2]+=1
                    filtering_statistics['Valid_read']+=1
                    self._write_fastq(writename,'R1',CB,writeCB,writeQual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(writename,'R1',RNA_read,rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                else:
                    if CB1read not in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) not in Barcode_correction_dict:
                        filtering_statistics['Invalid_Both_CB']+=1
                    else:
                        if CB1read not in Barcode_correction_dict:
                            filtering_statistics['Invalid_CB1']+=1
                        else:
                            filtering_statistics['Invalid_CB2']+=1
            CB.close()
            RNA_read.close()
            #os.system('gzip %s'%('%s/%s_STAR_CB.fastq'%(filtering_output_directory,sample)))
            #os.system('gzip %s'%('%s/%s_STAR_Reads.fastq'%(filtering_output_directory,sample)))
            for cell in Cell_statistics:
                Cell_statistics[cell][0]=0
            with open('%s/%s_filtering_statistics.tsv'%(filtering_output_directory,sample), "w", newline="") as csvfile:
                writer=csv.writer(csvfile,delimiter='\t')
                writer.writerow(['Category','Read Number','Percentage'])
                writer.writerow(['Total_read',filtering_statistics['Total_read'],filtering_statistics['Total_read']/filtering_statistics['Total_read']])
                writer.writerow(['Valid_read',filtering_statistics['Valid_read'],filtering_statistics['Valid_read']/filtering_statistics['Total_read']])
                writer.writerow(['Invalid_CB1',filtering_statistics['Invalid_CB1'],filtering_statistics['Invalid_CB1']/filtering_statistics['Total_read']])
                writer.writerow(['Invalid_CB2',filtering_statistics['Invalid_CB2'],filtering_statistics['Invalid_CB2']/filtering_statistics['Total_read']])
                writer.writerow(['Invalid_Both_CB',filtering_statistics['Invalid_Both_CB'],filtering_statistics['Invalid_Both_CB']/filtering_statistics['Total_read']])
            csvfile.close()
            with open('%s/%s_Cell_statistics.tsv'%(filtering_output_directory,sample), "w", newline="") as csvfile:
                writer=csv.writer(csvfile,delimiter='\t')
                writer.writerow(['Cellname','Number of unique UMI','Number of reads'])
                for cell in Cell_statistics.keys():
                    writer.writerow([cell,Cell_statistics[cell][1],Cell_statistics[cell][2]])
            csvfile.close()            

class inDrop_Individual_Sample_processing(inDrop_Data_processing):
    GEO=True
    def __init__(self,pathtolibraryindex,pathtocellbarcode1,pathtocellbarcode2umi,pathtorna,libraryindex,outputdir,GEO=True):
        try:
            super().__init__(pathtolibraryindex=pathtolibraryindex,
            pathtocellbarcode1=pathtocellbarcode1, 
            pathtocellbarcode2umi=pathtocellbarcode2umi,
            pathtorna=pathtorna,
            libraryindex=libraryindex,
            outputdir=outputdir)
            self.GEO=GEO
        except Exception as e:
            print(str(e))
    def _ParseGEOFastq(self,pathstofastqs):
        if pathstofastqs[0].endswith('.gz'):
            processes=[subprocess.Popen(['zcat',(fastq)],stdout=subprocess.PIPE) for fastq in pathstofastqs]
            totalreads = [r.stdout for r in processes]
        elif pathstofastqs[0].endswith('.bz2'):
            processes=[subprocess.Popen(['bzcat',(fastq)],stdout=subprocess.PIPE) for fastq in pathstofastqs]
            totalreads = [r.stdout for r in processes]
        elif pathstofastqs[0].endswith('.fastq'):
            processes=[subprocess.Popen(['cat',(fastq)],stdout=subprocess.PIPE) for fastq in pathstofastqs]
            totalreads = [r.stdout for r in processes]
        else:
            sys.exit('The format of the file %s is not recognized.'%(str(pathtofastq)))
        while True:
            names=[next(read).decode().split(' ')[1] for read in totalreads]
            Sequence=[next(read).decode() for read in totalreads]
            Blank=[next(read).decode() for read in totalreads]
            qualityscore= [next(read).decode() for read in totalreads]
            assert all(name==names[0] for name in names)
            if names:
                yield [names[0], Sequence, qualityscore]
            else:
                break
        for read in totalreads:
            read.close()    
    def Individual_sample_filter_kallisto(self):
        Barcode_correction_dict=build_barcode_neighborhoods(barcodelist=self.cellbarcodewhitelist)
        filtering_output_directory='%s/%s_kallisto'%(self.outputdir,str(list(self.libraryindex.keys())[0]))
        try:
            os.mkdir(filtering_output_directory)
        except:
            print('Nah, not a problem.')
        filtering_statistics={
        'Total_read':0,
        'Valid_read':0,
        'Invalid_CB1':0,
        'Invalid_CB2':0,
        'Invalid_Both_CB':0}
        sample=str(list(self.libraryindex.keys())[0])
        Cell_statistics={}#Structure: Cell_barcode:[[UMI_list],number of unique umis per cell, number of reads per cell]
        CB1=gzip.open('%s/%s_filtered_CB1.fastq.gz'%(filtering_output_directory,sample),'wt')
        CB2=gzip.open('%s/%s_filtered_CB2.fastq.gz'%(filtering_output_directory,sample),'wt')
        RNA_read=gzip.open('%s/%s_filtered_Reads.fastq.gz'%(filtering_output_directory,sample),'wt')
        if self.GEO is True:
            for read in self._ParseGEOFastq([self.pathtocellbarcode1,self.pathtocellbarcode2umi,self.pathtorna]):
                name=read[0].strip('\n')
                CB1read=read[1][0].strip('\n')
                CB2read=read[1][1].strip('\n')
                rnaread=read[1][2].strip('\n')
                CB1_qual=read[2][0].strip('\n')
                CB2_qual=read[2][1].strip('\n')
                rnaread_qual=read[2][2].strip('\n')
                filtering_statistics['Total_read']+=1
                umi=CB2read[8:]
                if CB1read in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) in Barcode_correction_dict:
                    writeCB1=Barcode_correction_dict[CB1read]
                    writeCB2=reverse_compliment(Barcode_correction_dict[reverse_compliment(CB2read[0:8])])
                    writename=''
                    if name.startswith('@'):
                        writename='%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    else:
                        writename='@%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    truecellname=writeCB1+writeCB2
                    if truecellname in Cell_statistics:
                        if umi not in Cell_statistics[truecellname][0]:
                            Cell_statistics[truecellname][0].append(umi)
                            Cell_statistics[truecellname][1]+=1
                            Cell_statistics[truecellname][2]+=1
                        else:
                            Cell_statistics[truecellname][2]+=1
                    else:
                        Cell_statistics[truecellname]=[[],0,0]
                        Cell_statistics[truecellname][0].append(umi)
                        Cell_statistics[truecellname][1]+=1
                        Cell_statistics[truecellname][2]+=1
                    filtering_statistics['Valid_read']+=1
                    self._write_fastq(writename,'R1',CB1,writeCB1,CB1_qual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(writename,'R1',CB2,writeCB2+umi,CB2_qual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(writename,'R1',RNA_read,rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                else:
                    if CB1read not in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) not in Barcode_correction_dict:
                        filtering_statistics['Invalid_Both_CB']+=1
                    else:
                        if CB1read not in Barcode_correction_dict:
                            filtering_statistics['Invalid_CB1']+=1
                        else:
                            filtering_statistics['Invalid_CB2']+=1
            CB1.close()
            CB2.close()
            RNA_read.close()
        else:
            for read in self._ParseFastq([self.pathtocellbarcode1,self.pathtocellbarcode2umi,self.pathtorna]):
                name=read[0].strip('\n')
                CB1read=read[1][0].strip('\n')
                CB2read=read[1][1].strip('\n')
                rnaread=read[1][2].strip('\n')
                CB1_qual=read[2][0].strip('\n')
                CB2_qual=read[2][1].strip('\n')
                rnaread_qual=read[2][2].strip('\n')
                filtering_statistics['Total_read']+=1
                umi=CB2read[8:]
                if CB1read in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) in Barcode_correction_dict:
                    writeCB1=Barcode_correction_dict[CB1read]
                    writeCB2=reverse_compliment(Barcode_correction_dict[reverse_compliment(CB2read[0:8])])
                    writename=''
                    if name.startswith('@'):
                        writename='%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    else:
                        writename='@%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    truecellname=writeCB1+writeCB2
                    if truecellname in Cell_statistics:
                        if umi not in Cell_statistics[truecellname][0]:
                            Cell_statistics[truecellname][0].append(umi)
                            Cell_statistics[truecellname][1]+=1
                            Cell_statistics[truecellname][2]+=1
                        else:
                            Cell_statistics[truecellname][2]+=1
                    else:
                        Cell_statistics[truecellname]=[[],0,0]
                        Cell_statistics[truecellname][0].append(umi)
                        Cell_statistics[truecellname][1]+=1
                        Cell_statistics[truecellname][2]+=1
                    filtering_statistics['Valid_read']+=1
                    self._write_fastq(writename,'R1',CB1,writeCB1,CB1_qual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(writename,'R1',CB2,writeCB2+umi,CB2_qual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(writename,'R1',RNA_read,rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                else:
                    if CB1read not in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) not in Barcode_correction_dict:
                        filtering_statistics['Invalid_Both_CB']+=1
                    else:
                        if CB1read not in Barcode_correction_dict:
                            filtering_statistics['Invalid_CB1']+=1
                        else:
                            filtering_statistics['Invalid_CB2']+=1
            CB1.close()
            CB2.close()
            RNA_read.close()
        for cell in Cell_statistics:
            Cell_statistics[cell][0]=0
        with open('%s/%s_filtering_statistics.tsv'%(filtering_output_directory,sample), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Category','Read Number','Percentage'])
            writer.writerow(['Total_read',filtering_statistics['Total_read'],filtering_statistics['Total_read']/filtering_statistics['Total_read']])
            writer.writerow(['Valid_read',filtering_statistics['Valid_read'],filtering_statistics['Valid_read']/filtering_statistics['Total_read']])
            writer.writerow(['Invalid_CB1',filtering_statistics['Invalid_CB1'],filtering_statistics['Invalid_CB1']/filtering_statistics['Total_read']])
            writer.writerow(['Invalid_CB2',filtering_statistics['Invalid_CB2'],filtering_statistics['Invalid_CB2']/filtering_statistics['Total_read']])
            writer.writerow(['Invalid_Both_CB',filtering_statistics['Invalid_Both_CB'],filtering_statistics['Invalid_Both_CB']/filtering_statistics['Total_read']])
        csvfile.close()
        with open('%s/%s_Cell_statistics.tsv'%(filtering_output_directory,sample), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Cellname','Number of unique UMI','Number of reads'])
            for cell in Cell_statistics.keys():
                writer.writerow([cell,Cell_statistics[cell][1],Cell_statistics[cell][2]])
        csvfile.close()

    def Individual_sample_filter_solo(self):
        Barcode_correction_dict=build_barcode_neighborhoods(barcodelist=self.cellbarcodewhitelist)
        filtering_output_directory='%s/%s_solo'%(self.outputdir,str(list(self.libraryindex.keys())[0]))
        try:
            os.mkdir(filtering_output_directory)
        except:
            print('Nah, not a problem.')
        filtering_statistics={
        'Total_read':0,
        'Valid_read':0,
        'Invalid_CB1':0,
        'Invalid_CB2':0,
        'Invalid_Both_CB':0}
        sample=str(list(self.libraryindex.keys())[0])
        Cell_statistics={}#Structure: Cell_barcode:[[UMI_list],number of unique umis per cell, number of reads per cell]
        CB=gzip.open('%s/%s_filtered_CB.fastq.gz'%(filtering_output_directory,sample),'wt')
        RNA_read=gzip.open('%s/%s_filtered_Reads.fastq.gz'%(filtering_output_directory,sample),'wt')
        if self.GEO is True:
            for read in self._ParseGEOFastq([self.pathtocellbarcode1,self.pathtocellbarcode2umi,self.pathtorna]):
                name=read[0].strip('\n')
                CB1read=read[1][0].strip('\n')
                CB2read=read[1][1].strip('\n')
                rnaread=read[1][2].strip('\n')
                CB1_qual=read[2][0].strip('\n')
                CB2_qual=read[2][1].strip('\n')
                rnaread_qual=read[2][2].strip('\n')
                filtering_statistics['Total_read']+=1
                umi=CB2read[8:]
                if CB1read in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) in Barcode_correction_dict:
                    writeCB1=Barcode_correction_dict[CB1read]
                    writeCB2=reverse_compliment(Barcode_correction_dict[reverse_compliment(CB2read[0:8])])
                    writeCB=writeCB1+writeCB2+umi
                    writeCBqual=CB1_qual+CB2_qual
                    writename=''
                    if name.startswith('@'):
                        writename='%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    else:
                        writename='@%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    truecellname=writeCB1+writeCB2
                    if truecellname in Cell_statistics:
                        if umi not in Cell_statistics[truecellname][0]:
                            Cell_statistics[truecellname][0].append(umi)
                            Cell_statistics[truecellname][1]+=1
                            Cell_statistics[truecellname][2]+=1
                        else:
                            Cell_statistics[truecellname][2]+=1
                    else:
                        Cell_statistics[truecellname]=[[],0,0]
                        Cell_statistics[truecellname][0].append(umi)
                        Cell_statistics[truecellname][1]+=1
                        Cell_statistics[truecellname][2]+=1
                    filtering_statistics['Valid_read']+=1
                    self._write_fastq(writename,'R1',CB,writeCB,writeCBqual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(writename,'R1',RNA_read,rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                else:
                    if CB1read not in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) not in Barcode_correction_dict:
                        filtering_statistics['Invalid_Both_CB']+=1
                    else:
                        if CB1read not in Barcode_correction_dict:
                            filtering_statistics['Invalid_CB1']+=1
                        else:
                            filtering_statistics['Invalid_CB2']+=1
            CB.close()
            RNA_read.close()
        else:
            for read in self._ParseFastq([self.pathtocellbarcode1,self.pathtocellbarcode2umi,self.pathtorna]):
                name=read[0].strip('\n')
                CB1read=read[1][0].strip('\n')
                CB2read=read[1][1].strip('\n')
                rnaread=read[1][2].strip('\n')
                CB1_qual=read[2][0].strip('\n')
                CB2_qual=read[2][1].strip('\n')
                rnaread_qual=read[2][2].strip('\n')
                filtering_statistics['Total_read']+=1
                umi=CB2read[8:]
                if CB1read in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) in Barcode_correction_dict:
                    writeCB1=Barcode_correction_dict[CB1read]
                    writeCB2=reverse_compliment(Barcode_correction_dict[reverse_compliment(CB2read[0:8])])
                    writeCB=writeCB1+writeCB2+umi
                    writeCBqual=CB1_qual+CB2_qual
                    writename=''
                    if name.startswith('@'):
                        writename='%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    else:
                        writename='@%s %s:%s:%s'%(name.split(' ')[0],writeCB1,writeCB2,umi)
                    truecellname=writeCB1+writeCB2
                    if truecellname in Cell_statistics:
                        if umi not in Cell_statistics[truecellname][0]:
                            Cell_statistics[truecellname][0].append(umi)
                            Cell_statistics[truecellname][1]+=1
                            Cell_statistics[truecellname][2]+=1
                        else:
                            Cell_statistics[truecellname][2]+=1
                    else:
                        Cell_statistics[truecellname]=[[],0,0]
                        Cell_statistics[truecellname][0].append(umi)
                        Cell_statistics[truecellname][1]+=1
                        Cell_statistics[truecellname][2]+=1
                    filtering_statistics['Valid_read']+=1
                    self._write_fastq(writename,'R1',CB,writeCB,writeCBqual)#ID,file_index,file,seq,quality_score
                    self._write_fastq(writename,'R1',RNA_read,rnaread,rnaread_qual)#ID,file_index,file,seq,quality_score
                else:
                    if CB1read not in Barcode_correction_dict and reverse_compliment(CB2read[0:8]) not in Barcode_correction_dict:
                        filtering_statistics['Invalid_Both_CB']+=1
                    else:
                        if CB1read not in Barcode_correction_dict:
                            filtering_statistics['Invalid_CB1']+=1
                        else:
                            filtering_statistics['Invalid_CB2']+=1
            CB.close()
            RNA_read.close()
        for cell in Cell_statistics:
            Cell_statistics[cell][0]=0
        with open('%s/%s_filtering_statistics.tsv'%(filtering_output_directory,sample), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Category','Read Number','Percentage'])
            writer.writerow(['Total_read',filtering_statistics['Total_read'],filtering_statistics['Total_read']/filtering_statistics['Total_read']])
            writer.writerow(['Valid_read',filtering_statistics['Valid_read'],filtering_statistics['Valid_read']/filtering_statistics['Total_read']])
            writer.writerow(['Invalid_CB1',filtering_statistics['Invalid_CB1'],filtering_statistics['Invalid_CB1']/filtering_statistics['Total_read']])
            writer.writerow(['Invalid_CB2',filtering_statistics['Invalid_CB2'],filtering_statistics['Invalid_CB2']/filtering_statistics['Total_read']])
            writer.writerow(['Invalid_Both_CB',filtering_statistics['Invalid_Both_CB'],filtering_statistics['Invalid_Both_CB']/filtering_statistics['Total_read']])
        csvfile.close()
        with open('%s/%s_Cell_statistics.tsv'%(filtering_output_directory,sample), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Cellname','Number of unique UMI','Number of reads'])
            for cell in Cell_statistics.keys():
                writer.writerow([cell,Cell_statistics[cell][1],Cell_statistics[cell][2]])
        csvfile.close()