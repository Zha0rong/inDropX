# InDropX

## What is InDropX?
InDropX is a pipeline built to convert the inDrop data format from CB1 + CB2/UMI + RNA to CB/UMI + RNA so that it can be analyzed by different recent developed bioinformatics analysis tools such as [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) or [kallistobus](https://www.kallistobus.tools).

It can demultiplex samples from the 4 ["Undetermined" fastq files from bcl2fastq](https://github.com/indrops/indrops) and correct cell barcodes using the whitelist provided by the original [InDrop pipeline](https://github.com/indrops/indrops).


## How do I use it?
Make sure Python is up and running in the system.
Install [PYYAML](https://pyyaml.org)
Specify location of fastq files, 8 base pairs long library index and output directory in the .yaml file.
Run:
```

python inDropX.py inDropX.parameters.yaml

```


The inDropX.parameters.yaml is the directory of inDropX.parameters.yaml, so if the supercomplicated.parameters.yaml is stored in excessively/sophisticated/directory, the command would be:

```

python inDropX.py excessively/sophisticated/directory/supercomplicated.parameters.yaml

```

### How do I specify location of fastq files, 8 base pairs long library index and output directory in the yaml file.
In the .yaml file three fields are required:
1. fastq_location

   In this field the entry should be at least one list with 4 fastq files in the following orders:
   [Biological Reads, Cell Barcode part 1, Cell Barcode part 2 and UMI, Library Index]. In the [indrops pipeline](https://github.com/indrops/indrops) the authors explained that if bcl2fastq is used to convert Illumina BCL files to fastq files: 
       > v3 : summer 2016 redesign requiring manual demultiplexing. R1 is the biological read. R2 carries the first half of the gel barcode, R3 carries the library index and R4 the second half of the gel barcode, the UMI and a fraction of the polyA tail.
2. Library_Index

   In this field the entry should be a dictionary with at least one key:value. The key should be the sample name and the alue is the 8 BP long library index sequence.
3. output_directory

   In this field the entry should be a directory where the output files needs to be stored. Storage space should be twice the size of the original fastq files.

## What does the output look like?

The output of the pipeline will be:
1. A tab delimited .tsv file called 'Read_statistics_for_demultiplexing.tsv', which shows the percentage of reads in every sample in the Multiplexed dataset.
2. One directory created for each sample, in which there are:
   1. [SampleName].filtering_statistics.tsv. This tab-delimited .tsv file shows how many reads are retained after cell barcode error correction.
   2. [SampleName].cell_statistics.tsv.  This tab-delimited .tsv file shows how many UMI and reads are detected in each cell in the sample.
   3. [SampleName].barcode1.fastq.gz. The unfiltered fastq file of part 1 of cell barcode. Intend for record keeping and upload to public repository.
   4. [SampleName].barcode2.fastq.gz. The unfiltered fastq file of part 3 of cell barcode and UMI. Intend for record keeping and upload to public repository.
   5. [SampleName].read.fastq.gz. The unfiltered fastq file of biological read. Intend for record keeping and upload to public repository.
   6. [SampleName].filtered.barcodes.umi.fastq.gz. The filtered and corrected fastq file of barcodes and UMI, first 16 bases of sequence are corrected cell barcode and next 6 bases are UMI.
   7. [SampleName].filtered.read.fastq.gz. The filtered and corrected fastq file of Biological Reads.

## How should I utilize the filtered results?

### Use STARsolo to align and quantify sample.

The pipeline is designed to anneal the two-parts cell barcodes and UMI together so it can be analyzed using [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) or [kallistobus](https://www.kallistobus.tools).
An example of [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) analysis command I use for the filtered results are the following:
```
STAR --soloType Droplet --outSAMtype BAM SortedByCoordinate \
--soloFeatures Gene GeneFull SJ\
--outSAMattributes NH HI AS nM CR CY UR UY CB UB \
--runThreadN $numberofthread --twopassMode Basic --sjdbGTFfile $gtf_location \
--genomeDir $reference_genome_location --soloUMIlen 6 \
--soloCBwhitelist whitelist.txt --soloCBmatchWLtype Exact \
--soloUMIdedup Exact --readFilesCommand zcat \
--soloUMIfiltering MultiGeneUMI --outFileNamePrefix $SampleName/$SampleName. \
--readFilesIn $SampleName.filtered.read.fastq.gz $SampleName.filtered.barcodes.umi.fastq.gz
```

### Use kallisto and bustools to align and quantify sample.

Kallisto and Bustools pipeline has been integrated together by the pipeline [kb](https://www.kallistobus.tools/kb_usage/kb_count/). 

Just like STARsolo mentioned previously, this pipeline allows user to specify the location and property of Cell barcodes and UMI using -x ("technology"). In the description of Kallisto, it says:
> Additionally kallisto bus will accept a string specifying a new technology in the format of bc:umi:seq where each of bc,umi and seq are a triplet of integers separated by a comma, 
> denoting the file index, start and stop of the sequence used. 
> For example to specify the 10xV2 technology we would use 0,0,16:0,16,26:1,0,0. 
> The first part bc is 0,0,16 indicating it is in the 0-th file (also known as the first file in plain english), 
> the barcode starts at the 0-th bp and ends at the 16-th bp in the sequence (i.e. 16bp barcode), 
> the UMI is similarly in the same file, right after the barcode in position 16-26 (a 10bp UMI), 
> finally the sequence is in a separate file, starts at 0 and ends at 0 (in this case stopping at 0 means there is no limit, we use the entire sequence).

In the case of the filtered.



