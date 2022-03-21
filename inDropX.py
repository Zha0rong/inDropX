#!/usr/bin/env python

import os
import bz2
import subprocess
from inDrop_Data_processing import *
import yaml
import argparse




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('YAML', type=argparse.FileType('r'), help='Specify parameters for InDropX Using YAML file.',default=None)
    args = parser.parse_args()
    parameters=yaml.safe_load(args.YAML)
    if parameters['Version']=='V3':
        pathtocellbarcode1=[]
        pathtocellbarcode2umi=[]
        pathtolibraryindex=[]
        pathtorna=[]
        allfiles=[]
        for i in range(len(parameters['fastq_location'])):
            pathtorna.append(parameters['fastq_location'][i][0])
            pathtocellbarcode1.append(parameters['fastq_location'][i][1])
            pathtocellbarcode2umi.append(parameters['fastq_location'][i][2])
            pathtolibraryindex.append(parameters['fastq_location'][i][3])
            allfiles.append(parameters['fastq_location'][i][0])
            allfiles.append(parameters['fastq_location'][i][1])
            allfiles.append(parameters['fastq_location'][i][2])
            allfiles.append(parameters['fastq_location'][i][3])

        if all(list(map(os.path.isfile,allfiles))) is False:
            sys.exit('Cannot find fastq file(s) at the directory specific in the YAML file.')
        CHECK_FOLDER = os.path.isdir(parameters['output_directory'])

        # If folder doesn't exist, then create it.
        if not CHECK_FOLDER:
            os.makedirs(parameters['output_directory'])
            print("created folder : ", parameters['output_directory'])

        else:
            print(parameters['output_directory'], "folder already exists.")
        inDrop_Data_processor=inDrop_Data_processing(pathtocellbarcode1=pathtocellbarcode1,
                            pathtocellbarcode2umi=pathtocellbarcode2umi,
                            pathtolibraryindex=pathtolibraryindex,
                            pathtorna=pathtorna,
                            libraryindex=parameters['Library_Index'], outputdir=parameters['output_directory'],version=parameters['Version'])
        inDrop_Data_processor.Demultiplexing_and_Correcting()
    elif parameters['Version']=='V2':
        pathtocellbarcodeumi=[]
        pathtorna=[]
        allfiles=[]
        for i in range(len(parameters['fastq_location'])):
            pathtorna.append(parameters['fastq_location'][i][0])
            pathtocellbarcodeumi.append(parameters['fastq_location'][i][1])
            allfiles.append(parameters['fastq_location'][i][0])
            allfiles.append(parameters['fastq_location'][i][1])

        if all(list(map(os.path.isfile,allfiles))) is False:
            sys.exit('Cannot find fastq file(s) at the directory specific in the YAML file.')
        CHECK_FOLDER = os.path.isdir(parameters['output_directory'])
        # If folder doesn't exist, then create it.
        if not CHECK_FOLDER:
            os.makedirs(parameters['output_directory'])
            print("created folder : ", parameters['output_directory'])

        else:
            print(parameters['output_directory'], "folder already exists.")
        inDrop_Data_processor=inDrop_Data_processing(pathtocellbarcode1=pathtocellbarcodeumi,
                            pathtocellbarcode2umi=[],
                            pathtolibraryindex=[],
                            pathtorna=pathtorna,
                            libraryindex=parameters['Library_Index'], outputdir=parameters['output_directory'],version=parameters['Version'])
        #inDrop_Data_processor.Demultiplexing_and_Correcting()


