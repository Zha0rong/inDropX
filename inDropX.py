import os
import bz2
import subprocess
from functions import read_function







if __name__ == '__main__':
    a=read_function.inDrop_Data_processing(pathtocellbarcode1='',
                            pathtocellbarcode2umi='',
                            pathtolibraryindex='',
                            pathtorna='',
                            libraryindex=sample_metadata, outputdir='')
    a.Demultiplex()
    a.correct_and_filter()

