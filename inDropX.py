import os
import bz2
import subprocess
import read_function





if __name__ == '__main__':
    a=read_function.inDrop_Data_processing(pathtocellbarcode1='',
                            pathtocellbarcode2umi='',
                            pathtolibraryindex='',
                            pathtorna='',
                            libraryindex='', outputdir='')
    a.Demultiplexing_and_Correcting()

