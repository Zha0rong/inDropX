import os
import bz2
import subprocess
import read_function
import project_recorder





if __name__ == '__main__':
    a=read_function.inDrop_Data_processing(pathtocellbarcode1='',
                            pathtocellbarcode2umi='',
                            pathtolibraryindex='',
                            pathtorna='',
                            libraryindex='', outputdir='')
    a.Demultiplex()
    a.correct_and_filter()

