import os
import bz2
import subprocess
from inDrop_Data_processing import *
import yaml





if __name__ == '__main__':
    inDrop_Data_processor=inDrop_Data_processing(pathtocellbarcode1='',
                            pathtocellbarcode2umi='',
                            pathtolibraryindex='',
                            pathtorna='',
                            libraryindex='', outputdir='')
    inDrop_Data_processor.Demultiplexing_and_Correcting()

