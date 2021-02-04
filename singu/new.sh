#!/bin/bash
export HDF5_USE_FILE_LOCKING=FALSE
#python /home/sstray/test/condor/singu/singtest.py -i $1
python /home/sstray/test/condor/singu/ridecalc.py
