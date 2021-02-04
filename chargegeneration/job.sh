#!/bin/bash
eval $(/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh)
source /home/sstray/scratch/RIDE/setup_RIDE_env.sh
/data/user/sstray/oscnext_metaproject/build/env-shell.sh python ./RIDE/processing/standard_simulation/RIDE_processing.py -i $1 -o $2 --cuts
mv expcharge.hdf5 expcharge$3.hdf5
mv realcharge.hdf5 realcharge$3.hdf5
