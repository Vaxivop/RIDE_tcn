#!/bin/bash
eval $(/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh)
python test.py -i $1
