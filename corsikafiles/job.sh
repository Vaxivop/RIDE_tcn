#!/bin/bash
/data/user/sstray/oscnext_metaproject/build/env-shell.sh python /home/sstray/test/condor/corsikafiles/corsica.py -i $1
#/opt/icetray/combo/build/env-shell.sh python /home/sstray/test/condor/corsikafiles/corsica.py -i $1
mv featurelist.hdf5 /data/user/sstray/no_theta/featurelist$2.hdf5
