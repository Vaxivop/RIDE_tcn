#!/bin/bash
/data/user/sstray/oscnext_metaproject/build/env-shell.sh python /home/sstray/test/condor/corsikafiles/tests/g83.py -i $1
mv group_data.hdf5 /data/user/sstray/group83/group83_data$2.hdf5
