#!/bin/bash
/data/user/sstray/oscnext_metaproject/build/env-shell.sh python /home/sstray/test/condor/corsikafiles/tests/g79.py -i $1
mv group_data.hdf5 /data/user/sstray/g79_zbin_2016_new/group79_data$2.hdf5
