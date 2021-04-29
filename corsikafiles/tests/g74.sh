#!/bin/bash
/data/user/sstray/oscnext_metaproject/build/env-shell.sh python /home/sstray/test/condor/corsikafiles/tests/g74.py -i $1
mv group_data.hdf5 /data/user/sstray/g74_zbin_2016_new/group74_data$2.hdf5
