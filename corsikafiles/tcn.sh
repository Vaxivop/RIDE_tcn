#!/bin/bash
/data/user/sstray/oscnext_metaproject/build/env-shell.sh python /home/sstray/test/condor/corsikafiles/tcncorsika.py -i $1
mv tcn_corsika_data.hdf5 /data/user/sstray/2016_nopri/tcn_corsika_data$2.hdf5
