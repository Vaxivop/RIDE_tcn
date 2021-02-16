#!/bin/bash
/home/jvmead/dom_eff/build/env-shell.sh python /home/jvmead/dom_eff/RIDE_tcn/corsikafiles/corsica.py -i $1
mv featurelist.hdf5 /home/jvmead/dom_eff/RIDE_tcn/test/featurelist$2.hdf5
