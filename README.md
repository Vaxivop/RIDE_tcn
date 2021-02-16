# Sofus' Relative Individual DOM Efficiency

Welcome. Here's a list of the important parts:

## How to run the main RIDE scripts

The main script used to produce hdf5 tables from a list of i3 files is ```corsica.py```. I will update it with more options later, but it is currently set to the following cuts:
* Any muon with endpoint not 100 meters inside the detector is cut
* Any muon with endpoint below -400 m is cut
* Any muon not within a zenith angle of 40 to 70 degrees is cut
* Any muon with less than 8 DOM hits is cut (SMT8 filter)
* Aditionally, all DOMs that are not between 20 and 100 m of the track are cut
* Any muon with zero DOMs with 20 m < r < 100 m is cut

The code includes a lot of options for TCN and different cuts, but that requires manual code modification and isn't included in arguments just yet.

To run the code you need to have the entire IceTray environment somewhere. You can run the code two ways: The first is mainly for testing that it works and is done via
```bash
python corsica.py -i /data/sim/IceCube/2012/filtered/level2/CORSIKA-in-ice/11499/52000-52999/Level2_IC86.2012_corsika.011499.052332.i3.bz2
```
The file in this instance is arbitrary. The script supports multiple files using
```bash
python corsica.py -i filepath1,filepath2,filepath3
```
(requires you to be in the IceTray environment) and so on, or by passing a .txt file containing a list of comma separated filespaths. If you want to pass hundreds or thousands of files, however, it is best to use a cluster. First you define how many scripts will process how many files in ```file_generation_corsika.py```. The current files, ```files_corsika.txt```, are set up for 100 scripts to run 100 files each. This will take a few hours on the cluster to process. You can modify this by changing the ```files_processed``` and/or ```filenum``` lines. Running the script tells you have many files and scripts will be run.

Next, you define where the files should be moved to. The submission script job.sh is currently set to use my own IceTray environment, but you can just substitute your own by editing the file. Because the /home/ directory doesn't handle large files well, the last line moves the hdf5 tables to a /data/user/ folder. Substitute the filepath with another folder you make in your own data folder. You can change the name too provided it still ends with ```$2.hdf5```.

After confirming the files and the script, you simply log onto the condor submitter (```ssh submitter```) and run ```condor_submit job.sub```. It is a bit more involved to run it with TCN at the moment, so don't try that just yet. You can see the progress using ```condor_q```. Check job.err or job.out if the job for some reason fails.


### How to run the single group script.

Exact same as above except the main script is now ```groupcorsika.py```. It uses group 83 as its default but you can change it by editing the script. The rest is the same. You still need to edit the job.sh file and so on, but it uses the same files.

## Converting HDF5 tables to the "RIDE file"

Edit the ```ridecalc.py``` script in the /singu/ folder so the line ```actualfiles = glob.glob('/filepath/*.hdf5')```points towards your folder. Optionally edit the line at the very bottom ```np.savetxt('RIDE_insertname.txt',RIDE)``` to whatever name you want.

NOTE: This script uses h5py which the default IceTray environment doesn't have access to. The best way around this that I know of is creating a singularity from dockerhub that has the hdf5 package. To do this, first build the singularity using
```bash
singularity build hfpy.sif docker://hdfgroup/h5py:2.7.0
```
and then run the script using
```bash
singularity run -B /filepath hfpy.sif python ridecalc.py
```
where /filepath is whatever folder you saved the hdf5 tables in. You only need to build the singularity once, but you have to run it every time you need to run the script.

You now have the RIDE.txt file. The best way to actually look through the data is opening up your favourite python editor and running
```python
import numpy as np
import pandas as pd
newcolumns = ('RIDE','mean_charge','total_charge','expected_hits','actual_hits','x','y','z','domstr','group','status')
RIDE = pd.DataFrame(np.loadtxt('RIDE.txt'),columns=newcolumns)
```
This gives you a nicely formatted pandas dataframe.

### Combining the files for a single group

If you don't want the RIDE.txt file but want all the event data in full, the procedure is the same as above except you need to run the ```combinefiles.py``` script instead. Be sure to edit the name for the input files first. Note that you *can* run this on more than a single group, but the script will most likely crash due to the filesize of the resulting hdf5 table unless you apply some cuts on the input hdf5 files.

## What are all these other files then?

A bunch of leftover garbage.
