import glob
from random import shuffle
#actualfiles = glob.glob('/data/sim/IceCube/2012/filtered/level2/CORSIKA-in-ice/11499/**/*.i3.bz2')
actual_files = glob.glob('/data/sim/IceCube/2016/generated/CORSIKA-in-ice/21269/IC86_2016_spe_templates_DOM_oversize1/level2/redo/eff100/0000000-0000999/Level2_eff100_IC86-2016_corsika_21269_00000.i3.zst')
#actualfiles = glob.glob('/data/sim/IceCube/2016/generated/CORSIKA-in-ice/21269/IC86_2016_spe_templates_DOM_oversize1/level2/redo/eff100/**/*.i3.zst')
shuffle(actualfiles)
print('Total  number of files:',len(actualfiles))
files_processed = 100000
actualfiles = actualfiles[:files_processed]
print('Number of files being processed:',len(actualfiles))
filenum = 100
print('Files processed per script:',filenum)
if len(actualfiles)%filenum != 0:
	print('Warning: File number is not a factor of number of files')
ls = [actualfiles[x:x+filenum] for x in range(0,len(actualfiles),filenum)]
print('Number of scripts:',len(ls))
import csv
with open('files_corsika.txt', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(ls)
