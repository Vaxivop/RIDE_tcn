import glob
import pickle
actualfiles = glob.glob('/data/user/sstray/forreal/*.hdf5')
print('Total  number of files:',len(actualfiles))
filenum = 1
print('Files processed per script:',filenum)
if len(actualfiles)%filenum != 0:
	print('Warning: File number is not a factor of number of files')
ls = [actualfiles[x:x+filenum] for x in range(0,len(actualfiles),filenum)]
print('Number of scripts:',len(ls))

import csv
with open('files_hdf5.txt', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(ls)
