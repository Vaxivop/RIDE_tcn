import glob
import pickle
actualfiles = glob.glob('/data/sim/IceCube/2012/filtered/level2/CORSIKA-in-ice/11499/77000-77999/*.i3.bz2')
ls = []
for i in range(len(actualfiles)/10):
	ls.append(actualfiles[int((i)*10):int((i+1)*10)])

import csv
with open('myfile.txt', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(ls)
