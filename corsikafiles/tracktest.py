def test(frame):
	T = frame['I3MCTree'].get_daughters(frame['I3MCTree'].get_primaries()[0])
	if len(T) > 1:
#		print('Cluster')
		return False
	if not (str(T[0].type)=='MuPlus' or str(T[0].type)=='MuMinus'):
#		print('Not a muon')
		return False
	mu = T[0]
	endp = np.array(mu.pos+mu.dir*mu.length)
	stopped_xy = border.contains_point((endp[0],endp[1]),radius=-trackradius)
	stopped_z = (endp[2] > (border_zminmax[0]+height)) * (endp[2] < (border_zminmax[1]-height))
	stoppedmuon = int(stopped_xy*stopped_z)

	x,y,z = endp[0],endp[1],endp[2]
	frame['muon_x'] = dataclasses.I3Double(x)
	frame['muon_y'] = dataclasses.I3Double(y)
	frame['muon_z'] = dataclasses.I3Double(z)
	frame['is_stopped'] = dataclasses.I3Double(stoppedmuon)

def savefeatures(frame,fin=None,names=None):
	for i,name in enumerate(names):
		fin[i].append(frame[name].value)


import os
from I3Tray import *
from icecube import dataio
from icecube import dataclasses
import numpy as np
import argparse
import matplotlib.path as mpath
bordercoords = np.array([(-256.1400146484375, -521.0800170898438), (-132.8000030517578, -501.45001220703125), (-9.13000011444092, -481.739990234375), (114.38999938964844, -461.989990234375), (237.77999877929688, -442.4200134277344), (361.0, -422.8299865722656), (405.8299865722656, -306.3800048828125), (443.6000061035156, -194.16000366210938), (500.42999267578125, -58.45000076293945), (544.0700073242188, 55.88999938964844), (576.3699951171875, 170.9199981689453), (505.2699890136719, 257.8800048828125), (429.760009765625, 351.0199890136719), (338.44000244140625, 463.7200012207031), (224.5800018310547, 432.3500061035156), (101.04000091552734, 412.7900085449219), (22.11000061035156, 509.5), (-101.05999755859375, 490.2200012207031), (-224.08999633789062, 470.8599853515625), (-347.8800048828125, 451.5199890136719), (-392.3800048828125, 334.239990234375), (-437.0400085449219, 217.8000030517578), (-481.6000061035156, 101.38999938964844), (-526.6300048828125, -15.60000038146973), (-570.9000244140625, -125.13999938964844), (-492.42999267578125, -230.16000366210938), (-413.4599914550781, -327.2699890136719), (-334.79998779296875, -424.5)])
border = mpath.Path(bordercoords)
border_zminmax = [-512.82,524.56]
trackradius = 100
height = 100

domfile = np.loadtxt('/home/sstray/test/condor/corsikafiles/dom_coords_spacing_HQE.txt',dtype='str')
testfilename = '/data/sim/IceCube/2012/filtered/level2/CORSIKA-in-ice/11499/77000-77999/Level2_IC86.2012_corsika.011499.077910.i3.bz2'

	
def test_my_little_function():
	my_features = ['muon_x','muon_y','muon_z','is_stopped']
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--infiles',required=True)
	args = parser.parse_args()
	filelist = args.infiles.split(',')
	for nums in range(len(filelist)):
		if nums == (len(filelist)-1):
			filelist[nums] = filelist[nums][:-1]
		featurestemp = [[] for i in range(len(my_features))]
		global event_run
		event_run = filelist[nums].split('/')[-1].split('.')[-3]

		tray = I3Tray()
		tray.AddModule('I3Reader','read_stuff',Filename=filelist[nums])
		tray.AddModule(test,'test')
		tray.AddModule(savefeatures,'asd',fin=featurestemp,names=my_features)
		tray.Execute()
		if nums == 0:
			featuresarray = np.column_stack(np.array(featurestemp))
		else:
			featuresarray = np.concatenate((featuresarray,np.column_stack(featurestemp)),axis=0)
	return featuresarray

#pr = cProfile.Profile()
#pr.enable()

thearray = test_my_little_function()

#pr.disable()
#s = StringIO.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr,stream=s).sort_stats(sortby)
#ps.print_stats()
#print(s.getvalue())


h5file = 1
if h5file == 1:
	import tables
	f = tables.open_file('muon_track_test.hdf5', 'w',filters=tables.Filters(complib='zlib', complevel=6))
	f.create_carray('/','array',obj=thearray)
	f.close()
