def test(frame):
	T = frame['I3MCTree'].get_daughters(frame['I3MCTree'].get_primaries()[0])
	if len(T) > 1:
		return False
	if not (str(T[0].type)=='MuPlus' or str(T[0].type)=='MuMinus'):
		return False
	dom = []
	string = []
	charge = []
	time = []
	pulse_width = []
	
	series = frame['SplitInIcePulses'].apply(frame)
	for i,pulse in enumerate(series):
		for hit in pulse[1]:
			string.append(pulse[0].string)
			dom.append(pulse[0].om)
			charge.append(hit.charge)
			time.append(hit.time)
			pulse_width.append(hit.width)

	domstr = np.column_stack([string,dom]).astype(str)
	domstr = np.char.zfill(domstr,2).astype(object)
	domstr = np.sum(domstr,axis=1).astype(str)
	domindex = np.array([np.where(domfile[:,3]==i)[0] for i in domstr])[:,0]
	
	x,y,z = domfile[domindex,5].astype(float),domfile[domindex,6].astype(float),domfile[domindex,7].astype(float)

	charge = np.array(charge).astype(float)
	time = np.array(time).astype(float)

	mu = T[0]	
	endp = np.array(mu.pos+mu.dir*mu.length)
        startp = np.array(mu.pos+mu.dir*(mu.length-200))
	stopped_xy = border.contains_point((endp[0],endp[1]),radius=-rad)
	stopped_z = (endp[2] > (border_zminmax[0]+height)) * (endp[2] < (border_zminmax[1]-height))
	stoppedmuon = int(stopped_xy*stopped_z)
	end_x,end_y,end_z = endp[0],endp[1],endp[2]
	start_x,start_y,start_z = startp[0],startp[1],startp[2]
	
	if len(x) == 0:
		print('No length')
		return False
	if frame['I3EventHeader'].sub_event_id != 0:
		return False
	tempid = event_run+str(frame['I3EventHeader'].event_id)
	eventid = np.ones(len(x)).astype(int)*int(tempid)

	frame['x'] = dataclasses.I3VectorDouble(x)
	frame['y'] = dataclasses.I3VectorDouble(y)
	frame['z'] = dataclasses.I3VectorDouble(z)
	frame['charge'] = dataclasses.I3VectorDouble(charge)
	frame['time'] = dataclasses.I3VectorDouble(time)
	frame['eventid'] = dataclasses.I3VectorDouble(eventid)

	frame['stopped_muon'] = dataclasses.I3Double(stoppedmuon)
	frame['end_x'] = dataclasses.I3Double(end_x)
	frame['end_y'] = dataclasses.I3Double(end_y)
	frame['end_z'] = dataclasses.I3Double(end_z)
	frame['start_x'] = dataclasses.I3Double(start_x)
	frame['start_y'] = dataclasses.I3Double(start_y)
	frame['start_z'] = dataclasses.I3Double(start_z)
def savefeatures(frame,fin=None,names=None):
	for i,name in enumerate(names):
		fin[i].append(frame[name])
def savetruth(frame,fin=None,names=None):
	for i,name in enumerate(names):
		if name == 'eventid':
			fin[i].append(frame[name][0])
		else:
			fin[i].append(frame[name].value)
def get_groups(features_in,truth_in,event_in):
	features = features_in
	eventnum = event_in
	unev = np.unique(eventnum)
	grouped_features = ([features[np.isin(eventnum,i)] for i in np.unique(eventnum)])
	if isinstance(truth_in,int) == False:
		truth = truth_in[np.isin(truth_in[:,-1],unev)]
		truth = np.array(truth)[truth[:,-1].argsort(),:]
		truth = np.array(truth[np.array([(len(i)<maxsize and len(i)>10) for i in grouped_features])])
	grouped_features = [np.vstack((np.array(i),np.zeros((maxsize,features.shape[1]))))[:maxsize] for i in grouped_features if len(i)<maxsize if len(i)>10]
	grouped_features = np.dstack(grouped_features)
	grouped_features = np.rollaxis(grouped_features,2,0)
	if isinstance(truth_in,int) == False:
		return grouped_features, truth
	return grouped_features
import os
from I3Tray import *
from icecube import dataio
from icecube import dataclasses
import numpy as np
import argparse

import matplotlib.path as mpath
bordercoords = np.array([(-256.1400146484375, -521.0800170898438), (-132.8000030517578, -501.45001220703125), (-9.13000011444092, -481.739990234375), (114.38999938964844, -461.989990234375), (237.77999877929688, -442.4200134277344), (361.0, -422.8299865722656), (405.8299865722656, -306.3800048828125), (443.6000061035156, -194.16000366210938), (500.42999267578125, -58.45000076293945), (544.0700073242188, 55.88999938964844), (576.3699951171875, 170.9199981689453), (505.2699890136719, 257.8800048828125), (429.760009765625, 351.0199890136719), (338.44000244140625, 463.7200012207031), (224.5800018310547, 432.3500061035156), (101.04000091552734, 412.7900085449219), (22.11000061035156, 509.5), (-101.05999755859375, 490.2200012207031), (-224.08999633789062, 470.8599853515625), (-347.8800048828125, 451.5199890136719), (-392.3800048828125, 334.239990234375), (-437.0400085449219, 217.8000030517578), (-481.6000061035156, 101.38999938964844), (-526.6300048828125, -15.60000038146973), (-570.9000244140625, -125.13999938964844), (-492.42999267578125, -230.16000366210938), (-413.4599914550781, -327.2699890136719), (-334.79998779296875, -424.5)])
border = mpath.Path(bordercoords)
from scipy.interpolate import interp1d
border_zminmax = [-512.82,524.56]
rad = 100
height = 100
maxsize = 100

domfile = np.loadtxt('/home/sstray/test/condor/corsikafiles/dom_coords_spacing_HQE.txt',dtype='str')
testfilename = '/data/sim/IceCube/2012/filtered/level2/CORSIKA-in-ice/11499/77000-77999/Level2_IC86.2012_corsika.011499.077910.i3.bz2'
outputname = 'tcn_corsika_data.hdf5'
meanstd = np.loadtxt('/home/sstray/test/condor/singu/meanstd.txt')
	
def test_my_little_function():
	my_features = ['x','y','z','charge','time','eventid']
	my_truth = ['stopped_muon','start_x','start_y','start_z','end_x','end_y','end_z','eventid']
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--infiles',required=True)
	args = parser.parse_args()
	filelist = args.infiles.split(',')
	for nums in range(len(filelist)):
#		if nums == (len(filelist)-1):
#			filelist[nums] = filelist[nums][:-1]
		featurestemp = [[] for i in range(len(my_features))]
		truthtemp = [[] for i in range(len(my_truth))]
		global event_run
		event_run = filelist[nums].split('/')[-1].split('.')[-3]
		
		tray = I3Tray()
		tray.AddModule('I3Reader','read_stuff',Filename=filelist[nums])
		tray.AddModule(test,'test')
		tray.AddModule(savefeatures,'asd',fin=featurestemp,names=my_features)
		tray.AddModule(savetruth,'asd2',fin=truthtemp, names=my_truth)
		tray.Execute()
		for i in range(len(featurestemp)):
			 featurestemp[i] = [item for sublist in featurestemp[i] for item in sublist]
		featuresarray = np.column_stack(featurestemp)
		trutharray = np.column_stack(truthtemp)
		trutharray = trutharray[np.unique(trutharray[:,-1],return_index=True)[1]]
		if len(trutharray) > 10:
			if nums == 0:
				print('Grouping')
				group,truth = get_groups(featuresarray[:,:-1],trutharray,featuresarray[:,-1])
			if nums != 0:
				print('Grouping and concatenating')
				tmp = get_groups(featuresarray[:,:-1],trutharray,featuresarray[:,-1])
				group = np.concatenate((group,np.copy(tmp[0])),axis=0)
				truth = np.concatenate((truth,np.copy(tmp[1])),axis=0)
		else:
			print(str(filelist[nums])+' is a bad file')
	return group, truth

#pr = cProfile.Profile()
#pr.enable()
group,truth = test_my_little_function()

#pr.disable()
#s = StringIO.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr,stream=s).sort_stats(sortby)
#ps.print_stats()
#print(s.getvalue())

print(group.shape)
print(truth.shape)
print(np.sum(truth[:,0]))
h5file = 1
if h5file == 1:
	import tables
	f = tables.open_file(outputname, 'w',filters=tables.Filters(complib='zlib', complevel=6))
	f.create_carray('/','features',obj=group)
	f.create_carray('/','truth',obj=truth)
	f.close()
