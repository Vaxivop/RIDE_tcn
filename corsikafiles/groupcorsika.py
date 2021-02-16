def test(frame):
	if frame['I3EventHeader'].sub_event_id != 0:
		return False
	T = frame['I3MCTree'].get_daughters(frame['I3MCTree'].get_primaries()[0])
	if len(T) > 1:
#		print('Cluster')
		return False
	if not (str(T[0].type)=='MuPlus' or str(T[0].type)=='MuMinus'):
#		print('Not a muon')
		return False

	mu = T[0]	
	endp = np.array(mu.pos+mu.dir*mu.length)
	startp = np.array(mu.pos)
	end_x,end_y,end_z = endp[0],endp[1],endp[2]
	start_x,start_y,start_z = startp[0],startp[1],startp[2]
	theta = mu.dir.zenith*180/np.pi
	
	stopped_xy = border.contains_point((endp[0],endp[1]),radius=-100)
	stopped_z = endp[2] > -400
	stopped_theta = (theta >= 40) * (theta <= 70)
	stoppedmuon = int(stopped_xy*stopped_z)
	if stoppedmuon == 0:
		return False

	dom = []
	string = []
	charge = []
	
	series = frame['SplitInIcePulses'].apply(frame)
	for i,pulse in enumerate(series):
		for hit in pulse[1]:
			string.append(pulse[0].string)
			dom.append(pulse[0].om)
			charge.append(hit.charge)
	if len(charge) < 8:
		return False
	domstr = np.column_stack([string,dom]).astype(str)
	domstr = np.char.zfill(domstr,2).astype(object)
	domstr = np.sum(domstr,axis=1).astype(str)
	domindex = np.array([np.where(domfile[:,3]==i)[0] for i in domstr])[:,0]
	
	x,y,z = domfile[:,5],domfile[:,6],domfile[:,7]
	chargep = np.array(charge).astype(float)
	charge = np.zeros(len(x))
	charge[domindex] = chargep
	charge = charge[grp].astype(float)

#	print('Accepted')
	
	coords = np.array([x,y,z]).T.astype(float)
	point_vector = coords-startp
	seg_vector = endp-startp
	t = np.dot(point_vector,seg_vector)/sum(seg_vector**2)
	close = startp + t[:,np.newaxis]*seg_vector
	close[t<0] = startp
	close[t>1] = endp
	domrad = np.linalg.norm(coords-close,axis=1)

	x = x[grp].astype(float)
	y = y[grp].astype(float)
	z = z[grp].astype(float)
	domrad = domrad[grp].astype(float)
	domstr = np.array(domfile[grp,3]).astype(float)
	HQE = domfile[grp,2].astype(str)
	HQE[HQE=='NQE'] = 0
	HQE[HQE=='HQE'] = 1
	HQE[HQE=='BAD'] = 2
	HQE=np.array(HQE).astype(float)

	tempid = event_run+str(frame['I3EventHeader'].event_id)
	eventid = np.ones(len(x)).astype(int)*int(tempid)

	frame['x'] = dataclasses.I3VectorDouble(x)
	frame['y'] = dataclasses.I3VectorDouble(y)
	frame['z'] = dataclasses.I3VectorDouble(z)
	frame['charge'] = dataclasses.I3VectorDouble(charge)
	frame['domrad'] = dataclasses.I3VectorDouble(domrad)
	frame['domstr'] = dataclasses.I3VectorDouble(domstr)
	frame['status'] = dataclasses.I3VectorDouble(HQE)
	frame['eventid'] = dataclasses.I3VectorDouble(eventid)

	frame['end_x'] = dataclasses.I3Double(end_x)
	frame['end_y'] = dataclasses.I3Double(end_y)
	frame['end_z'] = dataclasses.I3Double(end_z)
	frame['start_x'] = dataclasses.I3Double(start_x)
	frame['start_y'] = dataclasses.I3Double(start_y)
	frame['start_z'] = dataclasses.I3Double(start_z)
	frame['theta'] = dataclasses.I3Double(theta)
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
border_zminmax = [-512.82,524.56]


domfile = np.loadtxt('/home/sstray/test/condor/corsikafiles/dom_coords_spacing_HQE.txt',dtype='str')
outputname = 'group_data.hdf5'

strings_to_cut = [1,2,3,4,5,6,13,21,30,40,50,59,67,74,73,72,78,77,76,75,68,60,51,41,31,22,14,7]
dom_strings = np.array([i[:2] for i in domfile[:,3]]).astype(int)
good_doms = np.array(domfile[:,2]!='BAD')
not_outer = ~np.isin(dom_strings,strings_to_cut)
single_group = np.array(domfile[:,4]).astype(int) == 83
grp = single_group * not_outer * good_doms
maxsize = sum(grp)
print(maxsize)

def test_my_little_function():
	my_features = ['x','y','z','charge','domrad','domstr','status','eventid']
	my_truth = ['start_x','start_y','start_z','end_x','end_y','end_z','theta','eventid']
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--infiles',required=True)
	args = parser.parse_args()
	filelist = args.infiles.split(',')
	first_file_bad = 1
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
			if (nums == 0) or (first_file_bad == 1):
				print('Grouping')
				group,truth = get_groups(featuresarray[:,:-1],trutharray,featuresarray[:,-1])
				first_file_bad = 0
			else:
				print('Grouping and concatenating')
				tmp = get_groups(featuresarray[:,:-1],trutharray,featuresarray[:,-1])
				group = np.concatenate((group,np.copy(tmp[0])),axis=0)
				truth = np.concatenate((truth,np.copy(tmp[1])),axis=0)
		else:
			print('BAD FILE: '+str(filelist[nums]))
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
h5file = 1
if h5file == 1:
	import tables
	f = tables.open_file(outputname, 'w',filters=tables.Filters(complib='zlib', complevel=6))
	f.create_carray('/','features',obj=group)
	f.create_carray('/','truth',obj=truth)
	f.close()
