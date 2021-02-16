def test(frame):

	if frame['I3EventHeader'].sub_event_stream != 'InIceSplit':
		return False
	T = frame['I3MCTree'].get_daughters(frame['I3MCTree'].get_primaries()[0])
	if len(T) > 1:
		print('Cluster')
		return False
	if not (str(T[0].type)=='MuPlus' or str(T[0].type)=='MuMinus'):
		print('Not a muon')
		return False

	mu = T[0]
	if use_truth == 1:
		endp = np.array(mu.pos+mu.dir*mu.length)
		stopped_xy = border.contains_point((endp[0],endp[1]),radius=-trackradius)
		stopped_z = endp[2] > -400
		stopped_theta = (mu.dir.zenith >= 40*np.pi/180) * (mu.dir.zenith<=70*np.pi/180)
		stoppedmuon = int(stopped_xy*stopped_z*stopped_theta)
		if stoppedmuon == 0:
			print('Not stopped')
			return False
	dom = []
	string = []
	charge = []
	time = []
	series = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, "InIcePulses")

	for i,pulse in enumerate(series):
		for hit in pulse[1]:
			string.append(pulse[0].string)
			dom.append(pulse[0].om)
			charge.append(hit.charge)
			time.append(hit.time)
	if len(charge) < 8:# or len(charge) > 100:
		print('Not enough or too many DOM hits')
		return False
	domstr = np.column_stack([string,dom]).astype(str)
	domstr = np.char.zfill(domstr,2).astype(object)
	domstr = np.sum(domstr,axis=1).astype(str)
	domindex = np.array([np.where(domfile[:,3]==i)[0] for i in domstr])[:,0]
	x,y,z = domfile[:,5],domfile[:,6],domfile[:,7]

	xp = x[domindex].astype(float)
	yp = y[domindex].astype(float)
	zp = z[domindex].astype(float)
	chargep = np.array(charge).astype(float)
	time = np.array(time).astype(float)
	
	if use_truth == 0:
		features = np.column_stack((xp,yp,zp,chargep,time))
		for i in range(5):
			features[:,i] = (features[:,i]-meanstd[i,0])/meanstd[i,1]
		features = np.vstack((features,np.zeros((100,features.shape[1]))))[:100,:]
		features = np.rollaxis(features.reshape(features.shape[0],features.shape[1],1),2,0)
		prediction = model.predict(features)[0]
		if prediction < th1:
			print('Too low prediction')
			return False
	charge = np.zeros(len(x))
	charge[domindex] = chargep

	coords = np.array([x,y,z]).T.astype(float)
	mu = T[0]
	startp = np.array(mu.pos)
	endp = np.array(mu.pos+mu.dir*mu.length)
	point_vector = coords-startp
	seg_vector = endp-startp
	t = np.dot(point_vector,seg_vector)/sum(seg_vector**2)
	close = startp + t[:,np.newaxis]*seg_vector
	close[t<0] = startp
	close[t>1] = endp

	domrad = np.linalg.norm(coords-close,axis=1)
	domrad = (domrad < domradius) * (domrad > 20)

	if np.sum(domrad) == 0:
		print('Zero DOMs within track radius')
		return False

	print('Accepted')
	domstr = np.array(domfile[domrad,3]).astype(float)
	group = np.array(domfile[domrad,4]).astype(float)
	charge = charge[domrad].astype(float)
	x = x[domrad].astype(float)
	y = y[domrad].astype(float)
	z = z[domrad].astype(float)	
	HQE = domfile[domrad,2].astype(str)
	HQE[HQE=='NQE'] = 0
	HQE[HQE=='HQE'] = 1
	HQE[HQE=='BAD'] = 2
	HQE=np.array(HQE).astype(float)

	frame['x'] = dataclasses.I3VectorDouble(x)
	frame['y'] = dataclasses.I3VectorDouble(y)
	frame['z'] = dataclasses.I3VectorDouble(z)
	frame['charge'] = dataclasses.I3VectorDouble(charge)
	frame['domstr'] = dataclasses.I3VectorDouble(domstr)
	frame['group'] = dataclasses.I3VectorDouble(group)
	frame['HQE'] = dataclasses.I3VectorDouble(HQE)


def savefeatures(frame,fin=None,names=None):
	for i,name in enumerate(names):
		fin[i].append(frame[name])

use_truth = 1
#Tensorflow stuff
if use_truth == 0:
	from tensorflow.keras.models import load_model as lm
	model = lm('/home/sstray/test/condor/corsikafiles/corsika_model_test/classification_model')


import os
from I3Tray import *
from icecube import dataio
from icecube import dataclasses
from icecube import simclasses
from icecube import recclasses
import numpy as np
#import cProfile, pstats, StringIO
import argparse
import matplotlib.path as mpath
bordercoords = np.array([(-256.1400146484375, -521.0800170898438), (-132.8000030517578, -501.45001220703125), (-9.13000011444092, -481.739990234375), (114.38999938964844, -461.989990234375), (237.77999877929688, -442.4200134277344), (361.0, -422.8299865722656), (405.8299865722656, -306.3800048828125), (443.6000061035156, -194.16000366210938), (500.42999267578125, -58.45000076293945), (544.0700073242188, 55.88999938964844), (576.3699951171875, 170.9199981689453), (505.2699890136719, 257.8800048828125), (429.760009765625, 351.0199890136719), (338.44000244140625, 463.7200012207031), (224.5800018310547, 432.3500061035156), (101.04000091552734, 412.7900085449219), (22.11000061035156, 509.5), (-101.05999755859375, 490.2200012207031), (-224.08999633789062, 470.8599853515625), (-347.8800048828125, 451.5199890136719), (-392.3800048828125, 334.239990234375), (-437.0400085449219, 217.8000030517578), (-481.6000061035156, 101.38999938964844), (-526.6300048828125, -15.60000038146973), (-570.9000244140625, -125.13999938964844), (-492.42999267578125, -230.16000366210938), (-413.4599914550781, -327.2699890136719), (-334.79998779296875, -424.5)])
border = mpath.Path(bordercoords)
border_zminmax = [-512.82,524.56]
trackradius = 100
height = 100

domfile = np.loadtxt('/home/sstray/test/condor/corsikafiles/dom_coords_spacing_HQE.txt',dtype='str')
testfilename = '/data/sim/IceCube/2012/filtered/level2/CORSIKA-in-ice/11499/77000-77999/Level2_IC86.2012_corsika.011499.077910.i3.bz2'
outputname = 'featurelist.hdf5'
meanstd = np.loadtxt('/home/sstray/test/condor/corsikafiles/meanstd.txt')
domradius = 100
th1 = 0.9
th2 = 0.99

from icecube.sim_services import propagation
def test_my_little_function():
	my_features = ['x','y','z','charge','domstr','group','HQE']
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--infiles',required=True)
	args = parser.parse_args()
	filelist = args.infiles.split(',')
	for nums in range(len(filelist)):
		featurestemp = [[] for i in range(len(my_features))]
		global event_run
		event_run = filelist[nums].split('/')[-1].split('.')[-3]

                tray = I3Tray()
		tray.AddModule('I3Reader','read_stuff',Filename=filelist[nums])
		tray.AddModule('Delete',keys=["MMCTrackList"]) # <-- Add this
		tray.AddSegment(propagation.RecreateMCTree,"recreate",
                                RawMCTree="I3MCTree_preMuonProp",
				RNGState="I3MCTree_preMuonProp_RNGState",  # <-- Change this
				Paranoia=False)
                tray.AddModule(test,'test')
                tray.AddModule(savefeatures,'asd',fin=featurestemp,names=my_features)

                print("try: Excecute")
                try:
                        tray.Execute()
                except:
                        print("bad file: skipping...")
                        continue

		for i in range(len(featurestemp)):
			 featurestemp[i] = [item for sublist in featurestemp[i] for item in sublist]
		if nums == 0:
			featuresarray = np.column_stack(featurestemp)
		else:
			featuresarray = np.concatenate((featuresarray,np.column_stack(featurestemp)),axis=0)
		print(featuresarray.shape)
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
	f = tables.open_file(outputname, 'w',filters=tables.Filters(complib='zlib', complevel=6))
	f.create_carray('/','array',obj=thearray)
	f.close()
