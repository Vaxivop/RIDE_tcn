import os
from icecube import icetray, dataclasses
from icecube import dataio
from I3Tray import *
import numpy as np
import tables
import argparse

parser = argparse.ArgumentParser(description='script for reconstructing data and apply cuts on L2 sim files')
parser.add_argument('-i', '--infiles', help='data files to make the cuts on (corsika level2 i3 files).', required=False)
args = parser.parse_args()

file = '/data/sim/IceCube/2012/filtered/level2/CORSIKA-in-ice/11499/77000-77999/Level2_IC86.2012_corsika.011499.077780.i3.bz2'
domfile = np.genfromtxt('/home/sstray/test/condor/corsikafiles/dom_coords_spacing_HQE.txt',dtype='str')

def RIDE(frame):
	T = frame['I3MCTree'].get_daughters(frame['I3MCTree'].get_primaries()[0])
	if len(T) != 1:
		return False
	if not (str(T[0].type)=='MuPlus' or str(T[0].type)=='MuMinus'):
		return False
	muon = T[0]
	series = frame['SplitInIcePulses'].apply(frame)
	string = []
	dom = []
	charge = []
	time = []
	for i,pulse in enumerate(series):
		string.append(pulse[0].string)
		dom.append(pulse[0].om)
		charge.append(pulse[1][0].charge)
		time.append(pulse[1][0].time)
	domstr = np.array([string,dom]).T.astype(int).astype(str)
	domstr = np.char.zfill(domstr,2).astype(object)
	domstr = np.sum(domstr,axis=1).astype(str)
	dombool = np.isin(domfile[:,3],domstr)
	
	header = frame['I3EventHeader']
	events = (np.ones(len(domstr))*(header.event_id+header.sub_event_id)).astype(int)

	charge = np.array(charge)
	time = np.array(time)
	x = domfile[dombool,5].astype(float)
	y = domfile[dombool,6].astype(float)
	z = domfile[dombool,7].astype(float)
	group = domfile[dombool,4].astype(int)
	domstring = domfile[dombool,3].astype(int)
	HQE = domfile[dombool,2]
	HQE[HQE=='NQE'] = 0
	HQE[HQE=='HQE'] = 1
	HQE[HQE=='BAD'] = 2
	HQE = HQE.astype(int)
	coords = np.vstack((x,y,z)).T
	startp = np.array(muon.pos)
	endp = np.array((muon.pos)+muon.dir*muon.length)
	point_vector = coords-startp
	seg_vector = endp-startp
	t = np.dot(point_vector,seg_vector)/sum(seg_vector**2)
	close = startp + t[:,np.newaxis]*seg_vector
	close[t<0] = startp
	close[t>1] = endp
	domrad = np.linalg.norm(coords-close,axis=1)

	frame['events'] = dataclasses.I3VectorUInt64(events)
	frame['x'] = dataclasses.I3VectorDouble(x)
	frame['y'] = dataclasses.I3VectorDouble(y)
	frame['z'] = dataclasses.I3VectorDouble(z)
	frame['charge'] = dataclasses.I3VectorDouble(charge)
	frame['time'] = dataclasses.I3VectorDouble(time)
	frame['domstring'] = dataclasses.I3VectorUInt(domstring)
	frame['group'] = dataclasses.I3VectorUInt(group)
	frame['HQE'] = dataclasses.I3VectorUInt(HQE)
	frame['domrad'] = dataclasses.I3VectorDouble(domrad)
def featuresave(frame,col=None,features=None):
	for i, name in enumerate(col):
		features[i].append(frame[name])
def runtray():
	featurecol = ['events','x','y','z','charge','time','domstring','group','HQE','domrad']
	features = [[] for i in featurecol]
	tray = I3Tray()
	tray.AddModule('I3Reader', 'I3Reader', filename=file)
	tray.AddModule(RIDE,'RIDE')
	tray.AddModule(featuresave,'save_features',col=featurecol,features=features)
	tray.Execute()
	array = np.array([np.concatenate(i) for i in features]).T
runtray()

