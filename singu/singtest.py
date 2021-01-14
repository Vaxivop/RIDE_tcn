import h5py
import numpy as np
from tensorflow.keras.models import load_model as lm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--infiles',required=True)
parser.add_argument('-rad','--radius')
parser.add_argument('-th','--threshold')
args = parser.parse_args()

if args.radius:
	rad = args.radius
else:
	rad = 75
	print('No radius supplied. Will default to 75')
if args.threshold:
	th = args.threshold
else:
	th = 0.9
	print('No threshold supplied. Will default to 0.9')
meanstd = np.loadtxt('/home/sstray/test/condor/singu/meanstd.txt')
model = lm('/home/sstray/test/condor/singu/my_model/classification_model')

def get_groups(features_in,event_in):
  features = features_in
  eventnum = event_in
  grouped_features = ([features[np.isin(eventnum,[i])] for i in np.unique(eventnum)])
  grouped_features = [np.vstack((np.array(i),np.zeros((100,features.shape[1]))))[:100] for i in grouped_features if len(i)<100 if len(i)>10]
  grouped_features = np.dstack(grouped_features)
  grouped_features = np.rollaxis(grouped_features,2,0)
  return grouped_features


filelist = args.infiles.split(',')
print(filelist)
for nums in range(len(filelist)):
	hf = h5py.File(filelist[nums],'r')
	data = np.array(hf[list(hf.keys())[0]])
	hf.close()
	
	features = data[:,1:]
	print(features.shape)
	events = data[:,0]
	for i in range(5):
		features[:,i] = (features[:,i]-meanstd[i,0])/meanstd[i,1]
		print(np.mean(features[:,i]))
	events = events[features[:,-1] < rad]
	features = features[features[:,-1] < rad]
	print(features.shape)
	print(len(np.unique(events)))
	grouped = get_groups(features,events)
	print(grouped.shape)
	print(events)
	prediction = model.predict(grouped[:,:,:6])[:,0]
	grouped = grouped[prediction>th,:,:]
	print(grouped.shape)
	
	features = grouped.reshape(grouped.shape[0]*grouped.shape[1],grouped.shape[2])
	features = features[features[:,0]!=0]
	for i in range(5):
  		features[:,i] = features[:,i]*meanstd[i,1] + meanstd[i,0]
	print(features.shape)
	print(filelist[nums].split('/')[-1])
	hf2 = h5py.File(str(filelist[nums].split('/')[-1]),'w')
	hf2.create_dataset("array", data=features,compression="gzip",compression_opts=6)
	hf2.close()
