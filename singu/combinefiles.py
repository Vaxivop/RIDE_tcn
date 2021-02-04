import glob
import h5py
import numpy as np
files = glob.glob('/data/user/sstray/group83/*.hdf5')
out = h5py.File('output_83.hdf5','a')
for i,name in enumerate(files):
	temp = h5py.File(name,'r')
	features = temp['features']
	truth = temp['truth']
	if i == 0:
		fbig = out.create_dataset('features',data=features,maxshape=(None,features.shape[1],features.shape[2]),compression='gzip',compression_opts=6,chunks=True,dtype=np.float32)
		tbig = out.create_dataset('truth',data=truth,maxshape=(None,truth.shape[1]),compression='gzip',compression_opts=6,chunks=True,dtype=np.float32)
		print(fbig)
	else:
		fsize,tsize = len(fbig),len(tbig)
		print(fsize)
		fbig.resize((int(fsize+len(features)),features.shape[1],features.shape[2]))
		fbig[fsize:,:] = np.array(features,dtype=np.float32)
		print(fbig)
		tbig.resize((int(tsize+len(truth)),truth.shape[1]))
		tbig[tsize:] = np.array(truth,dtype=np.float32)
	temp.close()
print(tbig)
print(fbig)
out.close()
