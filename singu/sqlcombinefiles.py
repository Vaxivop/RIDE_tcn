import glob
import numpy as np
import sqlite3
import pandas as pd
import h5py
from random import shuffle
files = glob.glob('/data/user/sstray/tcn_data/*.hdf5')
print(files)
shuffle(files)
print(files)
conn = sqlite3.connect('corsikadb.db')
c = conn.cursor()
c.execute('''create table features([x] float,[y] float,[z] float,[charge] float,[time] float,[eventid] int)''')
c.execute('''create table truth([muon] int,[eventid] int)''')
for i,name in enumerate(files):
	print(i)
	temp = h5py.File(name,'r')
	features = temp['features'][:]
	truth = temp['truth'][:]
	temp.close()

	features=pd.DataFrame(features,columns=('x','y','z','charge','time','eventid'))
	features=features.apply(pd.to_numeric,downcast='integer',errors='ignore')
	truth=pd.DataFrame(truth,columns=('muon','eventid'))
	truth=truth.apply(pd.to_numeric,downcast='integer',errors='ignore')
	
	features.to_sql('features',conn,if_exists='append',index=False)
	truth.to_sql('truth',conn,if_exists='append',index=False)
