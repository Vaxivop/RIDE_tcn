import glob
import h5py
import numpy as np
actualfiles = glob.glob('/data/user/sstray/files_0.99_tcn/*.hdf5')
biglist = []
for i in actualfiles:
	biglist.append(np.array(h5py.File(i,'r')['array']))
print('Appended...')
bigarray = np.vstack(biglist)
print(len(bigarray))
x = bigarray[:,0]
y = bigarray[:,1]
z = bigarray[:,2]
charge = bigarray[:,3]
domstr = bigarray[:,4]
group = bigarray[:,5]
HQE = bigarray[:,6]

charge_dom = np.zeros(len(np.unique(domstr)))
group_dom = np.zeros_like(charge_dom)
domstr_dom = np.zeros_like(charge_dom)
x_dom = np.zeros_like(charge_dom)
y_dom = np.zeros_like(charge_dom)
z_dom = np.zeros_like(charge_dom)
totalhits_dom = np.zeros_like(charge_dom)
totalcharge_dom = np.zeros_like(charge_dom)
exphits_dom = np.zeros_like(charge_dom)
status_dom = np.zeros_like(charge_dom)
RIDE = np.zeros_like(charge_dom)
print('Starting RIDE calc...')
print(len(np.unique(domstr)))
for i, unds in enumerate(np.unique(domstr)):
	if i%100 == 0:
		print(i)
	temp = np.isin(domstr,unds)
	charge_dom[i] = np.mean(charge[temp])
	RIDE[i] = np.copy(charge_dom[i])
	group_dom[i] = group[temp][0]
	x_dom[i] = x[temp][0]
	y_dom[i] = y[temp][0]
	z_dom[i] = z[temp][0]
	domstr_dom[i] = domstr[temp][0]
#	print('Domstr '+str(domstr_dom[i])+' has a mean charge of '+str(charge_dom[i])+' and is in group '+str(group_dom[i]))
	status_dom[i] = HQE[temp][0]
	totalcharge_dom[i] = np.sum(charge[temp])
	exphits_dom[i] = len(charge[temp])
	totalhits_dom[i] = np.sum(charge[temp]>0)
d40b = group_dom == 40
dom40 = np.column_stack((RIDE[d40b],charge_dom[d40b],totalcharge_dom[d40b],exphits_dom[d40b],totalhits_dom[d40b],x_dom[d40b],y_dom[d40b],z_dom[d40b],domstr_dom[d40b],group_dom[d40b],status_dom[d40b]))
np.savetxt('dom40.txt',dom40)
print('Calculating RIDE...')
for ung in np.unique(group_dom):
	temp = np.isin(group_dom,ung)
	areNQE = status_dom[temp] == 0
	chargemed = np.median(RIDE[temp][areNQE])
	if chargemed !=0:
		RIDE[temp] = RIDE[temp]/chargemed
	else:
		RIDE[temp] = 0
RIDE = np.column_stack((RIDE,charge_dom,totalcharge_dom,exphits_dom,totalhits_dom,x_dom,y_dom,z_dom,domstr_dom,group_dom,status_dom))
np.savetxt('RIDE_tcn_0.99.txt',RIDE)
