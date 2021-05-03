#%run ride_df.py
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

newcolumns = ('RIDE','mean_charge','total_charge','expected_hits','actual_hits','x','y','z','domstr','group   ...:\ ','status')
all_doms = pd.DataFrame(np.loadtxt('RIDE_jvmead.txt'),columns=newcolumns)
doms = all_doms[ (all_doms.status!=2) ]
hqe = doms [ (doms.status==1) ]
nqe = doms [ (doms.status==0) ]

print(doms.head(5))

#####################

f = plt.figure()

plt.plot(nqe["z"],nqe["RIDE"],"o",markersize=1)
plt.plot(hqe["z"],hqe["RIDE"],"o",markersize=1)

plt.ylim((0,2.7))
plt.ylabel('RIDE')
plt.xlabel('z (m)')

plt.legend(["NQE","HQE"], loc='best')
plt.show()
f.savefig("ride_fz.pdf", bbox_inches='tight')
