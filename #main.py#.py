'''

import sncosmo
import pickle
data=sncosmo.load_example_data()
model=sncosmo.Model('salt2')
res, fitted_model = sncosmo.fit_lc(data, model,['z', 't0', 'x0', 'x1', 'c'],bounds={'z':(0.3, 0.7)})
#model._source=None
with open('test.pkl','wb') as handle:
    pickle.dump(sncosmo.Model('salt2'),handle,protocol=2)

'''
import sntd,timeit,glob
from multiprocessing import Pool
files=glob.glob('data/*.dat')
#filename2="example_photometric_data.dat"
#filename="refsdalS2_psfphot.dat"
#filename2="testRef2.dat"
#filename="myData.pkl"
#filename="test.rdb"
#tab=sncosmo.read_lc(filename,verbose=False,masked=True)
#temp1=sntd.read_data(filename)
curves=sntd.curveDict(telescopename='Hubble',object='Refsdal')
for f in files:
    temp2=sntd.read_data(f)
    print(temp2.table)
#sntd.fit_data(temp1,bounds={'z':(1.2,1.5)})
#sntd.fit_data(temp2,bounds={'z':(.3,.7)})
#print(timeit.timeit("sntd.fit_data(sntd.read_data('example_photometric_data.dat'),bounds={'z':(0.3, 0.7)})",setup="import sntd",number=1))
#print(timeit.timeit("sntd.fit_data(sntd.read_data('refsdalS2_psfphot.dat'),bounds={'z':(1.2, 1.5)})",setup="import sntd",number=1))
"""
sntd.write_data(temp,'myData.pkl')
temp=sntd.read_data('myData.pkl')
print(temp.F105.mags)
"""
'''
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np

x=np.linspace(8,12,100)
fig,ax=plt.subplots(1,1)
ax.plot(x,norm.pdf(x,10,.2))
#ax.plot(x,norm.pdf(x,0,.6))
plt.show()
'''