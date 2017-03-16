
import sntd,timeit
from multiprocessing import Pool
#filename2="example_photometric_data.dat"
filename="refsdalS2_psfphot.dat"
filename2="testRef2.dat"
#filename="myData.pkl"
#filename="test.rdb"
#tab=sncosmo.read_lc(filename,verbose=False,masked=True)
temp1=sntd.read_data(filename)
temp2=sntd.read_data(filename2)
#sntd.fit_data(temp1,bounds={'z':(1.2,1.5)})
print(timeit.timeit("sntd.fit_data(sntd.read_data('example_photometric_data.dat'))",setup="import sntd",number=1))
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