import sntd
#filename2="example_photometric_data.dat"
filename="refsdalS2_psfphot.dat"
filename2="testRef2.dat"
#filename="myData.pkl"
#filename="test.rdb"
#tab=sncosmo.read_lc(filename,verbose=False,masked=True)
temp1=sntd.read_data(filename)
temp2=sntd.read_data(filename2)
fits=sntd.fit_data([temp1,temp2])
"""
sntd.write_data(temp,'myData.pkl')
temp=sntd.read_data('myData.pkl')
print(temp.F105.mags)
"""