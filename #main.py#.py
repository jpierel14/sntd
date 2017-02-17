import sntd
#filename="example_photometric_data.dat"
#filename="refsdalS1_psfphot.dat"
filename="testRef2.dat"
#filename="myData.pkl"
#filename="test.rdb"
#tab=sncosmo.read_lc(filename,verbose=False,masked=True)
temp=sntd.read_data(filename)
print(temp.F105W.fluxes)
print(temp.F105W.table)
"""
sntd.write_data(temp,'myData.pkl')
temp=sntd.read_data('myData.pkl')
print(temp.F105.mags)
"""