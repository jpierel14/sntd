import sntd

filename="example_photometric_data.dat"
#filename="myData.pkl"
#filename="test.rdb"
#tab=sncosmo.read_lc(filename,verbose=False,masked=True)
temp=sntd.read_data(filename)
print(temp)
sntd.write_data(temp,'myData.pkl')
temp=sntd.read_data('myData.pkl')
print(temp.sdssr.jds)