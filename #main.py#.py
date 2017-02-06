import sntd
from astropy.table import Table
import sncosmo


filename="example_photometric_data.dat"
tab=sncosmo.read_lc(filename,verbose=False,masked=True)
temp=sntd.curve()
temp.table=tab
print(temp.table)
