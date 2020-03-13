from __future__ import print_function
import snsedextend

sedFile=snsedextend.example_sed
myLC=snsedextend.load_example_data()
print(myLC)

colorTable=snsedextend.curveToColor(myLC,colors=['U-B', 'r-J', 'r-H', 'r-K'], snType='Ic', zpsys='vega', bounds={'hostebv': (-1, 1), 't0': (53787.94, 53797.94)},constants={'mwr_v': 3.1, 'mwebv': '0.1267', 'z': '0.033529863', 'hostr_v': 3.1}, dust='CCM89Dust', effect_frames=['rest', 'obs'], effect_names=['host', 'mw'])

print(colorTable)

curveDict=snsedextend.fitColorCurve(colorTable)

print(curveDict.keys())

print(curveDict['U-B'])
