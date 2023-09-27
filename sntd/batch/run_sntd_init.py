import pickle,sys,sntd,os,traceback
from optparse import OptionParser
from copy import copy
import numpy as np
import sncosmo

njobs=njobsreplace
nlcs=nlcsreplace
parser = OptionParser()

(options,args)=parser.parse_args()

batchinitreplace

all_dat=pickle.load(open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
									  'sntd_data.pkl'),'rb'))
all_const=pickle.load(open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
										'sntd_constants.pkl'),'rb'))

inds=[int(nlcs/njobs)*int(sys.argv[1]),int(nlcs/njobs)*int(sys.argv[1])+int(nlcs/njobs)]
inds[1]=min(inds[-1],len(all_dat))
filename=os.path.join(os.path.abspath(os.path.dirname(__file__)),'sntd_fit%s.pkl'%sys.argv[1])
with open(filename.replace('.pkl','.DONE'),'w') as f:
	f.write('FALSE')

all_res=[]
for i in range(inds[0],inds[1]):
	if isinstance(all_dat[i],str):
		all_dat[i]=pickle.load(open(all_dat[i],'rb'))
	all_dat[i].constants={}
	if all_const is not None:
		for c in all_const.keys():
			if isinstance(all_const[c],(list,tuple,np.ndarray)):
				all_dat[i].constants[c]=all_const[c][i]
			else:
				all_dat[i].constants[c]=all_const[c]
	try:
		fitCurves=sntdcommandreplace
		all_res.append(copy(fitCurves))
	except Exception as e:
		print('Failed')
		print(traceback.format_exc())
		all_res.append(None)


try:
	pickle.dump(all_res,open(filename,'wb'))
except:
	for i in range(len(all_res)):
		if all_res[i] is not None:
			try:
				all_res[i].color.fits.model = None
			except:
				pass
			try:
				all_res[i].series.fits.model = None
			except:
				pass
			try:
				all_res[i].fits.model = None
			except:
				pass
	pickle.dump(all_res,open(filename,'wb'))
with open(filename.replace('.pkl','.DONE'),'w') as f:
	f.write('TRUE')

