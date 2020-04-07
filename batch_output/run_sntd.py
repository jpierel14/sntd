import pickle,sys,sntd,os,traceback
from optparse import OptionParser
from copy import copy
import numpy as np

nlcs_per=1
parser = OptionParser()

(options,args)=parser.parse_args()

print("Nothing to initialize...")

all_dat=pickle.load(open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                      'sntd_data.pkl'),'rb'))
all_const=pickle.load(open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                        'sntd_constants.pkl'),'rb'))
inds=[int(int(sys.argv[1])*nlcs_per),(int(sys.argv[1])+1)*int(nlcs_per)]
inds[1]=min(inds[-1],len(all_dat))

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
        fitCurves=sntd.fit_data(maxcall=5,minsnr=0,clip_data=False,verbose=True,n_per_node=1,warning_supress=True,color_curve=None,nMicroSamples=100,refImage="image_1",kernel="RBF",trial_fit=True,test_micro=False,guess_amplitude=True,wait_for_batch=True,batch_python_path=None,nbatch_jobs=2,batch_partition=None,par_or_batch="batch",fit_prior=None,min_points_per_band=3,color_bands=None,fitOrder=None,microlensing=None,flip=False,dust=None,cut_time=None,batch_init=None,effect_frames=[],effect_names=[],t0_guess={'image_1': 10, 'image_2': 70},method="parallel",constants=all_dat[i].constants,ignore=None,bounds={'t0': (-15, 15), 'x1': (-2, 2), 'c': (0, 1), 'td': (-15, 15), 'mu': (0.5, 2)},params=['x0', 'x1', 't0', 'c'],models="salt2-extended",bands=['bessellb', 'bessellr'],curves=all_dat[i],snType="Ia")
        all_res.append(copy(fitCurves))
    except Exception as e:
        print('Failed')
        print(traceback.format_exc())
        all_res.append(None)

pickle.dump(all_res,open(os.path.join(os.path.abspath(os.path.dirname(__file__)),'sntd_fit%s.pkl'%sys.argv[1]),'wb'))
