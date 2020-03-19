import pickle,sys,sntd,os
from optparse import OptionParser
from copy import copy

njobs=2
nlcs=2
parser = OptionParser()

(options,args)=parser.parse_args()


import test

print('hello')


all_dat=pickle.load(open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                      'sntd_data.pkl'),'rb'))
inds=[int(nlcs/njobs)*int(sys.argv[1]),int(nlcs/njobs)*int(sys.argv[1])+int(nlcs/njobs)]
inds[1]=min(inds[-1],len(all_dat))

all_res=[]
for i in range(inds[0],inds[1]):
    fitCurves=sntd.fit_data(modelcov=False,npoints=100,maxiter=None,verbose=True,warning_supress=True,color_curve=None,nMicroSamples=100,refImage="image_1",kernel="RBF",test_micro=True,guess_amplitude=True,wait_for_batch=False,batch_python_path=None,nbatch_jobs=2,batch_script=None,batch_partition="test",par_or_batch="batch",fit_prior=None,color_bands=['F110W', 'F160W'],fitOrder=['image_2', 'image_1'],microlensing=None,flip=False,dust=None,batch_init=None,effect_frames=[],effect_names=[],t0_guess=None,method="parallel",constants=all_dat[i].constants,ignore=None,bounds={'t0': (-20, 20), 'x1': (-3, 3), 'c': (-1, 1), 'mu': (0.5, 2), 'td': (-15, 15)},params=['x0', 't0', 'x1', 'c'],models="salt2-extended",bands=None,curves=all_dat[i],snType="Ia")
    fitCurves=sntd.fit_data(modelcov=False,npoints=100,maxiter=None,verbose=True,warning_supress=True,color_curve=None,nMicroSamples=100,refImage="image_1",kernel="RBF",test_micro=False,guess_amplitude=True,wait_for_batch=False,batch_python_path=None,nbatch_jobs=2,batch_script=None,batch_partition="test",par_or_batch="batch",fit_prior=fitCurves,color_bands=['F110W', 'F160W'],fitOrder=['image_2', 'image_1'],microlensing=None,flip=False,dust=None,batch_init=None,effect_frames=[],effect_names=[],t0_guess=None,method="color",constants=all_dat[i].constants,ignore=None,bounds={'t0': (-20, 20), 'x1': (-3, 3), 'c': (-1, 1), 'mu': (0.5, 2), 'td': (-15, 15)},params=['x0', 't0', 'x1', 'c'],models="salt2-extended",bands=fitCurves.micro_color_bands,curves=fitCurves,snType="Ia")

    all_res.append(copy(fitCurves))

pickle.dump(all_res,open(os.path.join(os.path.abspath(os.path.dirname(__file__)),'sntd_fit%s.pkl'%sys.argv[1]),'wb'))
