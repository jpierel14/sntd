import pickle,sys,sntd
from optparse import OptionParser
from copy import copy

njobs=10
nlcs=10
parser = OptionParser()

(options,args)=parser.parse_args()


all_dat=pickle.load(open('sntd_data.pkl','rb'))
inds=[int(nlcs/njobs)*int(sys.argv[1]),int(nlcs/njobs)*int(sys.argv[1])+int(nlcs/njobs)]

all_res=[]
for i in range(inds[0],inds[1]):
    fitCurves=sntd.fit_data(verbose=True,color_curve=None,nMicroSamples=100,refImage="image_1",kernel="RBF",batch_python_path=None,
                            nbatch_jobs=10,batch_script=None,batch_partition="test1",par_or_batch="batch",fit_prior=None,
                            fitOrder=['image_1', 'image_2'],microlensing=None,showPlots=False,seriesError=None,
                            guess_amplitude=True,flip=False,dust=None,effect_frames=[],effect_names=[],
                            refModel=None,t0_guess=None,method="parallel",constants={'z': 1.33},
                            ignore=None,bounds={'t0': (-20, 20), 'x1': (-3, 3), 'c': (-1, 1)},
                            params=['x0', 't0', 'x1', 'c'],models="salt2-extended",bands=['F110W', 'F160W'],
                            curves=all_dat[i],snType="Ia",kwargs={'modelcov': False, 'npoints': 100, 'maxiter': None},fitting_method="nest",modelcov=False,npoints=100,maxiter=None)
    all_res.append(copy(fitCurves))
pickle.dump(all_res,open('sntd_fit%s.pkl'%sys.argv[1],'wb'))