import pickle,sys,sntd,os,traceback
from optparse import OptionParser
from copy import copy
import numpy as np

nlcs_per=nlcsreplace
parser = OptionParser()

(options,args)=parser.parse_args()

batchinitreplace

all_dat=pickle.load(open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                      'sntd_data.pkl'),'rb'))
all_const=pickle.load(open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                        'sntd_constants.pkl'),'rb'))
inds=[int(int(sys.argv[1])*nlcs_per),(int(sys.argv[1])+1)*int(nlcs_per)]
inds[1]=min(inds[-1],len(all_dat))

all_input=[]
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
    all_input.append(all_dat[i])
try:
    fitCurves=sntdcommandreplace
except Exception as e:
    print('Failed')
    print(traceback.format_exc())
    fitCurves=traceback.format_exc()
    
for i in range(len(all_input)):
    filename=os.path.join(os.path.abspath(os.path.dirname(__file__)),'sntd_fit%s_%i.pkl'%(sys.argv[1],i))
    pickle.dump(fitCurves[i],open(filename,'wb'))
