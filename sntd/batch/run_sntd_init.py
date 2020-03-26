import pickle,sys,sntd,os
from optparse import OptionParser
from copy import copy

njobs=njobsreplace
nlcs=nlcsreplace
parser = OptionParser()

(options,args)=parser.parse_args()

batchinitreplace

all_dat=pickle.load(open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                      'sntd_data.pkl'),'rb'))
inds=[int(nlcs/njobs)*int(sys.argv[1]),int(nlcs/njobs)*int(sys.argv[1])+int(nlcs/njobs)]
inds[1]=min(inds[-1],len(all_dat))

all_res=[]
for i in range(inds[0],inds[1]):
    if isinstance(all_dat[i],str):
        all_dat[i]=pickle.load(open(all_dat[i],'rb'))
    try:
        fitCurves=sntdcommandreplace
        all_res.append(copy(fitCurves))
    except:
        all_res.append(None)

pickle.dump(all_res,open(os.path.join(os.path.abspath(os.path.dirname(__file__)),'sntd_fit%s.pkl'%sys.argv[1]),'wb'))
