import sntd
import matplotlib.pyplot as plt

#fitData=sntd.realizeMicro(nray=10,kappas=1,kappac=.3,gamma=.4)
myMISN = sntd.createMultiplyImagedSN('salt2', 'Ia', 1.33,z_lens=.53, bands=['F110W','F125W'],
                                     zp=[26.8,26.2], cadence=5., epochs=35.,skynoiseRange=(.001,.005),gain=70. ,
                                     time_delays=[10., 78.], magnifications=[7,3.5], objectName='My Type Ia SN',
                                     telescopename='HST',minsnr=5.0)#,microlensing_type='AchromaticMicrolensing',
                                     #microlensing_params=fitData)

myMISN.combine_curves(time_delays={'image_1':10,'image_2':78},magnifications={'image_1':7,'image_2':3.5})
fig=myMISN.plot_object(savefig=False,showModel=False,combined=True,showMicro=False)
plt.show()