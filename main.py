import sntd

myMISN = sntd.createMultiplyImagedSN('salt2', 'Ia', 1.33,z_lens=.53, bands=['F110W'],
                                     zp=[26.8], cadence=5., epochs=35.,skynoiseRange=(.001,.005),gain=70. , time_delays=[10., 78.],
                                     magnifications=[7,3.5], objectName='My Type Ia SN', telescopename='HST',minsnr=5.0)

