cimport cython
from cython.view cimport array as cvarray
from libc.math cimport pow
cimport numpy as np
import numpy as np
from libcpp cimport bool


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef float fast_chisq(np.float32_t [:] y,np.float32_t [:] yerr,np.float32_t [:] model_observations):    
		cdef Py_ssize_t i = 0
		cdef Py_ssize_t size = y.shape[0]
		cdef float chisquare = 0
		for i in range(size):
			temp=y[i]-model_observations[i]
			chisquare += temp*temp/(yerr[i]*yerr[i])
		return chisquare


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef np.ndarray[np.float32_t] salt2_source_flux_fast(np.float32_t [:,:] m0,np.float32_t [:,:] m1,np.float32_t [:] color,float x0,float x1,float c):
		cdef Py_ssize_t i,j = 0
		cdef Py_ssize_t size1 = m0.shape[0]
		cdef Py_ssize_t size2 = m0.shape[1]
		cdef np.float32_t [:,:] out=cvarray(shape=(size1,size2),itemsize=sizeof(np.float32_t),format="f")
		for i in range(size1):
			for j in range(size2):
				out[i][j]=x0*(m0[i][j] + x1 * m1[i][j]) *pow(10., (-0.4 * color[j] * c))

		return np.array(out)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef float salt2_flux_fast(np.float32_t [:] y,np.float32_t [:] yerr,np.float32_t [:] model_observations):
		cdef Py_ssize_t i = 0
		cdef Py_ssize_t size = y.shape[0]
		cdef float chisquare = 0
		for i in range(size):
			temp=y[i]-model_observations[i]
			chisquare += temp*temp/(yerr[i]*yerr[i])








@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef np.ndarray[np.float32_t] fast_bandflux_single(np.float32_t [:,:] flux,np.float32_t [:] wave,np.float32_t [:] trans,
													np.float32_t [:] dwave, np.float32_t [:] zp,
													np.float32_t [:] zpsys,np.int32_t [:,:] band_inds,
													float HC_ERG_AA):
		cdef Py_ssize_t i,j = 0
		cdef Py_ssize_t size = flux.shape[0]

		cdef np.float32_t [:] out=cvarray(shape=(size,),itemsize=sizeof(np.float32_t),format="f")
		cdef float total_flux = 0

		for i in range(size):
			total_flux=0
			for j in range(band_inds[i][0],band_inds[i][1]):
				total_flux+=wave[j]*trans[j]*flux[i][j]*dwave[j]


			out[i]=total_flux*pow(10,0.4*zp[i])/zpsys[i]/HC_ERG_AA
		return np.array(out)

