cimport cython
from cython.view cimport array as cvarray
cimport numpy as np


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