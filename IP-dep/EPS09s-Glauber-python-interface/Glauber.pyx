import numpy as np
cimport numpy as np
import h5py, vegas
from libc.math cimport *

cdef extern from "../EPS09s-interface-cpp/eps09s.h":
	cdef double Glauber_TA(int a, double s)

cpdef Glauber_TA_py(int A, double[:] s):
	
	return np.array([Glauber_TA(A, ss) for ss in s])
"""
def load_PbPb(path):
	db = 0.5
	brange = np.linspace(0,16,17)*db
	pTrange = np.linspace(0.1, 50, 20)

	# Pb Pb
	f = h5py.File(path, 'r')
	key_format = '2760-Pb-Pb-{}-{}'
	data3 = np.zeros([17, 17, 20])
	for i in range(17):
		b1 = db*i	
		for j in range(17):
			b2 = db*j
			bmin = np.min([b1, b2])
			bmax = np.max([b1, b2])
			key = key_format.format(bmin, bmax)
			data3[i,j] = f[key]['dsigma_dpT2_dy'].value[0]
	f.close()

# average over impact parameters for each centrality
cdef numerator(xvec):
	x, y, b = xvec
	s1 = np.sqrt((x+b/2)**2+y**2)
	s2 = np.sqrt((x-b/2)**2+y**2)
	y1 = interpn((brange,brange), data, [s1,s2], method='linear', bounds_error=False, fill_value=None)
	return TA(208, s1)*TA(208, s2)*y1

cdef denominator(xvec):
	x, y, b = xvec
	s1 = np.sqrt((x+b/2)**2+y**2)
	s2 = np.sqrt((x-b/2)**2+y**2)
	return TA(208, s1)*TA(208, s2)
"""






