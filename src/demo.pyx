# cython: language_level=3

# CuPy array creation

from libc.stdio cimport printf

cimport cupy as np

N = 16
# Define CuPy array of N bytes
cdef carr = np.ndarray[np.double_t, ndim=1] rr = np.zeros((N,), dtype=np.double)