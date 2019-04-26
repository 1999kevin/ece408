from numpy import *
from scipy.sparse import *
from scipy import *

import numpy

def map_bytes_to_dtype(bytes, int_type=False):
    # Should use a dict instead
    if (bytes == 4) and (not int_type):
        return "float32"
    elif (bytes == 8) and (not int_type):
        return "float64"
    elif (bytes == 4) and int_type:
        return "int32"
    elif (bytes == 8) and int_type:
        return "int64"
    else:
        raise ValueError("bytes must be 4 or 8 but {0} passed to function".format(bytes))

    
def import_full_c(filename):
    with open(filename, 'r') as fid:
        bytes = int(fromfile(fid, dtype=int32, count=1))
        m = int(fromfile(fid, dtype=int32, count=1))
        n = int(fromfile(fid, dtype=int32, count=1))

        x = reshape(fromfile(fid, dtype=map_bytes_to_dtype(bytes)), (m, n), order='F')
        
        return x

    
def import_full_r(filename):
    with open(filename, 'r') as fid:
        bytes = int(fromfile(fid, dtype=int32, count=1))
        m = int(fromfile(fid, dtype=int32, count=1))
        n = int(fromfile(fid, dtype=int32, count=1))

        x = fromfile(fid, dtype=map_bytes_to_dtype(bytes))
        x.shape = (m, n)
        return x

    
def import_vector(filename):
    with open(filename, 'r') as fid:
        bytes = int(fromfile(fid, dtype=int32, count=1))
        n = int(fromfile(fid, dtype=int32, count=1))

        x = fromfile(fid, dtype=map_bytes_to_dtype(bytes))
        assert(len(x) == n)
        return x

    
def import_ivector(filename):
    with open(filename, 'r') as fid:
        bytes = int(fromfile(fid, dtype=int32, count=1))
        n = int(fromfile(fid, dtype=int32, count=1))

        x = fromfile(fid, dtype=map_bytes_to_dtype(bytes, int_type=True))
        assert(len(x) == n)
        return x

    
def import_rcs(filename):
    with open(filename, 'r') as fid:
        bytes = int(fromfile(fid, dtype=int32, count=1))
        m = int(fromfile(fid, dtype=int32, count=1))
        n = int(fromfile(fid, dtype=int32, count=1))
        N = int(fromfile(fid, dtype=int32, count=1))

        v = fromfile(fid, dtype=map_bytes_to_dtype(bytes), count=N)
        j = fromfile(fid, dtype=int32, count=N)
        r = fromfile(fid, dtype=int32, count=m+1)

        return csr_matrix((v, j, r), shape=(m, n))


def import_diag(filename):
    with open(filename, 'r') as fid:
        bytes = int(fromfile(fid, dtype=int32, count=1))
        m = int(fromfile(fid, dtype=int32, count=1))
        n = int(fromfile(fid, dtype=int32, count=1))
        v = fromfile(fid, dtype=map_bytes_to_dtype(bytes), count=min(m,n))

        return spdiags(v, 0, m, n)

    
################################################################################


def export_full_r(filename, x, dtype=double):
    with open(filename, 'w') as fid:
        m, n = x.shape

        z = array([numpy.dtype(dtype).itemsize], dtype='int32')
        z.tofile(fid)
        
        z = array([m, n], dtype='int32')
        z.tofile(fid)

        if (x.dtype != dtype):
            x_dtype = array(x, dtype=dtype)
            x_dtype.tofile(fid)
        else:
            x.tofile(fid)
    
    
def export_vector(filename, x, dtype=double):
    with open(filename, 'w') as fid:
        n = len(x)

        z = array([numpy.dtype(dtype).itemsize], dtype='int32')
        z.tofile(fid)
        
        z = array([n], dtype='int32')
        z.tofile(fid)

        if (x.dtype != dtype):
            x_dtype = array(x, dtype=dtype)
            x_dtype.tofile(fid)
        else:
            x.tofile(fid)

            
def export_ivector(filename, x, dtype=int64):
    export_vector(filename, asarray(x, dtype=dtype), dtype=dtype)


def export_rcs(filename, A, dtype=double):
    with open(filename, 'w') as fid:
        A = csr_matrix(A, dtype=dtype)

        m, n = A.shape
        N = A.nnz

        v = A.data
        j = A.indices
        r = A.indptr
        
        z = array([numpy.dtype(dtype).itemsize, m, n, N], dtype='int32')
        z.tofile(fid)

        assert(j.dtype == 'int32')
        assert(r.dtype == 'int32')
        
        v.tofile(fid)
        j.tofile(fid)
        r.tofile(fid)
        

def export_diag(filename, A, dtype=double):
    with open(filename, 'w') as fid:
        m, n = A.shape
        v = asarray(A.diagonal(), dtype=dtype)
        
        z = array([numpy.dtype(dtype).itemsize, m, n], dtype=int32)
        z.tofile(fid)

        v.tofile(fid)


def export_r_filter_new(filename, H, dtype=double):
    with open(filename, 'w') as fid:
        header = array([numpy.dtype(dtype).itemsize, H.ndim], dtype="int32")
        n_log = array(H.shape, dtype="int32")
        header.tofile(fid)
        n_log.tofile(fid)

        H_dtype = asarray(H, dtype=dtype)
        H_dtype.tofile(fid)
