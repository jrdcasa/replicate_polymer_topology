import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "calc_dist.c":
  void c_setup_hist_bond(double deltaB, int maxBinB, double* bdist)
  void c_setup_hist_angle(double deltaA, int maxBinA, double* adist)
  int c_bondDist(int natoms, int nbonds, int dim2, int* bl, double* x, double* y, double* z, int* bhist)
  int c_angleDist(int natoms, int nangles, int dim3, int* al, double* x, double* y, double* z, int* ahist)
  int c_angleDist_omp(int natoms, int nangles, int dim3, int* al, double* x, double* y, double* z, int* ahist)
  void c_setup_hist_dih(double deltaD, int maxbinD, double* ddist)
  int c_dihDist(int natoms, int ndih,int dim4, int* dl, double* x, double* y, double* z, int* dihhist)
  int c_dihDistFlory(int natoms, int ndih,int dim4, int* dl, double* x, double* y, double* z, int* dihhist)
  int c_tacticity(int natoms, int ndih,int dim4, int* dl, double* x, double* y, double* z, int* dihhist)
  double dihedral(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4)

@cython.boundscheck(False)
@cython.wraparound(False)

####################################################################################  
def setup_hist_bondC(double deltaB, int maxBinB, np.ndarray[double,ndim=1] bdist):

  c_setup_hist_bond(deltaB,maxBinB, &bdist[0])
    
  return None

####################################################################################
def bondDistC(np.ndarray[int,ndim=2, mode="c"] bl, np.ndarray[double,ndim=1] x,\
              np.ndarray[double,ndim=1] y, np.ndarray[double,ndim=1] z, np.ndarray[int,ndim=1] bhist ):

  cdef int nbonds = len(bl)
  cdef int natoms = len(x)
  cdef int dim2 = 2
  cdef int iserror
  
  #The return value is 1 if the calculation is done without problems. If any problem is
  #detected the return value is 0.
  
  iserror = c_bondDist(natoms, nbonds, dim2, &bl[0,0], &x[0], &y[0], &z[0], &bhist[0])
  
  return iserror

####################################################################################  
def setup_hist_angleC(double deltaA, int maxBinA, np.ndarray[double, ndim=1] adist):

  c_setup_hist_angle(deltaA,maxBinA, &adist[0])
  
  return None

####################################################################################
def angleDistC(np.ndarray[int,ndim=2, mode="c"] al, np.ndarray[double,ndim=1] x,\
               np.ndarray[double,ndim=1] y, np.ndarray[double,ndim=1] z, np.ndarray[int,ndim=1] ahist):

  cdef int nangles = len(al)
  cdef int natoms = len(x)
  cdef int dim3 = 3
  cdef int iserror

  iserror = c_angleDist(natoms, nangles, dim3, &al[0,0], &x[0], &y[0], &z[0], &ahist[0])

  return iserror
  
####################################################################################
def angleDistC_omp(np.ndarray[int,ndim=2, mode="c"] al, np.ndarray[double,ndim=1] x,\
               np.ndarray[double,ndim=1] y, np.ndarray[double,ndim=1] z, np.ndarray[int,ndim=1] ahist):

  cdef int nangles = len(al)
  cdef int natoms = len(x)
  cdef int dim3 = 3
  cdef int iserror

  iserror = c_angleDist_omp(natoms, nangles, dim3, &al[0,0], &x[0], &y[0], &z[0], &ahist[0])

  return iserror  

####################################################################################  
def setup_hist_dihC(double deltaD, int maxBinD, np.ndarray[double, ndim=1] ddist):

  c_setup_hist_dih(deltaD, maxBinD, &ddist[0])
  
  return None

####################################################################################
def dihDistC(np.ndarray[int,ndim=2, mode="c"] dl, np.ndarray[double,ndim=1] x,\
             np.ndarray[double,ndim=1] y, np.ndarray[double,ndim=1] z, np.ndarray[int,ndim=1] dhist):

  cdef int ndih = len(dl)
  cdef int natoms = len(x)
  cdef int dim4 = 3
  cdef int iserror

  iserror = c_dihDist(natoms, ndih, dim4, &dl[0,0], &x[0], &y[0], &z[0], &dhist[0])

  return iserror

####################################################################################
def dihDistFloryC(np.ndarray[int, ndim=2, mode='c'] dl, np.ndarray[double,ndim=1] x,\
                  np.ndarray[double,ndim=1] y, np.ndarray[double,ndim=1] z, np.ndarray[int,ndim=1] dhist):


  cdef int ndih = len(dl)
  cdef int natoms = len(x)
  cdef int dim4 = 3
  cdef int iserror

  iserror = c_dihDistFlory(natoms, ndih, dim4, &dl[0,0], &x[0], &y[0], &z[0], &dhist[0])

  return iserror

####################################################################################
def dihedralC(double x1, double y1, double z1,
              double x2, double y2, double z2,
              double x3, double y3, double z3,
              double x4, double y4, double z4):

  cdef double a

  d = dihedral (x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
  return d

####################################################################################
def tactDistC(np.ndarray[int, ndim=2, mode='c'] dl, np.ndarray[double,ndim=1] x,\
              np.ndarray[double,ndim=1] y, np.ndarray[double,ndim=1] z, np.ndarray[int,ndim=1] dhist):


  cdef int ndih = len(dl)
  cdef int natoms = len(x)
  cdef int dim4 = 3
  cdef int iserror

  iserror = c_tacticity(natoms, ndih, dim4, &dl[0,0], &x[0], &y[0], &z[0], &dhist[0])

  return iserror
