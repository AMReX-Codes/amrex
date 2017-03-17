%module warpxC

%{
#include <WarpXConst.H>
#include <WarpXWrappers.h>

#define SWIG_FILE_WITH_INIT
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
%}
%include "numpy.i"
%init %{
import_array();
%}

%include "../Source/WarpXConst.H"

// For amrex_init(int argc, char *argv[]);
%include <argcargv.i>
%apply (int ARGC, char **ARGV) { (int argc, char *argv[]) }

%rename (addNParticles) wrapped_addNParticles;
%ignore addNParticles;

%include "../Source/WarpXWrappers.h"

%apply (int DIM1, double* IN_ARRAY1) {(int lenx, double* x),
                                      (int leny, double* y),
                                      (int lenz, double* z),
                                      (int lenvx, double* vx),
                                      (int lenvy, double* vy),
                                      (int lenvz, double* vz)}
%apply (int DIM1, int DIM2, double* IN_FARRAY2 ) {(int lena, int nattr, double* attr)} // Note Fortran ordering

%exception wrapped_addNParticles {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{
void wrapped_addNParticles(int speciesnumber, int lenx, double* x, int leny, double* y, int lenz, double* z,
                            int lenvx, double* vx, int lenvy, double* vy, int lenvz, double* vz,
                            int lena, int nattr, double* attr, int uniqueparticles) {
    if (lenx != leny || lenx != lenz || lenx != lenvx || lenx != lenvy || lenx != lenvz || lenx != lena ) {
        PyErr_Format(PyExc_ValueError,
                     "Lengths of arrays must be the same, (%d, %d, %d, %d, %d, %d, %d)",
                     lenx, leny, lenz, lenvx, lenvy, lenvz, lena);
        return;
    }
    addNParticles(speciesnumber, lenx, x, y, z, vx, vy, vz, nattr, attr, uniqueparticles);
}
%}

%include classwrapper.i
