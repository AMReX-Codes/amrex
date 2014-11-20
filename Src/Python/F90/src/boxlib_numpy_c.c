/*
 * BoxLib to NumPy routines.
 *
 * This is complied into its own Python extension module (instead of
 * being rolled into the libpyfboxlib.so shared library) so that the
 * NumPy API can be properly initialized (see import_array() below).
 * Slightly annoying, but it works.
 */

#include <stdlib.h>

#define PY_ARRAY_UNIQUE_SYMBOL _PYBOXLIB_ARRAY_

#include "Python.h"
#include "numpy/arrayobject.h"

// fortran prototypes
void multifab_as_numpy_f(void *cptr, int nbox, int *dims, void *ptr);

PyObject *
multifab_as_numpy (PyObject * self, PyObject * args)
{
  int ndim;
  long cptr;
  double *ptr;

  PyObject *arr = NULL;
  npy_intp dims[4];
  int idims[4];

  if (!PyArg_ParseTuple (args, "li", &cptr, &nbox))
    return NULL;

  multifab_as_numpy_f((void *) cptr, nbox, idims, &ptr);
  dims[0] = idims[0];
  dims[1] = idims[1];
  dims[2] = idims[2];
  dims[3] = idims[3];
  dims[ndim] = nc;

  arr = PyArray_NewFromDescr(&PyArray_Type,
                             PyArray_DescrFromType(NPY_DOUBLE), ndim+1, dims, NULL,
                             ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL);

  Py_INCREF(arr);
  return arr;
}

/* PyObject * */
/* plotfile_as_numpy (PyObject * self, PyObject * args) */
/* { */
/*   int level, nbox, n1, n2, n3, n4; */
/*   long cptr; */
/*   double *ptr; */

/*   PyObject *arr = NULL; */
/*   int ndim = 4; */
/*   npy_intp dims[4]; */

/*   if (!PyArg_ParseTuple (args, "lii", &cptr, &level, &nbox)) */
/*     return NULL; */

/*   plotfile_as_numpy_f((void *) cptr, &level, &nbox, &ptr, &n1, &n2, &n3, &n4); */

/*   dims[0] = n1; */
/*   dims[1] = n2; */
/*   dims[2] = n3; */
/*   dims[3] = n4; */

/*   arr = PyArray_NewFromDescr(&PyArray_Type, */
/*                              PyArray_DescrFromType(NPY_DOUBLE), ndim, dims, NULL, */
/*                              ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL); */

/*   Py_INCREF(arr); */
/*   return arr; */
/* } */

/* PyObject * */
/* lmultifab_as_numpy (PyObject * self, PyObject * args) */
/* { */
/*   int mfid, nbox, n1, n2, n3, n4; */
/*   double *ptr; */

/*   PyObject *arr = NULL; */
/*   int ndim = 4; */
/*   npy_intp dims[4]; */

/*   if (!PyArg_ParseTuple (args, "ii", &mfid, &nbox)) */
/*     return NULL; */

/*   lmultifab_as_numpy_f(&mfid, &nbox, &ptr, &n1, &n2, &n3, &n4); */

/*   dims[0] = n1; */
/*   dims[1] = n2; */
/*   dims[2] = n3; */
/*   dims[3] = n4; */

/*   arr = PyArray_NewFromDescr(&PyArray_Type, */
/*                              PyArray_DescrFromType(NPY_INT), ndim, dims, NULL, */
/*                              ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL); */


/*   Py_INCREF(arr); */
/*   return arr; */
/* } */

static PyMethodDef blnpy_methods[] = {
  {"multifab_as_numpy", multifab_as_numpy, METH_VARARGS,
   "Return NumPy array associated with a BoxLib multifab."},
  /* {"plotfile_as_numpy", plotfile_as_numpy, METH_VARARGS, */
  /*  "Return NumPy array associated with a BoxLib plotfile fab."}, */
  /* {"lmultifab_array", lmultifab_as_numpy, METH_VARARGS,  */
  /*  "Return NumPy array associated with a BoxLib lmultifab."}, */
  {NULL, NULL},
};

void
initblnpy(void)
{
  Py_InitModule("blnpy", blnpy_methods);
  import_array();
}
