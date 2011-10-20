/*
 * BoxLib to NumPy routines
 *
 * This file gets inserted into the pyboxlib.pyf file at compile time.
 */

// fortran prototypes
void multifab_as_numpy_f(int *mfid, int *nbox, void *ptr, int *nx, int *ny, int *nz, int *nc);
void lmultifab_as_numpy_f(int *mfid, int *nbox, void *ptr, int *nx, int *ny, int *nz, int *nc);

/*
divert(1)dnl
{"multifab_array", multifab_as_numpy, METH_VARARGS, "Return NumPy array associated with a BoxLib multifab."},
{"lmultifab_array", lmultifab_as_numpy, METH_VARARGS, "Return NumPy array associated with a BoxLib lmultifab."},
divert(0)dnl
*/

PyObject *
multifab_as_numpy (PyObject * self, PyObject * args)
{
  int mfid, nbox, nx, ny, nz, nc;
  double *ptr;

  PyArrayObject *arr = NULL;
  int ndim = 4;
  npy_intp dims[4];

  if (!PyArg_ParseTuple (args, "ii", &mfid, &nbox))
    return NULL;

  multifab_as_numpy_f(&mfid, &nbox, &ptr, &nx, &ny, &nz, &nc);

  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nz;
  dims[3] = nc;

  arr = PyArray_NewFromDescr(&PyArray_Type,
                             PyArray_DescrFromType(NPY_DOUBLE), ndim, dims, NULL,
                             ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL);


  Py_INCREF(arr);
  return arr;
}

PyObject *
lmultifab_as_numpy (PyObject * self, PyObject * args)
{
  int mfid, nbox, nx, ny, nz, nc;
  double *ptr;

  PyArrayObject *arr = NULL;
  int ndim = 4;
  npy_intp dims[4];

  if (!PyArg_ParseTuple (args, "ii", &mfid, &nbox))
    return NULL;

  lmultifab_as_numpy_f(&mfid, &nbox, &ptr, &nx, &ny, &nz, &nc);

  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nz;
  dims[3] = nc;

  arr = PyArray_NewFromDescr(&PyArray_Type,
                             PyArray_DescrFromType(NPY_INT), ndim, dims, NULL,
                             ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL);


  Py_INCREF(arr);
  return arr;
}
