/*
 * BoxLib to NumPy routines
 *
 * This file gets inserted into the pyboxlib.pyf file at compile time.
 */

// fortran prototypes
void multifab_as_numpy_f(int *mfid, int *nbox, void *ptr, int *nx, int *ny, int *nz, int *nc);
void lmultifab_as_numpy_f(int *mfid, int *nbox, void *ptr, int *nx, int *ny, int *nz, int *nc);

/*
{"multifab_array", multifab_as_numpy, METH_VARARGS, "Return NumPy array associated with a BoxLib multifab."},
{"lmultifab_array", lmultifab_as_numpy, METH_VARARGS, "Return NumPy array associated with a BoxLib lmultifab."},
*/

PyObject *
multifab_as_numpy (PyObject * self, PyObject * args)
{
  int mfid, nbox, n1, n2, n3, n4;
  double *ptr;

  PyObject *arr = NULL;
  int ndim = 4;
  npy_intp dims[4];

  if (!PyArg_ParseTuple (args, "ii", &mfid, &nbox))
    return NULL;

  multifab_as_numpy_f(&mfid, &nbox, &ptr, &n1, &n2, &n3, &n4);

  dims[0] = n1;
  dims[1] = n2;
  dims[2] = n3;
  dims[3] = n4;

  arr = PyArray_NewFromDescr(&PyArray_Type,
                             PyArray_DescrFromType(NPY_DOUBLE), ndim, dims, NULL,
                             ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL);

  Py_INCREF(arr);
  return arr;
}

PyObject *
lmultifab_as_numpy (PyObject * self, PyObject * args)
{
  int mfid, nbox, n1, n2, n3, n4;
  double *ptr;

  PyObject *arr = NULL;
  int ndim = 4;
  npy_intp dims[4];

  if (!PyArg_ParseTuple (args, "ii", &mfid, &nbox))
    return NULL;

  lmultifab_as_numpy_f(&mfid, &nbox, &ptr, &n1, &n2, &n3, &n4);

  dims[0] = n1;
  dims[1] = n2;
  dims[2] = n3;
  dims[3] = n4;

  arr = PyArray_NewFromDescr(&PyArray_Type,
                             PyArray_DescrFromType(NPY_INT), ndim, dims, NULL,
                             ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL);


  Py_INCREF(arr);
  return arr;
}
