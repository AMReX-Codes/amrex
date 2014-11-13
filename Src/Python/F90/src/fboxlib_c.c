/*
 * BoxLib Python wrappers.
 */

#include <stdlib.h>

#define PY_ARRAY_UNIQUE_SYMBOL _PYBOXLIB_ARRAY_

#include "Python.h"
#include "numpy/arrayobject.h"

////////////////////////////////////////////////////////////////////////////////

void pybl_open(void);

PyObject *
pyblw_open (PyObject * self, PyObject * args)
{
  pybl_open();
  Py_INCREF(Py_None);
  return Py_None;
}

////////////////////////////////////////////////////////////////////////////////

void pybl_close(void);

PyObject *
pyblw_close (PyObject * self, PyObject * args)
{
  pybl_close();
  Py_INCREF(Py_None);
  return Py_None;
}

////////////////////////////////////////////////////////////////////////////////

void pybl_mpi_rank(int*);

PyObject *
pyblw_mpi_rank (PyObject * self, PyObject * args)
{
  int rank;
  pybl_mpi_rank(&rank);
  return PyInt_FromLong(rank);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_mpi_size(int*);

PyObject *
pyblw_mpi_size (PyObject * self, PyObject * args)
{
  int size;
  pybl_mpi_size(&size);
  return PyInt_FromLong(size);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_mpi_reduce_max(double*);

PyObject *
pyblw_mpi_reduce_max (PyObject * self, PyObject * args)
{
  double r;
  if (!PyArg_ParseTuple (args, "d", &r))
    return NULL;
  pybl_mpi_reduce_max(&r);
  return PyFloat_FromDouble(r);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_boxarray_create_from_boxes(int*,int,int,void*);

PyObject *
pyblw_boxarray_create_from_boxes (PyObject * self, PyObject * args)
{
  PyArrayObject *boxes;
  int nboxes, dim;
  long ptr;
  if (!PyArg_ParseTuple (args, "O&ii", PyArray_Converter, &boxes, &nboxes, &dim))
    return NULL;
  pybl_boxarray_create_from_boxes((int*)PyArray_DATA(boxes),nboxes,dim,&ptr);
  return PyLong_FromLong(ptr);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_boxarray_print(void*);

PyObject *
pyblw_boxarray_print (PyObject * self, PyObject * args)
{
  long ptr;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  pybl_boxarray_print((void*)ptr);
  Py_INCREF(Py_None);
  return Py_None;
}

////////////////////////////////////////////////////////////////////////////////

void pybl_boxarray_nboxes(void*,int*);

PyObject *
pyblw_boxarray_nboxes (PyObject * self, PyObject * args)
{
  long ptr;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  int nboxes;
  pybl_boxarray_nboxes((void*)ptr,&nboxes);
  return PyInt_FromLong(nboxes);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_boxarray_dim(void*,int*);

PyObject *
pyblw_boxarray_dim (PyObject * self, PyObject * args)
{
  long ptr;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  int dim;
  pybl_boxarray_dim((void*)ptr,&dim);
  return PyInt_FromLong(dim);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_boxarray_maxsize(void*,int*,int);

PyObject *
pyblw_boxarray_maxsize (PyObject * self, PyObject * args)
{
  PyArrayObject *sizes;
  long ptr;
  if (!PyArg_ParseTuple (args, "lO&", &ptr, PyArray_Converter, &sizes))
    return NULL;
  pybl_boxarray_maxsize((void*)ptr,(int*)PyArray_DATA(sizes),PyArray_DIM(sizes,0));
  Py_INCREF(Py_None);
  return Py_None;
}

////////////////////////////////////////////////////////////////////////////////

void pybl_layout_create_from_boxarray(void*,int,int*,void*);

PyObject *
pyblw_layout_create_from_boxarray (PyObject * self, PyObject * args)
{
  PyArrayObject *pmask;
  long baptr, ptr;
  if (!PyArg_ParseTuple (args, "lO&", &baptr, PyArray_Converter, &pmask))
    return NULL;
  pybl_layout_create_from_boxarray((void*)baptr,PyArray_DIM(pmask, 0),(int*)PyArray_DATA(pmask),&ptr);
  return PyLong_FromLong(ptr);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_layout_print(void*);

PyObject *
pyblw_layout_print (PyObject * self, PyObject * args)
{
  long ptr;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  pybl_layout_print((void*)ptr);
  Py_INCREF(Py_None);
  return Py_None;
}

////////////////////////////////////////////////////////////////////////////////

void pybl_layout_nboxes(void*,int*);

PyObject *
pyblw_layout_nboxes (PyObject * self, PyObject * args)
{
  long ptr;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  int nboxes;
  pybl_layout_nboxes((void*)ptr,&nboxes);
  return PyInt_FromLong(nboxes);
}

////////////////////////////////////////////////////////////////////////////////

void multifab_as_numpy_f(void *cptr, int nbox, int *dims, void *ptr);

PyObject *
multifab_as_numpy (PyObject * self, PyObject * args)
{
  int nbox, ndim, nc;
  long cptr;
  double *ptr;

  PyObject *arr = NULL;
  npy_intp dims[4];
  int idims[4];

  if (!PyArg_ParseTuple (args, "liii", &cptr, &nbox, &ndim, &nc))
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




static PyMethodDef fcboxlib_methods[] = {
  {"open", pyblw_open, METH_VARARGS, "Oepn (initialize) BoxLib."},
  {"close", pyblw_close, METH_VARARGS, "Close BoxLib."},
  {"rank", pyblw_mpi_rank, METH_VARARGS, "Return MPI rank."},
  {"size", pyblw_mpi_size, METH_VARARGS, "Return MPI size."},
  {"reduce_max", pyblw_mpi_reduce_max, METH_VARARGS, "Parallel reduce (MAX)."},
  {"boxarray_create_from_boxes", pyblw_boxarray_create_from_boxes, METH_VARARGS, ""},
  {"boxarray_print", pyblw_boxarray_print, METH_VARARGS, ""},
  {"boxarray_nboxes", pyblw_boxarray_nboxes, METH_VARARGS, ""},
  {"boxarray_dim", pyblw_boxarray_dim, METH_VARARGS, ""},
  {"boxarray_maxsize", pyblw_boxarray_maxsize, METH_VARARGS, ""},
  {"layout_create_from_boxarray", pyblw_layout_create_from_boxarray, METH_VARARGS, ""},
  {"layout_print", pyblw_layout_print, METH_VARARGS, ""},
  {"layout_nboxes", pyblw_layout_nboxes, METH_VARARGS, ""},
  {"multifab_as_numpy", multifab_as_numpy, METH_VARARGS, "Return NumPy array associated with a BoxLib multifab."},
  {NULL, NULL},
};

void
initfcboxlib(void)
{
  Py_InitModule("fcboxlib", fcboxlib_methods);
  import_array();
}
