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

void pybl_layout_local(void*, int, int*);

PyObject *
pyblw_layout_local (PyObject * self, PyObject * args)
{
  long ptr;
  int n, local;
  if (!PyArg_ParseTuple (args, "li", &ptr, &n))
    return NULL;
  pybl_layout_local((void*)ptr,n,&local);
  if (local) {
    Py_RETURN_TRUE;
  }
  Py_RETURN_FALSE;
}

////////////////////////////////////////////////////////////////////////////////

void pybl_layout_get_box(void*, void*, int);

PyObject *
pyblw_layout_get_box (PyObject * self, PyObject * args)
{
  long ptr;
  int n, bx[7];
  if (!PyArg_ParseTuple (args, "li", &ptr, &n))
    return NULL;
  pybl_layout_get_box((void*)ptr,(void*)bx,n);

  npy_intp dims[1] = { bx[0] };
  PyObject *nplo = PyArray_EMPTY(1, dims, NPY_INT, 0);
  PyObject *nphi = PyArray_EMPTY(1, dims, NPY_INT, 0);

  for (int i=0; i<bx[0]; i++) {
    *(int*)PyArray_GETPTR1(nplo, i) = bx[1+i];
    *(int*)PyArray_GETPTR1(nphi, i) = bx[4+i];
  }

  return Py_BuildValue("(OO)", nplo, nphi);
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

void pybl_lmultifab_create_from_layout(void*,void*);

PyObject *
pyblw_lmultifab_create_from_layout (PyObject * self, PyObject * args)
{
  long laptr, ptr;
  if (!PyArg_ParseTuple (args, "l", &laptr))
    return NULL;
  pybl_lmultifab_create_from_layout((void*)laptr,&ptr);
  return PyLong_FromLong(ptr);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_lmultifab_info(void*,int*,int*,int*,int*);

PyObject *
pyblw_lmultifab_info (PyObject * self, PyObject * args)
{
  long ptr;
  int dim, nboxes, nc, ng;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  pybl_lmultifab_info((void*)ptr,&dim,&nboxes,&nc,&ng);
  return Py_BuildValue("(iiii)", dim, nboxes, nc, ng);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_multifab_create_from_layout(void*,int,int,void*);

PyObject *
pyblw_multifab_create_from_layout (PyObject * self, PyObject * args)
{
  int nc, ng;
  long laptr, ptr;
  if (!PyArg_ParseTuple (args, "lii", &laptr, &nc, &ng))
    return NULL;
  pybl_multifab_create_from_layout((void*)laptr,nc,ng,&ptr);
  return PyLong_FromLong(ptr);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_multifab_print(void*);

PyObject *
pyblw_multifab_print (PyObject * self, PyObject * args)
{
  long ptr;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  pybl_multifab_print((void*)ptr);
  Py_INCREF(Py_None);
  return Py_None;
}

////////////////////////////////////////////////////////////////////////////////

void pybl_multifab_info(void*,int*,int*,int*,int*);

PyObject *
pyblw_multifab_info (PyObject * self, PyObject * args)
{
  long ptr;
  int dim, nboxes, nc, ng;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  pybl_multifab_info((void*)ptr,&dim,&nboxes,&nc,&ng);
  return Py_BuildValue("(iiii)", dim, nboxes, nc, ng);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_multifab_layout(void*,void*);

PyObject *
pyblw_multifab_layout (PyObject * self, PyObject * args)
{
  long ptr, laptr;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  pybl_multifab_layout((void*)ptr,(void*)&laptr);
  return PyLong_FromLong(laptr);
}

////////////////////////////////////////////////////////////////////////////////

void pybl_multifab_fill_boundary(void*);

PyObject *
pyblw_multifab_fill_boundary (PyObject * self, PyObject * args)
{
  long ptr;
  if (!PyArg_ParseTuple (args, "l", &ptr))
    return NULL;
  pybl_multifab_fill_boundary((void*)ptr);
  Py_INCREF(Py_None);
  return Py_None;
}

////////////////////////////////////////////////////////////////////////////////

void pybl_multifab_write(void*, char*, int, char*, int);

PyObject *
pyblw_multifab_write (PyObject * self, PyObject * args)
{
  long ptr;
  char *dname, *hname;
  int dlen, hlen;
  if (!PyArg_ParseTuple (args, "lsisi", &ptr, &dname, &dlen, &hname, &hlen))
    return NULL;
  pybl_multifab_write((void*)ptr, dname, dlen, hname, hlen);
  Py_INCREF(Py_None);
  return Py_None;
}

////////////////////////////////////////////////////////////////////////////////

PyObject *tag_boxes_callback;

// this is called by fortran, and calls python. see tag_boxes.f90
void pyblw_tag_boxes_callback(void *mf, void *tb, double dx, int lev)
{
  PyObject *arglist, *result;
  arglist = Py_BuildValue("lldi", (long)mf, (long)tb, dx, lev);
  result = PyObject_CallObject(tag_boxes_callback, arglist);
  Py_DECREF(arglist);
  if (result==NULL)
    return;
  Py_DECREF(result);
}

void pybl_regrid(void**,void**,int*,int,void(*)(void*,void*,double,int));

PyObject *
pyblw_regrid (PyObject * self, PyObject * args)
{
  PyArrayObject *laptrs, *mfptrs;
  PyObject *tb_cb;
  int nlevs, max_levs;
  if (!PyArg_ParseTuple (args, "O&O&iiO", PyArray_Converter, &laptrs, PyArray_Converter, &mfptrs, &nlevs, &max_levs, &tb_cb))
    return NULL;
  if (!PyCallable_Check(tb_cb)) {
    PyErr_SetString(PyExc_TypeError, "parameter must be callable");
    return NULL;
  }
  tag_boxes_callback = tb_cb;
  pybl_regrid(PyArray_DATA(laptrs),PyArray_DATA(mfptrs),&nlevs,max_levs,pyblw_tag_boxes_callback);
  return PyInt_FromLong(nlevs);
}

////////////////////////////////////////////////////////////////////////////////

void multifab_as_numpy_f(void *cptr, int nbox, int *dims, void *ptr);

PyObject *
multifab_as_numpy (PyObject * self, PyObject * args)
{
  int nbox;
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

  arr = PyArray_NewFromDescr(&PyArray_Type,
                             PyArray_DescrFromType(NPY_DOUBLE), 4, dims, NULL,
                             ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL);

  Py_INCREF(arr);
  return arr;
}

////////////////////////////////////////////////////////////////////////////////

void lmultifab_as_numpy_f(void *cptr, int nbox, int *dims, void *ptr);

PyObject *
lmultifab_as_numpy (PyObject * self, PyObject * args)
{
  int nbox;
  long cptr;
  double *ptr;

  PyObject *arr = NULL;
  npy_intp dims[4];
  int idims[4];

  if (!PyArg_ParseTuple (args, "li", &cptr, &nbox))
    return NULL;

  lmultifab_as_numpy_f((void *) cptr, nbox, idims, &ptr);
  dims[0] = idims[0];
  dims[1] = idims[1];
  dims[2] = idims[2];
  dims[3] = idims[3];

  arr = PyArray_NewFromDescr(&PyArray_Type,
                             PyArray_DescrFromType(NPY_INT), 4, dims, NULL,
                             ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL);

  Py_INCREF(arr);
  return arr;
}

////////////////////////////////////////////////////////////////////////////////

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
  {"layout_local", pyblw_layout_local, METH_VARARGS, ""},
  {"layout_get_box", pyblw_layout_get_box, METH_VARARGS, ""},
  {"lmultifab_create_from_layout", pyblw_lmultifab_create_from_layout, METH_VARARGS, ""},
  {"lmultifab_info", pyblw_lmultifab_info, METH_VARARGS, ""},
  {"multifab_create_from_layout", pyblw_multifab_create_from_layout, METH_VARARGS, ""},
  {"multifab_print", pyblw_multifab_print, METH_VARARGS, ""},
  {"multifab_info", pyblw_multifab_info, METH_VARARGS, ""},
  {"multifab_layout", pyblw_multifab_layout, METH_VARARGS, ""},
  {"multifab_fill_boundary", pyblw_multifab_fill_boundary, METH_VARARGS, ""},
  {"multifab_write", pyblw_multifab_write, METH_VARARGS, ""},
  {"multifab_as_numpy", multifab_as_numpy, METH_VARARGS, "Return NumPy array associated with a BoxLib multifab."},
  {"lmultifab_as_numpy", lmultifab_as_numpy, METH_VARARGS, "Return NumPy array associated with a BoxLib multifab."},
  {"regrid", pyblw_regrid, METH_VARARGS, ""},
  {NULL, NULL},
};

void
initfcboxlib(void)
{
  Py_InitModule("fcboxlib", fcboxlib_methods);
  import_array();
}
