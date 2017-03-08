%{

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

#include <Array.H>
#include <IntVect.H>

#include <Particles.H>

#include <Box.H>
#include <FArrayBox.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <Geometry.H>

#include <WarpX.H>

%}

%inline %{
  std::ifstream & open_ifstream(const char *filename) {
    std::ifstream *infile = new std::ifstream(filename);
    return *infile;
  }
%}

%inline %{
  std::ofstream & open_ofstream(const char *filename) {
    std::ofstream *outfile = new std::ofstream(filename);
    return *outfile;
  }
%}

%inline %{
  void close_ofstream(std::ofstream& str) {
    str.close();
  }
%}

%inline %{
  void close_ifstream(std::ifstream & str) {
    str.close();
  }
%}

%include <std_string.i>
%rename(__str__) display;

typedef double Real;

%include "../../amrex/Src/C_BaseLib/Array.H"

%extend Array {
    T& __getitem__ (size_t i)
    {
        BL_ASSERT(i >= 0);
        BL_ASSERT(i < self->size());
        return self->std::vector<T>::operator[](i);
    }

    //
    // Same as above, except acts on const Array's.
    //
    const T& __getitem__ (size_t i) const
    {
        BL_ASSERT(i >= 0);
        BL_ASSERT(i < self->size());
        return self->std::vector<T>::operator[](i);
    }
}

%template(arrayReal) Array<Real>;
%template(arrayInt) Array<int>;
%template(arrayGeometry) Array<Geometry>;

%extend Array<Real> {
    PyObject *get_real() {
      PyObject *arr = 0;
      npy_intp dims[1];
      dims[0] = self->std::vector<Real>::size();
      arr = PyArray_NewFromDescr(&PyArray_Type,
                                 PyArray_DescrFromType(NPY_DOUBLE), 1, dims, NULL,
                                 self->dataPtr(), NPY_FORTRANORDER|NPY_ARRAY_WRITEABLE, NULL);
      Py_INCREF(arr);
      return arr;
    }
}
%extend Array<int> {
    PyObject *get_int() {
      PyObject *arr = 0;
      npy_intp dims[1];
      dims[0] = self->std::vector<int>::size();
      arr = PyArray_NewFromDescr(&PyArray_Type,
                                 PyArray_DescrFromType(NPY_INT), 1, dims, NULL,
                                 self->dataPtr(), NPY_FORTRANORDER|NPY_ARRAY_WRITEABLE, NULL);
      Py_INCREF(arr);
      return arr;
    }
}

// Note that IntVect.H cannot be directly included since swig cannot parse the line setting up "const int* getVect".
//%include "../../amrex/Src/C_BaseLib/IntVect.H"
class IntVect {
public:

    IntVect(int i, int j, int k);
    IntVect(const IntVect& rhs);

    IntVect& shift(int coord, int s);

    // static functions
    static const IntVect& TheZeroVector();
    static const IntVect& TheUnitVector();
    static const IntVect& TheNodeVector();
    static const IntVect& TheCellVector();
};

%extend IntVect {
    //void writeOn(std::ofstream *os){
    //    *os << *self;
    //}
    //void read(std::ifstream* ifs){
    //    *ifs >> *self;
    //}
    int __getitem__(int index){
        if (index < 0 || index >= BL_SPACEDIM) {
          // SWIG_SetErrorMsg(PyExc_IndexError,"Index out of bounds\n");
          return 0;
          }
        return (*self)[index];
    }
    int __len__() volatile { return BL_SPACEDIM; }
    void __setitem__(int index,int val){
        if (index < 0 || index >= BL_SPACEDIM) {
          //SWIG_SetErrorMsg(PyExc_IndexError,"Index out of bounds\n");
         } else {
          (*self).setVal(index,val);
        }
    }
    int __cmp__(const IntVect* other){
        if( (*self) == (*other) ) {
            return 0;
        }
        if( (*self)  <= (*other) ) {
            return -1;
        }
        return 1;
    }
    PyObject *get() {
      PyObject *arr = 0;
      npy_intp dims[1];
      dims[0] = BL_SPACEDIM;
      arr = PyArray_NewFromDescr(&PyArray_Type,
                                 PyArray_DescrFromType(NPY_INT), 1, dims, NULL,
                                 self->getVect(), NPY_FORTRANORDER|NPY_ARRAY_WRITEABLE, NULL);
      Py_INCREF(arr);
      return arr;
    }

    //std::string display() {
    //    std::ostringstream str;
    //        str << *self;
    //    return str.str();
    //}
}

//%include "../../amrex/Src/C_BaseLib/Box.H"
//%include "../../amrex/Src/C_BaseLib/FArrayBox.H"
//%include "../../amrex/Src/C_BaseLib/BoxArray.H"
//%include "../../amrex/Src/C_BaseLib/MultiFab.H"

//#if (BL_SPACEDIM > 2)
%ignore GetDLogA;
//#endif

%include "../../amrex/Src/C_BaseLib/Geometry.H"

%template(arrayBoxArray) Array<BoxArray>;

%include "../../amrex/Src/C_ParticleLib/Particles.H"

//%template("WarpXParticleBase") Particle<PIdx::nattribs,0>;
%template("WarpXParticleContainerBase") ParticleContainer<PIdx::nattribs,0,std::vector<Particle<PIdx::nattribs,0> > >;

%ignore GetChargeDensity;

%include "../Source/ParticleContainer.H"
%include "../Source/WarpXParticleContainer.H"

%extend WarpXParticleContainer {
    PyObject * getLocations() {
        Array<Real> result(0);
        self->GetParticleLocations(result);
        npy_intp dims[2] = {BL_SPACEDIM, self->TotalNumberOfParticles()};
        PyObject *arr = PyArray_NewFromDescr(&PyArray_Type,
                                             PyArray_DescrFromType(NPY_DOUBLE), 2, dims, NULL,
                                             result.dataPtr(), NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE, NULL);
        Py_INCREF(arr);
        return arr;
    }
    PyObject * getData(int start_comp, int num_comp) {
        Array<Real> result(0);
        self->GetParticleData(result, start_comp, num_comp);
        npy_intp dims[2] = {num_comp, self->TotalNumberOfParticles()};
        PyObject *arr = PyArray_NewFromDescr(&PyArray_Type,
                                             PyArray_DescrFromType(NPY_DOUBLE), 2, dims, NULL,
                                             result.dataPtr(), NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE, NULL);
        Py_INCREF(arr);
        return arr;
    }
    PyObject * getAllData() {
        int num_comp = PIdx::nattribs;
        Array<Real> result(0);
        self->GetParticleData(result, 0, num_comp);
        npy_intp dims[2] = {num_comp, self->TotalNumberOfParticles()};
        PyObject *arr = PyArray_NewFromDescr(&PyArray_Type,
                                             PyArray_DescrFromType(NPY_DOUBLE), 2, dims, NULL,
                                             result.dataPtr(), NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE, NULL);
        Py_INCREF(arr);
        return arr;
    }
};

%include "../../amrex/Src/C_AmrCoreLib/AmrCore.H"
%include "../Source/WarpX.H"
