%{

#include <WarpX.H>

using namespace amrex;

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

%include "../../amrex/Src/Base/AMReX_Array.H"

%extend amrex::Array {
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

// This is needed so that swig knows the amrex::Real is just a float or double
%include "../../amrex/Src/Base/AMReX_REAL.H"

%template(arrayReal) amrex::Array<amrex::Real>;
%template(arrayInt) amrex::Array<int>;
%template(arrayGeometry) amrex::Array<amrex::Geometry>;

%extend amrex::Array<amrex::Real> {
    PyObject *get_real() {
      // Get the data as a writeable numpy array, directly accessing the memory.
      PyObject *arr = 0;
      npy_intp dims[1];
      dims[0] = self->std::vector<amrex::Real>::size();
      arr = PyArray_NewFromDescr(&PyArray_Type,
                                 PyArray_DescrFromType(NPY_DOUBLE), 1, dims, NULL,
                                 self->dataPtr(), NPY_FORTRANORDER|NPY_ARRAY_WRITEABLE, NULL);
      Py_INCREF(arr);
      return arr;
    }
}
%extend amrex::Array<int> {
    PyObject *get_int() {
      // Get the data as a writeable numpy array, directly accessing the memory.
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

// This is needed by swig to define the macro D_DECL when including AMReX_IntVect.H
%include "../../amrex/Src/Base/AMReX_SPACE.H"

// This include can only be done with the modified AMReX_IntVect.H file that hides the ref-qualifiers from swig.
%include "../../amrex/Src/Base/AMReX_IntVect.H"
// Save this code, which is needed if AMReX_IntVect.H cannot be included.
// AMReX_IntVect.H uses ref-qualifiers in the lines setting up getVect that cannot be parsed by swig.
/*
class amrex::IntVect {
public:

    IntVect(int i, int j, int k);
    IntVect(const amrex::IntVect& rhs);

    IntVect& shift(int coord, int s);

    // static functions
    static const amrex::IntVect& TheZeroVector();
    static const amrex::IntVect& TheUnitVector();
    static const amrex::IntVect& TheNodeVector();
    static const amrex::IntVect& TheCellVector();
};
*/

%extend amrex::IntVect {
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
    int __cmp__(const amrex::IntVect* other){
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

//%include "../../amrex/Src/Base/AMReX_Box.H"
%include "../../amrex/Src/Base/AMReX_BaseFab.H"

%template(BaseFabReal) amrex::BaseFab<amrex::Real>;

%extend amrex::BaseFab<amrex::Real> {
    PyObject * get(size_t n=0) {
        PyObject *arr = 0;
        npy_intp dims[BL_SPACEDIM];
        amrex::IntVect size = self->box().size();
        dims[0] = size[0];
        dims[1] = size[1];
        dims[2] = size[2];
        arr = PyArray_NewFromDescr(&PyArray_Type,
                                   PyArray_DescrFromType(NPY_DOUBLE), BL_SPACEDIM, dims, NULL,
                                   self->dataPtr(n), NPY_FORTRANORDER|NPY_ARRAY_WRITEABLE, NULL);
        Py_INCREF(arr);
        return arr;
    }
}

//%include "../../amrex/Src/Base/AMReX_FArrayBox.H"
//%include "../../amrex/Src/Base/AMReX_BoxArray.H"

%include "../../amrex/Src/Base/AMReX_FabArray.H"

%extend amrex::FabArray {
    FAB& __getitem__ (size_t i)
    {
        BL_ASSERT(i >= 0);
        BL_ASSERT(i < self->size());
        return self->get(i);
    }

    //
    // Same as above, except acts on const Array's.
    //
    const FAB& __getitem__ (size_t i) const
    {
        BL_ASSERT(i >= 0);
        BL_ASSERT(i < self->size());
        return self->get(i);
    }
}

%template(FabArrayFArrayBox) amrex::FabArray< amrex::FArrayBox >;

%include "../../amrex/Src/Base/AMReX_MultiFab.H"

%template(arrayMultifab) amrex::Array< std::unique_ptr<amrex::MultiFab> >;
%template(arrayarrayMultifab) amrex::Array<amrex::Array< std::unique_ptr<amrex::MultiFab> > >;

#if (BL_SPACEDIM > 2)
// GetDLogA is only defined with BL_SPACEDIM <= 2
%ignore GetDLogA;
#endif

%include "../../amrex/Src/Base/AMReX_Geometry.H"

%template(arrayBoxArray) amrex::Array<amrex::BoxArray>;

%include "../../amrex/Src/Particle/AMReX_Particles.H"

// Becuase of an apparent problem in swig, the macro around the wrapping of tile_size gives an error during compilation
%ignore tile_size;

// Swig doesn't handle the unique_ptr return value causing an error during compilation
%ignore GetChargeDensity;

// Wrapping WarpXParIter fails since swig doesn't handle the alias SoA properly causing an error during compilation
%ignore WarpXParIter;

%template("WarpXParticleContainerBase") amrex::ParticleContainer<0,0,PIdx::nattribs>;

%include "../Source/ParticleContainer.H"
%include "../Source/WarpXParticleContainer.H"

%extend WarpXParticleContainer {
    int getNattribs() {
        return PIdx::nattribs;
    }
    PyObject * getLocations() {
        amrex::Array<amrex::Real> result(0);
        self->GetParticleLocations(result);
        npy_intp dims[2] = {BL_SPACEDIM, self->TotalNumberOfParticles()};
        PyObject *arr = PyArray_NewFromDescr(&PyArray_Type,
                                             PyArray_DescrFromType(NPY_DOUBLE), 2, dims, NULL,
                                             result.dataPtr(), NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE, NULL);
        Py_INCREF(arr);
        return arr;
    }
    PyObject * getData(int start_comp, int num_comp) {
        amrex::Array<amrex::Real> result(0);
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
        amrex::Array<amrex::Real> result(0);
        self->GetParticleData(result, 0, num_comp);
        npy_intp dims[2] = {num_comp, self->TotalNumberOfParticles()};
        PyObject *arr = PyArray_NewFromDescr(&PyArray_Type,
                                             PyArray_DescrFromType(NPY_DOUBLE), 2, dims, NULL,
                                             result.dataPtr(), NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE, NULL);
        Py_INCREF(arr);
        return arr;
    }
};

%include "../../amrex/Src/AmrCore/AMReX_AmrCore.H"
%include "../Source/WarpX.H"

