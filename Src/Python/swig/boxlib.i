#ifdef DIM1
%module bl1
#endif
#ifdef DIM2
%module bl2
#endif
 #ifdef DIM3
%module bl3
#endif
%{
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <IntVect.H>
#include <Box.H>
#include <FArrayBox.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <VisMF.H>
#include <Array.H>
#include <Geometry.H>
#include <ParallelDescriptor.H>

#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
    import_array();
%}

#ifdef DIM1
%constant int BL_DIM = 1;
#endif
#ifdef DIM2
%constant int BL_DIM = 2;
#endif
#ifdef DIM3
%constant int BL_DIM = 3;
#endif

%inline %{
  void StartParallel() {
    ParallelDescriptor::StartParallel();
  }
  int rank() {
    return ParallelDescriptor::MyProc();
  }
  int size() {
    return ParallelDescriptor::NProcs();
  }
  Real ReduceRealMax (Real lval) {
    Real rvar = lval;
    ParallelDescriptor::ReduceRealMax(rvar);
    return rvar;
  }
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

%include <numpy.i>

%include <std_string.i>
%rename(__str__) display;


typedef double Real;

template <class T>
class Array
    :
    public std::vector<T>
{
public:
    //
    // Constructs an Array<T> with no elements
    //
    Array ();
    //
    // Constructs an Array<T> of size len with the value of each
    // element defined by the default constructor for T.
    //
    explicit Array (size_t len);
    //
    // Constructs an Array<T> of size len with the value of each
    // elements given by initialvalue.
    //
    Array (size_t   len,
           const T& initialvalue);
    //
    // Constructs an Array<T> of size len in which the K'th
    // value is a copy of vec[K].
    //
    Array (const T* vec,
           size_t   len);

    int size () const { return self->std::vector<T>::size(); }
    //
    // Returns a reference to the K'th element in self Array<T>.
    // The element can be modified through self reference.  The
    // result can be used as an L-value.
    //
    %extend{ T& __getitem__ (size_t i)
    {
        BL_ASSERT(i < size()); return self->std::vector<T>::operator[](i);
    }
    }
    //
    // Same as above, except acts on const Array's.
    //
    %extend{ const T& __getitem__ (size_t i) const
    {
        BL_ASSERT(i < size()); return self->std::vector<T>::operator[](i);
    }
    }
    //
    // Different syntax for operator[].
    //
    T& get (size_t i) { return self->operator[](i); }
    //
    // Different syntax for const operator[].
    //
    const T& get (size_t i) const { return self->operator[](i); }
    //
    // Returns pointer to vector of data.  This function breaks object
    // encapsulation and should only be used for interfacing to
    // Fortran subroutines.
    //
    T* dataPtr () { return &self->operator[](0); }
    //
    // Same as above for constant arrays.
    //
    const T* dataPtr () const { return &self->operator[](0); }
    //
    // Changes the i'th element of self Array<T> to elem.
    //
    void set (size_t i, const T& elem) { get(i) = elem; }
private:
    //
    // This is disallowed.
    //
    Array<T>& operator= (int);
};



class IntVect {
public:

#ifdef DIM1
	IntVect(int i);
#endif
#ifdef DIM2
	IntVect(int i, int j);
#endif
#ifdef DIM3
	IntVect(int i, int j, int k);
#endif
	IntVect (const IntVect& rhs);

	IntVect& shift(int coord, int s);

	%extend {
		void writeOn(std::ofstream *os ){
			*os << *self;
		}
		void read(std::ifstream* ifs){
			*ifs >> *self;
		}
		int __getitem__(int index){
			return (*self)[index];
		}
		int __len__() volatile { return BL_SPACEDIM; }
		void __setitem__(int index,int val){
			(*self).setVal(index,val);
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

	        std::string display() {
	            std::ostringstream str;
                    str << *self;
	            return str.str();
	        }
	}
	// static functions
	static const IntVect& TheZeroVector();
	static const IntVect& TheUnitVector();
	static const IntVect& TheNodeVector();
	static const IntVect& TheCellVector();
};


class Box {
public:
	Box(const IntVect& small, const IntVect& big);

        Box (const IntVect& small,
	     const IntVect& big,
             const IntVect& typ);

	Box ( const Box& b);

        const IntVect smallEnd () const;
        const IntVect bigEnd () const;
        IntVect type () const;
        const IntVect size () const;
        bool contains (const IntVect& p) const;
	bool ok() const;
	bool contains (const Box& b) const;
        bool intersects (const Box& b) const;
        bool sameSize (const Box& b) const;
        bool sameType (const Box &b) const;
        bool cellCentered () const;
        long numPts () const;
        long volume () const;
        long index (const IntVect& v) const;
        Box& setSmall (const IntVect& sm);
        Box& setBig (const IntVect& bg);
        Box& shift (int dir,
                int nzones);
        Box& shiftHalf (int dir,
                    int num_halfs);
        Box& convert (const IntVect& typ);
        Box& surroundingNodes (int dir);
        Box& enclosedCells (int dir);
        Box& minBox (const Box &);
        Box chop (int dir,
                  int chop_pnt);
        Box& grow (int i);
        Box& grow (int idir, int n_cell);
        Box& growLo (int idir,
                     int n_cell=1);
        Box& growHi (int idir,
                     int n_cell=1);
        Box& refine (int refinement_ratio);
        Box& coarsen (int refinement_ratio);
        void next (IntVect &) const;

	%extend{
		void read(std::ifstream* ifs){
			*ifs >> *self;
		}
		void write(std::ofstream* os){
			*os << *self;
		}
		Box __and__(Box *right){
			Box result;
			result = *self & *right;
			return result;
		}
		void writeOn(std::ofstream *os){
			*os << *self;
		}

	        std::string display() {
	            std::ostringstream str;
                    str << *self;
	            return str.str();
	        }

		int __cmp__(const Box* other){
			if(*self == *other) return 0;
			else if( self->smallEnd().lexLT(other->smallEnd())){
				return 1;
			}
			else {
				return -1;
			}
		}
	}
};

class FArrayBox {
public:
	FArrayBox( const Box& b, int ncomp);
        FArrayBox();
	~FArrayBox();

	Box box() const;
	int nComp() const;
	void writeOn(std::ofstream& ofs) const;
	void readFrom(std::istream& ifs);

	double norm(int p, int comp, int numcomp ) const;

	void copy(FArrayBox&src, int srccomp, int destcomp, int numcomp);
	void copy(FArrayBox&src, Box& destbox);
	void copy(FArrayBox&src);

        void readFrom(std::ifstream& is);

	%extend{
	  std::string display() {
	    std::ostringstream str;
            str << *self;
	    return str.str();
	  }
        }

        void setVal (double     x,
                     const Box& bx,
                     int        nstart,
                     int        ncomp);
        void setVal (double     x,
                     const Box& bx,
                     int        N);
        void setVal (double x,
                     int    N);
	void setVal( double x );

        double min(int comp = 0) const;
        double min (const Box& subbox,
                    int        comp = 0) const;
        double max (int comp = 0) const;
        double max (const Box& subbox,
                    int        comp = 0) const;
        IntVect minIndex (int comp = 0) const;
        IntVect minIndex (const Box& subbox,
                          int        comp = 0) const;
        IntVect maxIndex (int comp = 0) const;
        IntVect maxIndex (const Box& subbox,
                          int        comp = 0) const;

	double sum (int comp,
		    int numcomp = 1) const;
	double sum (const Box& subbox,
		    int        comp,
		    int        numcomp = 1) const;
	%extend{
          PyObject *get_array() {
	    double   *ptr = self->dataPtr();
	    PyObject *arr = 0;
	    npy_intp dims[BL_SPACEDIM+1];

	    const IntVect size = self->box().size();
	    for (int i=0; i<BL_SPACEDIM; i++)
	      dims[i] = size[i];
	    dims[BL_SPACEDIM] = self->nComp();
	    arr = PyArray_NewFromDescr(&PyArray_Type,
				       PyArray_DescrFromType(NPY_DOUBLE), BL_SPACEDIM+1, dims, NULL,
				       ptr, NPY_FORTRAN|NPY_WRITEABLE, NULL);

	    Py_INCREF(arr);
	    return arr;
          }
        }

%newobject __add__;
%newobject __sub__;
%newobject __mul__;
%newobject __div__;

	%extend{
		double valIV( const IntVect* p, int n){
			return (*self)(*p,n);
		}
		void setValIV(double x, const IntVect* p, int n=1){
			(*self)(*p,n) = x;
		}
                void floor(double val) {
                  Box bx = self->box();
                  int nc = self->nComp();
                  for (IntVect iv=bx.smallEnd(), End=bx.bigEnd(); iv<=End; bx.next(iv)) {
                    for (int n=0; n<nc; ++n) {
                      (*self)(iv,n) = std::max((*self)(iv,n),val);
                    }
                  }
                }
		FArrayBox* __add__ (const FArrayBox* right ){
			Box bx = self->box();
			bx &= right->box();
			int nvar = self->nComp();
			int rightnvar = right->nComp();
			if (nvar > rightnvar ) nvar = rightnvar;
			FArrayBox* result = new FArrayBox(bx,nvar);
			result->copy(*self);
			(*result) += *right;
			return result;
		}
		FArrayBox* __sub__ (const FArrayBox* right ){
			Box bx = self->box();
			bx &= right->box();
			int nvar = self->nComp();
			int rightnvar = right->nComp();
			if (nvar > rightnvar ) nvar = rightnvar;
			FArrayBox* result = new FArrayBox(bx,nvar);
			result->copy(*self);
			(*result) -= *right;
			return result;
		}
		FArrayBox* __mul__ (const FArrayBox* right ){
			Box bx = self->box();
			bx &= right->box();
			int nvar = self->nComp();
			int rightnvar = right->nComp();
			if (nvar > rightnvar ) nvar = rightnvar;
			FArrayBox* result = new FArrayBox(bx,nvar);
			result->copy(*self);
			(*result) *= *right;
			return result;
		}
		FArrayBox* __mul__ (double right){
			Box bx = self->box();
			int nvar = self->nComp();
			FArrayBox* result = new FArrayBox(bx,nvar);
			result->copy(*self);
			result->mult(right);
			return result;
		}
		FArrayBox* __div__ (const FArrayBox* right ){
			Box bx = self->box();
			bx &= right->box();
			int nvar = self->nComp();
			int rightnvar = right->nComp();
			if (nvar > rightnvar ) nvar = rightnvar;
			FArrayBox* result = new FArrayBox(bx,nvar);
			result->copy(*self);
			(*result) /= *right;
			return result;
		}
        }
};

class BoxArray {
public:
	BoxArray();
	~BoxArray();
	BoxArray(BoxArray&ba);
        void resize (int len);
        void define (const BoxArray& bs);

	void writeOn(std::ostream& os);
        void readFrom (std::istream& is);

        int size () ;
	Box get( int index );
	void set( int i, Box& ibox);

	BoxArray& maxSize(int);
	BoxArray& maxSize(const IntVect& block_size);

        bool ok () const;
        bool isDisjoint () const;
        bool contains (const IntVect& v) const;
        bool contains (const Box& b) const;
        bool contains (const BoxArray& bl) const;
	Box minimalBox() ;
	void refine( int ratio ){
	    self->refine(ratio);
	}
	void coarsen(int ratio){
	    self->coarsen(ratio);
	}

%newobject __and__;
%newobject __or__;
%newobject complementIn;

	%extend{
		Box __getitem__(int index){
			return (*self)[index];
		}

                BoxArray complementIn (const Box& b) {
                    return BoxLib::complementIn(b,*self);
                }
		BoxArray __and__(const BoxArray& right){
                    BoxArray res = BoxLib::intersect(*self,right);
		    return res;
		}
		BoxArray __or__(const BoxArray& right){
                    BoxList bl(*self);
                    bl.join(BoxList(right));
		    return BoxArray(bl);
		}
		int __len__() { return (*self).size(); }
	        std::string display() {
	            std::ostringstream str;
                    str << *self;
	            return str.str();
	        }
		int __cmp__(const BoxArray* other){
		    if( (*self)==(*other)) return 0;
		    else return 1;
		}
	}
};

%include <std_list.i>
%include <std_pair.i>
%include <std_vector.i>

/* // From support.H */
/* // */
/* void SAXPYjpdf(FArrayBox& dest,double scale,const FArrayBox& src); */
/* std::vector<double> condMean(FArrayBox& src,bool cdir); */

/* // transfers between coarse and fine grids */
/* void tocoarse(int ratio,FArrayBox&fine,FArrayBox&crse,IntVect&iv); */
/* void tofine(int ratio,FArrayBox&crse,FArrayBox&fine,IntVect&iv,Real defval); */
/* void injectCoarse(FArrayBox&fine,int ratio, const IntVect&iv, */
/*                   FArrayBox&crse, const Box&tbox); */
/* // conditional vector merges */
/* void cvmgnfab(FArrayBox&res,FArrayBox&n,FArrayBox&p,FArrayBox&trg); */
/* void cvmgpfab(FArrayBox&res,FArrayBox&p,FArrayBox&n,FArrayBox&trg); */
/* // misc */
/* //void gradfillfab(FArrayBox&fab, const IntVect&iv); */

class MultiFab {
public:
        MultiFab (const BoxArray& bs,
                  int             ncomp,
                  int             ngrow);
	MultiFab();
	~MultiFab();
	%extend{
            void define (const BoxArray& bxs,
                         int             nvar,
                         int             ngrow,
                         FabAlloc        mem_mode = Fab_allocate) {
                self->FabArray<FArrayBox>::define(bxs,nvar,ngrow,mem_mode);
            }
        }

        bool ok () const;
        int nGrow () const;
        const BoxArray& boxArray () const;
        int size () const;
        int nComp () const;

	double min(int comp, int nghost=0) const;
	double max(int comp, int nghost=0) const;

	void copy(MultiFab &src){
	    self->copy(src);
	}

        void setVal (double val);
        void setVal (double val,
                     int        comp,
                     int        num_comp,
                     int        nghost = 0);
        void setVal (double val,
                     const Box& region,
                     int        comp,
                     int        num_comp,
                     int        nghost = 0);
        void setVal (double val,
                     int        nghost);
        void setVal (double val,
                     const Box& region,
                     int        nghost);
        void setBndry (double val);
        void setBndry (double val,
                       int        strt_comp,
                       int        ncomp);

	void FillBoundary (int scomp, int ncomp);
	void FillBoundary (bool local = false, bool cross = false);
	//void FillBoundary (int scomp, int ncomp, bool cross);

	%extend{
	    double sum(int comp = 0 ) const {
	        const BoxArray &ba = self->boxArray();
	        double retval = 0;
	        for(MFIter mfi(*self);mfi.isValid(); ++mfi){
		    int idx = mfi.index();
		    const Box &bx = ba[idx];
		    const FArrayBox &fab = (*self)[mfi];
		    retval += fab.sum(bx,comp,1);
	        }
		ParallelDescriptor::ReduceRealSum(retval);
	        return retval;
	    }
		void writeOut(char *name){
		    std::string str(name);
		    VisMF::Write(*self,str,VisMF::OneFilePerCPU);
		}
		//void copyToFab(FArrayBox &dest){
		//    self->copy(dest);
		//}
		void copyComp(MultiFab*src,int src_comp,int dest_comp,
				int num_comp){
		    self->copy(*src,src_comp,dest_comp,num_comp);
		}
		FArrayBox* __getitem__(int index) {
		    if ((*self).DistributionMap()[index]==ParallelDescriptor::MyProc() ) {
                       return &((*self)[index]);
                    }
                    return 0;
		}
	}
};

class RealBox {
public:
  RealBox (Real IN_ARRAY1[BL_SPACEDIM] /* lo */,
	   Real IN_ARRAY1[BL_SPACEDIM]) /* hi */;
};

class Geometry {
public:
    Geometry (const Box&     dom,
              const RealBox* rb     ,
              int            coord  ,
              int IN_ARRAY1[BL_SPACEDIM] );

    void FillPeriodicBoundary (MultiFab& mf,
                               bool      do_corners = false,
                               bool      local      = false);
};



