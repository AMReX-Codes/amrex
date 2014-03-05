%module boxlib
%{
#  include <iostream>
#  include <ostream>
#  include <fstream>
#  include <sstream>
#  include <stdio.h>
#  include <IntVect.H>
#  include <Box.H>
#  include <FArrayBox.H>
#  include <BoxArray.H>
#  include <MultiFab.H>
#  include <VisMF.H>
#  include <Array.H>
#  include <ParallelDescriptor.H>
#  include <support.H>

#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
    import_array();
%}

#define DIM2
#ifdef DIM2
%constant int BL_DIM = 2;
#else
%constant int BL_DIM = 3;
#endif

%inline %{  
  void StartParallel() {
    ParallelDescriptor::StartParallel();
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
#ifdef DIM2
	IntVect( int i, int j);
#else
	IntVect( int i, int j, int k);
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
		int __len__() { return BL_SPACEDIM; }
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
          PyObject * PyArr() {
             int n1, n2, n3, n4;
             double* ptr = self->dataPtr();
             PyObject *arr = 0;
             int ndim = 4;
             const IntVect size = self->box().size();
             n1 = size[0];
             n2 = (BL_SPACEDIM > 1 ? size[1] : 1);
             n3 = (BL_SPACEDIM > 2 ? size[2] : 1);
             n4 = self->nComp();
             npy_intp dims[4] = {n1, n2, n3, n4};
             arr = PyArray_NewFromDescr(&PyArray_Type,
                                        PyArray_DescrFromType(NPY_DOUBLE), ndim, dims, NULL,
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

// From support.H
//
void SAXPYjpdf(FArrayBox& dest,double scale,const FArrayBox& src);
std::vector<double> condMean(FArrayBox& src,bool cdir);

// transfers between coarse and fine grids
void tocoarse(int ratio,FArrayBox&fine,FArrayBox&crse,IntVect&iv);
void tofine(int ratio,FArrayBox&crse,FArrayBox&fine,IntVect&iv,Real defval);
void injectCoarse(FArrayBox&fine,int ratio, const IntVect&iv, 
                  FArrayBox&crse, const Box&tbox);
// conditional vector merges
void cvmgnfab(FArrayBox&res,FArrayBox&n,FArrayBox&p,FArrayBox&trg);
void cvmgpfab(FArrayBox&res,FArrayBox&p,FArrayBox&n,FArrayBox&trg);
// misc
void gradfillfab(FArrayBox&fab, const IntVect&iv);

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

	%extend{
	    double sum(int comp = 0 ) const {
	        const BoxArray &ba = self->boxArray();
	        double retval = 0;
	        for(MFIter mfi(*self);mfi.isValid(); ++mfi){
		    int idx = mfi.index();
		    const Box &bx = ba[idx];
		    const FArrayBox &fab = (*self)[idx];
		    retval += fab.sum(bx,comp,1);
	        }
		ParallelDescriptor::ReduceRealSum(retval);
	        return retval;
	    }
		void writeOut(char *name){
		    std::string str(name);
		    VisMF::Write(*self,str,VisMF::OneFilePerCPU);
		}
		void copyToFab(FArrayBox &dest){
		    self->copy(dest);
		}
		void copyComp(MultiFab*src,int src_comp,int dest_comp,
				int num_comp){
		    self->copy(*src,src_comp,dest_comp,num_comp);
		}
		// FArrayBox* __getitem__(int index) {
		//     if ((*self).DistributionMap()[index]==ParallelDescriptor::MyProc() ) {
                //        return : &((*self)[index]);
                //     }
                //     return 0;
		// }
	}
};

// Ignore definition of nested classes in ChemDriver header below
%nestedworkaround ChemDriver::Edge;
%nestedworkaround ChemDriver::Group;

%{
#include <ChemDriver.H>
#include <chemSupport.H>
%}

// We've fooled SWIG into thinking that Edge is a global class, so now we need
// to trick the C++ compiler into understanding this apparent global type.
%{
  typedef ChemDriver::Edge Edge;
%}

namespace std {
  %template(EdgeList) std::list<Edge>;
  %template(IntDoublePair) std::pair<int,double>;
  %template(IntDoubleVec) std::vector<std::pair<int,double> >;
  %template(DoubleVec) std::vector<double>;
};

//
// Data structure for generating chemical path diagrams,
//   an Edge represents the transfer of an atom from one chemical
//   species, sp1, to another, sp2, which occurs through a list
//   of reactions.  For each reaction, n of these atoms are
//   transferred (RWL is a list of pairs of rxn ID to atom transfer cnt)
//
class Edge
{
 public:
  //friend std::ostream& operator<< (std::ostream& os, const Edge& e);
  
  %extend{
    std::string display() {
      std::ostringstream str;
      str << *self;
      return str.str();
    }
  }

  Edge ();

  Edge (const std::string& n1, const std::string& n2, const Array<std::pair<int,Real> > rwl);
    
  Edge (const std::string& n1, const std::string& n2, int reac, Real weight );
    
  int equivSign (const Edge& rhs) const;

  void combine (const Edge& rhs, int sgn);

  bool touchesSp(const std::string& rhs) const;
    
  void reverse();

  %extend{
     std::vector<std::pair<int,double> > rwl () const
     { return self->rwl(); }
  }
  const std::string left() const;

  const std::string right() const;

  bool operator< (const Edge& rhs) const;
};



class ChemDriver
{
public:

  enum TRANSPORT
  {
    CD_EG, 
    CD_TRANLIB
  };

    ChemDriver ();	  

    ~ChemDriver ();

    static void SetTransport(const ChemDriver::TRANSPORT& tran_in);
    static ChemDriver::TRANSPORT Transport ();

    void solveTransient (FArrayBox&        Ynew,
                         FArrayBox&        Tnew,
                         const FArrayBox&  Yold,
                         const FArrayBox&  Told,
                         FArrayBox&        FuncCount,
                         const Box&        box,
                         int               sCompY,
                         int               sCompT,
                         Real              dt,
                         Real              Patm,
                         FArrayBox*        chemDiag=0,
                         bool              use_stiff_solver = true) const;

#ifdef LMC_SDC
    void solveTransient_sdc(FArrayBox&        rhoYnew,
                            FArrayBox&        rhoHnew,
                            FArrayBox&        Tnew,
                            const FArrayBox&  rhoYold,
                            const FArrayBox&  rhoHold,
                            const FArrayBox&  Told,
                            const FArrayBox&  const_src,
                            FArrayBox&        FuncCount,
                            const Box&        box,
                            int               sComprhoY,
                            int               sComprhoH,
                            int               sCompT,
                            Real              dt,
                            Real              Patm,
                            FArrayBox*        chemDiag,
                            bool              use_stiff_solver = true) const;
#endif

    void set_verbose_vode ();
    void set_max_vode_subcycles (int max_cyc);
    void set_species_Yscales (const std::string& scalesFile);
    //
    // Species info.
    //
    int numSpecies () const;
    int numElements () const;
    int numReactions () const;
    const Array<std::string>& speciesNames () const;
    const Array<std::string>& elementNames () const;
    Array<Real> speciesMolecWt () const;
    Array<Real> elementAtomicWt () const;
    int index (const std::string& specName) const;
    int indexElt (const std::string& eltName) const;
    Array<int> reactionsWithXonL (const std::string& lSpecName) const;
    Array<int> reactionsWithXonR (const std::string& rSpecName) const;
    Array<std::pair<std::string,int> > specCoeffsInReactions (int reacIdx) const;
    std::string reactionString (int reacIdx) const;
    const Array<int>& reactionMap() const;

    int numberOfElementXinSpeciesY (const std::string& eltX,
                                    const std::string& spcY) const;
    //
    // Thermo info.
    //
    void getRhoGivenPTY (FArrayBox&       Rho,
                         Real             Patm,
                         const FArrayBox& T,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompR) const;

    void getRhoGivenPvTY (FArrayBox&       Rho,
                          const FArrayBox& P,
                          const FArrayBox& T,
                          const FArrayBox& Y,
                          const Box&       box,
                          int              sCompP,
                          int              sCompT,
                          int              sCompY,
                          int              sCompR) const;

    void getPGivenRTY (FArrayBox&       p,
                       const FArrayBox& Rho,
                       const FArrayBox& T,
                       const FArrayBox& Y,
                       const Box&       box,
                       int              sCompR,
                       int              sCompT,
                       int              sCompY,
                       int              sCompP) const;

    void getTGivenPRY (FArrayBox&       T,
                       Real             p,
                       const FArrayBox& Rho,
                       const FArrayBox& Y,
                       const Box&       box,
                       int              sCompR,
                       int              sCompY,
                       int              sCompT) const;

    void getCpmixGivenTY (FArrayBox&       cpmix,
                          const FArrayBox& T,
                          const FArrayBox& Y,
                          const Box&       box,
                          int              sCompT,
                          int              sCompY,
                          int              sCompCp) const;

    void getCvmixGivenTY (FArrayBox&       cvmix,
                          const FArrayBox& T,
                          const FArrayBox& Y,
                          const Box&       box,
                          int              sCompT,
                          int              sCompY,
                          int              sCompCv) const;

    void getHmixGivenTY (FArrayBox&       hmix,
                         const FArrayBox& T,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompH) const;

    void getMwmixGivenY (FArrayBox&       mwmix,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompY,
                         int              sCompMw) const;

    void getCpGivenT (FArrayBox&       cp,
                      const FArrayBox& T,
                      const Box&       box,
                      int              sCompT,
                      int              sCompCp) const;

    void getHGivenT (FArrayBox&       h,
                     const FArrayBox& T,
                     const Box&       box,
                     int              sCompT,
                     int              sCompH) const;

    Real getRuniversal () const;
    Real getP1atm_MKS () const;
    //
    // Compute T that satisfies hmix=sum(h(T)), returns max Newton iterations
    // on any point over grid, return -1 if T jumped out of bounds during solve,
    // and -2 if solve failed anywhere.  Note that a temporary is not used, and
    // the solver kicks out when/if it fails, so it may exit after converting
    // only part of the temperature array.  Save a copy of T and check the return
    // code if you want to be extra careful...
    //
    int getTGivenHY (FArrayBox&       T,
                     const FArrayBox& H,
                     const FArrayBox& Y,
                     const Box&       box,
                     int              sCompH,
                     int              sCompY,
                     int              sCompT,
                     const Real&      errMAX = -1) const;

    Array<Real> massFracToMoleFrac (const Array<Real>& Y) const;
    Array<Real> moleFracToMassFrac (const Array<Real>& X) const;

    void molarProduction (FArrayBox&       Q,
                          const std::string&   specName,
                          const FArrayBox& C,
                          const FArrayBox& T,
                          const Box&       box,
                          int              sCompC,
                          int              sCompT,
                          int              sCompQ) const;
    //
    // Compute heat release (J/(m^3.s)) based on the temp, press and mass fractions
    //
    void heatRelease (FArrayBox&       Q,
                      const FArrayBox& Y,
                      const FArrayBox& T,
                      Real             Patm,
                      const Box&       box,
                      int              sCompY,
                      int              sCompT,
                      int              sCompQ) const;
    //
    // Compute dY/dt based on the input temp, press and mass fractions.
    //
    void reactionRateY (FArrayBox&       Ydot,
                        const FArrayBox& Y,
                        const FArrayBox& T,
                        Real             Patm,
                        const Box&       box,
                        int              sCompY,
                        int              sCompT,
                        int              sCompYdot) const;

#ifdef LMC_SDC
    //
    // Compute dRhoY/dt based on the input temp, press and mass densities.
    //
    void reactionRateRhoY(FArrayBox&       RhoYdot,
                          const FArrayBox& RhoY,
                          const FArrayBox& RhoH,
                          const FArrayBox& T,
                          Real             Patm,
                          const Box&       box,
                          int              sCompRhoY,
                          int              sCompRhoH,
                          int              sCompT,
                          int              sCompRhoYdot) const;
#endif

    void fwdRevReacRatesGivenXTP (FArrayBox&        FwdK,
                                  FArrayBox&        RevK,
                                  const Array<int>& rxnIDs,
                                  const FArrayBox&  X,
                                  const FArrayBox&  T,
                                  Real              Patm,
                                  const Box&        box,
                                  int               sCompX,
                                  int               sCompT,
                                  int               sCompFwdK,
                                  int               sCompRevK) const;
    
    void massFracToMoleFrac (FArrayBox&       X,
                             const FArrayBox& Y,
                             const Box&       box,
                             int              sCompY,
                             int              sCompX) const;

    void moleFracToMassFrac (FArrayBox&       Y,
                             const FArrayBox& X,
                             const Box&       box,
                             int              sCompX,
                             int              sCompY) const;

    void massFracToMolarConc (FArrayBox&       C,
                              const FArrayBox& Y,
                              const FArrayBox& T,
                              Real             Patm,
                              const Box&       box,
                              int              sCompY,
                              int              sCompT,
                              int              sCompC) const;

    void massFracToMolarConc (FArrayBox&       C,
                              const FArrayBox& Y,
                              const FArrayBox& T,
                              const FArrayBox& Rho,
                              const Box&       box,
                              int              sCompR,
                              int              sCompY,
                              int              sCompT,
                              int              sCompC) const;

    void molarConcToMoleFrac (FArrayBox&       X,
                              const FArrayBox& C,
                              const Box&       box,
                              int              sCompC,
                              int              sCompX) const;
    //
    // Normalize mass fractions to prevent negatives.
    //
    void normalizeMassFrac (FArrayBox&       Ynorm,
                            const FArrayBox& Y,
                            const std::string&   excessSpecies,
                            const Box&       box,
                            int              sCompY,
                            int              sCompYnorm) const;
    //
    // Chemical Diffusivities.
    //
    void getMixAveragedRhoDiff (FArrayBox&       rhoD,
                                const FArrayBox& Y,
                                const FArrayBox& T,
                                Real             Patm,
                                const Box&       box,
                                int              sCompY,
                                int              sCompT,
                                int              sCompRD) const;
    //
    // Viscosity.
    //
    void getMixShearVisc (FArrayBox&       eta,
                          const FArrayBox& Y,
                          const FArrayBox& T,
                          const Box&       box,
                          int              sCompY,
                          int              sCompT,
                          int              sCompE) const;
    
    void getElementMoles (FArrayBox&       C_elt,
                          const std::string&   name,
                          const FArrayBox& C,
                          const Box&       box,
                          int              sCompC,
                          int              sCompC_elt) const;
    //
    // Optically thin radiation model.
    //
    void getOTradLoss_TDF (FArrayBox&       Qloss,
                           const FArrayBox& T,
                           const FArrayBox& X,
                           const Real       Patm,
                           const Real       T_bg,
                           const Box&       box,
                           int              sCompX,
                           int              sCompT,
                           int              sCompQ) const;
    //
    // H - to - T solve parameter access.
    //
    Real getHtoTerrMAX () const;
    void setHtoTerrMAX (Real err);
    int getHtoTiterMAX () const;
    void setHtoTiterMAX (int err);    
    //
    // Handy functions.
    //
    static Array<int> encodeStringForFortran(const std::string& astr);
    static std::string decodeStringFromFortran(const int* coded, int length);

    //
    // Data structure for generating chemical path diagrams,
    //   an Edge represents the transfer of an atom from one chemical
    //   species, sp1, to another, sp2, which occurs through a list
    //   of reactions.  For each reaction, n of these atoms are
    //   transferred (RWL is a list of pairs of rxn ID to atom transfer cnt)
    //
    class Edge
    {
    public:
        friend std::ostream& operator<< (std::ostream& os, const Edge& e);
        
        Edge () {}

        Edge (const std::string& n1, const std::string& n2, const Array<std::pair<int,Real> > rwl);
        
        Edge (const std::string& n1, const std::string& n2, int reac, Real weight );
        
        int equivSign (const Edge& rhs) const;

        void combine (const Edge& rhs, int sgn);

        bool touchesSp(const std::string& rhs) const;
        
        void reverse();

        %extend{
           const IntDoubleArray& rwl () const;
        }
        const std::string left() const;

        const std::string right() const;

        bool operator< (const Edge& rhs) const;

        const Array<std::pair<int,Real> > RateWeightList() const {return RWL;}

    private:
        std::string sp1, sp2;
        Array<std::pair<int,Real> > RWL; // RateWeightList, each pair is (reac,coeff)
    };

    // 
    // Helper class for building chem edges.  A group is a list of constituent 
    // atoms, and this class allows a coupleof useful math operations on these
    // groups.  
    //
    class Group
    {
    public:
        friend std::ostream& operator<< (std::ostream& os, const Group& g);
        
        Group () {}

        Group (const std::map<std::string,int>& eltCnts);

        Group (const Group& rhs);

        Group operator- (const Group& rhs) const;

        Group operator* (int rhs) const;

        int operator[](const std::string& id) const;
            
        bool sameSign() const;

        bool contains(const std::string& id) const;
    
        Real awt(); // non-const because embedded lazy eval

        int size() const;
    
    private:
        void FillAtomicWeights ();
        std::map<std::string,int> mEltCnts;
        static std::map<std::string,Real> AtomicWeight;
    };

    // 
    // Compute edges from chem mechanism
    //
    std::list<Edge> getEdges (const std::string& trElt,
                              int PrintVerbose=0,
                              int HackSplitting=1) const;
    
protected:

    void getSpeciesNames ();
    void getElementNames ();
    void getStoichCoeffs ();

private:

    void initOnce ();

    Array<std::string> mSpeciesNames;
    Array<std::string> mElementNames;
    Real               mHtoTerrMAX;
    int                mHtoTiterMAX;
    Array<Real>        mTmpData;
    int mMaxreac, mMaxspec, mMaxelts, mMaxord, mMaxthrdb, mMaxtp, mMaxsp, mMaxspnml;
    Array<int>         mNu;
    Array<int> reaction_map;
};


%include <chemSupport.H>
