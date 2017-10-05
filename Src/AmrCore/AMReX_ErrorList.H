#ifndef _ErrorList_H_
#define _ErrorList_H_

#include <string>
#include <memory>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

extern "C"
{
    //
    // Type of extern "C" function called by ErrorRec to do tagging of cells for refinement.
    //
    typedef void (*ErrorFuncDefault)(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                                     const int* tagval, const int* clearval,
                                     amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                                     const int* lo, const int * hi, const int* nvar,
                                     const int* domain_lo, const int* domain_hi,
                                     const amrex::Real* dx, const amrex::Real* xlo,
                                     const amrex::Real* prob_lo, const amrex::Real* time,
                                     const int* level);

    typedef void (*ErrorFunc2Default)(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                                      const int* tagval, const int* clearval,
                                      amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                                      const int* lo, const int * hi, const int* nvar,
                                      const int* domain_lo, const int* domain_hi,
                                      const amrex::Real* dx, const int* level, const amrex::Real* avg);

    // Dimension agnostic version that always has three elements.
    // Note that this is only implemented for the ErrorFunc class, not ErrorFunc2.

    typedef void (*ErrorFunc3DDefault)(int* tag, const int* tlo, const int* thi,
                                       const int* tagval, const int* clearval,
         			       amrex::Real* data, const int* data_lo, const int* data_hi, 
                                       const int* lo, const int * hi, const int* nvar,
                                       const int* domain_lo, const int* domain_hi,
                                       const amrex::Real* dx, const amrex::Real* xlo,
                                       const amrex::Real* prob_lo, const amrex::Real* time,
                                       const int* level);  
  
}

namespace amrex {

//
// Error Record.
//
//
// ErrorRec is designed to tag cells for regridding based on the state
// data contained in AmrLevel and its derivatives.  It is conceivable that
// actual error tagging will be through derivation, so provision is made
// for this as well.
//

class ErrorRec
{
public:
    //
    // Error types.
    //
    enum ErrorType { Special=0, Standard, UseAverage };
    //
    // Class wrapper around ErrorFuncDefault.
    //
    class ErrorFunc
    {
    public:
        //
        // Bogus constructor.
        //
        ErrorFunc ();
        //
        // A Constructor.
        //
        ErrorFunc (ErrorFuncDefault inFunc);
        //
        // A Constructor.
        //
        ErrorFunc (ErrorFunc3DDefault inFunc);      
        //
        // Return a ptr to a clone of this object.
        // It is the responsibility of the caller to delete the result.
        //
        virtual ErrorFunc* clone () const;
        //
        // Destructor.
        //
        virtual ~ErrorFunc ();
        //
        // Tag cells using "regular" function.
        //
        virtual void operator () (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                                  const int* tagval, const int* clearval,
                                  Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const Real* xlo,
                                  const Real* prob_lo, const Real* time,
                                  const int* level) const;
        //
        // Tag cells using dimension-agnostic "regular" function.
        //
        virtual void operator () (int* tag, const int* tlo, const int* thi,
                                  const int* tagval, const int* clearval,
                                  Real* data, const int* data_lo, const int* data_hi,
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const Real* xlo,
                                  const Real* prob_lo, const Real* time,
                                  const int* level) const;      
    protected:

        ErrorFuncDefault   m_func;
        ErrorFunc3DDefault m_func3D;      
    };
    class ErrorFunc2
    {
    public:
        //
        // Bogus constructor.
        //
        ErrorFunc2 ();
        //
        // A Constructor.
        //
        ErrorFunc2 (ErrorFunc2Default inFunc);
        //
        // Return a ptr to a clone of this object.
        // It is the responsibility of the caller to delete the result.
        //
        virtual ErrorFunc2* clone () const;
        //
        // Destructor.
        //
        virtual ~ErrorFunc2 ();
        //
        // Tag cells cells using "v2" interface
        //
        virtual void operator () (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                                  const int* tagval, const int* clearval,
                                  Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const int* level, const Real* avg) const;
    protected:

        ErrorFunc2Default  m_func;
    };
    //
    // The constructors.
    //
    ErrorRec (const std::string&          nm,
              int                         ng,
              ErrorType                   etyp,
              const ErrorRec::ErrorFunc2& f2);
    ErrorRec (const std::string&         nm,
              int                        ng,
              ErrorType                  etyp,
              const ErrorRec::ErrorFunc& f);
    virtual ~ErrorRec();
    //
    // The name of the quantity to derive.
    //
    const std::string& name () const;
    //
    // The number of extra zones needed for derivation.
    //
    int nGrow () const;
    //
    // The type of the error tagging.
    //
    ErrorType errType () const;
    //
    // The extern "C" functions to do the error tagging.
    //
    virtual const ErrorRec::ErrorFunc&  errFunc () const;
    virtual const ErrorRec::ErrorFunc2& errFunc2() const;

private:
    //
    // Name of quantity to derive.
    //
    std::string derive_name;
    //
    // Number of extra zones.
    //
    int ngrow;
    //
    // The type of Error.
    //
    ErrorType err_type; 
    //
    // Functions to do error estimation.
    //
    ErrorFunc*  err_func;
    ErrorFunc2* err_func2;
};

//
// A List of ErrorRecs.
//
// Container class for ErrorRecs.
//

class ErrorList
{
public:
    ErrorList() {}
    //
    // The number of ErrorRecs in the list.
    //
    int size () const;
    //
    // Append a new ErrorRec to the list.
    //
    void add (const std::string&  name,
              int                 nextra, 
              ErrorRec::ErrorType typ,
              const ErrorRec::ErrorFunc&    func);
    void add (const std::string&  name,
              int                 nextra, 
              ErrorRec::ErrorType typ,
              const ErrorRec::ErrorFunc2&   func);
    //
    // The kth ErrorRec.
    //
    const ErrorRec& operator[] (int k) const;

    void clear (bool rs0 = false) { vec.clear(); if(rs0) { vec.resize(0); } }

private:

    Vector<std::unique_ptr<ErrorRec> > vec;
};

std::ostream& operator << (std::ostream& os, const ErrorList& elst);

}

#endif

