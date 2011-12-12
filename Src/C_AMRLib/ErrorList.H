#ifndef _ErrorList_H_
#define _ErrorList_H_

#include <string>
#include <vector>

#include <REAL.H>
#include <ArrayLim.H>

extern "C"
{
    //
    // Type of extern "C" function called by ErrorRec to do tagging of cells for refinement.
    //
    typedef void (*ErrorFunc)(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                              const int* tagval, const int* clearval,
                              Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                              const int* lo, const int * hi, const int* nvar,
                              const int* domain_lo, const int* domain_hi,
                              const Real* dx, const Real* xlo,
                              const Real* prob_lo, const Real* time,
                              const int* level);

    typedef void (*ErrorFunc2)(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                               const int* tagval, const int* clearval,
                               Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                               const int* lo, const int * hi, const int* nvar,
                               const int* domain_lo, const int* domain_hi,
                               const Real* dx, const int* level, const Real* avg);
}

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
    // The constructors.
    //
    ErrorRec (const std::string& nm,
              int                ng,
              ErrorType          etyp,
              ErrorFunc2         f2);
    ErrorRec (const std::string& nm,
              int                ng,
              ErrorType          etyp,
              ErrorFunc          f);
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
    ErrorFunc  errFunc () const;
    ErrorFunc2 errFunc2() const;

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
    ErrorFunc  err_func;
    ErrorFunc2 err_func2;
};

//
// A List of ErrorRecs.
//
// Container class for ErrorRecs.
//

class ErrorList
{
public:
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
              ErrorFunc           func);
    void add (const std::string&  name,
              int                 nextra, 
              ErrorRec::ErrorType typ,
              ErrorFunc2          func);
    //
    // The kth ErrorRec.
    //
    const ErrorRec& operator[] (int k) const;

    void clear () { vec.clear(); }

private:

    std::vector<ErrorRec> vec;
};

std::ostream& operator << (std::ostream& os, const ErrorList& elst);

#endif

