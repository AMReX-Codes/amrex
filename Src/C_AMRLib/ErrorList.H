//BL_COPYRIGHT_NOTICE

#ifndef _ErrorList_H_
#define _ErrorList_H_

//
// $Id: ErrorList.H,v 1.5 1998-03-23 21:02:56 lijewski Exp $
//

#include <REAL.H>
#include <List.H>
#include <aString.H>
#include <ArrayLim.H>


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
    enum ErrorType { Richardson=0, Special };
    //
    // Function called by ErrorRec to do tagging of cells for refinement.
    //
    typedef void (*ErrorFunc)(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                              const int* tagval, const int* clearval,
                              Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                              const int* lo, const int * hi, const int* nvar,
                              const int* domain_lo, const int* domain_hi,
                              const Real* dx, const Real* xlo,
                              const Real* prob_lo, const Real* time,
                              const int* level);

    ErrorRec (const aString& nm,
              int            ng,
              ErrorType      etyp,
              ErrorFunc      f)
        :
        derive_name(nm),
        ngrow(ng),
        err_func(f),
        err_type(etyp)
    {}

    ~ErrorRec () {}

    const aString& name () const { return derive_name; }

    int nGrow () const { return ngrow; }

    ErrorType errType () const { return err_type; }

    ErrorFunc errFunc () const { return err_func; }

private:
    //
    // Name of quantity to derive.
    //
    aString derive_name;
    //
    // Number of extra zones.
    //
    int ngrow;
    //
    // Function to do error estimation.
    //
    ErrorFunc err_func;
    //
    // The type of Error.
    //
    ErrorType err_type; 
};

class ErrorList
{
    friend ostream& operator<< (ostream& os, const ErrorList& el);

public:

    ErrorList () {}

    ~ErrorList () {}

    int length () const { return lst.length(); }

    void add (const aString&      name,
              int                 nextra, 
              ErrorRec::ErrorType typ,
              ErrorRec::ErrorFunc func)
    {
        //
        // Keep list in order of definition, append().
        //
        lst.append(ErrorRec(name, nextra, typ, func));
    }

    const ErrorRec* operator[] (int k) const;

protected:

    List<ErrorRec> lst;
};

#endif

