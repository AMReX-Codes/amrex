
#ifndef _ErrorList_H_
#define _ErrorList_H_

//
// $Id: ErrorList.H,v 1.9 2000-10-02 20:48:42 lijewski Exp $
//

#include <REAL.H>
#include <List.H>
#include <aString.H>
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
}

//
//@Man:
//@Memo: Error Record.
/*@Doc:

  ErrorRec is designed to tag cells for regridding based on the state
  data contained in AmrLevel and its derivatives.  It is conceivable that
  actual error tagging will be through derivation, so provision is made
  for this as well.
*/

class ErrorRec
{
public:
    //
    //@ManDoc: Error types.
    //
    enum ErrorType { Richardson=0, Special };
    //
    //@ManDoc: The constructor.
    //
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
    //
    //@ManDoc: The destructor.
    //
    ~ErrorRec () {}
    //
    //@ManDoc: The name of the quantity to derive.
    //
    const aString& name () const { return derive_name; }
    //
    //@ManDoc: The number of extra zones needed for derivation.
    //
    int nGrow () const { return ngrow; }
    //
    //@ManDoc: The type of the error tagging.
    //
    ErrorType errType () const { return err_type; }
    //
    //@ManDoc: The extern "C" function to do the error tagging.
    //
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

//
//@Man:
//@Memo: A List of ErrorRecs.
/*@Doc:

  Container class for ErrorRecs.
*/

class ErrorList
{
    friend ostream& operator<< (ostream& os, const ErrorList& el);

public:
    //
    //@ManDoc: The constructor.
    //
    ErrorList () {}
    //
    //@ManDoc: The destructor.
    //
    ~ErrorList () {}
    //
    //@ManDoc: The number of ErrorRecs in the list.
    //
    int length () const { return lst.length(); }
    //
    //@ManDoc: Append a new ErrorRec to the list.
    //
    void add (const aString&      name,
              int                 nextra, 
              ErrorRec::ErrorType typ,
              ErrorFunc           func)
    {
        //
        // Keep list in order of definition, append().
        //
        lst.append(ErrorRec(name, nextra, typ, func));
    }
    //
    //@ManDoc: The kth ErrorRec.
    //
    const ErrorRec* operator[] (int k) const;

protected:

    List<ErrorRec> lst;
};

#endif

