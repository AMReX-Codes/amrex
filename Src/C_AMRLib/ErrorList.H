//BL_COPYRIGHT_NOTICE

#ifndef _ErrorList_H_
#define _ErrorList_H_

//
// $Id: ErrorList.H,v 1.2 1997-11-24 18:52:29 lijewski Exp $
//

#include <REAL.H>
#include <List.H>
#include <aString.H>
#include <ArrayLim.H>

/// Function called by ErrorRec to do actual tagging of cells for refinement.

typedef void (*ErrorFunc)(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
			  const int* tagval, const int* clearval,
			  REAL* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const int* lo, const int * hi, const int* nvar,
			  const int* domain_lo, const int* domain_hi,
			  const REAL* dx, const REAL* xlo,
                          const REAL* prob_lo, const REAL* time,
			  const int* level);

enum ErrorType { Richardson=0, Special };

/// Contains necessary info and functions for tagging cells for regridding
class ErrorRec
{

  /** ErrorRec is designed to tag cells for regridding based on the state
    data contained in AmrLevel and its derivatives.  It is conceivable that
    actual error tagging will be through derivation, so provision is made for
    this as well.
    */

    friend class ErrorList;
    friend ostream& operator << (ostream&, const ErrorList& );

private:
 // name of quantity to derive
    aString derive_name; 
//  number of extra zones
    int	   ngrow;	
//  function to do error estimation
    ErrorFunc    err_func;
// 
    ErrorType    err_type;      
///
public:
//  constructor
    ErrorRec(const aString &nm, int ng, ErrorType et, ErrorFunc f);
//  destructor
    ~ErrorRec();      
///
    const aString& name() const { return derive_name; }
///
    int nGrow() const { return ngrow; }
///
    ErrorType errType() const { return err_type; }
///
    ErrorFunc errFunc() const { return err_func; }
};

/// Container class for ErrorRec
class ErrorList
{
///
    friend ostream& operator << (ostream&, const ErrorList& );

public:
///
    ErrorList();
///
    ~ErrorList();
///
    void clear();
///
    int length() const {return lst.length(); }
///
    void add(const aString &name, int nextra, 
	     ErrorType typ, ErrorFunc func);
///
    const ErrorRec* operator[] (int k) const;

protected:
///
    List<ErrorRec*> lst;
};

#endif

