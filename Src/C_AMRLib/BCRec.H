#ifndef _BCREC_H_
#define _BCREC_H_ 

//
// $Id: BCRec.H,v 1.1 1997-11-18 19:30:21 lijewski Exp $
//

#include <Misc.H>
#include <Box.H>
#include <BC_TYPES.H>

/// Necessary information and functions for computing boundary conditions
class BCRec
{

/* The BCRec class computes and stores  boundary condition
data for grids in the AMR hierarchy.  The taxonomy is defined
in Blob.C; boundary data can (usually) be proscribed in the "inputs" file.
*/

protected:
//  array of integer values describing boundary conditions
    int bc[2*BL_SPACEDIM];

public:
#ifdef __GNUG__
     bool operator== (const BCRec& rhs) const { return true; }
#endif 

//  sets to invalid bndry types
    BCRec(){};

//  constructor
    BCRec(D_DECL(int loX, int loY, int loZ),
	  D_DECL(int hiX, int hiY, int hiZ));

//  constructor
    BCRec(const int* lo, const int* hi);

/*  inherits bndry types from bc_domain when bx lies on edge of
    domain otherwise gets interior Dirchlet
*/
    BCRec(const BOX& bx, const BOX& domain, const BCRec &bc_domain);
    
//  Explicitly set bndry value
    void setLo(int dir, int bc_val) { bc[dir] = bc_val; }

//  Explicitly set bndry value
    void setHi(int dir, int bc_val) { bc[BL_SPACEDIM+dir] = bc_val; }

//  return bndry values (used in calls to FORTRAN)
    const int* vect() const { return bc; } 

/// access functions

//  return low-end boundary data
    const int* lo() const { return bc; }

//  return high-end boundary data
    const int* hi() const { return bc+BL_SPACEDIM; }

//  return low-end boundary data in direction <dir>
    int lo(int dir) const { return bc[dir]; }

//  return high-end boundary data in direction <dir>
    int hi(int dir) const { return bc[BL_SPACEDIM+dir]; }

/// io functions
//  pipe boundary data into ostream
    friend ostream& operator << (ostream&, const BCRec& );
};

///  function for setting a BC
void setBC(const BOX & bx, 
           const BOX & domain,
	   const BCRec & bc_dom, 
           BCRec & bcr);

///  function for setting array of BCs
void setBC(const BOX & bx, 
           const BOX & domain, 
	   int src_comp, 
           int dest_comp, 
           int ncomp,
	   const Array<BCRec> & bc_dom, 
           Array<BCRec> & bcr);
    
#endif
