//BL_COPYRIGHT_NOTICE

#ifndef _BCREC_H_
#define _BCREC_H_ 

//
// $Id: BCRec.H,v 1.6 1998-05-06 17:50:45 lijewski Exp $
//

#include <Misc.H>
#include <Box.H>
#include <BC_TYPES.H>

//
//@Man:
//@Memo: Boundary Condition Records.
/*@Doc:

  Necessary information and functions for computing boundary conditions.
*/

class BCRec
{
public:

#ifdef __GNUG__
    //
    // Here only to work around GNU bug.
    //
    bool operator== (const BCRec& rhs) const { return true; }
#endif 
    //
    //@ManDoc: The default constructor -- sets to invalid boundary types.
    //
    BCRec () {}
    //
    //@ManDoc: The constructor.
    //
    BCRec (D_DECL(int loX, int loY, int loZ),
           D_DECL(int hiX, int hiY, int hiZ));
    //
    //@ManDoc: Another constructor.
    //
    BCRec (const int* lo,
           const int* hi);

    /*@ManDoc: Yet another constructor.  Inherits bndry types from bc\_domain
               when bx lies on edge of domain otherwise gets interior Dirchlet.
    */
    BCRec (const Box&   bx,
           const Box&   domain,
           const BCRec& bc_domain);
    //
    //@ManDoc: Explicitly set lo bndry value.
    //
    void setLo (int dir,
                int bc_val)
    {
        bc[dir] = bc_val;
    }
    //
    //@ManDoc: Explicitly set hi bndry value.
    //
    void setHi (int dir,
                int bc_val)
    {
        bc[BL_SPACEDIM+dir] = bc_val;
    }
    //
    //@ManDoc: Return bndry values (used in calls to FORTRAN).
    //
    const int* vect () const { return bc; } 
    //
    //@ManDoc: Return low-end boundary data.
    //
    const int* lo () const { return bc; }
    //
    //@ManDoc: Return high-end boundary data.
    //
    const int* hi () const { return bc+BL_SPACEDIM; }
    //
    //@ManDoc: Return low-end boundary data in direction <dir>.
    //
    int lo (int dir) const { return bc[dir]; }
    //
    //@ManDoc: Return high-end boundary data in direction <dir>.
    //
    int hi (int dir) const { return bc[BL_SPACEDIM+dir]; }
    //
    //@ManDoc: ASCII write to ostream.
    //
    friend ostream& operator << (ostream&, const BCRec& );

protected:
    //
    // Array of integer values describing boundary conditions.
    //
    int bc[2*BL_SPACEDIM];
};

//
// Function for setting a BC.
//
void setBC (const Box&   bx, 
            const Box&   domain,
            const BCRec& bc_dom, 
            BCRec&       bcr);

//
// Function for setting array of BCs.
//
void setBC (const Box&          bx, 
            const Box&          domain, 
            int                 src_comp, 
            int                 dest_comp, 
            int                 ncomp,
            const Array<BCRec>& bc_dom, 
            Array<BCRec>&       bcr);

#endif /*_BCREC_H_*/
