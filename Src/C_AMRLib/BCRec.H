//BL_COPYRIGHT_NOTICE

#ifndef _BCREC_H_
#define _BCREC_H_ 

//
// $Id: BCRec.H,v 1.5 1998-03-30 21:57:32 lijewski Exp $
//

#include <Misc.H>
#include <Box.H>
#include <BC_TYPES.H>

//
// Necessary information and functions for computing boundary conditions.
//
class BCRec
{
public:
#ifdef __GNUG__
     bool operator== (const BCRec& rhs) const { return true; }
#endif 
    //
    // Sets to invalid bndry types.
    //
    BCRec () {}
    //
    // Constructor.
    //
    BCRec (D_DECL(int loX, int loY, int loZ),
           D_DECL(int hiX, int hiY, int hiZ));
    //
    // Constructor.
    //
    BCRec (const int* lo,
           const int* hi);
    //
    // Inherits bndry types from bc_domain when bx lies on edge of
    // domain otherwise gets interior Dirchlet.
    //
    BCRec (const Box&   bx,
           const Box&   domain,
           const BCRec& bc_domain);
    //
    // Explicitly set bndry value.
    //
    void setLo (int dir,
                int bc_val)
    {
        bc[dir] = bc_val;
    }
    //
    // Explicitly set bndry value.
    //
    void setHi (int dir,
                int bc_val)
    {
        bc[BL_SPACEDIM+dir] = bc_val;
    }
    //
    // Return bndry values (used in calls to FORTRAN).
    //
    const int* vect () const { return bc; } 
    //
    // Return low-end boundary data.
    //
    const int* lo () const { return bc; }
    //
    // Return high-end boundary data.
    //
    const int* hi () const { return bc+BL_SPACEDIM; }
    //
    // Return low-end boundary data in direction <dir>.
    //
    int lo (int dir) const { return bc[dir]; }
    //
    // Return high-end boundary data in direction <dir>.
    //
    int hi (int dir) const { return bc[BL_SPACEDIM+dir]; }
    //
    // Pipe boundary data into ostream.
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
