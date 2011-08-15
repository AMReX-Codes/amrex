
#ifndef _BCREC_H_
#define _BCREC_H_ 

#include <Box.H>
#include <BC_TYPES.H>
//
// Boundary Condition Records.
//
// Necessary information and functions for computing boundary conditions.
//

class BCRec
{
public:
    //
    // The default constructor -- sets to invalid boundary types.
    //
    BCRec ();
    //
    // The constructor.
    //
    BCRec (D_DECL(int loX, int loY, int loZ),
           D_DECL(int hiX, int hiY, int hiZ));
    //
    // Another constructor.
    //
    BCRec (const int* lo,
           const int* hi);
    //
    // Yet another constructor.  Inherits bndry types from bc_domain
    // when bx lies on edge of domain otherwise gets interior Dirchlet.
    //
    BCRec (const Box&   bx,
           const Box&   domain,
           const BCRec& bc_domain);
    //
    // Explicitly set lo bndry value.
    //
    void setLo (int dir,
                int bc_val);
    //
    // Explicitly set hi bndry value.
    //
    void setHi (int dir,
                int bc_val);
    //
    // Return bndry values (used in calls to FORTRAN).
    //
    const int* vect () const;
    //
    // Return low-end boundary data.
    //
    const int* lo () const;
    //
    // Return high-end boundary data.
    //
    const int* hi () const;
    //
    // Return low-end boundary data in direction <dir>.
    //
    int lo (int dir) const;
    //
    // Return high-end boundary data in direction <dir>.
    //
    int hi (int dir) const;
    //
    // ASCII write to ostream.
    //
    friend std::ostream& operator << (std::ostream&, const BCRec&);

protected:
    //
    // Array of integer values describing boundary conditions.
    //
    int bc[2*BL_SPACEDIM];
};

inline
void
BCRec::setLo (int dir,
              int bc_val)
{
    bc[dir] = bc_val;
}

inline
void
BCRec::setHi (int dir,
              int bc_val)
{
    bc[BL_SPACEDIM+dir] = bc_val;
}

inline
const int*
BCRec::vect () const
{
    return bc;
} 

inline
const int*
BCRec::lo () const
{
    return bc;
}

inline
const int*
BCRec::hi () const
{
    return bc+BL_SPACEDIM;
}

inline
int
BCRec::lo (int dir) const
{
    return bc[dir];
}

inline
int
BCRec::hi (int dir) const
{
    return bc[BL_SPACEDIM+dir];
}

namespace BoxLib
{
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
}

#endif /*_BCREC_H_*/
