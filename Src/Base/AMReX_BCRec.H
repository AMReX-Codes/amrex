
#ifndef AMREX_BCREC_H_
#define AMREX_BCREC_H_

#include <AMReX_Box.H>
#include <AMReX_BC_TYPES.H>

namespace amrex {
/**
* \brief Boundary Condition Records.
* Necessary information and functions for computing boundary conditions.
*/

// This class has standard layout.  And we should keep it so!

class BCRec
{
public:
    /**
    * \brief The default constructor, which does NOT set valid boundary types.
    */
    BCRec ();
    /**
    * \brief The constructor.
    */
    BCRec (AMREX_D_DECL(int loX, int loY, int loZ),
           AMREX_D_DECL(int hiX, int hiY, int hiZ));
    /**
    * \brief Another constructor.
    */
    BCRec (const int* lo,
           const int* hi);
    /*
    * \brief Yet another constructor.  Inherits bndry types from bc_domain
    * when bx lies on edge of domain otherwise gets interior Dirchlet.
    */
    BCRec (const Box&   bx,
           const Box&   domain,
           const BCRec& bc_domain);
    /*
    * \brief Explicitly set lo bndry value.
    */
    void setLo (int dir, int bc_val) { bc[dir] = bc_val; }
    /**
    * \brief Explicitly set hi bndry value.
    */
    void setHi (int dir, int bc_val) { bc[AMREX_SPACEDIM+dir] = bc_val; }
    /**
    * \brief  Return bndry values (used in calls to FORTRAN).
    */
    const int* vect () const& { return bc; }
    const int* vect () && = delete;

    const int* data () const& { return bc; }
    const int* data () && = delete;

    int vectSize () const { return (2 * AMREX_SPACEDIM); }
    void setVect (const int *bcvect) { for(int i(0); i < vectSize(); ++i) { bc[i] = bcvect[i]; } }
    /**
    * \brief Return low-end boundary data.
    */
    const int* lo () const& { return bc; }
    const int* lo () && = delete;
    /**
    * \brief Return high-end boundary data.
    */
    const int* hi () const& { return bc+AMREX_SPACEDIM; }
    const int* hi () && = delete;
    /**
    * \brief Return low-end boundary data in direction <dir>.
    */
    int lo (int dir) const { return bc[dir]; }
    /**
    * \brief Return high-end boundary data in direction <dir>.
    */
    int hi (int dir) const { return bc[AMREX_SPACEDIM+dir]; }
    /**
    * \brief ASCII write to ostream.
    */
    friend std::ostream& operator << (std::ostream&, const BCRec&);

private:
    /**
    * \brief Array of integer values describing boundary conditions.
    */
    int bc[2*AMREX_SPACEDIM];
};

    /**
    * \brief Function for setting a BC.
    */
    void setBC (const Box&   bx, 
                const Box&   domain,
                const BCRec& bc_dom, 
                BCRec&       bcr);
    /**
    * \brief Function for setting array of BCs.
    */
    void setBC (const Box&          bx, 
                const Box&          domain, 
                int                 src_comp, 
                int                 dest_comp, 
                int                 ncomp,
                const Vector<BCRec>& bc_dom, 
                Vector<BCRec>&       bcr);
}

#endif /*_BCREC_H_*/
