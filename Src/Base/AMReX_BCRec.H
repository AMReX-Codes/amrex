
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
    AMREX_GPU_HOST_DEVICE
    BCRec () noexcept
        : bc {AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus),
              AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus)}
        { }
    /**
    * \brief The constructor.
    */
    AMREX_GPU_HOST_DEVICE
    BCRec (AMREX_D_DECL(int loX, int loY, int loZ),
           AMREX_D_DECL(int hiX, int hiY, int hiZ)) noexcept
        : bc {AMREX_D_DECL(loX,loY,loZ),
              AMREX_D_DECL(hiX,hiY,hiZ)}
        {}
    /**
    * \brief Another constructor.
    */
    AMREX_GPU_HOST_DEVICE
    BCRec (const int* a_lo, const int* a_hi) noexcept
        : bc {AMREX_D_DECL(a_lo[0],a_lo[1],a_lo[2]),
              AMREX_D_DECL(a_hi[0],a_hi[1],a_hi[2])}
        {}
    /*
    * \brief Yet another constructor.  Inherits bndry types from bc_domain
    * when bx lies on edge of domain otherwise gets interior Dirchlet.
    */
    AMREX_GPU_HOST_DEVICE
    BCRec (const Box&   bx,
           const Box&   domain,
           const BCRec& bc_domain) noexcept
    {
        const int* bxlo = bx.loVect();
        const int* bxhi = bx.hiVect();
        const int* dlo  = domain.loVect();
        const int* dhi  = domain.hiVect();
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            int ilo = dir;
            int ihi = dir+AMREX_SPACEDIM;
            bc[ilo] = ( bxlo[dir]<=dlo[dir] ? bc_domain.bc[ilo] : BCType::int_dir );
            bc[ihi] = ( bxhi[dir]>=dhi[dir] ? bc_domain.bc[ihi] : BCType::int_dir );
        }
    }
    /*
    * \brief Explicitly set lo bndry value.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void setLo (int dir, int bc_val) noexcept { bc[dir] = bc_val; }
    /**
    * \brief Explicitly set hi bndry value.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void setHi (int dir, int bc_val) noexcept { bc[AMREX_SPACEDIM+dir] = bc_val; }
    /**
    * \brief  Return bndry values (used in calls to FORTRAN).
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    const int* vect () const& noexcept{ return bc; }
    const int* vect () && = delete;

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    const int* data () const& noexcept { return bc; }
    const int* data () && = delete;

    /**
    * \brief Return low-end boundary data.
    */
    AMREX_GPU_HOST_DEVICE  AMREX_INLINE
    const int* lo () const& noexcept { return bc; }
    const int* lo () && = delete;
    /**
    * \brief Return high-end boundary data.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    const int* hi () const& noexcept { return bc+AMREX_SPACEDIM; }
    const int* hi () && = delete;
    /**
    * \brief Return low-end boundary data in direction \<dir\>.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int lo (int dir) const noexcept { return bc[dir]; }
    /**
    * \brief Return high-end boundary data in direction \<dir\>.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int hi (int dir) const noexcept { return bc[AMREX_SPACEDIM+dir]; }
    /**
    * \brief Equal test.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator== (const BCRec& rhs) const noexcept {
        bool retval = true;
        for (int i = 0; i < 2*AMREX_SPACEDIM && retval; i++)
        {
            retval &= bc[i] == rhs.bc[i];
        }
        return retval;
    }
    /**
    * \brief Not equal test.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator!= (const BCRec& rhs) const noexcept { return !(*this == rhs); }
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
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
setBC (const Box& bx, const Box& domain, const BCRec& bc_dom, BCRec& bcr) noexcept
{
    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        bcr.setLo(dir, ( bxlo[dir]<=dlo[dir] ? bc_dom.lo(dir) : BCType::int_dir ));
        bcr.setHi(dir, ( bxhi[dir]>=dhi[dir] ? bc_dom.hi(dir) : BCType::int_dir ));
    }
}

/**
 * \brief Function for setting array of BCs.
 */
void
setBC (const Box&           bx,
       const Box&           domain,
       int                  src_comp,
       int                  dest_comp,
       int                  ncomp,
       const Vector<BCRec>& bc_dom,
       Vector<BCRec>&       bcr) noexcept;
}

#endif /*_BCREC_H_*/
