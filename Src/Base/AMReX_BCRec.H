
#ifndef AMREX_BCREC_H_
#define AMREX_BCREC_H_
#include <AMReX_Config.H>

#include <AMReX_Box.H>
#include <AMReX_BC_TYPES.H>

/**
 * \file AMReX_BCRec.H
 * \brief Boundary condition descriptors used throughout AMReX.
 */

namespace amrex {
/**
 * \class amrex::BCRec
 * \ingroup amrex_utilities
 * \brief Boundary condition record storing low/high face types per coordinate.
 *
 * `BCRec` has standard layout — this is a deliberate design constraint that
 * must be preserved for compatibility with Fortran code, host/device
 * transfers, and `std::memcpy`-based operations.
 */
class BCRec
{
public:
    /** \brief Default constructor leaves all face types set to `BCType::bogus`. */
    AMREX_GPU_HOST_DEVICE
    BCRec () noexcept

        = default;
    /**
     * \brief Construct a record with explicit low/high face types.
     *
     * The low-face values are provided first, followed by the high-face values.
     *
     * \param loX Low-face boundary type in the x direction.
     * \param loY Low-face boundary type in the y direction (used when \c AMREX_SPACEDIM \>= 2).
     * \param loZ Low-face boundary type in the z direction (used when \c AMREX_SPACEDIM \>= 3).
     * \param hiX High-face boundary type in the x direction.
     * \param hiY High-face boundary type in the y direction (used when \c AMREX_SPACEDIM \>= 2).
     * \param hiZ High-face boundary type in the z direction (used when \c AMREX_SPACEDIM \>= 3).
     */
    AMREX_GPU_HOST_DEVICE
    BCRec (AMREX_D_DECL(int loX, int loY, int loZ),
           AMREX_D_DECL(int hiX, int hiY, int hiZ)) noexcept
        : bc {AMREX_D_DECL(loX,loY,loZ),
              AMREX_D_DECL(hiX,hiY,hiZ)}
        {}
    /**
     * \brief Construct from arrays of low/high face types.
     *
     * \param a_lo Pointer to \c AMREX_SPACEDIM low-face entries.
     * \param a_hi Pointer to \c AMREX_SPACEDIM high-face entries.
     */
    AMREX_GPU_HOST_DEVICE
    BCRec (const int* a_lo, const int* a_hi) noexcept
        : bc {AMREX_D_DECL(a_lo[0],a_lo[1],a_lo[2]),
              AMREX_D_DECL(a_hi[0],a_hi[1],a_hi[2])}
        {}
    /**
     * \brief Initialize boundary faces relative to a domain record.
     *
     * Faces coincident with the domain inherit the boundary from \p bc_domain;
     * interior faces default to `BCType::int_dir`.
     *
     * \param bx Box whose faces inherit domain or interior boundaries.
     * \param domain Domain box describing the physical region.
     * \param bc_domain Boundary record providing domain face types.
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
    /**
     * \brief Set the low-face boundary for dimension \p dir.
     * \param dir Coordinate direction.
     * \param bc_val Boundary type code.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setLo (int dir, int bc_val) noexcept { bc[dir] = bc_val; }
    /**
     * \brief Set the high-face boundary for dimension \p dir.
     * \param dir Coordinate direction.
     * \param bc_val Boundary type code.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setHi (int dir, int bc_val) noexcept { bc[AMREX_SPACEDIM+dir] = bc_val; }
    /**
     * \brief Explicitly set the boundary value for a given face.
     * \param face Orientation describing the face.
     * \param bc_val Boundary type code.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set (Orientation face, int bc_val) noexcept {
        if (face.isLow()) {
            setLo(face.coordDir(), bc_val);
        } else {
            setHi(face.coordDir(), bc_val);
        }
    }
    /** \brief Return pointer to the internal boundary array (low faces, then high faces). Legacy name for Fortran callback compatibility. */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const int* vect () const& noexcept{ return bc; }
    const int* vect () && = delete;

    /** \brief Equivalent to vect(); returns the same pointer to the internal boundary array. */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const int* data () const& noexcept { return bc; }
    const int* data () && = delete;

    /** \brief Return pointer to the low-end face array. */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const int* lo () const& noexcept { return bc; }
    const int* lo () && = delete;
    /** \brief Return pointer to the high-end face array. */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const int* hi () const& noexcept { return bc+AMREX_SPACEDIM; }
    const int* hi () && = delete;
    /** \brief Return the low-end boundary code for dimension \p dir. */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int lo (int dir) const noexcept { return bc[dir]; }
    /** \brief Return the high-end boundary code for dimension \p dir. */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int hi (int dir) const noexcept { return bc[AMREX_SPACEDIM+dir]; }
    /** \brief Compare two records for identical boundary condition values. */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator== (const BCRec& rhs) const noexcept {
        bool retval = true;
        for (int i = 0; i < 2*AMREX_SPACEDIM && retval; i++)
        {
            retval &= bc[i] == rhs.bc[i];
        }
        return retval;
    }
    /** \brief Logical negation of operator==. */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator!= (const BCRec& rhs) const noexcept { return !(*this == rhs); }
    /**
    * \brief ASCII write to ostream.
    */
    friend std::ostream& operator << (std::ostream&, const BCRec&);

private:
    /**
     * \brief Array of integer values describing boundary conditions.
     */
    int bc[2*AMREX_SPACEDIM]{AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus),
              AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus)};
};

/**
 * \brief Reset a `BCRec` based on a domain record and box location.
 *
 * \param bx Box whose boundary faces are being initialized.
 * \param domain Domain box describing the physical region.
 * \param bc_dom Domain boundary record to inherit from.
 * \param bcr Record to populate.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
 * \brief Populate a span of `BCRec` objects based on a domain record.
 *
 * \param bx Box whose boundary faces are being initialized.
 * \param domain Domain box describing the physical region.
 * \param src_comp Starting component index in \p bc_dom.
 * \param dest_comp Starting component index in \p bcr.
 * \param ncomp Number of components to fill.
 * \param bc_dom Domain boundary records to inherit from.
 * \param bcr Output boundary records to populate.
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
