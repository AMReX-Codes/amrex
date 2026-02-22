
#ifndef AMREX_RealBox_H_
#define AMREX_RealBox_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_RealBox.H
 * \brief A box in real (floating-point) space.
 */

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_RealVect.H>
#include <AMReX_SPACE.H>

#include <iosfwd>
#include <array>

//NOLINTBEGIN(clang-analyzer-security.ArrayBound)

namespace amrex {

/**
 * \class amrex::RealBox
 * \ingroup amrex_real_space
 * \brief A box in real (floating-point) space.
 *
 * A RealBox is valid (ok()) if and only if its volume is non-negative.
 */
class RealBox
{
public:
    //! The default constructor.  Builds invalid RealBox.
    constexpr RealBox () noexcept = default;

    /**
     * \brief Construct region given diagonal points.
     * \param a_lo Pointer to the lower corner coordinates.
     * \param a_hi Pointer to the upper corner coordinates.
     */
    AMREX_GPU_HOST_DEVICE
    RealBox (const Real* a_lo, const Real* a_hi) noexcept
        : xlo{AMREX_D_DECL(a_lo[0],a_lo[1],a_lo[2])}, xhi{AMREX_D_DECL(a_hi[0],a_hi[1],a_hi[2])} {}

    /**
     * \brief Construct region given lo and hi corners as std::array.
     * \param a_lo Lower corner coordinates.
     * \param a_hi Upper corner coordinates.
     */
    RealBox (const std::array<Real,AMREX_SPACEDIM>& a_lo,
             const std::array<Real,AMREX_SPACEDIM>& a_hi) noexcept;

    /**
    * \brief Construct region given index box, cell spacing
    * and physical location of index (0,0,0).
    *
    * \param bx   Integer AMR box.
    * \param dx   Cell spacing along each dimension.
    * \param base Physical location of cell (0,0,0).
    */
    RealBox (const Box&  bx, const Real* dx, const Real* base) noexcept;
    //! Explicit dimension specific constructors.
    AMREX_GPU_HOST_DEVICE
    RealBox (AMREX_D_DECL(Real x0, Real y0, Real z0),
             AMREX_D_DECL(Real x1, Real y1, Real z1)) noexcept
        : xlo{AMREX_D_DECL(x0,y0,z0)}, xhi{AMREX_D_DECL(x1,y1,z1)} {}
    //! Return pointer to the lo-side coordinate array.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    const Real* lo () const& noexcept { return xlo; }
    /// \cond DOXYGEN_IGNORE
    AMREX_GPU_HOST_DEVICE
    const Real* lo () && = delete;
    /// \endcond
    //! Return pointer to the hi-side coordinate array.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    const Real* hi () const& noexcept { return xhi; }
    /// \cond DOXYGEN_IGNORE
    AMREX_GPU_HOST_DEVICE
    const Real* hi () && = delete;
    /// \endcond
    //! Return lo-side coordinate in direction \p dir.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Real lo (int dir) const noexcept { return xlo[dir]; }
    //! Return hi-side coordinate in direction \p dir.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Real hi (int dir) const noexcept { return xhi[dir]; }
    //! Returns length in specified direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Real length (int dir) const noexcept { return xhi[dir]-xlo[dir]; }
    //! Sets lo side to the coordinates in \p a_lo.
    void setLo (const Real* a_lo) noexcept { AMREX_D_EXPR(xlo[0] = a_lo[0], xlo[1] = a_lo[1], xlo[2] = a_lo[2]); }
    //! Sets lo side to the coordinates in \p a_lo.
    void setLo (const Vector<Real>& a_lo) noexcept { AMREX_D_EXPR(xlo[0] = a_lo[0], xlo[1] = a_lo[1], xlo[2] = a_lo[2]); }
    //! Sets lo side to the coordinates in \p a_lo.
    void setLo (const RealVect& a_lo) noexcept { AMREX_D_EXPR(xlo[0] = a_lo[0], xlo[1] = a_lo[1], xlo[2] = a_lo[2]); }
    //! Sets lo side in direction \p dir to the value \p a_lo.
    void setLo (int dir, Real a_lo) noexcept { BL_ASSERT(dir >= 0 && dir < AMREX_SPACEDIM); xlo[dir] = a_lo; }
    //! Sets hi side to the coordinates in \p a_hi.
    void setHi (const Real* a_hi) noexcept { AMREX_D_EXPR(xhi[0] = a_hi[0], xhi[1] = a_hi[1], xhi[2] = a_hi[2]); }
    //! Sets hi side to the coordinates in \p a_hi.
    void setHi (const Vector<Real>& a_hi) noexcept { AMREX_D_EXPR(xhi[0] = a_hi[0], xhi[1] = a_hi[1], xhi[2] = a_hi[2]); }
    //! Sets hi side to the coordinates in \p a_hi.
    void setHi (const RealVect& a_hi) noexcept { AMREX_D_EXPR(xhi[0] = a_hi[0], xhi[1] = a_hi[1], xhi[2] = a_hi[2]); }
    //! Sets hi side in direction \p dir to the value \p a_hi.
    void setHi (int dir, Real a_hi) noexcept { BL_ASSERT(dir >= 0 && dir < AMREX_SPACEDIM); xhi[dir] = a_hi; }
    //! Is the RealBox OK; i.e. does it have non-negative volume?
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool ok () const noexcept {
        return (length(0) >= 0.0)
#if (AMREX_SPACEDIM > 1)
            && (length(1) >= 0.0)
#endif
#if (AMREX_SPACEDIM > 2)
            && (length(2) >= 0.0)
#endif
            ;
    }

    //! Returns the volume of the RealBox. If this RealBox is invalid,
    //! it's volume is considered to be zero.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Real volume () const noexcept {
        if (ok()) { return AMREX_D_TERM(length(0), *length(1), *length(2)); }
        return 0.0;
    }

    /**
     * \brief Return true if the specified point lies inside the box.
     * \param point Cartesian coordinates to test.
     * \param eps   Optional tolerance applied symmetrically to each face.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool contains (const Real* point,
                   Real eps = 0.0) const noexcept {
        return  AMREX_D_TERM((xlo[0]-eps < point[0]) && (point[0] < xhi[0]+eps),
                          && (xlo[1]-eps < point[1]) && (point[1] < xhi[1]+eps),
                          && (xlo[2]-eps < point[2]) && (point[2] < xhi[2]+eps));
    }

    /**
     * \brief Return true if the specified point lies inside the box.
     * \param point Cartesian coordinates to test.
     * \param eps   Optional tolerance applied symmetrically to each face.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool contains (XDim3 point, Real eps = 0.0) const noexcept {
        return  AMREX_D_TERM((xlo[0]-eps < point.x) && (point.x < xhi[0]+eps),
                          && (xlo[1]-eps < point.y) && (point.y < xhi[1]+eps),
                          && (xlo[2]-eps < point.z) && (point.z < xhi[2]+eps));
    }

    /**
     * \brief Return true if the specified RealVect lies inside the box.
     * \param rv  Point expressed as a RealVect.
     * \param eps Optional tolerance applied symmetrically to each face.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool contains (const RealVect& rv,
                   Real eps=0.0) const noexcept { return contains(rv.dataPtr(), eps); }

    /**
     * \brief Return true if the specified RealBox is fully contained in this box.
     * \param rb  Candidate RealBox to test.
     * \param eps Optional tolerance applied symmetrically to each face.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool contains (const RealBox& rb,
                   Real eps = 0.0) const noexcept {
        return contains(rb.xlo, eps) && contains(rb.xhi, eps);
    }

    /**
     * \brief Return true if the specified RealBox intersects this box.
     * \param bx Other RealBox to test.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    bool intersects (const RealBox& bx) const noexcept {
        return  ! (AMREX_D_TERM((xlo[0] > bx.xhi[0]) || (xhi[0] < bx.xlo[0]),
                             || (xlo[1] > bx.xhi[1]) || (xhi[1] < bx.xlo[1]),
                             || (xlo[2] > bx.xhi[2]) || (xhi[2] < bx.xlo[2])));
    }

private:
    //
    // The data.
    //
    Real xlo[AMREX_SPACEDIM] = {AMREX_D_DECL(Real(0.),Real(0.),Real(0.))};
    Real xhi[AMREX_SPACEDIM] = {AMREX_D_DECL(Real(-1.),Real(-1.),Real(-1.))};
};

//! Write a RealBox to an output stream in ASCII format.
std::ostream& operator<< (std::ostream&, const RealBox&);

//! Read a RealBox from an input stream in ASCII format.
std::istream& operator>> (std::istream&, RealBox&);

/**
 * \brief Return true if two RealBoxes are equal to within tolerance \p eps.
 * \param box1 First RealBox.
 * \param box2 Second RealBox.
 * \param eps  Absolute tolerance applied to each coordinate independently.
 */
bool AlmostEqual (const RealBox& box1,
                  const RealBox& box2,
                  Real eps = 0.0) noexcept;

}

//NOLINTEND(clang-analyzer-security.ArrayBound)

#endif /*_RealBox_H_*/
