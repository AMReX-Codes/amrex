
#ifndef AMREX_RealBox_H_
#define AMREX_RealBox_H_

#include <iosfwd>
#include <array>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_RealVect.H>

namespace amrex {

//!A Box with real dimensions.  A RealBox is OK iff volume >= 0.

class RealBox // : public Gpu::Managed
{
public:
    //! The default constructor.  Builds invalid RealBox.
    AMREX_GPU_HOST_DEVICE
    RealBox ();

    //! Construct region given diagonal points.
    AMREX_GPU_HOST_DEVICE
    RealBox (const Real* lo,
             const Real* hi);

    RealBox (const std::array<Real,AMREX_SPACEDIM>& lo,
             const std::array<Real,AMREX_SPACEDIM>& hi);

    /**
    * \brief Construct region given index box, cell spacing
    * and physical location of index (0,0,0).
    */
    AMREX_GPU_HOST_DEVICE
    RealBox (const Box&  bx,
             const Real* dx,
             const Real* base);
    //! Explicit dimension specific constructors.
    AMREX_GPU_HOST_DEVICE
    RealBox (AMREX_D_DECL(Real x0, Real y0, Real z0),
             AMREX_D_DECL(Real x1, Real y1, Real z1));
    //! Returns lo side.
    AMREX_GPU_HOST_DEVICE
    const Real* lo () const& { return xlo; }
    AMREX_GPU_HOST_DEVICE
    const Real* lo () && = delete;
    //! Returns hide side.
    AMREX_GPU_HOST_DEVICE
    const Real* hi () const& { return xhi; }
    AMREX_GPU_HOST_DEVICE
    const Real* hi () && = delete;
    //! Returns length in specified direction.
    AMREX_GPU_HOST_DEVICE
    Real lo (int dir) const { return xlo[dir]; }
    //! Returns hi side in specified direction.
    AMREX_GPU_HOST_DEVICE
    Real hi (int dir) const { return xhi[dir]; }
    //! Returns length in specified direction.
    AMREX_GPU_HOST_DEVICE
    Real length (int dir) const { return xhi[dir]-xlo[dir]; }
    //! Sets lo side.
    void setLo (const Real* a_lo) { AMREX_D_EXPR(xlo[0] = a_lo[0], xlo[1] = a_lo[1], xlo[2] = a_lo[2]); }
    //! Sets lo side.
    void setLo (const Vector<Real>& a_lo) { AMREX_D_EXPR(xlo[0] = a_lo[0], xlo[1] = a_lo[1], xlo[2] = a_lo[2]); }
    //! Sets lo side in specified direction.
    void setLo (int dir, Real a_lo) { BL_ASSERT(dir >= 0 && dir < AMREX_SPACEDIM); xlo[dir] = a_lo; }
    //! Sets hi side.
    void setHi (const Real* a_hi) { AMREX_D_EXPR(xhi[0] = a_hi[0], xhi[1] = a_hi[1], xhi[2] = a_hi[2]); }
    //! Sets hi side.
    void setHi (const Vector<Real>& a_hi) { AMREX_D_EXPR(xhi[0] = a_hi[0], xhi[1] = a_hi[1], xhi[2] = a_hi[2]); }
    //! Sets hi side in specified direction.
    void setHi (int dir, Real a_hi) { BL_ASSERT(dir >= 0 && dir < AMREX_SPACEDIM); xhi[dir] = a_hi; }
    //! Is the RealBox OK; i.e. does it have non-negative volume?
    AMREX_GPU_HOST_DEVICE
    bool ok () const;
    //! Returns the volume of the RealBox. If this RealBox is invalid, 
    //! it's volume is considered to be zero.
    AMREX_GPU_HOST_DEVICE
    Real volume () const;
    //! Is the specified point contained in the RealBox?
    AMREX_GPU_HOST_DEVICE
    bool contains (const Real* point, Real eps=0.0) const;
    AMREX_GPU_HOST_DEVICE
    bool contains (const RealVect& rv, Real eps=0.0) const { return contains(rv.dataPtr(), eps); }
    //! Is the specified RealBox contained in this RealBox?
    AMREX_GPU_HOST_DEVICE
    bool contains (const RealBox& bx, Real eps=0.0) const;

    //! Does the specified RealBox intersect with this RealBox?
    AMREX_GPU_HOST_DEVICE
    bool intersects (const RealBox& bx) const;
    
private:
    //
    // The data.
    //
    Real xlo[AMREX_SPACEDIM];
    Real xhi[AMREX_SPACEDIM];
};

//
// Nice ASCII output.
//
std::ostream& operator<< (std::ostream&, const RealBox&);
//
// Nice ASCII input.
//
std::istream& operator>> (std::istream&, RealBox&);

}

#endif /*_RealBox_H_*/
