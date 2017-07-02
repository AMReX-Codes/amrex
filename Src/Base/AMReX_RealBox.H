
#ifndef _RealBox_H_
#define _RealBox_H_ 

#include <iosfwd>
#include <array>

#include <AMReX_Array.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>

namespace amrex {

//!A Box with real dimensions.  A RealBox is OK iff volume > 0.

class RealBox
{
public:
    //! The default constructor.  Builds invalid RealBox.
    RealBox ();

    //! Construct region given diagonal points.
    RealBox (const Real* lo,
             const Real* hi);

    RealBox (const std::array<Real,BL_SPACEDIM>& lo,
             const std::array<Real,BL_SPACEDIM>& hi);

    /**
    * \brief Construct region given index box, cell spacing
    * and physical location of index (0,0,0).
    */
    RealBox (const Box&  bx,
             const Real* dx,
             const Real* base);
    //! Explicit dimension specific constructors.
    RealBox (AMREX_D_DECL(Real x0, Real y0, Real z0),
             AMREX_D_DECL(Real x1, Real y1, Real z1));

    void initialize_device_memory();
    void copy_device_memory();
    void copy_xlo();
    void copy_xhi();

    //! Set epsilon used by contains() function when comparing floating point numbers.
    static void setEpsilon (Real epsilon) { eps = epsilon; }
    //! Returns epsilon.
    static Real epsilon () { return eps; }
    //! Returns lo side.
    const Real* lo () const& { return xlo; }

    const Real* loF () const& {
#ifdef CUDA
	return (Real*) Device::get_host_pointer(xlo_d.get());
#else
	return xlo;
#endif
    }
    const Real* lo () && = delete;
    //! Returns hide side.
    const Real* hi () const& { return xhi; }

    const Real* hiF () const& {
#ifdef CUDA
	return xhi_d.get();
#else
	return xhi;
#endif
    }
    const Real* hi () && = delete;
    //! Returns length in specified direction.
    Real lo (int dir) const { return xlo[dir]; }
    //! Returns hi side in specified direction.
    Real hi (int dir) const { return xhi[dir]; }
    //! Returns length in specified direction.
    Real length (int dir) const { return xhi[dir]-xlo[dir]; }
    //! Sets lo side.
    void setLo (const Real* lo) { initialize_device_memory(); AMREX_D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]); copy_device_memory(); }
    //! Sets lo side.
    void setLo (const Array<Real>& lo) { initialize_device_memory(); AMREX_D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]); copy_device_memory(); }
    //! Sets lo side in specified direction.
    void setLo (int dir, Real lo) { BL_ASSERT(dir >= 0 && dir < BL_SPACEDIM); initialize_device_memory(); xlo[dir] = lo; copy_device_memory(); }
    //! Sets hi side.
    void setHi (const Real* hi) { initialize_device_memory(); AMREX_D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]); copy_device_memory(); }
    //! Sets hi side.
    void setHi (const Array<Real>& hi) { initialize_device_memory(); AMREX_D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]); copy_device_memory(); }
    //! Sets hi side in specified direction.
    void setHi (int dir, Real hi) { BL_ASSERT(dir >= 0 && dir < BL_SPACEDIM); initialize_device_memory(); xhi[dir] = hi; copy_device_memory(); }
    //! Is the RealBox OK; i.e. does it have positive volume?
    bool ok () const;
    //! Is the specified point contained in the RealBox?
    bool contains (const Real* point) const;
    //! Is the specified RealBox contained in this RealBox?
    bool contains (const RealBox& bx) const;

private:
    //
    // Static data.
    //
    static Real eps;
    //
    // The data.
    //
    Real xlo[BL_SPACEDIM];
    Real xhi[BL_SPACEDIM];

#ifdef CUDA
    std::shared_ptr<Real> xlo_d;
    std::shared_ptr<Real> xhi_d;
#endif

};

//
// Nice ASCII output.
//
std::ostream& operator<< (std::ostream&, const RealBox&);
//
// Nice ASCII input.
//
std::istream& operator>> (std::istream&, RealBox&);

#ifdef CUDA
    // \cond CODEGEN
    class RealBox_init
    {
    public:
        RealBox_init ();
        ~RealBox_init ();
    private:
        static int m_cnt;
    };
    // \endcond
#endif

}

#ifdef CUDA
static amrex::RealBox_init file_scope_RealBox_init_object;
#endif

#endif /*_RealBox_H_*/
