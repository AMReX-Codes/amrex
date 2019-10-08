#ifndef AMREX_GEOMETRY_H_
#define AMREX_GEOMETRY_H_

#include <iosfwd>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Array.H>
#include <AMReX_CoordSys.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>
#include <AMReX_Periodicity.H>

namespace amrex {
/**
* \brief Rectangular problem domain geometry.
*
* This class describes problem domain and coordinate system for
* RECTANGULAR problem domains.  Since the problem domain is RECTANGULAR,
* periodicity is meaningful.
*/

class MultiFab;
class DistributionMapping;
class BoxArray;

struct GeometryData
{
 //! Returns the cellsize for each coordinate direction.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real* CellSize () const noexcept { return dx; }
    //
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real CellSize (int dir) const noexcept { return dx[dir]; }
 //! Returns the lo end of the problem domain in each dimension.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real* ProbLo () const noexcept { return prob_domain.lo(); }
    //
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ProbLo (int dir) const noexcept { return prob_domain.lo(dir); }
 //! Returns the hi end of the problem domain in each dimension.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real* ProbHi () const noexcept { return prob_domain.hi(); }
    //
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ProbHi (int dir) const noexcept { return prob_domain.hi(dir); }
 //! Returns our rectangular domain.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Box& Domain () const noexcept { return domain; }
 //! Returns whether the domain is periodic in the given direction.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int isPeriodic (const int i) const noexcept { return is_periodic[i]; }
 //! Coordinates type
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int Coord () const noexcept { return coord; }

public:
    RealBox prob_domain;
    Box domain;
    Real dx[AMREX_SPACEDIM];
    int is_periodic[AMREX_SPACEDIM];
    int coord;
};

class Geometry
    :
    public CoordSys
{
public:
    //! The default constructor.
    Geometry () noexcept;
    //! Constructor taking the rectangular domain.
    explicit Geometry (const Box&     dom,
		       const RealBox* rb     =  nullptr,
		       int            coord  = -1,
		       int const*     is_per =  nullptr) noexcept;
    Geometry (const Box& dom, const RealBox& rb, int coord,
              Array<int,AMREX_SPACEDIM> const& is_per) noexcept;

    ~Geometry () = default;
    Geometry (const Geometry& rhs) = default;
    Geometry (Geometry&& rhs) noexcept = default;
    Geometry& operator= (const Geometry& rhs) = default;
    Geometry& operator= (Geometry&& rhs) noexcept = default;

    //! Returns non-static copy of geometry's stored data.
    GeometryData data() const noexcept {
        return { prob_domain, domain, {AMREX_D_DECL(dx[0],dx[1],dx[2])},
            {AMREX_D_DECL(is_periodic[0], is_periodic[1], is_periodic[2])},
                static_cast<int>(c_sys) };
    }

    //! Read static values from ParmParse database.
    static void Setup (const RealBox* rb = nullptr, int coord = -1, int const* is_per = nullptr) noexcept;
    static void ResetDefaultProbDomain (const RealBox& rb) noexcept;
    static void ResetDefaultPeriodicity (const Array<int,AMREX_SPACEDIM>& is_per) noexcept;
    static void ResetDefaultCoord (int coord) noexcept;

    //! Set the rectangular domain after using default constructor.
    void define (const Box& dom, const RealBox* rb = 0, int coord = -1, int const* is_per = nullptr) noexcept;
    void define (const Box& dom, const RealBox& rb, int coord, Array<int,AMREX_SPACEDIM> const& is_per) noexcept;
    //! Returns the problem domain.
    const RealBox& ProbDomain () const noexcept { return prob_domain; }
    //! Sets the problem domain.
    void ProbDomain (const RealBox& rb) noexcept
    {
        prob_domain = rb;
        SetOffset(rb.lo());
    }
    //! Returns the lo end of the problem domain in each dimension.
    const Real* ProbLo () const noexcept { return prob_domain.lo(); }
    //! Returns the hi end of the problem domain in each dimension.
    const Real* ProbHi () const noexcept { return prob_domain.hi(); }
    //! Returns the lo end of the problem domain in specified direction.
    Real ProbLo (int dir) const noexcept { return prob_domain.lo(dir); }
    //! Returns the hi end of the problem domain in specified direction.
    Real ProbHi (int dir) const noexcept { return prob_domain.hi(dir); }

    GpuArray<Real,AMREX_SPACEDIM> ProbLoArray () const noexcept {
        return {AMREX_D_DECL(prob_domain.lo(0),prob_domain.lo(1),prob_domain.lo(2))};
    }

    GpuArray<Real,AMREX_SPACEDIM> ProbHiArray () const noexcept {
        return {AMREX_D_DECL(prob_domain.hi(0),prob_domain.hi(1),prob_domain.hi(2))};
    }

    //! Returns the overall size of the domain by multiplying the ProbLength's together
    Real ProbSize () const noexcept
    {
        return AMREX_D_TERM(prob_domain.length(0),*prob_domain.length(1),*prob_domain.length(2));
    }
    //! Returns length of problem domain in specified dimension.
    Real ProbLength (int dir) const noexcept { return prob_domain.length(dir); }
    //! Returns our rectangular domain.
    const Box& Domain () const noexcept { return domain; }
    //! Sets our rectangular domain.
    void Domain (const Box& bx) noexcept { domain = bx; }
    //! Define a multifab of areas and volumes with given grow factor.
    void GetVolume (MultiFab&       vol,
                    const BoxArray& grds,
                    const DistributionMapping& dm,
                    int             grow) const;
    //! Fill the pre-built multifab with volume
    void GetVolume (MultiFab&       vol) const;

    void GetVolume (FArrayBox&      vol,
                    const BoxArray& grds,
                    int             idx,
                    int             grow) const;
    /**
    * \brief Compute d(log(A))/dr at cell centers in given region and
    *           stuff the results into the passed MultiFab.
    */
    void GetDLogA (MultiFab&       dloga,
                   const BoxArray& grds,
                   const DistributionMapping& dm,
                   int             dir,
                   int             grow) const;
    /**
    * \brief Compute area of cell faces in given region and stuff
    * stuff the results into the passed MultiFab.
    */
    void GetFaceArea (MultiFab&       area,
                      const BoxArray& grds,
		      const DistributionMapping& dm,
                      int             dir,
                      int             grow) const;
    //! Fill the pre-built multifab with area
    void GetFaceArea (MultiFab&       area,
                      int             dir) const;

    void GetFaceArea (FArrayBox&      area,
                      const BoxArray& grds,
                      int             idx,
                      int             dir,
                      int             grow) const;
    //! Is the domain periodic in the specified direction?
    bool isPeriodic (int dir) const noexcept { return is_periodic[dir]; }
    //! Is domain periodic in any direction?
    bool isAnyPeriodic () const noexcept
    {
        return AMREX_D_TERM(isPeriodic(0),||isPeriodic(1),||isPeriodic(2));
    }
    //! Is domain periodic in all directions?
    bool isAllPeriodic ()  const noexcept
    {
        return AMREX_D_TERM(isPeriodic(0),&&isPeriodic(1),&&isPeriodic(2));
    }
    Array<int,AMREX_SPACEDIM> isPeriodic () const noexcept {
        return {AMREX_D_DECL(static_cast<int>(is_periodic[0]),
                             static_cast<int>(is_periodic[1]),
                             static_cast<int>(is_periodic[2]))};
    }
    GpuArray<int,AMREX_SPACEDIM> isPeriodicArray () const noexcept {
        return {AMREX_D_DECL(static_cast<int>(is_periodic[0]),
                             static_cast<int>(is_periodic[1]),
                             static_cast<int>(is_periodic[2]))};
    }
    //! What's period in specified direction?
    int period (int dir) const noexcept { BL_ASSERT(is_periodic[dir]); return domain.length(dir); }

    Periodicity periodicity () const noexcept {
	return Periodicity(IntVect(AMREX_D_DECL(domain.length(0) * is_periodic[0],
					  domain.length(1) * is_periodic[1],
					  domain.length(2) * is_periodic[2])));
    }

    Periodicity periodicity (const Box& b) const noexcept {
        AMREX_ASSERT(b.cellCentered());
	return Periodicity(IntVect(AMREX_D_DECL(b.length(0) * is_periodic[0],
                                                b.length(1) * is_periodic[1],
                                                b.length(2) * is_periodic[2])));
    }

    /**
    * \brief Compute Array of shifts which will translate src so that it will
    * intersect target with non-zero intersection.  the array will be
    * resized internally, so anything previously there will be gone
    * DO NOT return non-periodic shifts, even if the box's do
    * intersect without shifting.  The logic is that you will only do
    * this as a special case if there is some periodicity.
    */
    void periodicShift (const Box&      target,
                        const Box&      src,
                        Vector<IntVect>& out) const noexcept;

    //! Return domain box with non-periodic directions grown by ngrow.
    Box growNonPeriodicDomain (int ngrow) const noexcept;
    //! Return domain box with periodic directions grown by ngrow.
    Box growPeriodicDomain (int ngrow) const noexcept;

    //! Set periodicity flags and return the old flags.
    //! Note that, unlike Periodicity class, the flags are just boolean.
    //!
    Array<int,AMREX_SPACEDIM>
    setPeriodicity (Array<int,AMREX_SPACEDIM> const& period) noexcept {
        Array<int,AMREX_SPACEDIM> r{AMREX_D_DECL(is_periodic[0],
                                                 is_periodic[1],
                                                 is_periodic[2])};
        AMREX_D_TERM(is_periodic[0] = period[0];,
                     is_periodic[1] = period[1];,
                     is_periodic[2] = period[2];);
        return r;
    }

    void coarsen (IntVect const& rr) {
        domain.coarsen(rr);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            dx[i] = (ProbHi(i) - ProbLo(i)) / domain.length(i);
            inv_dx[i] = 1./dx[i];
        }
    }

    void refine (IntVect const& rr) {
        domain.refine(rr);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            dx[i] = (ProbHi(i) - ProbLo(i)) / domain.length(i);
            inv_dx[i] = 1./dx[i];
        }
    }

private:
    void read_params ();

    // is_periodic and RealBox used to be static
    bool    is_periodic[AMREX_SPACEDIM] = {AMREX_D_DECL(false,false,false)};
    RealBox prob_domain;
    //
    Box     domain;
};


//! Nice ASCII output.
std::ostream& operator<< (std::ostream&, const Geometry&);
//! Nice ASCII input.
std::istream& operator>> (std::istream&, Geometry&);

inline
Geometry
coarsen (Geometry const& fine, IntVect const& rr) {
    Geometry r{fine};
    r.coarsen(rr);
    return r;
}

inline
Geometry
coarsen (Geometry const& fine, int rr) { return coarsen(fine, IntVect(rr)); }

inline
Geometry
refine (Geometry const& crse, IntVect const& rr) {
    Geometry r{crse};
    r.refine(rr);
    return r;
}

inline
Geometry
refine (Geometry const& crse, int rr) { return refine(crse, IntVect(rr)); }

inline
const Geometry&
DefaultGeometry () {
    return *(AMReX::top()->getDefaultGeometry());
}

}

#endif /*_GEOMETRY_H_*/
