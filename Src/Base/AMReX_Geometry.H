#ifndef AMREX_GEOMETRY_H_
#define AMREX_GEOMETRY_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_CoordSys.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>
#include <AMReX_Periodicity.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <iosfwd>
#include <map>

namespace amrex {
/**
* \class Geometry
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
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real* CellSize () const noexcept { return dx; }
    //
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real CellSize (int dir) const noexcept { return dx[dir]; }
 //! Returns the lo end of the problem domain in each dimension.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real* ProbLo () const noexcept { return prob_domain.lo(); }
    //
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ProbLo (int dir) const noexcept { return prob_domain.lo(dir); }
 //! Returns the hi end of the problem domain in each dimension.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real* ProbHi () const noexcept { return prob_domain.hi(); }
    //
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ProbHi (int dir) const noexcept { return prob_domain.hi(dir); }
 //! Returns our rectangular domain.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Box& Domain () const noexcept { return domain; }
//! Returns whether the domain is periodic in the given direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int isPeriodic (const int i) const noexcept { return is_periodic[i]; }
 //! Coordinates type
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int Coord () const noexcept { return coord; }

public:
    RealBox prob_domain;
    Box domain;
    Real dx[AMREX_SPACEDIM];
    int is_periodic[AMREX_SPACEDIM]; /**< For each dimension, 0 if the domain is non-periodic and 1 if it is. */
    int coord;
};

class Geometry
    :
    public CoordSys
{
public:
    /** \brief The default constructor.
     *
     *   Leaves object in an unusable state. A "define" method must be called before use.
     */
    Geometry () noexcept;

    /**
     * Constructs a fully-defined geometry object that describes the problem domain and coordinate system.
     *
     * \param dom cell-centered `Box` specifying the index space bounds of the domain.
     * \param rb A pointer to a `RealBox` specifying the real space bounds of the domain. Optional. If not specified,
     *        will be ParmParsed from "geometry.prob_lo" and "geometry.prob_hi".
     * \param coord Specifies the coordinate system type. Optional. Can be one of:
     *        CoordSys::cartesian
     *        CoordSys::RZ
     *        CoordSys::SPHERICAL
     *        Default is cartesian.
     * \param is_per pointer to memory space specifying periodicity in each coordinate direction.
     *        Optional. Default is non-periodic.
     */
    explicit Geometry (const Box&     dom,
                       const RealBox* rb     =  nullptr,
                       int            coord  = -1,
                       int const*     is_per =  nullptr) noexcept;
    /**
     * Constructs a fully-defined geometry object that describes the problem domain and coordinate system.
     *
     * \param dom cell-centered `Box` specifying the index space bounds of the domain.
     * \param rb A pointer to a `RealBox` specifying the real space bounds of the domain. Optional. If not specified,
     *        will be ParmParsed from "geometry.prob_lo" and "geometry.prob_hi".
     * \param coord Specifies the coordinate system type. Optional. Can be one of:
     *        CoordSys::cartesian
     *        CoordSys::RZ
     *        CoordSys::SPHERICAL
     *        Default is cartesian.
     * \param is_per amrex::Array specifying periodicity in each coordinate direction. Optional.
     *        Default is non-periodic.
     */
    Geometry (const Box& dom, const RealBox& rb, int coord,
              Array<int,AMREX_SPACEDIM> const& is_per) noexcept;

    ~Geometry () = default;
    Geometry (const Geometry& rhs) = default;
    Geometry (Geometry&& rhs) noexcept = default;
    Geometry& operator= (const Geometry& rhs) = default;
    Geometry& operator= (Geometry&& rhs) noexcept = default;

    //! Returns non-static copy of geometry's stored data.
    [[nodiscard]] GeometryData data() const noexcept {
        return { prob_domain, domain, {AMREX_D_DECL(dx[0],dx[1],dx[2])},
            {AMREX_D_DECL(is_periodic[0], is_periodic[1], is_periodic[2])},
                static_cast<int>(c_sys) };
    }

    //! Read static values from ParmParse database.
    static void Setup (const RealBox* rb = nullptr, int coord = -1, int const* isper = nullptr) noexcept;
    static void ResetDefaultProbDomain (const RealBox& rb) noexcept;
    static void ResetDefaultPeriodicity (const Array<int,AMREX_SPACEDIM>& is_per) noexcept;
    static void ResetDefaultCoord (int coord) noexcept;

    /**
     * Defines a geometry object that describes the problem domain and coordinate system for rectangular problem
     * domains. Useful for when the Geometry object has been constructed with the default constructor.
     *
     * \param dom cell-centered `Box` specifying the index space bounds of the domain.
     * \param rb A pointer to a `RealBox` specifying the real space bounds of the domain. Optional. If not specified,
     *        will be ParmParsed from "geometry.prob_lo" and "geometry.prob_hi".
     * \param coord Specifies the coordinate system type. Optional. Can be one of:
     *        CoordSys::cartesian
     *        CoordSys::RZ
     *        CoordSys::SPHERICAL
     *        Default is cartesian.
     * \param is_per pointer to memory space specifying periodicity in each coordinate direction.
     *        Optional. Default is non-periodic.
     */
    void define (const Box& dom, const RealBox* rb = nullptr, int coord = -1, int const* is_per = nullptr) noexcept;

    /**
     * Defines a geometry object that describes the problem domain and coordinate system for rectangular problem
     * domains. Useful for when the Geometry object has been constructed with the default constructor.
     *
     * \param dom cell-centered `Box` specifying the index space bounds of the domain.
     * \param rb A pointer to a `RealBox` specifying the real space bounds of the domain. Optional. If not specified,
     *        will be ParmParsed from "geometry.prob_lo" and "geometry.prob_hi".
     * \param coord Specifies the coordinate system type. Optional. Can be one of:
     *        CoordSys::cartesian
     *        CoordSys::RZ
     *        CoordSys::SPHERICAL
     *        Default is cartesian.
     * \param is_per amrex::Array specifying periodicity in each coordinate direction. Optional.
     *        Default is non-periodic.
     */
    void define (const Box& dom, const RealBox& rb, int coord, Array<int,AMREX_SPACEDIM> const& is_per) noexcept;

    //! Returns the problem domain.
    [[nodiscard]] const RealBox& ProbDomain () const noexcept { return prob_domain; }
    //! Sets the problem domain.
    void ProbDomain (const RealBox& rb) noexcept
    {
        prob_domain = rb;
        computeRoundoffDomain();
    }
    //! Returns the lo end of the problem domain in each dimension.
    [[nodiscard]] const Real* ProbLo () const noexcept { return prob_domain.lo(); }
    //! Returns the hi end of the problem domain in each dimension.
    [[nodiscard]] const Real* ProbHi () const noexcept { return prob_domain.hi(); }
    //! Returns the lo end of the problem domain in specified direction.
    [[nodiscard]] Real ProbLo (int dir) const noexcept { return prob_domain.lo(dir); }
    //! Returns the hi end of the problem domain in specified direction.
    [[nodiscard]] Real ProbHi (int dir) const noexcept { return prob_domain.hi(dir); }

    [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> ProbLoArray () const noexcept {
        return {{AMREX_D_DECL(prob_domain.lo(0),prob_domain.lo(1),prob_domain.lo(2))}};
    }

    [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> ProbHiArray () const noexcept {
        return {{AMREX_D_DECL(prob_domain.hi(0),prob_domain.hi(1),prob_domain.hi(2))}};
    }

    [[nodiscard]] GpuArray<ParticleReal,AMREX_SPACEDIM> ProbLoArrayInParticleReal () const noexcept {
        return roundoff_lo;
    }
    [[nodiscard]] GpuArray<ParticleReal,AMREX_SPACEDIM> ProbHiArrayInParticleReal () const noexcept {
        return roundoff_hi;
    }

    //! Returns the overall size of the domain by multiplying the ProbLength's together
    [[nodiscard]] Real ProbSize () const noexcept
    {
        return AMREX_D_TERM(prob_domain.length(0),*prob_domain.length(1),*prob_domain.length(2));
    }
    //! Returns length of problem domain in specified dimension.
    [[nodiscard]] Real ProbLength (int dir) const noexcept { return prob_domain.length(dir); }
    //! Returns our rectangular domain.
    [[nodiscard]] const Box& Domain () const noexcept { return domain; }
    //! Sets our rectangular domain.
    void Domain (const Box& bx) noexcept
    {
        domain = bx;
        computeRoundoffDomain();
    }
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

    //! Return the volume of the specified cell.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static Real Volume (const IntVect& point, const GeometryData& geomdata)
    {
        const auto *dx = geomdata.CellSize();

        Real vol;

#if AMREX_SPACEDIM == 1

        auto coord = geomdata.Coord();

        if (coord == CoordSys::cartesian) {
            // Cartesian

            vol = dx[0];
        }
        else {
            // Spherical

            Real rl = geomdata.ProbLo()[0] + static_cast<Real>(point[0]) * dx[0];
            Real rr = rl + dx[0];

            constexpr Real pi = Real(3.1415926535897932);
            vol = (4.0_rt / 3.0_rt) * pi * dx[0] * (rl * rl + rl * rr + rr * rr);
        }

#elif AMREX_SPACEDIM == 2

        auto coord = geomdata.Coord();

        if (coord == CoordSys::cartesian) {
            // Cartesian

            vol = dx[0] * dx[1];
        }
        else {
            // Cylindrical

            Real r_l = geomdata.ProbLo()[0] + static_cast<Real>(point[0]) * dx[0];
            Real r_r = geomdata.ProbLo()[0] + static_cast<Real>(point[0]+1) * dx[0];

            constexpr Real pi = Real(3.1415926535897932);
            vol = pi * (r_l + r_r) * dx[0] * dx[1];
        }

#else

        amrex::ignore_unused(point);

        // Cartesian

        vol = dx[0] * dx[1] * dx[2];

#endif

        return vol;
    }

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
    [[nodiscard]] bool isPeriodic (int dir) const noexcept { return is_periodic[dir]; }
    //! Is domain periodic in any direction?
    [[nodiscard]] bool isAnyPeriodic () const noexcept
    {
        return AMREX_D_TERM(isPeriodic(0),||isPeriodic(1),||isPeriodic(2));
    }
    //! Is domain periodic in all directions?
    [[nodiscard]] bool isAllPeriodic ()  const noexcept
    {
        return AMREX_D_TERM(isPeriodic(0),&&isPeriodic(1),&&isPeriodic(2));
    }
    [[nodiscard]] Array<int,AMREX_SPACEDIM> isPeriodic () const noexcept {
        return {{AMREX_D_DECL(static_cast<int>(is_periodic[0]),
                              static_cast<int>(is_periodic[1]),
                              static_cast<int>(is_periodic[2]))}};
    }
    [[nodiscard]] GpuArray<int,AMREX_SPACEDIM> isPeriodicArray () const noexcept {
        return {{AMREX_D_DECL(static_cast<int>(is_periodic[0]),
                              static_cast<int>(is_periodic[1]),
                              static_cast<int>(is_periodic[2]))}};
    }
    //! What's period in specified direction?
    [[nodiscard]] int period (int dir) const noexcept { BL_ASSERT(is_periodic[dir]); return domain.length(dir); }

    [[nodiscard]] Periodicity periodicity () const noexcept {
        return Periodicity(IntVect(AMREX_D_DECL(domain.length(0) * is_periodic[0],
                                          domain.length(1) * is_periodic[1],
                                          domain.length(2) * is_periodic[2])));
    }

    [[nodiscard]] Periodicity periodicity (const Box& b) const noexcept {
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
    [[nodiscard]] Box growNonPeriodicDomain (IntVect const& ngrow) const noexcept;
    //! Return domain box with non-periodic directions grown by ngrow.
    [[nodiscard]] Box growNonPeriodicDomain (int ngrow) const noexcept;
    //! Return domain box with periodic directions grown by ngrow.
    [[nodiscard]] Box growPeriodicDomain (IntVect const& ngrow) const noexcept;
    //! Return domain box with periodic directions grown by ngrow.
    [[nodiscard]] Box growPeriodicDomain (int ngrow) const noexcept;

    //! Set periodicity flags and return the old flags.
    //! Note that, unlike Periodicity class, the flags are just boolean.
    //!
    Array<int,AMREX_SPACEDIM>
    setPeriodicity (Array<int,AMREX_SPACEDIM> const& period) noexcept {
        Array<int,AMREX_SPACEDIM> r{{AMREX_D_DECL(is_periodic[0],
                                                  is_periodic[1],
                                                  is_periodic[2])}};
        AMREX_D_TERM(is_periodic[0] = period[0];,
                     is_periodic[1] = period[1];,
                     is_periodic[2] = period[2];);
        return r;
    }

    void coarsen (IntVect const& rr) {
        domain.coarsen(rr);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            dx[i] = (ProbHi(i) - ProbLo(i)) / static_cast<Real>(domain.length(i));
            inv_dx[i] = Real(1.)/dx[i];
        }
    }

    void refine (IntVect const& rr) {
        domain.refine(rr);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            dx[i] = (ProbHi(i) - ProbLo(i)) / static_cast<Real>(domain.length(i));
            inv_dx[i] = Real(1.)/dx[i];
        }
    }

    /**
    * \brief Returns true if a point is outside the roundoff domain.
    *        All particles with positions inside the roundoff domain
    *        are sure to be mapped to cells inside the Domain() box. Note that
    *        the same need not be true for all points inside ProbDomain().
    */
    [[nodiscard]] bool outsideRoundoffDomain (AMREX_D_DECL(ParticleReal x, ParticleReal y, ParticleReal z)) const;

    /**
    * \brief Returns true if a point is inside the roundoff domain.
    *        All particles with positions inside the roundoff domain
    *        are sure to be mapped to cells inside the Domain() box. Note that
    *        the same need not be true for all points inside ProbDomain().
    */
    [[nodiscard]] bool insideRoundoffDomain (AMREX_D_DECL(ParticleReal x, ParticleReal y, ParticleReal z)) const;

    /**
    * \brief Compute the roundoff domain. Public because it contains an
    *        extended host / device lambda.
    */
    void computeRoundoffDomain ();

private:
    void read_params ();

    // is_periodic and RealBox used to be static
    bool    is_periodic[AMREX_SPACEDIM] = {AMREX_D_DECL(false,false,false)};
    RealBox prob_domain;

    // Due to round-off errors, not all floating point numbers for which plo >= x < phi
    // will map to a cell that is inside "domain". "roundoff_{lo,hi}" store
    // positions that are very close to that in prob_domain, and for which all doubles and floats
    // between them will map to a cell inside domain.
    GpuArray<ParticleReal , AMREX_SPACEDIM> roundoff_lo, roundoff_hi;

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
