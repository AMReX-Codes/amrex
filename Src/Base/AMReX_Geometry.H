#ifndef AMREX_GEOMETRY_H_
#define AMREX_GEOMETRY_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Geometry.H
 * \brief Problem-domain geometry: maps between index space and physical space.
 */

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
#include <numbers>

namespace amrex {

class MultiFab;
class DistributionMapping;
class BoxArray;

/**
 * \struct amrex::GeometryData
 * \ingroup amrex_real_space
 * \brief GPU-safe snapshot of Geometry data suitable for capture in GPU kernels.
 */
struct GeometryData
{
 //! Returns the cellsize for each coordinate direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real* CellSize () const noexcept { return dx; }
    //! Returns the cellsize for direction \p dir.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real CellSize (int dir) const noexcept { return dx[dir]; }
 //! Returns the lo end of the problem domain in each dimension.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real* ProbLo () const noexcept { return prob_domain.lo(); }
    //! Returns the lo end of the problem domain in direction \p dir.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ProbLo (int dir) const noexcept { return prob_domain.lo(dir); }
 //! Returns the hi end of the problem domain in each dimension.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real* ProbHi () const noexcept { return prob_domain.hi(); }
    //! Returns the hi end of the problem domain in direction \p dir.
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
    RealBox prob_domain; //!< Physical bounds of the problem domain.
    Box domain; //!< Cell-centered index-space domain.
    Real dx[AMREX_SPACEDIM]; //!< Cell size in each dimension.
    int is_periodic[AMREX_SPACEDIM]; //!< For each dimension, 0 if the domain is non-periodic and 1 if it is.
    int coord; //!< Coordinate-system type: 0 for Cartesian, 1 for RZ, and 2 for spherical.
};

/**
 * \class Geometry
 * \ingroup amrex_real_space
 * \brief Rectangular problem domain geometry.
 *
 * This class describes problem domain and coordinate system for
 * RECTANGULAR problem domains.  Since the problem domain is RECTANGULAR,
 * periodicity is meaningful.
 */
class Geometry
    :
    public CoordSys
{
public:
    /** \brief Construct a copy of the default Geometry, if one is available.
     *
     * If there is no AMReX instance, or the default Geometry has not been set
     * up yet, the object is left in an unusable state and a define method must
     * be called before use.
     */
    Geometry () noexcept;

    /**
     * \brief Construct a Geometry for the specified index-space domain.
     *
     * \param dom Cell-centered Box specifying the index-space domain.
     * \param rb Physical bounds of the domain. If nullptr, the bounds of the
     *           default Geometry are used (see Setup()).
     * \param coord Coordinate-system type: 0 for Cartesian
     *              (CoordSys::cartesian), 1 for RZ (CoordSys::RZ), or 2 for
     *              spherical (CoordSys::SPHERICAL). Any other value selects
     *              the coordinate system of the default Geometry.
     * \param is_per Periodicity flags for each coordinate direction. If
     *               nullptr, the flags of the default Geometry are used.
     */
    explicit Geometry (const Box&     dom,
                       const RealBox* rb     =  nullptr,
                       int            coord  = -1,
                       int const*     is_per =  nullptr);
    /**
     * \brief Construct a Geometry with explicitly supplied physical bounds and
     * periodicity.
     *
     * \param dom Cell-centered Box specifying the index-space domain.
     * \param rb Physical bounds of the domain.
     * \param coord Coordinate-system type: 0 for Cartesian
     *              (CoordSys::cartesian), 1 for RZ (CoordSys::RZ), or 2 for
     *              spherical (CoordSys::SPHERICAL). Any other value selects
     *              the coordinate system of the default Geometry.
     * \param is_per Periodicity flags for each coordinate direction.
     */
    Geometry (const Box& dom, const RealBox& rb, int coord,
              Array<int,AMREX_SPACEDIM> const& is_per);

    ~Geometry () = default;
    Geometry (const Geometry& rhs) = default;
    Geometry (Geometry&& rhs) noexcept = default;
    Geometry& operator= (const Geometry& rhs) = default;
    Geometry& operator= (Geometry&& rhs) noexcept = default;

    //! Return a GPU-safe snapshot of the Geometry data.
    [[nodiscard]] GeometryData data() const noexcept {
        return {
            .prob_domain = prob_domain,
            .domain = domain,
            .dx = {AMREX_D_DECL(dx[0],dx[1],dx[2])},
            .is_periodic = {AMREX_D_DECL(is_periodic[0], is_periodic[1], is_periodic[2])},
            .coord = static_cast<int>(c_sys)
        };
    }

    /**
     * \brief Set up the default Geometry, unless it has already been set up.
     *
     * A Geometry constructed or defined without explicit physical bounds,
     * coordinate system, or periodicity takes those from the default Geometry,
     * which this function initializes. Anything not passed in here is read
     * from the ParmParse prefix `geometry`.
     *
     * \param rb Problem-domain bounds. If nullptr, the bounds are read from
     *           `geometry.prob_lo` (default zero) together with either
     *           `geometry.prob_hi` or `geometry.prob_extent`, exactly one of
     *           which must be provided.
     * \param coord Coordinate-system identifier. If outside the range [0,2],
     *              it is read from `geometry.coord_sys` and defaults to
     *              Cartesian.
     * \param isper Periodicity flags. If nullptr, they are read from
     *              `geometry.is_periodic` and default to non-periodic.
     */
    static void Setup (const RealBox* rb = nullptr, int coord = -1, int const* isper = nullptr);
    //! Set the problem domain of the default Geometry to \p rb.
    static void ResetDefaultProbDomain (const RealBox& rb) noexcept;
    //! Set the periodicity flags of the default Geometry to \p is_per.
    static void ResetDefaultPeriodicity (const Array<int,AMREX_SPACEDIM>& is_per) noexcept;
    //! Set the coordinate-system type of the default Geometry to \p coord.
    static void ResetDefaultCoord (int coord) noexcept;

    /**
     * \brief Define this Geometry for the specified index-space domain.
     *
     * \param dom Cell-centered Box specifying the index-space domain.
     * \param rb Physical bounds of the domain. If nullptr, the bounds of the
     *           default Geometry are used (see Setup()).
     * \param coord Coordinate-system type: 0 for Cartesian
     *              (CoordSys::cartesian), 1 for RZ (CoordSys::RZ), or 2 for
     *              spherical (CoordSys::SPHERICAL). Any other value selects
     *              the coordinate system of the default Geometry.
     * \param is_per Periodicity flags for each coordinate direction. If
     *               nullptr, the flags of the default Geometry are used.
     */
    void define (const Box& dom, const RealBox* rb = nullptr, int coord = -1, int const* is_per = nullptr);

    /**
     * \brief Define this Geometry with explicitly supplied physical bounds and
     * periodicity.
     *
     * \param dom Cell-centered Box specifying the index-space domain.
     * \param rb Physical bounds of the domain.
     * \param coord Coordinate-system type: 0 for Cartesian
     *              (CoordSys::cartesian), 1 for RZ (CoordSys::RZ), or 2 for
     *              spherical (CoordSys::SPHERICAL). Any other value selects
     *              the coordinate system of the default Geometry.
     * \param is_per Periodicity flags for each coordinate direction.
     */
    void define (const Box& dom, const RealBox& rb, int coord, Array<int,AMREX_SPACEDIM> const& is_per);

    //! Returns the problem domain.
    [[nodiscard]] const RealBox& ProbDomain () const noexcept { return prob_domain; }
    //! Sets the problem domain.
    void ProbDomain (const RealBox& rb)
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

    //! Return the lo end of the problem domain in a GpuArray for device code.
    [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> ProbLoArray () const noexcept {
        return {{AMREX_D_DECL(prob_domain.lo(0),prob_domain.lo(1),prob_domain.lo(2))}};
    }

    //! Return the hi end of the problem domain in a GpuArray for device code.
    [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> ProbHiArray () const noexcept {
        return {{AMREX_D_DECL(prob_domain.hi(0),prob_domain.hi(1),prob_domain.hi(2))}};
    }

    //! Return the lo end of the roundoff domain in ParticleReal precision.
    [[nodiscard]] GpuArray<ParticleReal,AMREX_SPACEDIM> ProbLoArrayInParticleReal () const noexcept {
        return roundoff_lo;
    }

    //! Return the hi end of the roundoff domain in ParticleReal precision.
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
    void Domain (const Box& bx)
    {
        AMREX_ASSERT(bx.cellCentered());
        domain = bx;
        computeRoundoffDomain();
    }
    /**
     * \brief Define \p vol on \p grds with \p dm and fill it with cell volumes.
     *
     * The result has one component and \p grow ghost cells in every direction.
     */
    void GetVolume (MultiFab&       vol,
                    const BoxArray& grds,
                    const DistributionMapping& dm,
                    int             grow) const;
    //! Fill all valid and ghost cells of the pre-defined \p vol with cell volumes.
    void GetVolume (MultiFab&       vol) const;

    /**
     * \brief Resize \p vol to one component on entry \p idx of \p grds grown
     * by \p grow in every direction, and fill it with cell volumes.
     */
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

            constexpr Real pi = std::numbers::pi_v<Real>;
            vol = (4.0_rt / 3.0_rt) * pi * dx[0] * (rl * rl + rl * rr + rr * rr);
        }

#elif AMREX_SPACEDIM == 2

        auto coord = geomdata.Coord();

        if (coord == CoordSys::cartesian) {
            // Cartesian

            vol = dx[0] * dx[1];
        }
        else if (coord == CoordSys::RZ) {
            // Cylindrical

            Real r_l = geomdata.ProbLo()[0] + static_cast<Real>(point[0]) * dx[0];
            Real r_r = geomdata.ProbLo()[0] + static_cast<Real>(point[0]+1) * dx[0];

            constexpr Real pi = std::numbers::pi_v<Real>;
            vol = pi * (r_l + r_r) * dx[0] * dx[1];
        }
        else {
            // Spherical

            Real r_l = geomdata.ProbLo()[0] + static_cast<Real>(point[0]) * dx[0];
            Real r_r = geomdata.ProbLo()[0] + static_cast<Real>(point[0]+1) * dx[0];

            Real theta_l = geomdata.ProbLo()[1] + static_cast<Real>(point[1]) * dx[1];
            Real theta_r = geomdata.ProbLo()[1] + static_cast<Real>(point[1]+1) * dx[1];

            constexpr Real twoThirdsPi = static_cast<Real>((2.0 / 3.0) * std::numbers::pi_v<Real>);
            vol = twoThirdsPi * (std::cos(theta_l) - std::cos(theta_r)) * dx[0] *
                (r_r*r_r + r_r*r_l + r_l*r_l);
        }

#else

        amrex::ignore_unused(point);

        // Cartesian

        vol = dx[0] * dx[1] * dx[2];

#endif

        return vol;
    }

    /**
    * \brief Define \p dloga on \p grds with \p dm and fill it with
    * d(log(A))/dr at cell centers.
    *
    * The result has one component and \p grow ghost cells in every direction.
    * \p dir is the coordinate direction of the face area A.
    *
    * \note This function is only available in 1D and 2D.
    */
    void GetDLogA (MultiFab&       dloga,
                   const BoxArray& grds,
                   const DistributionMapping& dm,
                   int             dir,
                   int             grow) const;
    /**
    * \brief Define \p area on the faces of \p grds normal to \p dir, using
    * \p dm, and fill it with face areas.
    *
    * The result has one component and \p grow ghost cells in every direction.
    */
    void GetFaceArea (MultiFab&       area,
                      const BoxArray& grds,
                      const DistributionMapping& dm,
                      int             dir,
                      int             grow) const;
    /**
    * \brief Fill all valid and ghost cells of the pre-defined \p area with face
    * areas normal to \p dir.
    */
    void GetFaceArea (MultiFab&       area,
                      int             dir) const;

    /**
     * \brief Resize \p area to one component on the faces normal to \p dir of
     * entry \p idx of \p grds grown by \p grow in every direction, and fill it
     * with face areas.
     */
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
    //! Return the periodicity flags for all directions as an Array.
    [[nodiscard]] Array<int,AMREX_SPACEDIM> isPeriodic () const noexcept {
        return {{AMREX_D_DECL(static_cast<int>(is_periodic[0]),
                              static_cast<int>(is_periodic[1]),
                              static_cast<int>(is_periodic[2]))}};
    }
    //! Return the periodicity flags in a GpuArray for device code.
    [[nodiscard]] GpuArray<int,AMREX_SPACEDIM> isPeriodicArray () const noexcept {
        return {{AMREX_D_DECL(static_cast<int>(is_periodic[0]),
                              static_cast<int>(is_periodic[1]),
                              static_cast<int>(is_periodic[2]))}};
    }
    //! What's period in specified direction?
    [[nodiscard]] int period (int dir) const noexcept { BL_ASSERT(is_periodic[dir]); return domain.length(dir); }

    //! Return the Periodicity based on the length of the domain.
    [[nodiscard]] Periodicity periodicity () const noexcept {
        return Periodicity(IntVect(AMREX_D_DECL(domain.length(0) * is_periodic[0],
                                          domain.length(1) * is_periodic[1],
                                          domain.length(2) * is_periodic[2])));
    }

    //! Return the Periodicity based on the length of the cell-centered box \p b.
    [[nodiscard]] Periodicity periodicity (const Box& b) const noexcept {
        AMREX_ASSERT(b.cellCentered());
        return Periodicity(IntVect(AMREX_D_DECL(b.length(0) * is_periodic[0],
                                                b.length(1) * is_periodic[1],
                                                b.length(2) * is_periodic[2])));
    }

    /**
    * \brief Compute Array of shifts which will translate \p src so that it will
    * intersect \p target with non-zero intersection. The \p out array will be
    * resized internally, so anything previously there will be gone.
    * DO NOT return non-periodic shifts, even if the boxes do
    * intersect without shifting.  The logic is that you will only do
    * this as a special case if there is some periodicity. Each returned IntVect
    * is a nonzero index-space shift whose component is zero in a non-periodic
    * direction and zero or plus/minus the domain length in a periodic
    * direction.
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

    /**
     * \brief Set the periodicity flags and return the old flags.
     *
     * Note that, unlike the Periodicity class, the flags are just boolean. The
     * caller can use the returned array to restore the previous state after a
     * temporary override.
     */
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

    /**
     * \brief Coarsen the Geometry by \p rr.
     *
     * This does not refresh roundoff_lo/hi. Call computeRoundoffDomain() before
     * using particle roundoff bounds from the modified Geometry.
     */
    void coarsen (IntVect const& rr) {
        domain.coarsen(rr);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            dx[i] = (ProbHi(i) - ProbLo(i)) / static_cast<Real>(domain.length(i));
            inv_dx[i] = Real(1.)/dx[i];
        }
    }

    /**
     * \brief Refine the Geometry by \p rr.
     *
     * This does not refresh roundoff_lo/hi. Call computeRoundoffDomain() before
     * using particle roundoff bounds from the modified Geometry.
     */
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

    //! Compute the roundoff domain.
    void computeRoundoffDomain ();

    //! Returns roundoff domain's lower end
    [[nodiscard]] GpuArray<ParticleReal, AMREX_SPACEDIM> const&
    RoundOffLo () const { return roundoff_lo; }

    //! Returns roundoff domain's higher end
    [[nodiscard]] GpuArray<ParticleReal, AMREX_SPACEDIM> const&
    RoundOffHi () const { return roundoff_hi; }

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

    friend std::istream& operator>> (std::istream&, Geometry&);
};


/**
 * \ingroup amrex_real_space
 * \brief Write a Geometry to an output stream in ASCII format.
 */
std::ostream& operator<< (std::ostream&, const Geometry&);
/**
 * \ingroup amrex_real_space
 * \brief Read a Geometry from an input stream in ASCII format.
 */
std::istream& operator>> (std::istream&, Geometry&);

/**
 * \brief Return a Geometry coarsened by \p rr.
 *
 * The returned Geometry's roundoff_lo/hi are not refreshed. Call
 * computeRoundoffDomain() before using particle roundoff bounds from it.
 */
inline
Geometry
coarsen (Geometry const& fine, IntVect const& rr) {
    Geometry r{fine};
    r.coarsen(rr);
    return r;
}

/**
 * \brief Return a Geometry coarsened by \p rr.
 *
 * The returned Geometry's roundoff_lo/hi are not refreshed. Call
 * computeRoundoffDomain() before using particle roundoff bounds from it.
 */
inline
Geometry
coarsen (Geometry const& fine, int rr) { return coarsen(fine, IntVect(rr)); }

/**
 * \brief Return a Geometry refined by \p rr.
 *
 * The returned Geometry's roundoff_lo/hi are not refreshed. Call
 * computeRoundoffDomain() before using particle roundoff bounds from it.
 */
inline
Geometry
refine (Geometry const& crse, IntVect const& rr) {
    Geometry r{crse};
    r.refine(rr);
    return r;
}

/**
 * \brief Return a Geometry refined by \p rr.
 *
 * The returned Geometry's roundoff_lo/hi are not refreshed. Call
 * computeRoundoffDomain() before using particle roundoff bounds from it.
 */
inline
Geometry
refine (Geometry const& crse, int rr) { return refine(crse, IntVect(rr)); }

/**
 * \ingroup amrex_real_space
 * \brief Return the default Geometry held by the top-level AMReX instance.
 */
inline
const Geometry&
DefaultGeometry () {
    return *(AMReX::top()->getDefaultGeometry());
}

}

#endif /*_GEOMETRY_H_*/
