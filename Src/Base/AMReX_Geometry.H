#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <iosfwd>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef BL_USE_UPCXX
#include <AMReX_BLPgas.H>
#endif

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

class Geometry
    :
    public CoordSys
{
public:
    //! The default constructor.
    Geometry ();
    //! Constructor taking the rectangular domain.
    explicit Geometry (const Box&     dom,
		       const RealBox* rb     =  nullptr,
		       int            coord  = -1,
		       int*           is_per =  nullptr);

    ~Geometry () = default;
    Geometry (const Geometry& rhs) = default;
    Geometry (Geometry&& rhs) noexcept = default;
    Geometry& operator= (const Geometry& rhs) = default;
    Geometry& operator= (Geometry&& rhs) = default;

    static void Finalize ();
    //! Read static values from ParmParse database.
    static void Setup (const RealBox* rb = nullptr, int coord = -1, int* is_per = nullptr);
    //! Set the rectangular domain after using default constructor.
    void define (const Box& dom, const RealBox* rb = 0, int coord = -1, int* is_per = 0);
    //! Returns the problem domain.
    static const RealBox& ProbDomain () { return prob_domain; }
    //! Sets the problem domain.
    static void ProbDomain (const RealBox& rb) 
    { 
        prob_domain = rb; 
        SetOffset(rb.lo());
    }
    //! Returns the lo end of the problem domain in each dimension.
    static const Real* ProbLo ()  { return prob_domain.lo(); }
    //! Returns the hi end of the problem domain in each dimension.
    static const Real* ProbHi () { return prob_domain.hi(); }
    //! Returns the lo end of the problem domain in specified direction.
    static Real ProbLo (int dir) { return prob_domain.lo(dir); }
    //! Returns the hi end of the problem domain in specified direction.
    static Real ProbHi (int dir) { return prob_domain.hi(dir); }
    //! Returns the overall size of the domain by multiplying the ProbLength's together
    static Real ProbSize ()
    {
        return AMREX_D_TERM(prob_domain.length(0),*prob_domain.length(1),*prob_domain.length(2));
    }
    //! Returns length of problem domain in specified dimension.
    static Real ProbLength (int dir) { return prob_domain.length(dir); }
    //! Returns our rectangular domain.
    const Box& Domain () const { return domain; }
    //! Sets our rectangular domain.
    void Domain (const Box& bx) { domain = bx; }
    //! Define a multifab of areas and volumes with given grow factor.
    void GetVolume (MultiFab&       vol,
                    const BoxArray& grds,
                    const DistributionMapping& dm,
                    int             grow) const;
    //! Fill the pre-built multifab with volume
    void GetVolume (MultiFab&       vol) const;

    void GetVolume (FArrayBox&       vol,
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
    static bool isPeriodic (int dir) { return is_periodic[dir] != 0; }
    //! Is domain periodic in any direction?
    static bool isAnyPeriodic ()
    {
        return AMREX_D_TERM(isPeriodic(0),||isPeriodic(1),||isPeriodic(2));
    }
    //! Is domain periodic in all directions?
    static bool isAllPeriodic ()
    {
        return AMREX_D_TERM(isPeriodic(0),&&isPeriodic(1),&&isPeriodic(2));
    }
    //! What's period in specified direction?
    int period (int dir) const { BL_ASSERT(is_periodic[dir]); return domain.length(dir); }

    Periodicity periodicity () const {
	return Periodicity(IntVect(AMREX_D_DECL(domain.length(0) * is_periodic[0],
					  domain.length(1) * is_periodic[1],
					  domain.length(2) * is_periodic[2])));
    }

    static Periodicity periodicity (const Box& b) {
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
                        Vector<IntVect>& out) const;

#ifdef BL_USE_MPI
    //! Send Geometry to sidecar procs.
    static void SendGeometryToSidecar (Geometry *geom, int whichSidecar);
    //
    static void BroadcastGeometry (Geometry &geom, int fromProc, MPI_Comm comm);
    static void BroadcastGeometry (Geometry &geom, int fromProc, MPI_Comm comm, bool bcastSource);
#endif

private:
    //! Helper functions.
    void read_params ();
    //
    // Static data.
    //
    static int     spherical_origin_fix;
    static bool    is_periodic[AMREX_SPACEDIM]; // 0 means not periodic
    static RealBox prob_domain;
    // Non-static data.
    //
    Box domain;
};


//! Nice ASCII output.
std::ostream& operator<< (std::ostream&, const Geometry&);
//! Nice ASCII input.
std::istream& operator>> (std::istream&, Geometry&);

}

#endif /*_GEOMETRY_H_*/
