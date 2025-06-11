#ifndef AMREX_EDGE_FLUX_REGISTER_H_
#define AMREX_EDGE_FLUX_REGISTER_H_
#include <AMReX_Config.H>

#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>

namespace amrex {

/**
 * Edge Flux Register for Constrained Transport
 *
 * This Flux Register is useful for solving system like dB/dt + curl E = 0
 * on a staggered mesh. (Here d is of course partial derivation.) B is a
 * vector on cell faces, and E is a vector on cell edges. In 2D, E has only
 * one component, Ez, and it is on the nodes of a 2d mesh.
 *
 * At the beginning of a coarse step, `reset()` is called.  In MFIter for
 * the coarse level advance, `CrseAdd` is called with coarse flux (i.e., E).
 * The flux is not scaled.  In MFIter for the fine level advance, `FineAdd`
 * is called.  After the fine level finishes its time steps, `Reflux` is
 * called to update the coarse level B on the coarse/fine boundary. The user
 * is also expected to call this version of average_down_faces from
 * AMReX_MultiFabUtil.H to synchronize the coarse level data with the fine
 * level.
 *
 * \vertbatim
    template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int>>
    void average_down_faces (const Array<const MF*,AMREX_SPACEDIM>& fine,
                             const Array<MF*,AMREX_SPACEDIM>& crse,
                             const IntVect& ratio, const Geometry& crse_geom)
 * \endverbatim
 *
 * Note that both CrseAdd and FineAdd are async for GPU builds. That means
 * it's the user's responsibility to keep the FArrayBox arguments alive or
 * call Gpu::streamSynchronize() when necessary.
 *
 * Because staggered grids are used, tiling could be very confusing. To avoid
 * confusion, this class assumes that tiling is not enabled for the MFIter
 * loop containing calls to CrseAdd and FineAdd.
 *
 * If the equation has an extra factor due to the choice of units, the
 * factor can be absorbed into dt. If we have `v x B` instead of E, the sign
 * can also been absorbed into dt. Note that whatever the choice of sign is,
 * the dt arguments passed to CrseAdd and FineAdd should have the same sign.
 *
 * We try to keep the interface simple by not providing overloads that
 * specify the component index. If the user's data does not start with
 * component 0, it can be worked around by creating alias FArrayBox and
 * MultiFab.
 */
class EdgeFluxRegister
{
public:

    EdgeFluxRegister () = default;

    EdgeFluxRegister (const BoxArray& fba, const BoxArray& cba,
                      const DistributionMapping& fdm, const DistributionMapping& cdm,
                      const Geometry& fgeom, const Geometry& cgeom,
                      int nvar = 1);

    void define (const BoxArray& fba, const BoxArray& cba,
                 const DistributionMapping& fdm, const DistributionMapping& cdm,
                 const Geometry& fgeom, const Geometry& cgeom,
                 int nvar = 1);

    void reset ();

#if (AMREX_SPACEDIM == 3)

    void CrseAdd (MFIter const& mfi, const Array<FArrayBox const*,3>& E_crse, Real dt_crse);
    void FineAdd (MFIter const& mfi, const Array<FArrayBox const*,3>& E_fine, Real dt_fine);

#else /* 2D */

    void CrseAdd (MFIter const& mfi, FArrayBox const& E_crse, Real dt_crse);
    void FineAdd (MFIter const& mfi, FArrayBox const& E_fine, Real dt_fine);

#endif

    void Reflux (Array<MultiFab*,AMREX_SPACEDIM> const& B_crse) const;

private:

    Geometry m_fine_geom;
    Geometry m_crse_geom;

    IntVect m_ratio;
    int m_ncomp;

#if (AMREX_SPACEDIM == 3)

    Array<MultiFab,AMREX_SPACEDIM> m_E_crse; // on original grids

    // There are AMREX_SPACEDIM*2 faces. For each face, we need to store two
    // component. For example, at the x-faces, we need to store Ey and Ez.
    Array<Array<MultiFab,2>,AMREX_SPACEDIM*2> m_E_fine;

    // Mask on the coarse level indicating overlap with m_E_fine
    Array<iMultiFab,AMREX_SPACEDIM> m_fine_mask;

#else

    MultiFab m_E_crse;
    Array<MultiFab,AMREX_SPACEDIM*2> m_E_fine;
    iMultiFab m_fine_mask;

#endif

    LayoutData<int> m_has_cf; // Flag on the coarse level indicating c/f interface
};

}

#endif
