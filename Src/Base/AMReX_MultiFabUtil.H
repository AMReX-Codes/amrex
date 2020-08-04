#ifndef AMREX_MultiFabUtil_H_
#define AMREX_MultiFabUtil_H_

#include <AMReX_BaseFabUtility.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFabUtil_C.H>

#include <AMReX_MultiFabUtilI.H>

namespace amrex
{
    //! Average nodal-based MultiFab onto cell-centered MultiFab.
    void average_node_to_cellcenter (MultiFab& cc, int dcomp,
                                     const MultiFab& nd, int scomp,
                                     int ncomp, int ngrow = 0);

    //! Average edge-based MultiFab onto cell-centered MultiFab. This fills in
    //! ngrow ghost cells in the cell-centered MultiFab. Both cell centered and
    //! edge centered MultiFabs need to have ngrow ghost values
    void average_edge_to_cellcenter (MultiFab& cc, int dcomp,
                                     const Vector<const MultiFab*>& edge,
                                     int ngrow = 0);

    //! Average face-based MultiFab onto cell-centered MultiFab.
    void average_face_to_cellcenter (MultiFab& cc, int dcomp,
                                     const Vector<const MultiFab*>& fc,
                                     int ngrow = 0);
    void average_face_to_cellcenter (MultiFab& cc, int dcomp,
                                     const Array<const MultiFab*,AMREX_SPACEDIM>& fc,
                                     int ngrow = 0);
    void average_face_to_cellcenter (MultiFab& cc,
                                     const Vector<const MultiFab*>& fc,
                                     const Geometry& geom);
    void average_face_to_cellcenter (MultiFab& cc,
                                     const Array<const MultiFab*,AMREX_SPACEDIM>& fc,
                                     const Geometry& geom);

    //! Average cell-centered MultiFab onto face-based MultiFab.
    void average_cellcenter_to_face (const Vector<MultiFab*>& fc,
                                     const MultiFab& cc,
                                     const Geometry& geom);
    void average_cellcenter_to_face (const Array<MultiFab*,AMREX_SPACEDIM>& fc,
                                     const MultiFab& cc,
                                     const Geometry& geom);

    //! Average fine face-based MultiFab onto crse face-based MultiFab.
    void average_down_faces (const Vector<const MultiFab*>& fine,
                             const Vector<MultiFab*>& crse,
                             const IntVect& ratio,
                             int ngcrse = 0);
    void average_down_faces (const Vector<const MultiFab*>& fine,
                             const Vector<MultiFab*>& crse,
                             int ratio,
                             int ngcrse = 0);
    void average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                             const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                             const IntVect& ratio,
                             int ngcrse = 0);
    void average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                             const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                             int ratio,
                             int ngcrse = 0);
    //! This version does average down for one direction.
    //! It uses the IndexType of MultiFabs to determine the direction.
    //! It is expected that one direction is nodal and the rest are cell-centered.
    void average_down_faces (const MultiFab& fine, MultiFab& crse,
                             const IntVect& ratio, int ngcrse=0);

    //  This version takes periodicity into account.
    void average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                             const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                             const IntVect& ratio, const Geometry& crse_geom);
    //  This version takes periodicity into account.
    void average_down_faces (const MultiFab& fine, MultiFab& crse,
                             const IntVect& ratio, const Geometry& crse_geom);

    //! Average fine edge-based MultiFab onto crse edge-based MultiFab.
    void average_down_edges (const Vector<const MultiFab*>& fine,
                             const Vector<MultiFab*>& crse,
                             const IntVect& ratio,
                             int ngcrse = 0);
    void average_down_edges (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                             const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                             const IntVect& ratio,
                             int ngcrse = 0);
    //! This version does average down for one direction.
    //! It uses the IndexType of MultiFabs to determine the direction.
    //! It is expected that one direction is cell-centered and the rest are nodal.
    void average_down_edges (const MultiFab& fine, MultiFab& crse,
                             const IntVect& ratio, int ngcrse=0);

    //! Average fine node-based MultiFab onto crse node-centered MultiFab.
    template <typename FAB>
    void average_down_nodal (const FabArray<FAB>& S_fine,
                             FabArray<FAB>& S_crse,
                             const IntVect& ratio,
                             int ngcrse = 0);

    //! Average fine cell-based MultiFab onto crse cell-centered MultiFab using
    //! volume-weighting. This routine DOES NOT assume that the crse BoxArray is
    //! a coarsened version of the fine BoxArray.
    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
                       const Geometry& fine_geom, const Geometry& crse_geom,
                       int scomp, int ncomp, const IntVect& ratio);
    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
                       const Geometry& fine_geom, const Geometry& crse_geom,
                       int scomp, int ncomp, const int ratio);

    //! Average MultiFab onto crse MultiFab without volume weighting. This
    //! routine DOES NOT assume that the crse BoxArray is a coarsened version of
    //! the fine BoxArray. Work for both cell-centered and nodal MultiFabs.
    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
                       int scomp, int ncomp, const IntVect& ratio);
    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
                       int scomp, int ncomp, int ratio);

    //! Add a coarsened version of the data contained in the S_fine MultiFab to
    //! S_crse, including ghost cells.
    void sum_fine_to_coarse (const MultiFab& S_Fine, MultiFab& S_crse,
                             int scomp, int ncomp,
                             const IntVect& ratio,
                             const Geometry& cgeom, const Geometry& fgeom);

    //! Output state data for a single zone
    void print_state (const MultiFab& mf, const IntVect& cell, const int n=-1,
                      const IntVect& ng = IntVect::TheZeroVector());

    //! Write each fab individually
    void writeFabs (const MultiFab& mf, const std::string& name);
    void writeFabs (const MultiFab& mf, int comp, int ncomp, const std::string& name);

    //! Extract a slice from the given cell-centered MultiFab at coordinate
    //! "coord" along direction "dir".
    std::unique_ptr<MultiFab> get_slice_data(int dir, Real coord,
                                             const MultiFab& cc,
                                             const Geometry& geom, int fstart, int ncomp,
                                             bool interpolate=false);

    //! Return an iMultiFab that has the same BoxArray and DistributionMapping
    //! as the coarse MultiFab cmf. Cells covered by the coarsened fine grids
    //! are set to fine_value, whereas other cells are set to crse_value.
    template <typename FAB>
    iMultiFab makeFineMask (const FabArray<FAB>& cmf, const BoxArray& fba, const IntVect& ratio,
                            int crse_value = 0, int fine_value = 1);
    iMultiFab makeFineMask (const BoxArray& cba, const DistributionMapping& cdm,
                            const BoxArray& fba, const IntVect& ratio,
                            int crse_value = 0, int fine_value = 1);
    template <typename FAB>
    iMultiFab makeFineMask (const FabArray<FAB>& cmf, const BoxArray& fba, const IntVect& ratio,
                            Periodicity const& period, int crse_value, int fine_value);
    iMultiFab makeFineMask (const BoxArray& cba, const DistributionMapping& cdm,
                            const IntVect& cnghost, const BoxArray& fba, const IntVect& ratio,
                            Periodicity const& period, int crse_value, int fine_value);
    template <typename FAB>
    iMultiFab makeFineMask (const FabArray<FAB>& cmf, const FabArray<FAB>& fmf,
                            const IntVect& cnghost, const IntVect& ratio,
                            Periodicity const& period, int crse_value, int fine_value);
    template <typename FAB>
    iMultiFab makeFineMask (const FabArray<FAB>& cmf, const FabArray<FAB>& fmf,
                            const IntVect& cnghost, const IntVect& ratio,
                            Periodicity const& period, int crse_value, int fine_value,
                            LayoutData<int>& has_cf);

    MultiFab makeFineMask (const BoxArray& cba, const DistributionMapping& cdm,
                           const BoxArray& fba, const IntVect& ratio,
                           Real crse_value, Real fine_value);

    //! Computes divergence of face-data stored in the umac MultiFab.
    void computeDivergence (MultiFab& divu, const Array<MultiFab const*,AMREX_SPACEDIM>& umac,
                            const Geometry& geom);

    //! Computes gradient of face-data stored in the umac MultiFab.
    void computeGradient (MultiFab& grad, const Array<MultiFab const*,AMREX_SPACEDIM>& umac,
                          const Geometry& geom);

    //! Convert iMultiFab to MultiFab
    MultiFab ToMultiFab (const iMultiFab& imf);
    //! Convert iMultiFab to Long
    FabArray<BaseFab<Long> > ToLongMultiFab (const iMultiFab& imf);

    //! Periodic shift MultiFab
    MultiFab periodicShift (MultiFab const& mf, IntVect const& offset,
                            Periodicity const& period);

    //! example: auto mf = amrex::cast<MultiFab>(imf);
    template <typename T, typename U>
    T cast (U const& mf_in)
    {
        T mf_out(mf_in.boxArray(), mf_in.DistributionMap(), mf_in.nComp(), mf_in.nGrowVect());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf_in); mfi.isValid(); ++mfi)
        {
            const Long n = mfi.fabbox().numPts() * mf_in.nComp();
            auto pdst = mf_out[mfi].dataPtr();
            auto psrc = mf_in [mfi].dataPtr();
            AMREX_HOST_DEVICE_PARALLEL_FOR_1D ( n, i,
            {
                pdst[i] = static_cast<typename U::value_type>(psrc[i]);
            });
        }
        return mf_out;
    }
}

namespace amrex {

template <typename FAB>
iMultiFab
makeFineMask (const FabArray<FAB>& cmf, const BoxArray& fba, const IntVect& ratio,
              int crse_value, int fine_value)
{
    return makeFineMask(cmf.boxArray(), cmf.DistributionMap(), cmf.nGrowVect(),
                        fba, ratio, Periodicity::NonPeriodic(), crse_value, fine_value);
}

template <typename FAB>
iMultiFab
makeFineMask (const FabArray<FAB>& cmf, const BoxArray& fba, const IntVect& ratio,
              Periodicity const& period, int crse_value, int fine_value)
{
    return makeFineMask(cmf.boxArray(), cmf.DistributionMap(), cmf.nGrowVect(),
                        fba, ratio, period, crse_value, fine_value);
}

template <typename FAB>
iMultiFab
makeFineMask (const FabArray<FAB>& cmf, const FabArray<FAB>& fmf,
              const IntVect& cnghost, const IntVect& ratio,
              Periodicity const& period, int crse_value, int fine_value)
{
    iMultiFab mask(cmf.boxArray(), cmf.DistributionMap(), 1, cnghost);
    mask.setVal(crse_value);

    iMultiFab foo(amrex::coarsen(fmf.boxArray(),ratio), fmf.DistributionMap(),
                  1, 0, MFInfo().SetAlloc(false));
    const FabArrayBase::CPC& cpc = mask.getCPC(cnghost,foo,IntVect::TheZeroVector(),period);
    mask.setVal(fine_value, cpc, 0, 1);

    return mask;
}

template <typename FAB>
iMultiFab
makeFineMask (const FabArray<FAB>& cmf, const FabArray<FAB>& fmf,
              const IntVect& cnghost, const IntVect& ratio,
              Periodicity const& period, int crse_value, int fine_value,
              LayoutData<int>& has_cf)
{
    iMultiFab mask(cmf.boxArray(), cmf.DistributionMap(), 1, cnghost);
    mask.setVal(crse_value);

    iMultiFab foo(amrex::coarsen(fmf.boxArray(),ratio), fmf.DistributionMap(),
                  1, 0, MFInfo().SetAlloc(false));
    const FabArrayBase::CPC& cpc = mask.getCPC(cnghost,foo,IntVect::TheZeroVector(),period);
    mask.setVal(fine_value, cpc, 0, 1);

    has_cf = mask.RecvLayoutMask(cpc);

    return mask;
}

//! Average fine node-based MultiFab onto crse node-based MultiFab.
//! This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
template <typename FAB>
void average_down_nodal (const FabArray<FAB>& fine, FabArray<FAB>& crse,
                         const IntVect& ratio, int ngcrse)
{
    AMREX_ASSERT(fine.is_nodal());
    AMREX_ASSERT(crse.is_nodal());
    AMREX_ASSERT(crse.nComp() == fine.nComp());

    int ncomp = crse.nComp();
    using value_type = typename FAB::value_type;

    if (isMFIterSafe(fine, crse))
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(ngcrse);
            Array4<value_type> const& crsearr = crse.array(mfi);
            Array4<value_type const> const& finearr = fine.const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                amrex_avgdown_nodes(tbx,crsearr,finearr,0,0,ncomp,ratio);
            });
        }
    }
    else
    {
        FabArray<FAB> ctmp(amrex::coarsen(fine.boxArray(),ratio), fine.DistributionMap(),
                           ncomp, ngcrse);
        average_down_nodal(fine, ctmp, ratio, ngcrse);
        crse.ParallelCopy(ctmp,0,0,ncomp,ngcrse,ngcrse);
    }
}

}

#endif
