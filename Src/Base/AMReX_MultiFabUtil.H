#ifndef AMREX_MultiFabUtil_H_
#define AMREX_MultiFabUtil_H_
#include <AMReX_Config.H>

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

    template <typename CMF, typename FMF,
              std::enable_if_t<IsFabArray_v<CMF> && IsFabArray_v<FMF>, int> = 0>
    void average_face_to_cellcenter (CMF& cc, int dcomp,
                                     const Array<const FMF*,AMREX_SPACEDIM>& fc,
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
                                     const Geometry& geom,
                                     int ncomp = 1,
                                     bool use_harmonic_averaging = false);
    void average_cellcenter_to_face (const Array<MultiFab*,AMREX_SPACEDIM>& fc,
                                     const MultiFab& cc,
                                     const Geometry& geom,
                                     int ncomp = 1,
                                     bool use_harmonic_averaging = false);

    //! Average fine face-based MultiFab onto crse face-based MultiFab.
    template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int> = 0>
    void average_down_faces (const Vector<const MF*>& fine,
                             const Vector<MF*>& crse,
                             const IntVect& ratio,
                             int ngcrse = 0);
    template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int> = 0>
    void average_down_faces (const Vector<const MF*>& fine,
                             const Vector<MF*>& crse,
                             int ratio,
                             int ngcrse = 0);
    template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int> = 0>
    void average_down_faces (const Array<const MF*,AMREX_SPACEDIM>& fine,
                             const Array<MF*,AMREX_SPACEDIM>& crse,
                             const IntVect& ratio,
                             int ngcrse = 0);
    template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int> = 0>
    void average_down_faces (const Array<const MF*,AMREX_SPACEDIM>& fine,
                             const Array<MF*,AMREX_SPACEDIM>& crse,
                             int ratio,
                             int ngcrse = 0);
    //! This version does average down for one direction.
    //! It uses the IndexType of MultiFabs to determine the direction.
    //! It is expected that one direction is nodal and the rest are cell-centered.
    template <typename FAB>
    void average_down_faces (const FabArray<FAB>& fine, FabArray<FAB>& crse,
                             const IntVect& ratio, int ngcrse=0);

    //  This version takes periodicity into account.
    template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int> = 0>
    void average_down_faces (const Array<const MF*,AMREX_SPACEDIM>& fine,
                             const Array<MF*,AMREX_SPACEDIM>& crse,
                             const IntVect& ratio, const Geometry& crse_geom);
    //  This version takes periodicity into account.
    template <typename FAB>
    void average_down_faces (const FabArray<FAB>& fine, FabArray<FAB>& crse,
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
                             int ngcrse = 0,
                             bool mfiter_is_definitely_safe=false);

    //! Average fine cell-based MultiFab onto crse cell-centered MultiFab using
    //! volume-weighting. This routine DOES NOT assume that the crse BoxArray is
    //! a coarsened version of the fine BoxArray.
    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
                       const Geometry& fgeom, const Geometry& cgeom,
                       int scomp, int ncomp, const IntVect& ratio);
    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
                       const Geometry& fgeom, const Geometry& cgeom,
                       int scomp, int ncomp, int rr);

    //! Average MultiFab onto crse MultiFab without volume weighting. This
    //! routine DOES NOT assume that the crse BoxArray is a coarsened version of
    //! the fine BoxArray. Work for both cell-centered and nodal MultiFabs.
    template<typename FAB>
    void average_down (const FabArray<FAB>& S_fine, FabArray<FAB>& S_crse,
                       int scomp, int ncomp, const IntVect& ratio);
    template<typename FAB>
    void average_down (const FabArray<FAB>& S_fine, FabArray<FAB>& S_crse,
                       int scomp, int ncomp, int rr);

    //! Add a coarsened version of the data contained in the S_fine MultiFab to
    //! S_crse, including ghost cells.
    void sum_fine_to_coarse (const MultiFab& S_Fine, MultiFab& S_crse,
                             int scomp, int ncomp,
                             const IntVect& ratio,
                             const Geometry& cgeom, const Geometry& fgeom);

    //! Output state data for a single zone
    void print_state (const MultiFab& mf, const IntVect& cell, int n=-1,
                      const IntVect& ng = IntVect::TheZeroVector());

    //! Write each fab individually
    void writeFabs (const MultiFab& mf, const std::string& name);
    void writeFabs (const MultiFab& mf, int comp, int ncomp, const std::string& name);

    //! Extract a slice from the given cell-centered MultiFab at coordinate
    //! "coord" along direction "dir".
    std::unique_ptr<MultiFab> get_slice_data(int dir, Real coord,
                                             const MultiFab& cc,
                                             const Geometry& geom, int start_comp, int ncomp,
                                             bool interpolate=false);

    /**
     * \brief Get data in a cell of MultiFab/FabArray
     *
     * This returns a Vector containing the data in a given cell, if it's
     * available on a process. The returned Vector is empty if a process
     * does not have the given cell.
     */
    template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int> FOO = 0>
    Vector<typename MF::value_type> get_cell_data (MF const& mf, IntVect const& cell);

    /**
     * \brief Get data in a line of MultiFab/FabArray
     *
     * This returns a MultiFab/FabArray containing the data in a line
     * specified by a direction and a cell.
     */
    template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int> FOO = 0>
    MF get_line_data (MF const& mf, int dir, IntVect const& cell);

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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf_in); mfi.isValid(); ++mfi)
        {
            const Long n = mfi.fabbox().numPts() * mf_in.nComp();
            auto pdst = mf_out[mfi].dataPtr();
            auto psrc = mf_in [mfi].dataPtr();
            AMREX_HOST_DEVICE_PARALLEL_FOR_1D ( n, i,
            {
                pdst[i] = static_cast<typename U::value_type>(psrc[i]); // NOLINT(bugprone-signed-char-misuse)
            });
        }
        return mf_out;
    }

    /**
     * \brief Reduce FabArray/MultiFab data to a plane.
     *
     * This function takes a FabArray/MultiFab and reduces its data to a
     * plane.  The return data are stored in a BaseFab with only one cell in
     * the normal direction of the plane.  The index range of the BaseFab in
     * the other directions is the same as the provided domain Box.  If data
     * do not exist along a certain line, the value is set to the minimum,
     * maximum and zero, for reduce max, min and sum, respectively.  The
     * reduction is local and the user may need to do MPI communication
     * afterwards if needed.
     *
     * In the example code below, the sum along each line at (i,j) in the
     * z-direction is computed and stored at (i,j,0) of the returned
     * BaseFab.
     \verbatim
         int dir = 2; // z-direction
         auto const& domain_box = geom.Domain();
         auto const& ma = mf.const_arrays();
         auto rr = ReduceToPlane<ReduceOpSum,Real>(dir, domain_box, mf,
             [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> Real
             {
                 return ma[box_no](i,j,k); // data at (i,j,k) of Box box_no
             });
     \endverbatim
     *
     * Below is another example.  This finds the maximum value in the
     * x-direction and stores the maximum value and the i-index.  An MPI
     * reduce is then called to further reduce the data to the root process
     * 0.
     \verbatim
         int dir = 0; // x-direction
         auto const& domain_box = geom.Domain().surroundingNodes(); // nodal data
         auto const& ma = mf.const_arrays();
         auto rr = ReduceToPlane<ReduceOpMax,KeyValuePair<Real,int>>
             (dir, domain_box, mf,
              [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
                  -> KeyValuePair<Real,int>
             {
                 return {ma[box_no](i,j,k), i};
             });
         ParallelReduce::Max(rr.dataPtr(), rr.size(), root,
                             ParallelDescriptor::Communicator());
         // Process root now has the final results.
     \endverbatim
     *
     * \tparam Op  reduce operator (e.g., ReduceOpSum, ReduceOpMin and ReduceOpMax)
     * \tparam T   data type of reduction result
     * \tparam FAB FabArray/MultiFab type
     * \tparam F   callable type like a lambda function
     *
     * \param direction normal direction of the plane (e.g., 0, 1 and 2)
     * \param domain    domain Box
     * \param mf        a FabArray/MultiFab object specifying the iteration space
     * \param f         a callable object returning T.  It takes four ints,
     *                  where the first int is the local box index and the others
     *                  are spatial indices for x, y, and z-directions.
     *
     * \return reduction result (BaseFab<T>)
     */
    template <typename Op, typename T, typename FAB, typename F,
              std::enable_if_t<IsBaseFab<FAB>::value
#ifndef AMREX_USE_CUDA
                               && IsCallableR<T,F,int,int,int,int>::value
#endif
                               , int> FOO = 0>
    BaseFab<T>
    ReduceToPlane (int direction, Box const& domain, FabArray<FAB> const& mf, F&& f);

    /**
     * \brief Sum MultiFab data to line
     *
     * Return a HostVector that contains the sum of the given MultiFab data in the plane
     * with the given normal direction.  The size of the vector is
     * domain.length(direction) x ncomp.  The vector is actually a 2D array, where the
     * element for component icomp at spatial index k is at [icomp*ncomp+k].
     *
     * \param mf MultiFab data for summing
     * \param icomp starting component
     * \param ncomp number of components
     * \param domain the domain
     * \param direction the direction of the line
     * \param local If false, reduce across MPI processes.
     */
    Gpu::HostVector<Real> sumToLine (MultiFab const& mf, int icomp, int ncomp,
                                     Box const& domain, int direction, bool local = false);

    /** \brief Volume weighted sum for a vector of MultiFabs
     *
     * Return a volume weighted sum of MultiFabs of AMR data.  The sum is
     * perform on a single component of the data.  If the MultiFabs are
     * built with EB Factories, the cut cell volume fraction will be
     * included in the weight.
     */
    Real volumeWeightedSum (Vector<MultiFab const*> const& mf, int icomp,
                            Vector<Geometry> const& geom,
                            Vector<IntVect> const& ratio,
                            bool local = false);

    /**
     * \brief Fourth-order interpolation from fine to coarse level.
     *
     * This is for high-order "average-down" of finite-difference data.  If
     * ghost cell data are used, it's the caller's responsibility to fill
     * the ghost cells before calling this function.
     *
     * \param cmf   coarse data
     * \param scomp starting component
     * \param ncomp number of component
     * \param fmf   fine data
     * \param ratio refinement ratio.
     */
    void FourthOrderInterpFromFineToCoarse (MultiFab& cmf, int scomp, int ncomp,
                                            MultiFab const& fmf,
                                            IntVect const& ratio);
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
                         const IntVect& ratio, int ngcrse, bool mfiter_is_definitely_safe)
{
    AMREX_ASSERT(fine.is_nodal());
    AMREX_ASSERT(crse.is_nodal());
    AMREX_ASSERT(crse.nComp() == fine.nComp());

    int ncomp = crse.nComp();
    using value_type = typename FAB::value_type;

    if (mfiter_is_definitely_safe || isMFIterSafe(fine, crse))
    {
#ifdef AMREX_USE_OMP
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

// *************************************************************************************************************

// Average fine cell-based MultiFab onto crse cell-centered MultiFab.
// We do NOT assume that the coarse layout is a coarsened version of the fine layout.
// This version does NOT use volume-weighting
template<typename FAB>
void average_down (const FabArray<FAB>& S_fine, FabArray<FAB>& S_crse, int scomp, int ncomp, int rr)
{
    average_down(S_fine,S_crse,scomp,ncomp,rr*IntVect::TheUnitVector());
}

template<typename FAB>
void average_down (const FabArray<FAB>& S_fine, FabArray<FAB>& S_crse,
                   int scomp, int ncomp, const IntVect& ratio)
{
    BL_PROFILE("amrex::average_down");
    AMREX_ASSERT(S_crse.nComp() == S_fine.nComp());
    AMREX_ASSERT((S_crse.is_cell_centered() && S_fine.is_cell_centered()) ||
                 (S_crse.is_nodal()         && S_fine.is_nodal()));

    using value_type = typename FAB::value_type;

    bool is_cell_centered = S_crse.is_cell_centered();

    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

    if (crse_S_fine_BA == S_crse.boxArray() && S_fine.DistributionMap() == S_crse.DistributionMap())
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && S_crse.isFusingCandidate()) {
            auto const& crsema = S_crse.arrays();
            auto const& finema = S_fine.const_arrays();
            if (is_cell_centered) {
                ParallelFor(S_crse, IntVect(0), ncomp,
                            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
                            {
                                amrex_avgdown(i,j,k,n,crsema[box_no],finema[box_no],scomp,scomp,ratio);
                            });
            } else {
                ParallelFor(S_crse, IntVect(0), ncomp,
                            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
                            {
                                amrex_avgdown_nodes(i,j,k,n,crsema[box_no],finema[box_no],scomp,scomp,ratio);
                            });
            }
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(S_crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& bx = mfi.tilebox();
                Array4<value_type> const& crsearr = S_crse.array(mfi);
                Array4<value_type const> const& finearr = S_fine.const_array(mfi);

                if (is_cell_centered) {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                                                      {
                                                          amrex_avgdown(i,j,k,n,crsearr,finearr,scomp,scomp,ratio);
                                                      });
                } else {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                                                      {
                                                          amrex_avgdown_nodes(i,j,k,n,crsearr,finearr,scomp,scomp,ratio);
                                                      });
                }
            }
        }
    }
    else
    {
        FabArray<FAB> crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, 0, MFInfo(),DefaultFabFactory<FAB>());

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && crse_S_fine.isFusingCandidate()) {
            auto const& crsema = crse_S_fine.arrays();
            auto const& finema = S_fine.const_arrays();
            if (is_cell_centered) {
                ParallelFor(crse_S_fine, IntVect(0), ncomp,
                            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
                            {
                                amrex_avgdown(i,j,k,n,crsema[box_no],finema[box_no],0,scomp,ratio);
                            });
            } else {
                ParallelFor(crse_S_fine, IntVect(0), ncomp,
                            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
                            {
                                amrex_avgdown_nodes(i,j,k,n,crsema[box_no],finema[box_no],0,scomp,ratio);
                            });
            }
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(crse_S_fine,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& bx = mfi.tilebox();
                Array4<value_type> const& crsearr = crse_S_fine.array(mfi);
                Array4<value_type const> const& finearr = S_fine.const_array(mfi);

                //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
                //        because the crse fab is a temporary which was made starting at comp 0, it is
                //        not part of the actual crse multifab which came in.

                if (is_cell_centered) {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                                                      {
                                                          amrex_avgdown(i,j,k,n,crsearr,finearr,0,scomp,ratio);
                                                      });
                } else {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                                                      {
                                                          amrex_avgdown_nodes(i,j,k,n,crsearr,finearr,0,scomp,ratio);
                                                      });
                }
            }
        }

        S_crse.ParallelCopy(crse_S_fine,0,scomp,ncomp);
    }
}





   /**
    * \brief Returns part of a norm based on two MultiFabs
    * The MultiFabs MUST have the same underlying BoxArray.
    * The function f is applied elementwise as f(x(i,j,k,n),y(i,j,k,n))
    * inside the summation (subject to a valid mask entry pf(mask(i,j,k,n)
    */

template <typename F>
Real
NormHelper (const MultiFab& x, int xcomp,
            const MultiFab& y, int ycomp,
            F && f,
            int numcomp, IntVect nghost, bool local)
{
    BL_ASSERT(x.boxArray() == y.boxArray());
    BL_ASSERT(x.DistributionMap() == y.DistributionMap());
    BL_ASSERT(x.nGrowVect().allGE(nghost) && y.nGrowVect().allGE(nghost));

    Real sm = Real(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& xma = x.const_arrays();
        auto const& yma = y.const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, x, nghost,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            auto const& xfab = xma[box_no];
            auto const& yfab = yma[box_no];
            for (int n = 0; n < numcomp; ++n) {
                t += f(xfab(i,j,k,xcomp+n) , yfab(i,j,k,ycomp+n));
            }
            return t;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            Array4<Real const> const& xfab = x.const_array(mfi);
            Array4<Real const> const& yfab = y.const_array(mfi);
            AMREX_LOOP_4D(bx, numcomp, i, j, k, n,
            {
                sm += f(xfab(i,j,k,xcomp+n) , yfab(i,j,k,ycomp+n));
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}

   /**
    * \brief Returns part of a norm based on three MultiFabs
    * The MultiFabs MUST have the same underlying BoxArray.
    * The Predicate pf is used to test the mask
    * The function f is applied elementwise as f(x(i,j,k,n),y(i,j,k,n))
    * inside the summation (subject to a valid mask entry pf(mask(i,j,k,n)
    */

template <typename MMF, typename Pred, typename F>
Real
NormHelper (const MMF& mask,
               const MultiFab& x, int xcomp,
               const MultiFab& y, int ycomp,
               Pred && pf,
               F && f,
               int numcomp, IntVect nghost, bool local)
{
    BL_ASSERT(x.boxArray() == y.boxArray());
    BL_ASSERT(x.boxArray() == mask.boxArray());
    BL_ASSERT(x.DistributionMap() == y.DistributionMap());
    BL_ASSERT(x.DistributionMap() == mask.DistributionMap());
    BL_ASSERT(x.nGrowVect().allGE(nghost) && y.nGrowVect().allGE(nghost));
    BL_ASSERT(mask.nGrowVect().allGE(nghost));

    Real sm = Real(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& xma = x.const_arrays();
        auto const& yma = y.const_arrays();
        auto const& mma = mask.const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, x, nghost,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            if (pf(mma[box_no](i,j,k))) {
                auto const& xfab = xma[box_no];
                auto const& yfab = yma[box_no];
                for (int n = 0; n < numcomp; ++n) {
                    t += f(xfab(i,j,k,xcomp+n) , yfab(i,j,k,ycomp+n));
                }
            }
            return t;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            Array4<Real const> const& xfab = x.const_array(mfi);
            Array4<Real const> const& yfab = y.const_array(mfi);
            auto const& mfab = mask.const_array(mfi);
            AMREX_LOOP_4D(bx, numcomp, i, j, k, n,
            {
                if (pf(mfab(i,j,k))) {
                    sm += f(xfab(i,j,k,xcomp+n) , yfab(i,j,k,ycomp+n));
                }
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}

template <typename CMF, typename FMF,
          std::enable_if_t<IsFabArray_v<CMF> && IsFabArray_v<FMF>, int> FOO>
void average_face_to_cellcenter (CMF& cc, int dcomp,
                                 const Array<const FMF*,AMREX_SPACEDIM>& fc,
                                 int ngrow)
{
    AMREX_ASSERT(cc.nComp() >= dcomp + AMREX_SPACEDIM);
    AMREX_ASSERT(fc[0]->nComp() == 1);

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && cc.isFusingCandidate()) {
        auto const& ccma = cc.arrays();
        AMREX_D_TERM(auto const& fxma = fc[0]->const_arrays();,
                     auto const& fyma = fc[1]->const_arrays();,
                     auto const& fzma = fc[2]->const_arrays(););
        ParallelFor(cc, IntVect(ngrow),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
#if (AMREX_SPACEDIM == 1)
            GeometryData gd{};
            gd.coord = 0;
#endif
            amrex_avg_fc_to_cc(i,j,k, ccma[box_no], AMREX_D_DECL(fxma[box_no],
                                                                 fyma[box_no],
                                                                 fzma[box_no]),
                               dcomp
#if (AMREX_SPACEDIM == 1)
                               , gd
#endif
                );
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.growntilebox(ngrow);
            auto const& ccarr = cc.array(mfi);
            AMREX_D_TERM(auto const& fxarr = fc[0]->const_array(mfi);,
                         auto const& fyarr = fc[1]->const_array(mfi);,
                         auto const& fzarr = fc[2]->const_array(mfi););

#if (AMREX_SPACEDIM == 1)
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D( bx, i, j, k,
            {
                GeometryData gd;
                gd.coord = 0;
                amrex_avg_fc_to_cc(i,j,k, ccarr, fxarr, dcomp, gd);
            });
#else
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D( bx, i, j, k,
            {
                amrex_avg_fc_to_cc(i,j,k, ccarr, AMREX_D_DECL(fxarr,fyarr,fzarr), dcomp);
            });
#endif
        }
    }
}

template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int>>
void average_down_faces (const Vector<const MF*>& fine,
                         const Vector<MF*>& crse,
                         const IntVect& ratio, int ngcrse)
{
    AMREX_ASSERT(fine.size() == AMREX_SPACEDIM && crse.size() == AMREX_SPACEDIM);
    average_down_faces(Array<const MF*,AMREX_SPACEDIM>
                               {{AMREX_D_DECL(fine[0],fine[1],fine[2])}},
                       Array<MF*,AMREX_SPACEDIM>
                               {{AMREX_D_DECL(crse[0],crse[1],crse[2])}},
                       ratio, ngcrse);
}

template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int>>
void average_down_faces (const Vector<const MF*>& fine,
                         const Vector<MF*>& crse, int ratio, int ngcrse)
{
    average_down_faces(fine,crse,IntVect{ratio},ngcrse);
}

template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int>>
void average_down_faces (const Array<const MF*,AMREX_SPACEDIM>& fine,
                         const Array<MF*,AMREX_SPACEDIM>& crse,
                         int ratio, int ngcrse)
{
    average_down_faces(fine,crse,IntVect{ratio},ngcrse);
}

template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int>>
void average_down_faces (const Array<const MF*,AMREX_SPACEDIM>& fine,
                         const Array<MF*,AMREX_SPACEDIM>& crse,
                         const IntVect& ratio, int ngcrse)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        average_down_faces(*fine[idim], *crse[idim], ratio, ngcrse);
    }
}

template <typename FAB>
void average_down_faces (const FabArray<FAB>& fine, FabArray<FAB>& crse,
                         const IntVect& ratio, int ngcrse)
{
    BL_PROFILE("average_down_faces");

    AMREX_ASSERT(crse.nComp() == fine.nComp());
    AMREX_ASSERT(fine.ixType() == crse.ixType());
    const auto type = fine.ixType();
    int dir;
    for (dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (type.nodeCentered(dir)) break;
    }
    auto tmptype = type;
    tmptype.unset(dir);
    if (dir >= AMREX_SPACEDIM || !tmptype.cellCentered()) {
        amrex::Abort("average_down_faces: not face index type");
    }
    const int ncomp = crse.nComp();
    if (isMFIterSafe(fine, crse))
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && crse.isFusingCandidate()) {
            auto const& crsema = crse.arrays();
            auto const& finema = fine.const_arrays();
            ParallelFor(crse, IntVect(ngcrse), ncomp,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                amrex_avgdown_faces(i,j,k,n, crsema[box_no], finema[box_no], 0, 0, ratio, dir);
            });
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(ngcrse);
                auto const& crsearr = crse.array(mfi);
                auto const& finearr = fine.const_array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                {
                    amrex_avgdown_faces(i,j,k,n, crsearr, finearr, 0, 0, ratio, dir);
                });
            }
        }
    }
    else
    {
        FabArray<FAB> ctmp(amrex::coarsen(fine.boxArray(),ratio), fine.DistributionMap(),
                      ncomp, ngcrse, MFInfo(), DefaultFabFactory<FAB>());
        average_down_faces(fine, ctmp, ratio, ngcrse);
        crse.ParallelCopy(ctmp,0,0,ncomp,ngcrse,ngcrse);
    }
}

template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int>>
void average_down_faces (const Array<const MF*,AMREX_SPACEDIM>& fine,
                         const Array<MF*,AMREX_SPACEDIM>& crse,
                         const IntVect& ratio, const Geometry& crse_geom)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        average_down_faces(*fine[idim], *crse[idim], ratio, crse_geom);
    }
}

template <typename FAB>
void average_down_faces (const FabArray<FAB>& fine, FabArray<FAB>& crse,
                         const IntVect& ratio, const Geometry& crse_geom)
{
    FabArray<FAB> ctmp(amrex::coarsen(fine.boxArray(),ratio), fine.DistributionMap(),
                  crse.nComp(), 0);
    average_down_faces(fine, ctmp, ratio, 0);
    crse.ParallelCopy(ctmp,0,0,crse.nComp(),0,0,crse_geom.periodicity());
}

template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int> FOO>
Vector<typename MF::value_type> get_cell_data (MF const& mf, IntVect const& cell)
{
    using T = typename MF::value_type;
    const int ncomp = mf.nComp();
    Gpu::DeviceVector<T> dv(ncomp);
    auto dp = dv.data();
    bool found = false;
    auto loc = cell.dim3();
    for (MFIter mfi(mf); mfi.isValid() && !found; ++mfi)
    {
        Box const& box = mfi.validbox();
        if (box.contains(cell)) {
            found = true;
            auto const& fab = mf.const_array(mfi);
            amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept
            {
                for (int n = 0; n < ncomp; ++n) {
                    dp[n] = fab(loc.x,loc.y,loc.z,n);
                }
            });
        }
    }
    Vector<T> hv;
    if (found) {
        hv.resize(ncomp);
        Gpu::copy(Gpu::deviceToHost, dv.begin(), dv.end(), hv.begin());
    }
    return hv;
}

template <typename MF, std::enable_if_t<IsFabArray<MF>::value,int> FOO>
MF get_line_data (MF const& mf, int dir, IntVect const& cell)
{
    BoxArray const& ba = mf.boxArray();
    DistributionMapping const& dm = mf.DistributionMap();
    const auto nboxes = static_cast<int>(ba.size());

    BoxList bl(ba.ixType());
    Vector<int> procmap;
    Vector<int> index_map;
    for (int i = 0; i < nboxes; ++i) {
        Box const& b = ba[i];
        IntVect lo = cell;
        lo[dir] = b.smallEnd(dir);
        if (b.contains(lo)) {
            IntVect hi = lo;
            hi[dir] = b.bigEnd(dir);
            Box b1d(lo,hi,b.ixType());
            bl.push_back(b1d);
            procmap.push_back(dm[i]);
            index_map.push_back(i);
        }
    }

    if (bl.isEmpty()) {
        return MF();
    } else {
        BoxArray rba(std::move(bl));
        DistributionMapping rdm(std::move(procmap));
        MF rmf(rba, rdm, mf.nComp(), IntVect(0),
               MFInfo().SetArena(mf.arena()));
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rmf); mfi.isValid(); ++mfi) {
            Box const& b = mfi.validbox();
            auto const& dfab = rmf.array(mfi);
            auto const& sfab = mf.const_array(index_map[mfi.index()]);
            amrex::ParallelFor(b, mf.nComp(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dfab(i,j,k,n) = sfab(i,j,k,n);
            });
        }
        return rmf;
    }
}

template <typename Op, typename T, typename FAB, typename F,
          std::enable_if_t<IsBaseFab<FAB>::value
#ifndef AMREX_USE_CUDA
                           && IsCallableR<T,F,int,int,int,int>::value
#endif
                           , int> FOO>
BaseFab<T>
ReduceToPlane (int direction, Box const& domain, FabArray<FAB> const& mf, F&& f)
{
    Box domain2d = domain;
    domain2d.setRange(direction, 0);

    T initval;
    Op().init(initval);

    BaseFab<T> r(domain2d);
    r.template setVal<RunOn::Device>(initval);
    auto const& ar = r.array();

    for (MFIter mfi(mf,MFItInfo().UseDefaultStream().DisableDeviceSync());
         mfi.isValid(); ++mfi)
    {
        Box bx = mfi.validbox() & domain;
        if (bx.ok()) {
            int box_no = mfi.LocalIndex();
#if defined(AMREX_USE_GPU)
            Box b2d = bx;
            b2d.setRange(direction,0);
            const auto blo = amrex::lbound(bx);
            const auto len = amrex::length(bx);
            constexpr int nthreads = 128;
            auto nblocks = static_cast<int>(b2d.numPts());
#ifdef AMREX_USE_SYCL
            constexpr std::size_t shared_mem_bytes = sizeof(T)*Gpu::Device::warp_size;
            amrex::launch<nthreads>(nblocks, shared_mem_bytes, Gpu::gpuStream(),
                                    [=] AMREX_GPU_DEVICE (Gpu::Handler const& h)
            {
                int bid = h.blockIdx();
                int tid = h.threadIdx();
#else
            amrex::launch<nthreads>(nblocks, Gpu::gpuStream(),
                                    [=] AMREX_GPU_DEVICE ()
            {
                int bid = blockIdx.x;
                int tid = threadIdx.x;
#endif
                T tmp;
                Op().init(tmp);
                T* p;
                if (direction == 0) {
                    int k = bid /   len.y;
                    int j = bid - k*len.y;
                    k += blo.z;
                    j += blo.y;
                    for (int i = blo.x + tid; i < blo.x+len.x; i += nthreads) {
                        Op().local_update(tmp, f(box_no,i,j,k));
                    }
                    p = ar.ptr(0,j,k);
                } else if (direction == 1) {
                    int k = bid /   len.x;
                    int i = bid - k*len.x;
                    k += blo.z;
                    i += blo.x;
                    for (int j = blo.y + tid; j < blo.y+len.y; j += nthreads) {
                        Op().local_update(tmp, f(box_no,i,j,k));
                    }
                    p = ar.ptr(i,0,k);
                } else {
                    int j = bid /   len.x;
                    int i = bid - j*len.x;
                    j += blo.y;
                    i += blo.x;
                    for (int k = blo.z + tid; k < blo.z+len.z; k += nthreads) {
                        Op().local_update(tmp, f(box_no,i,j,k));
                    }
                    p = ar.ptr(i,j,0);
                }
#ifdef AMREX_USE_SYCL
                Op().template parallel_update<T>(*p, tmp, h);
#else
                Op().template parallel_update<T,nthreads>(*p, tmp);
#endif
            });
#else
            // CPU
            if (direction == 0) {
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    Op().local_update(ar(0,j,k), f(box_no,i,j,k));
                });
            } else if (direction == 1) {
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    Op().local_update(ar(i,0,k), f(box_no,i,j,k));
                });
            } else {
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    Op().local_update(ar(i,j,0), f(box_no,i,j,k));
                });
            }
#endif
        }
    }
    Gpu::streamSynchronize();

    return r;
}

}

#endif
