#ifndef AMREX_FILLPATCHER_H_
#define AMREX_FILLPATCHER_H_
#include <AMReX_Config.H>

#include <AMReX_FillPatchUtil.H>

namespace amrex {

/**
 * \brief FillPatcher is for filling a fine level MultiFab/FabArray.
 *
 * This class is not as general as the FillPatchTwoLevels functions.  It
 * fills the fine ghost cells not overlapping any fine level valid cells
 * with interpolation of the coarse data.  Then it fills the fine ghost
 * cells overlapping fine level valid cells with the fine level data.  If
 * the valid cells of the destination need to be filled, it will be done as
 * well.  Finally, it will fill the physical bounbary using the user
 * provided functor.  The `fill` member function can be used to do the
 * operations just described.  Alternatively, one can also use the
 * `fillCoarseFineBounary` to fill the ghost cells at the coarse/fine
 * boundary only.  Then one can manually call FillBoundary to fill the other
 * ghost cells, and use the physical BC functor to handle the physical
 * boundeary.
 *
 * The communication of the coarse data needed for spatial interpolation is
 * optimized at the cost of being error-prone.  One must follow the
 * following guidelines.
 *
 * (1) This class is for filling data during time stepping, not during
 * regrid.  The fine level data passed as input must have the same BoxArray
 * and DistributionMapping as the destination.  It's OK they are the same
 * MultiFab.  For AmrLevel based codes, AmrLevel::FillPatcherFill wil try to
 * use FillPatcher if it can, and AmrLevel::FillPatch will use the fillpatch
 * functions.
 *
 * (2) When to build?  It is recommended that one uses `std::unique_ptr` to
 * store the FillPatcher object, and build it only when it is needed and
 * it's a nullptr.  For AmrLevel based codes, the AmrLevel class will build
 * it for you as needed when you call the AmrLevel::FillPatcherFill
 * function.
 *
 * (3) When to destroy?  Usually, we do time steppig on a coarse level
 * first.  Then we recursively do time stepping on fine levels.  After the
 * finer level finishes, we do reflux and averge the fine data down to the
 * coarse level.  After that we should destroy the FillPatcher object
 * associated with these two levels, because the coarse data stored in the
 * object has become outdated.  For AmrCore based codes, you could use
 * Tests/Amr/Advection_AmrCore as an example.  For AmrLevel based codes, you
 * should do this in the post_timestep virtual function (see
 * Tests/Amr/Advection_AmrLevel for an example).
 *
 * (4) The source MultiFabs/FabArrays (i.e., the crse_data and fine_data
 * arguments of the fill function) need to have exactly the same number of
 * components as the ncomp argument of the constructor, even though it's
 * allowed to fill only some of the components with the fill function.
 *
 * (5) This only works for cell-centered and nodal data.
 */

template <class MF = MultiFab>
class FillPatcher
{
public:

    /**
     * \brief Constructor of FillPatcher
     *
     * \param fba    fine level BoxArray
     * \param fdm    fine level DistributionMapping
     * \param fgeom  fine level Geometry
     * \param cba    coarse level BoxArray
     * \param cdm    coarse level DistributionMapping
     * \param cgeom  coarse level Geometry
     * \param nghost max number of ghost cells to be filled at coarse/fine boundary
     * \param ncomp  the number of components
     * \param interp for spatial interpolation
     * \param eb_index_space optional argument for specifying EB IndexSpace
     */
    FillPatcher (BoxArray const& fba, DistributionMapping const& fdm,
                 Geometry const& fgeom,
                 BoxArray const& cba, DistributionMapping const& cdm,
                 Geometry const& cgeom,
                 IntVect const& nghost, int ncomp, InterpBase* interp,
#ifdef AMREX_USE_EB
                 EB2::IndexSpace const* eb_index_space = EB2::TopIndexSpaceIfPresent());
#else
                 EB2::IndexSpace const* eb_index_space = nullptr);
#endif

    /**
     * \brief Function to fill data
     *
     * \param mf          destination MultiFab/FabArray
     * \param nghost      number of ghost cells to fill. This must be <= what's
     *                    provided to the constructor
     * \param time        time associated with the destination
     * \param crse_data   coarse level data
     * \param crse_time   time associated with the coarse data
     * \param fine_data   fine level data
     * \param fine_time   time associated with the fine data
     * \param scomp       starting component of the source
     * \param dcomp       starting component of the destination
     * \param ncomp       the number of components to fill
     * \param cbc         for filling coarse level physical BC
     * \param cbccomp     starting component of the coarse level BC functor
     * \param fbc         for filling fine level physical BC
     * \param fbccomp     starting component of the fine level BC functor
     * \param bcs         BCRec specifying physical boundary types
     * \parame bcscomp    starting component of the BCRec Vector.
     * \param pre_interp  optional pre-interpolation hook for modifying the coarse data
     * \param post_interp optional post-interpolation hook for modifying the fine data
     */
    template <typename BC,
              typename PreInterpHook=NullInterpHook<MF>,
              typename PostInterpHook=NullInterpHook<MF> >
    void fill (MF& mf, IntVect const& nghost, Real time,
               Vector<MF*> const& crse_data, Vector<Real> const& crse_time,
               Vector<MF*> const& fine_data, Vector<Real> const& fine_time,
               int scomp, int dcomp, int ncomp,
               BC& cbc, int cbccomp, BC& fbc, int fbccomp,
               Vector<BCRec> const& bcs, int bcscomp,
               PreInterpHook const& pre_interp = {},
               PostInterpHook const& post_interp = {});

    /**
     * \brief Function to fill data at coarse/fine boundary only
     *
     * \param mf          destination MultiFab/FabArray
     * \param nghost      number of ghost cells to fill. This must be <= what's
     *                    provided to the constructor
     * \param time        time associated with the destination
     * \param crse_data   coarse level data
     * \param crse_time   time associated with the coarse data
     * \param scomp       starting component of the source
     * \param dcomp       starting component of the destination
     * \param ncomp       the number of components to fill
     * \param cbc         for filling coarse level physical BC
     * \param cbccomp     starting component of the coarse level BC functor
     * \param bcs         BCRec specifying physical boundary types
     * \param bcscomp     starting component of the BCRec Vector.
     * \param pre_interp  optional pre-interpolation hook for modifying the coarse data
     * \param post_interp optional post-interpolation hook for modifying the fine data
     */
    template <typename BC,
              typename PreInterpHook=NullInterpHook<MF>,
              typename PostInterpHook=NullInterpHook<MF> >
    void fillCoarseFineBoundary (MF& mf, IntVect const& nghost, Real time,
                                 Vector<MF*> const& crse_data,
                                 Vector<Real> const& crse_time,
                                 int scomp, int dcomp, int ncomp,
                                 BC& cbc, int cbccomp,
                                 Vector<BCRec> const& bcs, int bcscomp,
                                 PreInterpHook const& pre_interp = {},
                                 PostInterpHook const& post_interp = {});

private:

    BoxArray m_fba;
    BoxArray m_cba;
    DistributionMapping m_fdm;
    DistributionMapping m_cdm;
    Geometry m_fgeom;
    Geometry m_cgeom;
    IntVect m_nghost;
    int m_ncomp;
    InterpBase* m_interp;
    EB2::IndexSpace const* m_eb_index_space = nullptr;
    Vector<std::pair<Real,std::unique_ptr<MF>>> m_cf_crse_data;
    std::unique_ptr<MF> m_cf_fine_data;
};

template <class MF>
FillPatcher<MF>::FillPatcher (BoxArray const& fba, DistributionMapping const& fdm,
                              Geometry const& fgeom,
                              BoxArray const& cba, DistributionMapping const& cdm,
                              Geometry const& cgeom,
                              IntVect const& nghost, int ncomp, InterpBase* interp,
                              EB2::IndexSpace const* eb_index_space)
    : m_fba(fba),
      m_cba(cba),
      m_fdm(fdm),
      m_cdm(cdm),
      m_fgeom(fgeom),
      m_cgeom(cgeom),
      m_nghost(nghost),
      m_ncomp(ncomp),
      m_interp(interp),
      m_eb_index_space(eb_index_space)
{
    static_assert(IsFabArray<MF>::value,
                  "FillPatcher<MF>: MF must be FabArray type");
    AMREX_ALWAYS_ASSERT(m_fba.ixType().cellCentered() || m_fba.ixType().nodeCentered());
}

template <class MF>
template <typename BC, typename PreInterpHook, typename PostInterpHook>
void
FillPatcher<MF>::fill (MF& mf, IntVect const& nghost, Real time,
                       Vector<MF*> const& cmf, Vector<Real> const& ct,
                       Vector<MF*> const& fmf, Vector<Real> const& ft,
                       int scomp, int dcomp, int ncomp,
                       BC& cbc, int cbccomp,
                       BC& fbc, int fbccomp,
                       Vector<BCRec> const& bcs, int bcscomp,
                       PreInterpHook const& pre_interp,
                       PostInterpHook const& post_interp)
{
    BL_PROFILE("FillPatcher::fill()");

    AMREX_ALWAYS_ASSERT(m_fba == fmf[0]->boxArray() &&
                        m_fdm == fmf[0]->DistributionMap());

    fillCoarseFineBoundary(mf, nghost, time, cmf, ct, scomp, dcomp, ncomp,
                           cbc, cbccomp, bcs, bcscomp, pre_interp, post_interp);

    FillPatchSingleLevel(mf, nghost, time, fmf, ft, scomp, dcomp, ncomp,
                         m_fgeom, fbc, fbccomp);
}

template <class MF>
template <typename BC, typename PreInterpHook, typename PostInterpHook>
void
FillPatcher<MF>::fillCoarseFineBoundary (MF& mf, IntVect const& nghost, Real time,
                                         Vector<MF*> const& cmf,
                                         Vector<Real> const& ct,
                                         int scomp, int dcomp, int ncomp,
                                         BC& cbc, int cbccomp,
                                         Vector<BCRec> const& bcs, int bcscomp,
                                         PreInterpHook const& pre_interp,
                                         PostInterpHook const& post_interp)
{
    BL_PROFILE("FillPatcher::fillCFB");

    AMREX_ALWAYS_ASSERT(nghost.allLE(m_nghost) &&
                        m_fba == mf.boxArray() &&
                        m_fdm == mf.DistributionMap() &&
                        m_cba == cmf[0]->boxArray() &&
                        m_cdm == cmf[0]->DistributionMap() &&
                        m_ncomp >= ncomp &&
                        m_ncomp == cmf[0]->nComp());

    IntVect ratio;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        ratio[idim] = m_fgeom.Domain().length(idim) / m_cgeom.Domain().length(idim);
    }
    AMREX_ASSERT(m_fgeom.Domain() == amrex::refine(m_cgeom.Domain(),ratio));

    const InterpolaterBoxCoarsener& coarsener = m_interp->BoxCoarsener(ratio);
    const FabArrayBase::FPinfo& fpc = FabArrayBase::TheFPinfo(mf, mf,
                                                              m_nghost,
                                                              coarsener,
                                                              m_fgeom,
                                                              m_cgeom,
                                                              m_eb_index_space);

    if ( ! fpc.ba_crse_patch.empty())
    {
        if (m_cf_fine_data == nullptr) {
            m_cf_fine_data = std::make_unique<MF>
                (make_mf_fine_patch<MF>(fpc, m_ncomp));
        }

        int ncmfs = cmf.size();
        for (int icmf = 0; icmf < ncmfs; ++icmf) {
            Real t = ct[icmf];
            auto it = std::find_if(m_cf_crse_data.begin(), m_cf_crse_data.end(),
                                   [=] (auto const& x) {
                                       return amrex::almostEqual(x.first,t,5);
                                   });

            if (it == std::end(m_cf_crse_data)) {
                MF mf_crse_patch = make_mf_crse_patch<MF>(fpc, m_ncomp);
                mf_crse_patch.ParallelCopy(*cmf[icmf], m_cgeom.periodicity());

                std::pair<Real,std::unique_ptr<MF>> tmp;
                tmp.first = t;
                tmp.second = std::make_unique<MF>(std::move(mf_crse_patch));
                m_cf_crse_data.push_back(std::move(tmp));
            }
        }

        MF mf_crse_patch;
        if (m_cf_crse_data.size() > 0 &&
            amrex::almostEqual(time, m_cf_crse_data[0].first,5))
        {
            mf_crse_patch = MF(*m_cf_crse_data[0].second, amrex::make_alias,
                               scomp, ncomp);
        }
        else if (m_cf_crse_data.size() > 1 &&
                 amrex::almostEqual(time, m_cf_crse_data[1].first,5))
        {
            mf_crse_patch = MF(*m_cf_crse_data[1].second, amrex::make_alias,
                               scomp, ncomp);
        }
        else if (m_cf_crse_data.size() == 2)
        {
            mf_crse_patch = make_mf_crse_patch<MF>(fpc, ncomp);
            int const ng_space_interp = 8; // Need to be big enough
            Box domain = m_cgeom.growPeriodicDomain(ng_space_interp);
            domain.convert(mf.ixType());
            Real t0 = m_cf_crse_data[0].first;
            Real t1 = m_cf_crse_data[1].first;
            Real alpha = (t1-time)/(t1-t0);
            Real beta = (time-t0)/(t1-t0);
            AMREX_ASSERT(alpha >= 0._rt && beta >= 0._rt);
            auto const& a = mf_crse_patch.arrays();
            auto const& a0 = m_cf_crse_data[0].second->const_arrays();
            auto const& a1 = m_cf_crse_data[1].second->const_arrays();
            amrex::ParallelFor(mf_crse_patch, IntVect(0), ncomp,
                               [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k, int n) noexcept
                               {
                                   if (domain.contains(i,j,k)) {
                                       a[bi](i,j,k,n)
                                           = alpha*a0[bi](i,j,k,scomp+n)
                                           +  beta*a1[bi](i,j,k,scomp+n);
                                   }
                               });
            Gpu::streamSynchronize();
        }
        else
        {
            amrex::Abort("FillPatcher: High order interpolation in time not supported.  Or FillPatcher was not properly deleted.");
        }

        cbc(mf_crse_patch, 0, ncomp, nghost, time, cbccomp);

        pre_interp(mf_crse_patch, 0, ncomp);

        FillPatchInterp(*m_cf_fine_data, scomp, mf_crse_patch, 0,
                        ncomp, IntVect(0), m_cgeom, m_fgeom,
                        amrex::grow(amrex::convert(m_fgeom.Domain(),
                                                   mf.ixType()),nghost),
                        ratio, m_interp, bcs, bcscomp);

        post_interp(*m_cf_fine_data, scomp, ncomp);

        mf.ParallelCopy(*m_cf_fine_data, scomp, dcomp, ncomp, IntVect{0}, nghost);
    }
}

}

#endif
