#ifndef AMREX_FILLPATCHER_H_
#define AMREX_FILLPATCHER_H_
#include <AMReX_Config.H>

#include <AMReX_FillPatchUtil.H>
#include <utility>

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
 *
 * This class also provides support for RungeKutta::RK3 and RungeKutta::RK4.
 * The storeRKCoarseData function can be used to store coarse AMR level
 * data that are needed for filling fine level data's ghost cells in this
 * class.  The `fillRK` function can be used to fill ghost cells for fine
 * AMR levels.  This operation at the coarse/fine boundary is non-trivial
 * for RK orders higher than 2.  Note that it is expected that time stepping
 * on the coarse level is perform before any fine level time stepping, and
 * it's the user's reponsibility to properly create and destroy this object.
 * See AmrLevel::RK for an example of using the RungeKutta functions and
 * FillPatcher together.
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
                 BoxArray const& cba, DistributionMapping const& cdm, // NOLINT
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
               Vector<MF*> const& cmf, Vector<Real> const& ct,
               Vector<MF*> const& fmf, Vector<Real> const& ft,
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
                                 Vector<MF*> const& cmf,
                                 Vector<Real> const& ct,
                                 int scomp, int dcomp, int ncomp,
                                 BC& cbc, int cbccomp,
                                 Vector<BCRec> const& bcs, int bcscomp,
                                 PreInterpHook const& pre_interp = {},
                                 PostInterpHook const& post_interp = {});

    /**
     * \brief Store coarse AMR level data for RK3 and RK4
     *
     * \tparam order RK order.  Must be 3 or 4.
     * \param time   time at the beginning of the step
     * \param dt     time step
     * \param S_old  data at time
     * \param RK_k   right-hand side at RK stages
     */
    template <std::size_t order>
    void storeRKCoarseData (Real time, Real dt, MF const& S_old,
                            Array<MF,order> const& RK_k);

    /**
     * \brief Fill ghost cells of fine AMR level for RK3 and RK4
     *
     * \param stage     RK stage number starting from 1
     * \param iteration iteration number on fine level during a coarse time
     *                  step.  For an AMR simulation with subcycling and a
     *                  refinement ratio of 2, the number is either 1 or 2,
     *                  denoting the first and second substep, respectively.
     * \param ncycle    number of subcyling steps.  It's usually 2 or 4.
     *                  Without subcycling, this will be 1.
     * \param cbc       filling physical boundary on coarse level
     * \param fbc       filling physical boundary on fine level
     * \param bcs       physical BC types
     */
    template <typename BC>
    void fillRK (int stage, int iteration, int ncycle, MF& mf, Real time,
                 BC& cbc, BC& fbc, Vector<BCRec> const& bcs);

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
    MF m_sfine;
    IntVect m_ratio;
    Vector<std::pair<Real,std::unique_ptr<MF>>> m_cf_crse_data;
    std::unique_ptr<MF> m_cf_crse_data_tmp;
    std::unique_ptr<MF> m_cf_fine_data;
    Real m_dt_coarse = std::numeric_limits<Real>::lowest();

    FabArrayBase::FPinfo const& getFPinfo ();
};

template <class MF>
FillPatcher<MF>::FillPatcher (BoxArray const& fba, DistributionMapping const& fdm,
                              Geometry const& fgeom,
                              BoxArray const& cba, DistributionMapping const& cdm, // NOLINT
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
      m_eb_index_space(eb_index_space),
      m_sfine(fba, fdm, 1, nghost, MFInfo().SetAlloc(false))
{
    static_assert(IsFabArray<MF>::value,
                  "FillPatcher<MF>: MF must be FabArray type");
    AMREX_ALWAYS_ASSERT(m_fba.ixType().cellCentered() || m_fba.ixType().nodeCentered());

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_ratio[idim] = m_fgeom.Domain().length(idim) / m_cgeom.Domain().length(idim);
    }
    AMREX_ASSERT(m_fgeom.Domain() == amrex::refine(m_cgeom.Domain(),m_ratio));
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
FabArrayBase::FPinfo const&
FillPatcher<MF>::getFPinfo ()
{
    const InterpolaterBoxCoarsener& coarsener = m_interp->BoxCoarsener(m_ratio);
    return FabArrayBase::TheFPinfo(m_sfine, m_sfine, m_nghost, coarsener,
                                   m_fgeom, m_cgeom, m_eb_index_space);
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

    auto const& fpc = getFPinfo();

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

        if (m_cf_crse_data_tmp == nullptr) {
            m_cf_crse_data_tmp = std::make_unique<MF>
                (make_mf_crse_patch<MF>(fpc, m_ncomp));
        }

        int const ng_space_interp = 8; // Need to be big enough
        Box domain = m_cgeom.growPeriodicDomain(ng_space_interp);
        domain.convert(mf.ixType());

        int idata = -1;
        if (m_cf_crse_data.size() == 1) {
            idata = 0;
        } else if (m_cf_crse_data.size() == 2) {
            Real const teps = std::abs(m_cf_crse_data[1].first -
                                       m_cf_crse_data[0].first) * 1.e-3_rt;
            if (time > m_cf_crse_data[0].first - teps &&
                time < m_cf_crse_data[0].first + teps) {
                idata = 0;
            } else if (time > m_cf_crse_data[1].first - teps &&
                       time < m_cf_crse_data[1].first + teps) {
                idata = 1;
            } else {
                idata = 2;
            }
        }

        if (idata == 0 || idata == 1) {
            auto const& dst = m_cf_crse_data_tmp->arrays();
            auto const& src = m_cf_crse_data[idata].second->const_arrays();
            amrex::ParallelFor(*m_cf_crse_data_tmp, IntVect(0), ncomp,
                               [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k, int n) noexcept
                               {
                                   if (domain.contains(i,j,k)) {
                                       dst[bi](i,j,k,n) = src[bi](i,j,k,n+scomp);
                                   }
                               });
        } else if (idata == 2) {
            Real t0 = m_cf_crse_data[0].first;
            Real t1 = m_cf_crse_data[1].first;
            Real alpha = (t1-time)/(t1-t0);
            Real beta = (time-t0)/(t1-t0);
            auto const& a = m_cf_crse_data_tmp->arrays();
            auto const& a0 = m_cf_crse_data[0].second->const_arrays();
            auto const& a1 = m_cf_crse_data[1].second->const_arrays();
            amrex::ParallelFor(*m_cf_crse_data_tmp, IntVect(0), ncomp,
                               [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k, int n) noexcept
                               {
                                   if (domain.contains(i,j,k)) {
                                       a[bi](i,j,k,n)
                                           = alpha*a0[bi](i,j,k,scomp+n)
                                           +  beta*a1[bi](i,j,k,scomp+n);
                                   }
                               });
        }
        else
        {
            amrex::Abort("FillPatcher: High order interpolation in time not supported.  Or FillPatcher was not properly deleted.");
        }
        Gpu::streamSynchronize();

        cbc(*m_cf_crse_data_tmp, 0, ncomp, nghost, time, cbccomp);

        detail::call_interp_hook(pre_interp, *m_cf_crse_data_tmp, 0, ncomp);

        FillPatchInterp(*m_cf_fine_data, scomp, *m_cf_crse_data_tmp, 0,
                        ncomp, IntVect(0), m_cgeom, m_fgeom,
                        amrex::grow(amrex::convert(m_fgeom.Domain(),
                                                   mf.ixType()),nghost),
                        m_ratio, m_interp, bcs, bcscomp);

        detail::call_interp_hook(post_interp, *m_cf_fine_data, scomp, ncomp);

        mf.ParallelCopy(*m_cf_fine_data, scomp, dcomp, ncomp, IntVect{0}, nghost);
    }
}

template <typename MF>
template <std::size_t order>
void FillPatcher<MF>::storeRKCoarseData (Real /*time*/, Real dt, MF const& S_old,
                                         Array<MF,order> const& RK_k)
{
    m_dt_coarse = dt;
    m_cf_crse_data.resize(order+1);

    auto const& fpc = getFPinfo();

    for (auto& tmf : m_cf_crse_data) {
        tmf.first = std::numeric_limits<Real>::lowest(); // because we dont' need it
        tmf.second = std::make_unique<MF>(make_mf_crse_patch<MF>(fpc, m_ncomp));
    }
    m_cf_crse_data[0].second->ParallelCopy(S_old, m_cgeom.periodicity());
    for (std::size_t i = 0; i < order; ++i) {
        m_cf_crse_data[i+1].second->ParallelCopy(RK_k[i], m_cgeom.periodicity());
    }
}

template <typename MF>
template <typename BC>
void FillPatcher<MF>::fillRK (int stage, int iteration, int ncycle,
                              MF& mf, Real time, BC& cbc, BC& fbc,
                              Vector<BCRec> const& bcs)
{
    int rk_order = m_cf_crse_data.size()-1;
    if (rk_order != 3 && rk_order != 4) {
        amrex::Abort("FillPatcher: unsupported RK order "+std::to_string(rk_order));
        return;
    }
    AMREX_ASSERT(stage > 0 && stage <= rk_order);

    auto const& fpc = getFPinfo();
    if (m_cf_crse_data_tmp == nullptr) {
        m_cf_crse_data_tmp = std::make_unique<MF>
            (make_mf_crse_patch<MF>(fpc, m_ncomp));
    }

    auto const& u = m_cf_crse_data_tmp->arrays();
    auto const& u0 = m_cf_crse_data[0].second->const_arrays();
    auto const& k1 = m_cf_crse_data[1].second->const_arrays();
    auto const& k2 = m_cf_crse_data[2].second->const_arrays();
    auto const& k3 = m_cf_crse_data[3].second->const_arrays();

    Real dtc = m_dt_coarse;
    Real r = Real(1) / Real(ncycle);
    Real xsi = Real(iteration-1) / Real(ncycle);

    int const ng_space_interp = 8; // Need to be big enough
    Box cdomain = m_cgeom.growPeriodicDomain(ng_space_interp);
    cdomain.convert(m_cf_crse_data_tmp->ixType());

    if (rk_order == 3) {
        // coefficients for U
        Real b1 = xsi - Real(5./6.)*xsi*xsi;
        Real b2 = Real(1./6.)*xsi*xsi;
        Real b3 = Real(2./3)*xsi*xsi;
        // coefficients for Ut
        Real c1 = Real(1.) - Real(5./3.)*xsi;
        Real c2 = Real(1./3.)*xsi;
        Real c3 = Real(4./3.)*xsi;
        // coefficients for Utt
        constexpr Real d1 = Real(-5./3.);
        constexpr Real d2 = Real(1./3.);
        constexpr Real d3 = Real(4./3.);
        if (stage == 1) {
            amrex::ParallelFor(*m_cf_crse_data_tmp, IntVect(0), m_ncomp,
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k, int n) noexcept
            {
                if (cdomain.contains(i,j,k)) {
                    Real kk1 = k1[bi](i,j,k,n);
                    Real kk2 = k2[bi](i,j,k,n);
                    Real kk3 = k3[bi](i,j,k,n);
                    Real uu  = b1*kk1 + b2*kk2 + b3*kk3;
                    u[bi](i,j,k,n) = u0[bi](i,j,k,n) + dtc*uu;
                }
            });
        } else if (stage == 2) {
            amrex::ParallelFor(*m_cf_crse_data_tmp, IntVect(0), m_ncomp,
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k, int n) noexcept
            {
                if (cdomain.contains(i,j,k)) {
                    Real kk1 = k1[bi](i,j,k,n);
                    Real kk2 = k2[bi](i,j,k,n);
                    Real kk3 = k3[bi](i,j,k,n);
                    Real uu = b1*kk1 + b2*kk2 + b3*kk3;
                    Real ut = c1*kk1 + c2*kk2 + c3*kk3;
                    u[bi](i,j,k,n) = u0[bi](i,j,k,n) + dtc*(uu + r*ut);
                }
            });
        } else if (stage == 3) {
            amrex::ParallelFor(*m_cf_crse_data_tmp, IntVect(0), m_ncomp,
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k, int n) noexcept
            {
                if (cdomain.contains(i,j,k)) {
                    Real kk1 = k1[bi](i,j,k,n);
                    Real kk2 = k2[bi](i,j,k,n);
                    Real kk3 = k3[bi](i,j,k,n);
                    Real uu  = b1*kk1 + b2*kk2 + b3*kk3;
                    Real ut  = c1*kk1 + c2*kk2 + c3*kk3;
                    Real utt = d1*kk1 + d2*kk2 + d3*kk3;
                    u[bi](i,j,k,n) = u0[bi](i,j,k,n) + dtc*
                        (uu + Real(0.5)*r*ut + Real(0.25)*r*r*utt);
                }
            });
        }
    } else if (rk_order == 4) {
        auto const& k4 = m_cf_crse_data[4].second->const_arrays();
        Real xsi2 = xsi*xsi;
        Real xsi3 = xsi2*xsi;
        // coefficients for U
        Real b1 = xsi - Real(1.5)*xsi2 + Real(2./3.)*xsi3;
        Real b2 = xsi2 - Real(2./3.)*xsi3;
        Real b3 = b2;
        Real b4 = Real(-0.5)*xsi2 + Real(2./3.)*xsi3;
        // coefficients for Ut
        Real c1 = Real(1.) - Real(3.)*xsi + Real(2.)*xsi2;
        Real c2 = Real(2.)*xsi - Real(2.)*xsi2;
        Real c3 = c2;
        Real c4 = -xsi + Real(2.)*xsi2;
        // coefficients for Utt
        Real d1 = Real(-3.) + Real(4.)*xsi;
        Real d2 = Real( 2.) - Real(4.)*xsi;
        Real d3 = d2;
        Real d4 = Real(-1.) + Real(4.)*xsi;
        // coefficients for Uttt
        constexpr Real e1 = Real( 4.);
        constexpr Real e2 = Real(-4.);
        constexpr Real e3 = Real(-4.);
        constexpr Real e4 = Real( 4.);
        if (stage == 1) {
            amrex::ParallelFor(*m_cf_crse_data_tmp, IntVect(0), m_ncomp,
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k, int n) noexcept
            {
                if (cdomain.contains(i,j,k)) {
                    Real kk1 = k1[bi](i,j,k,n);
                    Real kk2 = k2[bi](i,j,k,n);
                    Real kk3 = k3[bi](i,j,k,n);
                    Real kk4 = k4[bi](i,j,k,n);
                    Real uu  = b1*kk1 + b2*kk2 + b3*kk3 + b4*kk4;
                    u[bi](i,j,k,n) = u0[bi](i,j,k,n) + dtc*uu;
                }
            });
        } else if (stage == 2) {
            amrex::ParallelFor(*m_cf_crse_data_tmp, IntVect(0), m_ncomp,
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k, int n) noexcept
            {
                if (cdomain.contains(i,j,k)) {
                    Real kk1 = k1[bi](i,j,k,n);
                    Real kk2 = k2[bi](i,j,k,n);
                    Real kk3 = k3[bi](i,j,k,n);
                    Real kk4 = k4[bi](i,j,k,n);
                    Real uu = b1*kk1 + b2*kk2 + b3*kk3 + b4*kk4;
                    Real ut = c1*kk1 + c2*kk2 + c3*kk3 + c4*kk4;
                    u[bi](i,j,k,n) = u0[bi](i,j,k,n) + dtc*(uu + Real(0.5)*r*ut);
                }
            });
        } else if (stage == 3 || stage == 4) {
            Real r2 = r*r;
            Real r3 = r2*r;
            Real at = (stage == 3) ? Real(0.5)*r : r;
            Real att = (stage == 3) ? Real(0.25)*r2 : Real(0.5)*r2;
            Real attt = (stage == 3) ? Real(0.0625)*r3 : Real(0.125)*r3;
            Real akk = (stage == 3) ? Real(-4.) : Real(4.);
            amrex::ParallelFor(*m_cf_crse_data_tmp, IntVect(0), m_ncomp,
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k, int n) noexcept
            {
                if (cdomain.contains(i,j,k)) {
                    Real kk1 = k1[bi](i,j,k,n);
                    Real kk2 = k2[bi](i,j,k,n);
                    Real kk3 = k3[bi](i,j,k,n);
                    Real kk4 = k4[bi](i,j,k,n);
                    Real uu   = b1*kk1 + b2*kk2 + b3*kk3 + b4*kk4;
                    Real ut   = c1*kk1 + c2*kk2 + c3*kk3 + c4*kk4;
                    Real utt  = d1*kk1 + d2*kk2 + d3*kk3 + d4*kk4;
                    Real uttt = e1*kk1 + e2*kk2 + e3*kk3 + e4*kk4;
                    u[bi](i,j,k,n) = u0[bi](i,j,k,n) + dtc *
                        (uu + at*ut + att*utt + attt*(uttt+akk*(kk3-kk2)));
                }
            });
        }
    }
    Gpu::streamSynchronize();

    cbc(*m_cf_crse_data_tmp, 0, m_ncomp, m_nghost, time, 0);

    if (m_cf_fine_data == nullptr) {
        m_cf_fine_data = std::make_unique<MF>(make_mf_fine_patch<MF>(fpc, m_ncomp));
    }

    FillPatchInterp(*m_cf_fine_data, 0, *m_cf_crse_data_tmp, 0,
                    m_ncomp, IntVect(0), m_cgeom, m_fgeom,
                    amrex::grow(amrex::convert(m_fgeom.Domain(),
                                               mf.ixType()),m_nghost),
                    m_ratio, m_interp, bcs, 0);

    // xxxxx We can optimize away this ParallelCopy by making a special fpinfo.
    mf.ParallelCopy(*m_cf_fine_data, 0, 0, m_ncomp, IntVect(0), m_nghost);

    mf.FillBoundary(m_fgeom.periodicity());
    fbc(mf, 0, m_ncomp, m_nghost, time, 0);
}

}

#endif
