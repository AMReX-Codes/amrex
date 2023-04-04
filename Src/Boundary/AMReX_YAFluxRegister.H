#ifndef AMREX_YAFLUXREGISTER_H_
#define AMREX_YAFLUXREGISTER_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_YAFluxRegister_K.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex {

/**
  YAFluxRegister is yet another Flux Register for refluxing.

  At the beginning of a coarse step, `reset()` is called.  In MFIter
  for the coarse level advance, `CrseAdd` is called with coarse flux.
  The flux is not scaled.  In MFIter for the fine level advance,
  `FineAdd` is called.  After the fine level finished its time steps,
  `Reflux` is called to update the coarse cells next to the
  coarse/fine boundary.
*/

template <typename MF>
class YAFluxRegisterT
{
public:

    using T = typename MF::value_type;
    using FAB = typename MF::fab_type;

    YAFluxRegisterT () = default;

    YAFluxRegisterT (const BoxArray& fba, const BoxArray& cba,
                     const DistributionMapping& fdm, const DistributionMapping& cdm,
                     const Geometry& fgeom, const Geometry& cgeom,
                     const IntVect& ref_ratio, int fine_lev, int nvar);

    void define (const BoxArray& fba, const BoxArray& cba,
                 const DistributionMapping& fdm, const DistributionMapping& cdm,
                 const Geometry& fgeom, const Geometry& cgeom,
                 const IntVect& ref_ratio, int fine_lev, int nvar);

    void reset ();

    void CrseAdd (const MFIter& mfi,
                  const std::array<FAB const*, AMREX_SPACEDIM>& flux,
                  const Real* dx, Real dt, RunOn runon) noexcept;

    void FineAdd (const MFIter& mfi,
                  const std::array<FAB const*, AMREX_SPACEDIM>& flux,
                  const Real* dx, Real dt, RunOn runon) noexcept;

    void Reflux (MF& state, int dc = 0);

    bool CrseHasWork (const MFIter& mfi) const noexcept {
        return m_crse_fab_flag[mfi.LocalIndex()] != crse_cell;
    }

    bool FineHasWork (const MFIter& mfi) const noexcept {
        return !(m_cfp_fab[mfi.LocalIndex()].empty());
    }

    MF& getFineData ();

    MF& getCrseData ();

    enum CellType : int {
        // must be same as in AMReX_YAFluxRegiser_K.H
        crse_cell = 0, crse_fine_boundary_cell, fine_cell
    };

protected:

    MF m_crse_data;
    iMultiFab m_crse_flag;
    Vector<int> m_crse_fab_flag;

    MF m_cfpatch;                   //!< This is built on crse/fine patches
    MF m_cfp_mask;
    Vector<Vector<FAB*> > m_cfp_fab;  //!< The size of this is (# of local fine grids (# of crse/fine patches for that grid))
    Vector<int> m_cfp_localindex;

    Geometry m_fine_geom;
    Geometry m_crse_geom;

    IntVect m_ratio;
    int m_fine_level;
    int m_ncomp;
};

template <typename MF>
YAFluxRegisterT<MF>::YAFluxRegisterT (const BoxArray& fba, const BoxArray& cba,
                                      const DistributionMapping& fdm, const DistributionMapping& cdm,
                                      const Geometry& fgeom, const Geometry& cgeom,
                                      const IntVect& ref_ratio, int fine_lev, int nvar)
{
    define(fba, cba, fdm, cdm, fgeom, cgeom, ref_ratio, fine_lev, nvar);
}

template <typename MF>
void
YAFluxRegisterT<MF>::define (const BoxArray& fba, const BoxArray& cba,
                             const DistributionMapping& fdm, const DistributionMapping& cdm,
                             const Geometry& fgeom, const Geometry& cgeom,
                             const IntVect& ref_ratio, int fine_lev, int nvar)
{
    m_fine_geom = fgeom;
    m_crse_geom = cgeom;
    m_ratio = ref_ratio;
    m_fine_level = fine_lev;
    m_ncomp = nvar;

    m_crse_data.define(cba, cdm, nvar, 0);

    m_crse_flag.define(cba, cdm, 1, 1);

    const auto& cperiod = m_crse_geom.periodicity();
    const std::vector<IntVect>& pshifts = cperiod.shiftIntVect();

    BoxArray cfba = fba;
    cfba.coarsen(ref_ratio);

    Box cdomain = m_crse_geom.Domain();
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (m_crse_geom.isPeriodic(idim)) {
            cdomain.grow(idim,1);
        }
    }

    m_crse_fab_flag.resize(m_crse_flag.local_size(), crse_cell);

    m_crse_flag.setVal(crse_cell);
    {
        iMultiFab foo(cfba, fdm, 1, 1, MFInfo().SetAlloc(false));
        const FabArrayBase::CPC& cpc1 = m_crse_flag.getCPC(IntVect(1), foo, IntVect(1), cperiod);
        m_crse_flag.setVal(crse_fine_boundary_cell, cpc1, 0, 1);
        const FabArrayBase::CPC& cpc0 = m_crse_flag.getCPC(IntVect(1), foo, IntVect(0), cperiod);
        m_crse_flag.setVal(fine_cell, cpc0, 0, 1);
        auto recv_layout_mask = m_crse_flag.RecvLayoutMask(cpc0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_crse_flag); mfi.isValid(); ++mfi) {
            if (recv_layout_mask[mfi]) {
                m_crse_fab_flag[mfi.LocalIndex()] = fine_cell;
            }
        }
    }

    BoxList cfp_bl;
    Vector<int> cfp_procmap;
    int nlocal = 0;
    const int myproc = ParallelDescriptor::MyProc();
    const auto n_cfba = static_cast<int>(cfba.size());
    cfba.uniqify();

#ifdef AMREX_USE_OMP

    const int nthreads = omp_get_max_threads();
    Vector<BoxList> bl_priv(nthreads, BoxList());
    Vector<Vector<int> > procmap_priv(nthreads);
    Vector<Vector<int> > localindex_priv(nthreads);
#pragma omp parallel
    {
        BoxList bl_tmp;
        const int tid = omp_get_thread_num();
        BoxList& bl = bl_priv[tid];
        Vector<int>& pmp = procmap_priv[tid];
        Vector<int>& lid = localindex_priv[tid];
#pragma omp for
        for (int i = 0; i < n_cfba; ++i)
        {
            Box bx = amrex::grow(cfba[i], 1);
            bx &= cdomain;

            cfba.complementIn(bl_tmp, bx);
            const auto ntmp = static_cast<int>(bl_tmp.size());
            bl.join(bl_tmp);

            int proc = fdm[i];
            for (int j = 0; j < ntmp; ++j) {
                pmp.push_back(proc);
            }

            if (proc == myproc) {
                lid.push_back(ntmp);
            }
        }
    }

    for (auto const& bl : bl_priv) {
        cfp_bl.join(bl);
    }

    for (auto const& pmp : procmap_priv) {
        cfp_procmap.insert(std::end(cfp_procmap), std::begin(pmp), std::end(pmp));
    }

    for (auto& lid : localindex_priv) {
        for (int nl : lid) {
            for (int j = 0; j < nl; ++j) {
                m_cfp_localindex.push_back(nlocal);
            }
            ++nlocal;
        }
    }

#else

    BoxList bl_tmp;
    for (int i = 0; i < n_cfba; ++i)
    {
        Box bx = amrex::grow(cfba[i], 1);
        bx &= cdomain;

        cfba.complementIn(bl_tmp, bx);
        const auto ntmp = static_cast<int>(bl_tmp.size());
        cfp_bl.join(bl_tmp);

        int proc = fdm[i];
        for (int j = 0; j < ntmp; ++j) {
            cfp_procmap.push_back(proc);
        }

        if (proc == myproc) {
            for (int j = 0; j < ntmp; ++j) {
                m_cfp_localindex.push_back(nlocal);  // This Array store local index in fine ba/dm.
            }                                        // Its size is local size of cfp.
            ++nlocal;
        }
    }

#endif

    // It's safe even if cfp_bl is empty.

    BoxArray cfp_ba(std::move(cfp_bl));
    DistributionMapping cfp_dm(std::move(cfp_procmap));
    m_cfpatch.define(cfp_ba, cfp_dm, nvar, 0);

    m_cfp_fab.resize(nlocal);
    for (MFIter mfi(m_cfpatch); mfi.isValid(); ++mfi)
    {
        const int li = mfi.LocalIndex();
        const int flgi = m_cfp_localindex[li];
        FAB& fab = m_cfpatch[mfi];
        m_cfp_fab[flgi].push_back(&fab);
    }

    bool is_periodic = m_fine_geom.isAnyPeriodic();
    if (is_periodic) {
        m_cfp_mask.define(cfp_ba, cfp_dm, 1, 0);
        m_cfp_mask.setVal(T(1.0));

        Vector<Array4BoxTag<T> > tags;

        bool run_on_gpu = Gpu::inLaunchRegion();

        const Box& domainbox = m_crse_geom.Domain();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (!run_on_gpu)
#endif
        {
            std::vector< std::pair<int,Box> > isects;

            for (MFIter mfi(m_cfp_mask); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.fabbox();
                if (!domainbox.contains(bx))  // part of the box is outside periodic boundary
                {
                    FAB& fab = m_cfp_mask[mfi];
                    auto const& arr = m_cfp_mask.array(mfi);
                    for (const auto& iv : pshifts)
                    {
                        if (iv != IntVect::TheZeroVector())
                        {
                            cfba.intersections(bx+iv, isects);
                            for (const auto& is : isects)
                            {
                                const Box& ibx = is.second - iv;
                                if (run_on_gpu) {
                                    tags.push_back({arr,ibx});
                                } else {
                                    fab.template setVal<RunOn::Host>(T(0.0), ibx);
                                }
                            }
                        }
                    }
                }
            }
        }

#ifdef AMREX_USE_GPU
        amrex::ParallelFor(tags, 1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, Array4BoxTag<T> const& tag) noexcept
        {
            tag.dfab(i,j,k,n) = T(0);
        });
#endif
    }
}

template <typename MF>
void
YAFluxRegisterT<MF>::reset ()
{
    m_crse_data.setVal(T(0.0));
    m_cfpatch.setVal(T(0.0));
}

template <typename MF>
void
YAFluxRegisterT<MF>::CrseAdd (const MFIter& mfi,
                              const std::array<FAB const*, AMREX_SPACEDIM>& flux,
                              const Real* dx, Real dt, RunOn runon) noexcept
{
    BL_ASSERT(m_crse_data.nComp() == flux[0]->nComp());

    if (m_crse_fab_flag[mfi.LocalIndex()] == crse_cell) {
        return;  // this coarse fab is not close to fine fabs.
    }

    const Box& bx = mfi.tilebox();
    const int nc = m_crse_data.nComp();
    AMREX_D_TERM(const T dtdx = static_cast<T>(dt/dx[0]);,
                 const T dtdy = static_cast<T>(dt/dx[1]);,
                 const T dtdz = static_cast<T>(dt/dx[2]););
    AMREX_D_TERM(FAB const* fx = flux[0];,
                 FAB const* fy = flux[1];,
                 FAB const* fz = flux[2];);

    auto fab = m_crse_data.array(mfi);
    auto const flag = m_crse_flag.const_array(mfi);

    AMREX_D_TERM(Array4<T const> fxarr = fx->const_array();,
                 Array4<T const> fyarr = fy->const_array();,
                 Array4<T const> fzarr = fz->const_array(););

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, bx, tbx,
    {
        yafluxreg_crseadd(tbx, fab, flag, AMREX_D_DECL(fxarr,fyarr,fzarr),
                          AMREX_D_DECL(dtdx,dtdy,dtdz),nc);
    });
}

template <typename MF>
void
YAFluxRegisterT<MF>::FineAdd (const MFIter& mfi,
                              const std::array<FAB const*, AMREX_SPACEDIM>& a_flux,
                              const Real* dx, Real dt, RunOn runon) noexcept
{
    BL_ASSERT(m_cfpatch.nComp() == a_flux[0]->nComp());

    const int li = mfi.LocalIndex();
    Vector<FAB*>& cfp_fabs = m_cfp_fab[li];
    if (cfp_fabs.empty()) return;

    const Box& tbx = mfi.tilebox();
    const Box& bx = amrex::coarsen(tbx, m_ratio);
    const Box& fbx = amrex::refine(bx, m_ratio);
    const int nc = m_cfpatch.nComp();

    const T ratio = static_cast<T>(AMREX_D_TERM(m_ratio[0],*m_ratio[1],*m_ratio[2]));
    std::array<T,AMREX_SPACEDIM> dtdx{{AMREX_D_DECL(static_cast<T>(dt/(dx[0]*ratio)),
                                                    static_cast<T>(dt/(dx[1]*ratio)),
                                                    static_cast<T>(dt/(dx[2]*ratio)))}};
    const Dim3 rr = m_ratio.dim3();

    std::array<FAB const*,AMREX_SPACEDIM> flux{{AMREX_D_DECL(a_flux[0],a_flux[1],a_flux[2])}};
    bool use_gpu = (runon == RunOn::Gpu) && Gpu::inLaunchRegion();
    amrex::ignore_unused(use_gpu);
    std::array<FAB,AMREX_SPACEDIM> ftmp;
    if (fbx != tbx) {
        AMREX_ASSERT(!use_gpu);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const Box& b = amrex::surroundingNodes(fbx,idim);
            ftmp[idim].resize(b,nc);
            ftmp[idim].template setVal<RunOn::Host>(T(0.0));
            ftmp[idim].template copy<RunOn::Host>(*a_flux[idim]);
            flux[idim] = &ftmp[idim];
        }
    }

    AMREX_ASSERT(bx.cellCentered());

    for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
    {
        const Box& lobx = amrex::adjCellLo(bx, idim);
        const Box& hibx = amrex::adjCellHi(bx, idim);
        FAB const* f = flux[idim];
        for (FAB* cfp : cfp_fabs)
        {
            {
                const Box& lobx_is = lobx & cfp->box();
                const int side = 0;
                if (lobx_is.ok())
                {
                    auto d = cfp->array();
                    auto dtdxs = dtdx[idim];
                    int dirside = idim*2+side;
                    Array4<T const> farr = f->const_array();
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(runon, lobx_is, tmpbox,
                    {
                        yafluxreg_fineadd(tmpbox, d, farr, dtdxs, nc, dirside, rr);
                    });
                }
            }
            {
                const Box& hibx_is = hibx & cfp->box();
                const int side = 1;
                if (hibx_is.ok())
                {
                    auto d = cfp->array();
                    auto dtdxs = dtdx[idim];
                    int dirside = idim*2+side;
                    Array4<T const> farr = f->const_array();
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(runon, hibx_is, tmpbox,
                    {
                        yafluxreg_fineadd(tmpbox, d, farr, dtdxs, nc, dirside, rr);
                    });
                }
            }
        }
    }
}

template <typename MF>
void
YAFluxRegisterT<MF>::Reflux (MF& state, int dc)
{
    if (!m_cfp_mask.empty())
    {
        const int ncomp = m_ncomp;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_cfpatch); mfi.isValid(); ++mfi)
        {
            const Box& bx = m_cfpatch[mfi].box();
            auto const maskfab = m_cfp_mask.array(mfi);
            auto       cfptfab = m_cfpatch.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                cfptfab(i,j,k,n) *= maskfab(i,j,k);
            });
        }
    }

    m_crse_data.ParallelCopy(m_cfpatch, m_crse_geom.periodicity(), FabArrayBase::ADD);

    BL_ASSERT(state.nComp() >= dc + m_ncomp);
    amrex::Add(state, m_crse_data, 0, dc, m_ncomp, 0);
}

template <typename MF>
MF&
YAFluxRegisterT<MF>::getFineData ()
{
    return m_cfpatch;
}

template <typename MF>
MF&
YAFluxRegisterT<MF>::getCrseData ()
{
    return m_crse_data;
}

using YAFluxRegister = YAFluxRegisterT<MultiFab>;

}

#endif
