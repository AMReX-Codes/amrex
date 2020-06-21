
#include <AMReX_YAFluxRegister.H>
#include <AMReX_YAFluxRegister_K.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

YAFluxRegister::YAFluxRegister (const BoxArray& fba, const BoxArray& cba,
                                const DistributionMapping& fdm, const DistributionMapping& cdm,
                                const Geometry& fgeom, const Geometry& cgeom,
                                const IntVect& ref_ratio, int fine_lev, int nvar)
{
    define(fba, cba, fdm, cdm, fgeom, cgeom, ref_ratio, fine_lev, nvar);
}

void
YAFluxRegister::define (const BoxArray& fba, const BoxArray& cba,
                        const DistributionMapping& fdm, const DistributionMapping& cdm,
                        const Geometry& fgeom, const Geometry& cgeom,
                        const IntVect& ref_ratio, int fine_lev, int nvar)
{
    m_fine_geom = fgeom;
    m_crse_geom = cgeom;
    m_ratio = ref_ratio;
    m_fine_level = fine_lev;
    m_ncomp = nvar;

    m_crse_data.define(cba, cdm, nvar, 0, MFInfo(), FArrayBoxFactory());

    m_crse_flag.define(cba, cdm, 1, 1, MFInfo(), DefaultFabFactory<IArrayBox>());

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
#ifdef _OPENMP
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
    const int n_cfba = cfba.size();
    cfba.uniqify();

#ifdef _OPENMP
    
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
            const int ntmp = bl_tmp.size();
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
        const int ntmp = bl_tmp.size();
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
    m_cfpatch.define(cfp_ba, cfp_dm, nvar, 0, MFInfo(), FArrayBoxFactory());

    m_cfp_fab.resize(nlocal);
    for (MFIter mfi(m_cfpatch); mfi.isValid(); ++mfi)
    {
        const int li = mfi.LocalIndex();
        const int flgi = m_cfp_localindex[li];
        FArrayBox& fab = m_cfpatch[mfi];
        m_cfp_fab[flgi].push_back(&fab);
    }

    bool is_periodic = m_fine_geom.isAnyPeriodic();
    if (is_periodic) {
        m_cfp_mask.define(cfp_ba, cfp_dm, 1, 0, MFInfo(), FArrayBoxFactory());
        m_cfp_mask.setVal(1.0);

        Vector<Array4BoxTag<Real> > tags;

        bool run_on_gpu = Gpu::inLaunchRegion();

        const Box& domainbox = m_crse_geom.Domain();

#ifdef _OPENMP
#pragma omp parallel if (!run_on_gpu)
#endif
        {
            std::vector< std::pair<int,Box> > isects;

            for (MFIter mfi(m_cfp_mask); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.fabbox();
                if (!domainbox.contains(bx))  // part of the box is outside periodic boundary
                {
                    FArrayBox& fab = m_cfp_mask[mfi];
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
                                    fab.setVal<RunOn::Host>(0.0, ibx);
                                }
                            }
                        }
                    }
                }
            }
        }

#ifdef AMREX_USE_GPU
        amrex::ParallelFor(tags, 1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, Array4<Real> const& a) noexcept
        {
            a(i,j,k,n) = 0.0;
        });
#endif
    }
}


void
YAFluxRegister::reset ()
{
    m_crse_data.setVal(0.0);
    m_cfpatch.setVal(0.0);
}


void
YAFluxRegister::CrseAdd (const MFIter& mfi,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                         const Real* dx, Real dt, RunOn runon) noexcept
{
    BL_ASSERT(m_crse_data.nComp() == flux[0]->nComp());

    if (m_crse_fab_flag[mfi.LocalIndex()] == crse_cell) {
        return;  // this coarse fab is not close to fine fabs.
    }

    const Box& bx = mfi.tilebox();
    const int nc = m_crse_data.nComp();
    AMREX_D_TERM(const Real dtdx = dt/dx[0];,
                 const Real dtdy = dt/dx[1];,
                 const Real dtdz = dt/dx[2];);
    AMREX_D_TERM(FArrayBox const* fx = flux[0];,
                 FArrayBox const* fy = flux[1];,
                 FArrayBox const* fz = flux[2];);

    auto fab = m_crse_data.array(mfi);
    auto const flag = m_crse_flag.array(mfi);

    AMREX_D_TERM(Array4<Real const> fxarr = fx->const_array();,
                 Array4<Real const> fyarr = fy->const_array();,
                 Array4<Real const> fzarr = fz->const_array(););

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, bx, tbx,
    {
        yafluxreg_crseadd(tbx, fab, flag, AMREX_D_DECL(fxarr,fyarr,fzarr),
                          AMREX_D_DECL(dtdx,dtdy,dtdz),nc);
    });
}


void
YAFluxRegister::FineAdd (const MFIter& mfi,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& a_flux,
                         const Real* dx, Real dt, RunOn runon) noexcept
{
    BL_ASSERT(m_cfpatch.nComp() == a_flux[0]->nComp());

    const int li = mfi.LocalIndex();
    Vector<FArrayBox*>& cfp_fabs = m_cfp_fab[li];
    if (cfp_fabs.empty()) return;

    const Box& tbx = mfi.tilebox();
    const Box& bx = amrex::coarsen(tbx, m_ratio);
    const Box& fbx = amrex::refine(bx, m_ratio);
    const int nc = m_cfpatch.nComp();

    const Real ratio = static_cast<Real>(AMREX_D_TERM(m_ratio[0],*m_ratio[1],*m_ratio[2]));
    std::array<Real,AMREX_SPACEDIM> dtdx{AMREX_D_DECL(dt/(dx[0]*ratio),
                                                      dt/(dx[1]*ratio),
                                                      dt/(dx[2]*ratio))};
    const Dim3 rr = m_ratio.dim3();

    std::array<FArrayBox const*,AMREX_SPACEDIM> flux{AMREX_D_DECL(a_flux[0],a_flux[1],a_flux[2])};
    bool use_gpu = (runon == RunOn::Gpu) && Gpu::inLaunchRegion();
    amrex::ignore_unused(use_gpu);
    std::array<FArrayBox,AMREX_SPACEDIM> ftmp;
    if (fbx != tbx) {
        AMREX_ASSERT(!use_gpu);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const Box& b = amrex::surroundingNodes(fbx,idim);
            ftmp[idim].resize(b,nc);
            ftmp[idim].setVal<RunOn::Host>(0.0);
            ftmp[idim].copy<RunOn::Host>(*a_flux[idim]);
            flux[idim] = &ftmp[idim];
        }
    }
    
    AMREX_ASSERT(bx.cellCentered());

    for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
    {
        const Box& lobx = amrex::adjCellLo(bx, idim);
        const Box& hibx = amrex::adjCellHi(bx, idim);
        FArrayBox const* f = flux[idim];
        for (FArrayBox* cfp : cfp_fabs)
        {
            {
                const Box& lobx_is = lobx & cfp->box();
                const int side = 0;
                if (lobx_is.ok())
                {
                    auto d = cfp->array();
                    Real dtdxs = dtdx[idim];
                    int dirside = idim*2+side;
                    Array4<Real const> farr = f->const_array();
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
                    Real dtdxs = dtdx[idim];
                    int dirside = idim*2+side;
                    Array4<Real const> farr = f->const_array();
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(runon, hibx_is, tmpbox,
                    {
                        yafluxreg_fineadd(tmpbox, d, farr, dtdxs, nc, dirside, rr);
                    });
                }
            }
        }
    }
}


void
YAFluxRegister::Reflux (MultiFab& state, int dc)
{
    if (!m_cfp_mask.empty())
    {
        const int ncomp = m_ncomp;
#ifdef _OPENMP
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
    MultiFab::Add(state, m_crse_data, 0, dc, m_ncomp, 0);
}

}
