#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_MLEBNodeFDLap_K.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLinOp_K.H>
#include <AMReX_MLNodeTensorLap_K.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParReduce.H>
#include <AMReX_SingleBoxCGSolver.H>

#ifdef AMREX_USE_GPU
#include <AMReX_Scan.H>
#endif

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

namespace amrex {

namespace {

// Fill ghost cells outside the domain by reflection.  A direction in which the
// data are nodal is mirrored about the boundary node, a cell-centered direction
// about the boundary face.  Periodic directions are skipped.  When the ghost
// cell is outside the domain in direction flip_dir, an EB position p is mirrored
// to mlebndfdlap_hm(p), or to mlebndfdlap_pmax() if p is 0; -1 disables this.
void fill_domain_ghost (MultiFab& mf, Geometry const& geom, int flip_dir)
{
    mf.FillBoundary(geom.periodicity());
    Box const domain = amrex::convert(geom.Domain(), mf.ixType());
    auto const dlo = amrex::lbound(domain);
    auto const dhi = amrex::ubound(domain);
    GpuArray<bool,AMREX_SPACEDIM> is_periodic;
    GpuArray<int,AMREX_SPACEDIM> nodal;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        is_periodic[idim] = geom.isPeriodic(idim);
        nodal[idim] = domain.type(idim) == IndexType::NODE ? 1 : 0;
    }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        Box const& gbx = mfi.fabbox();
        if (domain.contains(gbx)) { continue; }
        auto const& a = mf.array(mfi);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int idx[3];
            idx[0] = i; idx[1] = j; idx[2] = k;
            int lo[3];
            lo[0] = dlo.x; lo[1] = dlo.y; lo[2] = dlo.z;
            int hi[3];
            hi[0] = dhi.x; hi[1] = dhi.y; hi[2] = dhi.z;
            bool outside = false;
            bool flip = false;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (!is_periodic[d]) {
                    if (idx[d] < lo[d]) {
                        idx[d] = 2*lo[d] - idx[d] - (1-nodal[d]);
                        outside = true;
                        flip = flip || (d == flip_dir);
                    } else if (idx[d] > hi[d]) {
                        idx[d] = 2*hi[d] - idx[d] + (1-nodal[d]);
                        outside = true;
                        flip = flip || (d == flip_dir);
                    }
                }
            }
            if (outside) {
                Real const v = a(idx[0],idx[1],idx[2]);
                a(i,j,k) = !flip ? v : ((v == Real(0.0)) ? mlebndfdlap_pmax() : mlebndfdlap_hm(v));
            }
        });
    }
}

template <typename S, typename P>
void fapply_box (Box const& box, Array4<Real> const& y, Array4<Real const> const& x,
                 Array4<int const> const& dmsk, S const& sig,
                 [[maybe_unused]] bool has_eb, [[maybe_unused]] Array4<Real const> const& levset,
                 [[maybe_unused]] GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& ebp,
                 [[maybe_unused]] P const& phieb,
                 GpuArray<Real,AMREX_SPACEDIM> const& b,
                 [[maybe_unused]] bool rz, [[maybe_unused]] Real dx0, [[maybe_unused]] Real dx1,
                 [[maybe_unused]] Real xlo, [[maybe_unused]] Real alpha)
{
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
    if (has_eb) {
#if (AMREX_SPACEDIM == 2)
        if (rz) {
            AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
            {
                mlebndfdlap_adotx_rz_eb(i,j,k,y,x,levset,dmsk,ebp[0],ebp[1],sig,phieb,
                                        dx0,dx1,xlo,alpha);
            });
        } else
#endif
        {
            AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
            {
                mlebndfdlap_adotx_eb(i,j,k,y,x,levset,dmsk,AMREX_D_DECL(ebp[0],ebp[1],ebp[2]),
                                     sig,phieb,AMREX_D_DECL(b[0],b[1],b[2]));
            });
        }
    } else
#endif
    {
#if (AMREX_SPACEDIM == 2)
        if (rz) {
            AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
            {
                mlebndfdlap_adotx_rz(i,j,k,y,x,dmsk,sig,dx0,dx1,xlo,alpha);
            });
        } else
#endif
        {
            AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
            {
                mlebndfdlap_adotx(i,j,k,y,x,dmsk,sig,AMREX_D_DECL(b[0],b[1],b[2]));
            });
        }
    }
}

template <typename S>
void fsmooth_box (Box const& box, Array4<Real> const& sol, Array4<Real const> const& rhs,
                  Array4<int const> const& dmsk, S const& sig,
                  [[maybe_unused]] bool has_eb, [[maybe_unused]] Array4<Real const> const& levset,
                  [[maybe_unused]] GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& ebp,
                  GpuArray<Real,AMREX_SPACEDIM> const& b,
                  [[maybe_unused]] bool rz, [[maybe_unused]] Real dx0, [[maybe_unused]] Real dx1,
                  [[maybe_unused]] Real xlo, [[maybe_unused]] Real alpha, int redblack)
{
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
    if (has_eb) {
#if (AMREX_SPACEDIM == 2)
        if (rz) {
            AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
            {
                mlebndfdlap_gsrb_rz_eb(i,j,k,sol,rhs,levset,dmsk,ebp[0],ebp[1],sig,
                                       dx0,dx1,xlo,redblack,alpha);
            });
        } else
#endif
        {
            AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
            {
                mlebndfdlap_gsrb_eb(i,j,k,sol,rhs,levset,dmsk,AMREX_D_DECL(ebp[0],ebp[1],ebp[2]),
                                    sig,AMREX_D_DECL(b[0],b[1],b[2]),redblack);
            });
        }
    } else
#endif
    {
#if (AMREX_SPACEDIM == 2)
        if (rz) {
            AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
            {
                mlebndfdlap_gsrb_rz(i,j,k,sol,rhs,dmsk,sig,dx0,dx1,xlo,redblack,alpha);
            });
        } else
#endif
        {
            AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
            {
                mlebndfdlap_gsrb(i,j,k,sol,rhs,dmsk,sig,AMREX_D_DECL(b[0],b[1],b[2]),redblack);
            });
        }
    }
}

}

#ifdef AMREX_USE_EB
MLEBNodeFDLaplacian::MLEBNodeFDLaplacian (
    const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info,
    const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}
#endif

MLEBNodeFDLaplacian::MLEBNodeFDLaplacian (
    const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info)
{
    define(a_geom, a_grids, a_dmap, a_info);
}

void
MLEBNodeFDLaplacian::setSigma (Array<Real,AMREX_SPACEDIM> const& a_sigma) noexcept
{
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_sigma[i] = a_sigma[i];
    }
}

void
MLEBNodeFDLaplacian::setSigma (int amrlev, MultiFab const& a_sigma)
{
    m_needs_update = true;
    m_has_sigma_mf = true;
    m_sigma_mf[amrlev] = std::make_unique<MultiFab>
        (this->m_grids[amrlev][0], this->m_dmap[amrlev][0], 1, 1, MFInfo{},
         *(this->m_factory[amrlev][0]));
    MultiFab::Copy(*m_sigma_mf[amrlev], a_sigma, 0, 0, 1, 0);
#ifdef AMREX_USE_EB
    amrex::EB_set_covered(*m_sigma_mf[amrlev], Real(0.0));
#endif
}

void
MLEBNodeFDLaplacian::setRZ (bool flag) // NOLINT
{
#if (AMREX_SPACEDIM == 2)
    m_rz = flag;
#else
    amrex::ignore_unused(flag, m_rz);
#endif
}

void
MLEBNodeFDLaplacian::setAlpha (Real a_alpha) // NOLINT
{
#if (AMREX_SPACEDIM == 2)
    m_rz_alpha = a_alpha;
#else
    amrex::ignore_unused(a_alpha);
#endif
}

#ifdef AMREX_USE_EB

void
MLEBNodeFDLaplacian::setEBDirichlet (Real a_phi_eb)
{
    m_s_phi_eb = a_phi_eb;
}

void
MLEBNodeFDLaplacian::define (const Vector<Geometry>& a_geom,
                             const Vector<BoxArray>& a_grids,
                             const Vector<DistributionMapping>& a_dmap,
                             const LPInfo& a_info,
                             const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    static_assert(AMREX_SPACEDIM > 1, "MLEBNodeFDLaplacian: 1D not supported");

    BL_PROFILE("MLEBNodeFDLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    if (a_grids.size() > 1) {
        amrex::Abort("MLEBNodeFDLaplacian: multi-level not supported");
    }

    Vector<FabFactory<FArrayBox> const*> _factory;
    for (const auto *x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }

    // The EB data on coarse MG levels are built here, so the EB index space
    // does not limit coarsening.
    int eb_limit_coarsening = false;
    m_coarsening_strategy = CoarseningStrategy::Sigma; // This will fill nodes outside Neumann BC
    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, _factory, eb_limit_coarsening);

    build_eb_data();

    m_sigma_mf.resize(this->m_num_amr_levels);
    m_sigma_edge.resize(this->m_num_amr_levels);
    for (int ilev = 0; ilev < this->m_num_amr_levels; ++ilev) {
        m_sigma_edge[ilev].resize(this->m_num_mg_levels[ilev]);
    }
}

void
MLEBNodeFDLaplacian::build_eb_data ()
{
    BL_PROFILE("MLEBNodeFDLaplacian::build_eb_data()");

    m_levset.resize(m_num_amr_levels);
    m_eb_pos.resize(m_num_amr_levels);
    m_has_eb.resize(m_num_amr_levels);

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        auto const* factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        if (!factory || factory->isAllRegular()) { continue; }

        int nmglevs = m_num_mg_levels[amrlev];
        m_levset[amrlev].resize(nmglevs);
        m_eb_pos[amrlev].resize(nmglevs);
        m_has_eb[amrlev].resize(nmglevs);
        for (int mglev = 0; mglev < nmglevs; ++mglev) {
            BoxArray const& ba = m_grids[amrlev][mglev];
            DistributionMapping const& dm = m_dmap[amrlev][mglev];
            m_levset[amrlev][mglev].define(amrex::convert(ba,IntVect(1)), dm, 1, 1);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_eb_pos[amrlev][mglev][idim].define
                    (amrex::convert(ba,IntVect::TheEdgeVector(idim)), dm, 1, 1);
            }
            m_has_eb[amrlev][mglev].define(ba, dm);
        }

        // Level 0: copy from the factory.
        {
            auto const& levset_f = factory->getLevelSet();
            auto const& edgecent = factory->getEdgeCent();
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(levset_f.nGrow() >= 1 &&
                                             edgecent[0]->nGrow() >= 1,
                "MLEBNodeFDLaplacian: the EB factory needs at least one ghost cell");
            MultiFab::Copy(m_levset[amrlev][0], levset_f, 0, 0, 1, 1);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto& ebp = m_eb_pos[amrlev][0][idim];
                auto const off = IntVect::TheDimensionVector(idim).dim3();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(ebp); mfi.isValid(); ++mfi) {
                    if (edgecent[idim]->ok(mfi)) {
                        Box const& bx = mfi.fabbox();
                        Array4<Real> const& ebpa = ebp.array(mfi);
                        Array4<Real const> const& eca = edgecent[idim]->const_array(mfi);
                        Array4<Real const> const& lsa = levset_f.const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                        {
                            mlebndfdlap_eb_pos_from_cent(i,j,k,ebpa,eca,lsa,off);
                        });
                    } else {
                        ebp[mfi].setVal<RunOn::Device>(Real(1.0)); // regular or covered
                    }
                }
            }
        }

        // Coarse levels: from the next finer level.
        for (int mglev = 1; mglev < nmglevs; ++mglev)
        {
            IntVect const ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[mglev-1];
            Dim3 const rr = ratio.dim3(1);
            auto const& flevset = m_levset[amrlev][mglev-1];
            auto& clevset = m_levset[amrlev][mglev];
            auto const& febp = m_eb_pos[amrlev][mglev-1];
            auto& cebp = m_eb_pos[amrlev][mglev];

            bool const need_parallel_copy = !amrex::isMFIterSafe(clevset, flevset);
            MultiFab clevset_tmp;
            Array<MultiFab,AMREX_SPACEDIM> cebp_tmp;
            MultiFab* pclevset = &clevset;
            Array<MultiFab*,AMREX_SPACEDIM> pcebp = GetArrOfPtrs(cebp);
            if (need_parallel_copy) {
                BoxArray const& cba = amrex::coarsen(m_grids[amrlev][mglev-1], ratio);
                clevset_tmp.define(amrex::convert(cba,IntVect(1)), flevset.DistributionMap(),
                                   1, 0, MFInfo().SetArena(The_Async_Arena()));
                pclevset = &clevset_tmp;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    cebp_tmp[idim].define(amrex::convert(cba,IntVect::TheEdgeVector(idim)),
                                         flevset.DistributionMap(), 1, 0,
                                         MFInfo().SetArena(The_Async_Arena()));
                    pcebp[idim] = &cebp_tmp[idim];
                }
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*pclevset, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& cls = pclevset->array(mfi);
                Array4<Real const> const& fls = flevset.const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    cls(i,j,k) = fls(i*rr.x, j*rr.y, k*rr.z);
                });
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    Box const& ebx = mfi.tilebox(IntVect::TheEdgeVector(idim));
                    Array4<Real> const& cebpa = pcebp[idim]->array(mfi);
                    Array4<Real const> const& febpa = febp[idim].const_array(mfi);
                    auto const off = IntVect::TheDimensionVector(idim).dim3();
                    AMREX_HOST_DEVICE_FOR_3D(ebx, i, j, k,
                    {
                        mlebndfdlap_coarsen_eb_pos(i,j,k,cebpa,febpa,fls,off,rr);
                    });
                }
            }

            if (need_parallel_copy) {
                clevset.ParallelCopy(clevset_tmp);
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    cebp[idim].ParallelCopy(cebp_tmp[idim]);
                }
            }

            fill_domain_ghost(clevset, m_geom[amrlev][mglev], -1);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                fill_domain_ghost(cebp[idim], m_geom[amrlev][mglev], idim);
            }
        }

        // Whether each box, including its ghost nodes, is mixed: it has open
        // nodes and also covered nodes or cut edges.  An all-covered box is
        // not mixed.  All the reductions are launched before the first
        // value() call syncs.
        for (int mglev = 0; mglev < nmglevs; ++mglev) {
            auto const& levset = m_levset[amrlev][mglev];
            auto const& ebp = m_eb_pos[amrlev][mglev];
            using ROps = ReduceOps<ReduceOpLogicalOr, ReduceOpLogicalOr>;
            using RData = ReduceData<int, int>; // open, covered or cut
            using ReduceTuple = RData::Type;
            int const nboxes = levset.local_size();
            Vector<std::unique_ptr<ROps>> rops(nboxes);
            Vector<std::unique_ptr<RData>> rdata(nboxes);
            for (MFIter mfi(levset, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi) {
                int const li = mfi.LocalIndex();
                rops[li] = std::make_unique<ROps>();
                rdata[li] = std::make_unique<RData>(*rops[li]);
                auto const& lsa = levset.const_array(mfi);
                rops[li]->eval(mfi.fabbox(), *rdata[li],
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                {
                    bool const covered = lsa(i,j,k) >= Real(0.0);
                    return { !covered, covered };
                });
                // The level set is injected, so a coarse edge can be cut
                // between two open nodes when only the fine midpoint node
                // is covered.  Such a box has no covered node but needs
                // the EB stencil.
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    auto const& ea = ebp[idim].const_array(mfi);
                    rops[li]->eval(ebp[idim][mfi].box(), *rdata[li],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                    {
                        return { false, ea(i,j,k) < Real(1.0) };
                    });
                }
            }
            for (MFIter mfi(levset, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi) {
                int const li = mfi.LocalIndex();
                auto const r = rdata[li]->value(*rops[li]);
                m_has_eb[amrlev][mglev][mfi] = amrex::get<0>(r) && amrex::get<1>(r);
            }
        }
    }
}

void
MLEBNodeFDLaplacian::limit_coarsening ()
{
    BL_PROFILE("MLEBNodeFDLaplacian::limit_coarsening()");

    int nmglevs = m_num_mg_levels[0];
    if (m_levset.empty() || m_levset[0].empty() || nmglevs <= 1) { return; }

    // Stop coarsening at the last MG level that still has enough unknowns
    // for a meaningful bottom solve and still has covered nodes.
    constexpr int min_open_nodes = 18;
    // Levels with more cells than this are assumed to have enough unknowns
    // and are not examined.
    constexpr Long max_npts_to_check = 65536;

    // nopen and covered share one buffer for a single MPI reduction.
    Vector<int> buf(2*nmglevs, 0);
    int* nopen = buf.data();
    int* covered = buf.data() + nmglevs;

    // Unknowns: nodes that are neither Dirichlet nor covered, each counted
    // once by its owner.  Only the ntest coarsest levels are examined.
    int ntest = 0;
    for (int mglev = nmglevs-1; mglev > 0; --mglev) {
        if (m_grids[0][mglev].numPts() > max_npts_to_check) { break; }
        auto const& dmask = *m_dirichlet_mask[0][mglev];
        auto omask = makeOwnerMask(m_grids[0][mglev], m_dmap[0][mglev], m_geom[0][mglev]);
        auto const& dma = dmask.const_arrays();
        auto const& oma = omask->const_arrays();
        nopen[mglev] = ParReduce(TypeList<ReduceOpSum>{}, TypeList<int>{},
                                 dmask, IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> GpuTuple<int>
        {
            return { (dma[box_no](i,j,k) == 0 && oma[box_no](i,j,k)) ? 1 : 0 };
        });
        ++ntest;
    }

    // A level without covered nodes has lost the EB Dirichlet condition,
    // which can make it singular.  The level set is injected, so a node
    // covered on a level is covered on all finer levels.  Hence the local
    // search can stop at the first level with covered nodes.
    for (int mglev = nmglevs-1; mglev > 0; --mglev) {
        auto const& levset = m_levset[0][mglev];
        auto const& ma = levset.const_arrays();
        covered[mglev] = ParReduce(TypeList<ReduceOpLogicalOr>{}, TypeList<int>{},
                                   levset, IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> GpuTuple<int>
        {
            return { (ma[box_no](i,j,k) >= Real(0.0)) ? 1 : 0 };
        });
        if (covered[mglev]) {
            for (int lev = 1; lev < mglev; ++lev) { covered[lev] = 1; }
            break;
        }
    }

    ParallelAllReduce::Sum(buf.data(), int(buf.size()), ParallelContext::CommunicatorSub());

    int last_good = 0;
    for (int mglev = nmglevs-1; mglev > 0; --mglev) {
        if (mglev < nmglevs-ntest || nopen[mglev] >= min_open_nodes) {
            last_good = mglev;
            break;
        }
    }
    int last_covered = 0;
    for (int mglev = nmglevs-1; mglev > 0; --mglev) {
        if (covered[mglev]) {
            last_covered = mglev;
            break;
        }
    }

    int const new_nmglevs = std::min(last_good, last_covered) + 1;
    if (new_nmglevs < nmglevs) {
        resizeMultiGrid(new_nmglevs);
        nmglevs = m_num_mg_levels[0];
        m_levset[0].resize(nmglevs);
        m_eb_pos[0].resize(nmglevs);
        m_has_eb[0].resize(nmglevs);
        m_sigma_edge[0].resize(nmglevs);
    }
}

#endif

void
MLEBNodeFDLaplacian::define (const Vector<Geometry>& a_geom,
                             const Vector<BoxArray>& a_grids,
                             const Vector<DistributionMapping>& a_dmap,
                             const LPInfo& a_info)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM>1, "MLEBNodeFDLaplacian: 1D not supported");

    BL_PROFILE("MLEBNodeFDLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    if (a_grids.size() > 1) {
        amrex::Abort("MLEBNodeFDLaplacian: multi-level not supported");
    }

    m_coarsening_strategy = CoarseningStrategy::Sigma; // This will fill nodes outside Neumann BC
    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info);

#ifdef AMREX_USE_EB
    // No EB factory here, but the per-level vectors must exist.
    m_levset.resize(this->m_num_amr_levels);
    m_eb_pos.resize(this->m_num_amr_levels);
    m_has_eb.resize(this->m_num_amr_levels);
#endif

    m_sigma_mf.resize(this->m_num_amr_levels);
    m_sigma_edge.resize(this->m_num_amr_levels);
    for (int ilev = 0; ilev < this->m_num_amr_levels; ++ilev) {
        m_sigma_edge[ilev].resize(this->m_num_mg_levels[ilev]);
    }
}

void
MLEBNodeFDLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    IntVect const ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[cmglev-1];
#if (AMREX_SPACEDIM == 1)
    int semicoarsening_dir = 0;
#else
    // Direction NOT coarsened by this MG step. Derived from the level's
    // ratio, because info.semicoarsening_direction is -1 when the direction
    // is chosen automatically.
    int semicoarsening_dir = 2;
    if (ratio[1] == 1) {
        semicoarsening_dir = 1;
    } else if (ratio[0] == 1) {
        semicoarsening_dir = 0;
    }
#endif

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), ratio);
        cfine.define(ba, fine.DistributionMap(), 1, 0, MFInfo().SetArena(The_Async_Arena()));
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][cmglev-1];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> cfab = pcrse->array(mfi);
        Array4<Real const> const& ffab = fine.const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
        if (ratio == 2) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction(i,j,k,cfab,ffab,mfab);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_semi_restriction(i,j,k,cfab,ffab,mfab, semicoarsening_dir);
            });
        }
    }

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }
}

void
MLEBNodeFDLaplacian::interpolation (int amrlev, int fmglev, MultiFab& fine,
                                    const MultiFab& crse) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::interpolation()");

    IntVect const ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[fmglev];
#if (AMREX_SPACEDIM == 1)
    int semicoarsening_dir = 0;
#else
    // Direction NOT coarsened by this MG step. Derived from the level's
    // ratio, because info.semicoarsening_direction is -1 when the direction
    // is chosen automatically.
    int semicoarsening_dir = 2;
    if (ratio[1] == 1) {
        semicoarsening_dir = 1;
    } else if (ratio[0] == 1) {
        semicoarsening_dir = 0;
    }
#endif

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    const MultiFab* cmf = &crse;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), ratio);
        cfine.define(ba, fine.DistributionMap(), 1, 0, MFInfo().SetArena(The_Async_Arena()));
        cfine.ParallelCopy(crse);
        cmf = &cfine;
    }

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][fmglev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& ffab = fine.array(mfi);
        Array4<Real const> const& cfab = cmf->const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
        if (ratio == 2) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndtslap_interpadd(i,j,k,ffab,cfab,mfab);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndtslap_semi_interpadd(i,j,k,ffab,cfab,mfab,semicoarsening_dir);
            });
        }
    }
}

void
MLEBNodeFDLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLEBNodeFDLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();

#ifdef AMREX_USE_EB
    // If neither setEBDirichlet overload was called, m_s_phi_eb still holds
    // the "use the m_phi_eb array" sentinel while m_phi_eb is empty. Default
    // to homogeneous Dirichlet on the EB instead.
    if (m_s_phi_eb == std::numeric_limits<Real>::lowest() && m_phi_eb.empty()) {
        m_s_phi_eb = Real(0.0);
    }

    // Set covered nodes to Dirichlet, but with a negative value.
    // compGrad relies on the negative value to detect EB.
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        if (m_levset[amrlev].empty()) { continue; }
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            auto const& levset_ar = m_levset[amrlev][mglev].const_arrays();
            auto& dmask_mf = *m_dirichlet_mask[amrlev][mglev];
            auto const& dmask_ar = dmask_mf.arrays();
            amrex::ParallelFor(dmask_mf,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                if (levset_ar[box_no](i,j,k) >= Real(0.0)) {
                    dmask_ar[box_no](i,j,k) = -1;
                }
            });
        }
    }

    limit_coarsening();
#endif

    {
        int amrlev = 0;
        int mglev = m_num_mg_levels[amrlev]-1;
        auto const& dotmasks = m_bottom_dot_mask.arrays();
        auto const& dirmasks = m_dirichlet_mask[amrlev][mglev]->const_arrays();
        amrex::ParallelFor(m_bottom_dot_mask,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            if (dirmasks[box_no](i,j,k)) {
                dotmasks[box_no](i,j,k) = Real(0.);
            }
        });
    }

    AMREX_ASSERT(!isBottomSingular());

    Gpu::streamSynchronize();

#if (AMREX_SPACEDIM == 2)
    if (m_rz) {
        if (m_geom[0][0].ProbLo(0) == 0._rt) {
            // With the alpha/r^2 term the solution must vanish on the axis, so
            // the caller has to declare Dirichlet there. That in turn lets
            // buildMasks mark the axis nodes, which keeps the residual,
            // restriction, interpolation and dot products consistent with what
            // the operator kernels do.
            if (m_rz_alpha != 0._rt) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_lobc[0][0] == BCType::Dirichlet,
                    "The lo-x BC must be Dirichlet for 2d RZ with a non-zero alpha");
            } else {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_lobc[0][0] == BCType::Neumann,
                                                 "The lo-x BC must be Neumann for 2d RZ");
            }
        }
        if (m_sigma[0] == 0._rt) {
            m_sigma[0] = 1._rt; // For backward compatibility
        }
    }
#endif

    if (m_has_sigma_mf) {
        update_sigma();
    }
}

#ifdef AMREX_USE_EB
bool
MLEBNodeFDLaplacian::scaleRHS (int amrlev, MultiFab* rhs) const
{
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());

    if (!factory) {return false; }

    if (rhs && !m_levset[amrlev].empty()) {
        auto const& dmask = *m_dirichlet_mask[amrlev][0];
        auto const& ebp = m_eb_pos[amrlev][0];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            if (m_has_eb[amrlev][0][mfi]) {
                const Box& box = mfi.tilebox();
                Array4<Real> const& rhsarr = rhs->array(mfi);
                Array4<int const> const& dmarr = dmask.const_array(mfi);
                AMREX_D_TERM(Array4<Real const> const& ebpx = ebp[0].const_array(mfi);,
                             Array4<Real const> const& ebpy = ebp[1].const_array(mfi);,
                             Array4<Real const> const& ebpz = ebp[2].const_array(mfi));
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_scale_rhs(i,j,k,rhsarr,dmarr,AMREX_D_DECL(ebpx,ebpy,ebpz));
                });
            }
        }
    }

    return true;
}
#endif

void
MLEBNodeFDLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::Fapply()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    bool rz = false;
    Real sig0 = Real(1.0), dx0 = Real(0.0), dx1 = Real(0.0), xlo = Real(0.0), alpha = Real(0.0);
#if (AMREX_SPACEDIM == 2)
    rz = m_rz;
    sig0 = m_sigma[0];
    dx0 = m_geom[amrlev][mglev].CellSize(0);
    dx1 = m_geom[amrlev][mglev].CellSize(1)/std::sqrt(m_sigma[1]);
    xlo = m_geom[amrlev][mglev].ProbLo(0);
    alpha = m_rz_alpha;
#endif
    GpuArray<Real,AMREX_SPACEDIM> const b
        {AMREX_D_DECL(m_sigma[0]*dxinv[0]*dxinv[0],
                      m_sigma[1]*dxinv[1]*dxinv[1],
                      m_sigma[2]*dxinv[2]*dxinv[2])};

    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];

    Real phieb = std::numeric_limits<Real>::lowest();
#ifdef AMREX_USE_EB
    phieb = (m_in_solution_mode && !this->m_precond_mode) ? m_s_phi_eb : Real(0.0);
    bool const has_eb_level = !m_levset[amrlev].empty();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        Array4<Real const> const& xarr = in.const_array(mfi);
        Array4<Real> const& yarr = out.array(mfi);
        Array4<int const> const& dmarr = dmask.const_array(mfi);

        bool has_eb = false;
        Array4<Real const> levset;
        GpuArray<Array4<Real const>,AMREX_SPACEDIM> ebp;
        Array4<Real const> phiebarr;
#ifdef AMREX_USE_EB
        if (has_eb_level && m_has_eb[amrlev][mglev][mfi]) {
            has_eb = true;
            levset = m_levset[amrlev][mglev].const_array(mfi);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                ebp[idim] = m_eb_pos[amrlev][mglev][idim].const_array(mfi);
            }
            if (phieb == std::numeric_limits<Real>::lowest()) {
                phiebarr = m_phi_eb[amrlev].const_array(mfi);
            }
        }
#endif

        if (m_has_sigma_mf) {
            EBNodeFDEdgeSigma sig;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                sig.s[idim] = m_sigma_edge[amrlev][mglev][idim].const_array(mfi);
            }
            if (phiebarr) {
                fapply_box(box, yarr, xarr, dmarr, sig, has_eb, levset, ebp, phiebarr,
                           b, rz, dx0, dx1, xlo, alpha);
            } else {
                fapply_box(box, yarr, xarr, dmarr, sig, has_eb, levset, ebp, phieb,
                           b, rz, dx0, dx1, xlo, alpha);
            }
        } else if (rz) {
            EBNodeFDRZConstSigma const sig{sig0};
            if (phiebarr) {
                fapply_box(box, yarr, xarr, dmarr, sig, has_eb, levset, ebp, phiebarr,
                           b, rz, dx0, dx1, xlo, alpha);
            } else {
                fapply_box(box, yarr, xarr, dmarr, sig, has_eb, levset, ebp, phieb,
                           b, rz, dx0, dx1, xlo, alpha);
            }
        } else {
            EBNodeFDConstSigma const sig{};
            if (phiebarr) {
                fapply_box(box, yarr, xarr, dmarr, sig, has_eb, levset, ebp, phiebarr,
                           b, rz, dx0, dx1, xlo, alpha);
            } else {
                fapply_box(box, yarr, xarr, dmarr, sig, has_eb, levset, ebp, phieb,
                           b, rz, dx0, dx1, xlo, alpha);
            }
        }
    }
}

void
MLEBNodeFDLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::Fsmooth()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    bool rz = false;
    Real sig0 = Real(1.0), dx0 = Real(0.0), dx1 = Real(0.0), xlo = Real(0.0), alpha = Real(0.0);
#if (AMREX_SPACEDIM == 2)
    rz = m_rz;
    sig0 = m_sigma[0];
    dx0 = m_geom[amrlev][mglev].CellSize(0);
    dx1 = m_geom[amrlev][mglev].CellSize(1)/std::sqrt(m_sigma[1]);
    xlo = m_geom[amrlev][mglev].ProbLo(0);
    alpha = m_rz_alpha;
#endif
    GpuArray<Real,AMREX_SPACEDIM> const b
        {AMREX_D_DECL(m_sigma[0]*dxinv[0]*dxinv[0],
                      m_sigma[1]*dxinv[1]*dxinv[1],
                      m_sigma[2]*dxinv[2]*dxinv[2])};

    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];

#ifdef AMREX_USE_EB
    bool const has_eb_level = !m_levset[amrlev].empty();
#endif

    for (int redblack = 0; redblack < 2; ++redblack) {
        if (redblack > 0) {
            applyBC(amrlev, mglev, sol, BCMode::Homogeneous, StateMode::Correction);
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sol,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            Array4<Real> const& solarr = sol.array(mfi);
            Array4<Real const> const& rhsarr = rhs.const_array(mfi);
            Array4<int const> const& dmskarr = dmask.const_array(mfi);

            bool has_eb = false;
            Array4<Real const> levset;
            GpuArray<Array4<Real const>,AMREX_SPACEDIM> ebp;
#ifdef AMREX_USE_EB
            if (has_eb_level && m_has_eb[amrlev][mglev][mfi]) {
                has_eb = true;
                levset = m_levset[amrlev][mglev].const_array(mfi);
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    ebp[idim] = m_eb_pos[amrlev][mglev][idim].const_array(mfi);
                }
            }
#endif

            if (m_has_sigma_mf) {
                EBNodeFDEdgeSigma sig;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    sig.s[idim] = m_sigma_edge[amrlev][mglev][idim].const_array(mfi);
                }
                fsmooth_box(box, solarr, rhsarr, dmskarr, sig, has_eb, levset, ebp,
                            b, rz, dx0, dx1, xlo, alpha, redblack);
            } else if (rz) {
                EBNodeFDRZConstSigma const sig{sig0};
                fsmooth_box(box, solarr, rhsarr, dmskarr, sig, has_eb, levset, ebp,
                            b, rz, dx0, dx1, xlo, alpha, redblack);
            } else {
                EBNodeFDConstSigma const sig{};
                fsmooth_box(box, solarr, rhsarr, dmskarr, sig, has_eb, levset, ebp,
                            b, rz, dx0, dx1, xlo, alpha, redblack);
            }
        }
    }

    nodalSync(amrlev, mglev, sol);
}

void
MLEBNodeFDLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    amrex::ignore_unused(amrlev, mglev, mf);
}

void
MLEBNodeFDLaplacian::fixUpResidualMask (int /*amrlev*/, iMultiFab& /*resmsk*/)
{
    amrex::Abort("MLEBNodeFDLaplacian::fixUpResidualMask: TODO");
}

void
MLEBNodeFDLaplacian::compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                               MultiFab& sol, Location /*loc*/) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::compGrad()");
    if (amrex::isMFIterSafe(*grad[0], sol)) {
        this->compGrad_doit(amrlev, grad, sol);
    } else {
        Array<MultiFab,AMREX_SPACEDIM> grad_tmp;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            grad_tmp[idim].define(amrex::convert(sol.boxArray(), grad[idim]->ixType()),
                                  sol.DistributionMap(), 1, 0,
                                  MFInfo{}.SetArena(The_Async_Arena()));
        }
        this->compGrad_doit(amrlev, GetArrOfPtrs(grad_tmp), sol);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            grad[idim]->ParallelCopy(grad_tmp[idim], 0, 0, 1);
        }
    }
}

void
MLEBNodeFDLaplacian::compGrad_doit (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                                    MultiFab& sol) const
{
    AMREX_ASSERT(AMREX_D_TERM(grad[0]->ixType() == IndexType(IntVect::TheEdgeVector(0)),
                           && grad[1]->ixType() == IndexType(IntVect::TheEdgeVector(1)),
                           && grad[2]->ixType() == IndexType(IntVect::TheEdgeVector(2))));

    const int mglev = 0;
    AMREX_D_TERM(const auto dxi = m_geom[amrlev][mglev].InvCellSize(0);,
                 const auto dyi = m_geom[amrlev][mglev].InvCellSize(1);,
                 const auto dzi = m_geom[amrlev][mglev].InvCellSize(2);)

#ifdef AMREX_USE_EB
    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];
    const auto phieb = m_s_phi_eb;
    bool const has_eb = !m_levset[amrlev].empty();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*grad[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM(const Box& xbox = mfi.tilebox(IntVect::TheEdgeVector(0));,
                     const Box& ybox = mfi.tilebox(IntVect::TheEdgeVector(1));,
                     const Box& zbox = mfi.tilebox(IntVect::TheEdgeVector(2));)
        Array4<Real const> const& p = sol.const_array(mfi);
        AMREX_D_TERM(Array4<Real> const& gpx = grad[0]->array(mfi);,
                     Array4<Real> const& gpy = grad[1]->array(mfi);,
                     Array4<Real> const& gpz = grad[2]->array(mfi);)
#ifdef AMREX_USE_EB
        if (has_eb) {
            Array4<int const> const& dmarr = dmask.const_array(mfi);
            AMREX_D_TERM(Array4<Real const> const& ebpx = m_eb_pos[amrlev][mglev][0].const_array(mfi);,
                         Array4<Real const> const& ebpy = m_eb_pos[amrlev][mglev][1].const_array(mfi);,
                         Array4<Real const> const& ebpz = m_eb_pos[amrlev][mglev][2].const_array(mfi);)
            if (phieb == std::numeric_limits<Real>::lowest()) {
                auto const& phiebarr = m_phi_eb[amrlev].const_array(mfi);
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(
                    xbox, txbox,
                    {
                        mlebndfdlap_grad_x(txbox, gpx, p, dmarr, ebpx, phiebarr, dxi);
                    }
                    , ybox, tybox,
                    {
                        mlebndfdlap_grad_y(tybox, gpy, p, dmarr, ebpy, phiebarr, dyi);
                    }
                    , zbox, tzbox,
                    {
                        mlebndfdlap_grad_z(tzbox, gpz, p, dmarr, ebpz, phiebarr, dzi);
                    });
            } else {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(
                    xbox, txbox,
                    {
                        mlebndfdlap_grad_x(txbox, gpx, p, dmarr, ebpx, phieb, dxi);
                    }
                    , ybox, tybox,
                    {
                        mlebndfdlap_grad_y(tybox, gpy, p, dmarr, ebpy, phieb, dyi);
                    }
                    , zbox, tzbox,
                    {
                        mlebndfdlap_grad_z(tzbox, gpz, p, dmarr, ebpz, phieb, dzi);
                    });
            }
        } else
#endif
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(
                    xbox, txbox,
                    {
                        mlebndfdlap_grad_x(txbox, gpx, p, dxi);
                    }
                    , ybox, tybox,
                    {
                        mlebndfdlap_grad_y(tybox, gpy, p, dyi);
                    }
                    , zbox, tzbox,
                    {
                        mlebndfdlap_grad_z(tzbox, gpz, p, dzi);
                    });
        }
    }
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
void
MLEBNodeFDLaplacian::fillIJMatrix (MFIter const& mfi,
                                   Array4<HypreNodeLap::AtomicInt const> const& gid,
                                   Array4<int const> const& lid,
                                   HypreNodeLap::Int* ncols,
                                   HypreNodeLap::Int* cols,
                                   Real* mat) const
{
    const int amrlev = 0;
    const int mglev  = m_num_mg_levels[amrlev]-1;

#if (AMREX_SPACEDIM == 2)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_rz,
        "MLEBNodeFDLaplacian::fillIJMatrix: RZ is not supported yet");
#endif

    const Geometry& geom = m_geom[amrlev][mglev];
    const auto dxinv = geom.InvCellSizeArray();
    const GpuArray<Real,AMREX_SPACEDIM> bcoef
        {AMREX_D_DECL(m_sigma[0]*dxinv[0]*dxinv[0],
                      m_sigma[1]*dxinv[1]*dxinv[1],
                      m_sigma[2]*dxinv[2]*dxinv[2])};

    const Box& nddom = amrex::surroundingNodes(geom.Domain());
    const auto ndlo = amrex::lbound(nddom);
    const auto ndhi = amrex::ubound(nddom);

    const auto lobc = LoBC();
    const auto hibc = HiBC();
    GpuArray<bool,AMREX_SPACEDIM> reflect_lo{};
    GpuArray<bool,AMREX_SPACEDIM> reflect_hi{};
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        reflect_lo[idim] = (lobc[idim] == LinOpBCType::Neumann ||
                            lobc[idim] == LinOpBCType::inflow);
        reflect_hi[idim] = (hibc[idim] == LinOpBCType::Neumann ||
                            hibc[idim] == LinOpBCType::inflow);
    }

    const bool has_sig = m_has_sigma_mf;
    GpuArray<Array4<Real const>,AMREX_SPACEDIM> sig{};
    if (has_sig) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            sig[idim] = m_sigma_edge[amrlev][mglev][idim].const_array(mfi);
        }
    }

    bool has_eb = false;
    Array4<Real const> levset{};
    GpuArray<Array4<Real const>,AMREX_SPACEDIM> ebp{};
#ifdef AMREX_USE_EB
    if (!m_levset[amrlev].empty() && m_has_eb[amrlev][mglev][mfi]) {
        has_eb = true;
        levset = m_levset[amrlev][mglev].const_array(mfi);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            ebp[idim] = m_eb_pos[amrlev][mglev][idim].const_array(mfi);
        }
    }
#endif

    const Box& ndbx = mfi.validbox();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE
        (ndbx.numPts()*(2*AMREX_SPACEDIM+1) <
         static_cast<Long>(std::numeric_limits<int>::max()),
         "The Box is too big.  We could use Long here, but it would be much slower.");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        const auto blo = amrex::lbound(ndbx);
        const auto blen = amrex::length(ndbx);
        const int npts = static_cast<int>(ndbx.numPts());
        Scan::PrefixSum<int>
            (npts,
             [=] AMREX_GPU_DEVICE (int offset) noexcept -> int
             {
                 int const k = offset / (blen.x*blen.y);
                 int const j = (offset - k*blen.x*blen.y) / blen.x;
                 int const i = offset - k*blen.x*blen.y - j*blen.x;
                 if (lid(i+blo.x,j+blo.y,k+blo.z) < 0) { return 0; }
                 return mlebndfdlap_ijmat_row(i+blo.x, j+blo.y, k+blo.z, gid,
                                              has_eb, has_sig, bcoef, sig,
                                              levset, ebp, ndlo, ndhi,
                                              reflect_lo, reflect_hi).n;
             },
             [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
             {
                 int const k = offset / (blen.x*blen.y);
                 int const j = (offset - k*blen.x*blen.y) / blen.x;
                 int const i = offset - k*blen.x*blen.y - j*blen.x;
                 int const row_lid = lid(i+blo.x,j+blo.y,k+blo.z);
                 if (row_lid < 0) { return; }
                 auto const& row = mlebndfdlap_ijmat_row
                     (i+blo.x, j+blo.y, k+blo.z, gid, has_eb, has_sig, bcoef,
                      sig, levset, ebp, ndlo, ndhi, reflect_lo, reflect_hi);
                 ncols[row_lid] = row.n;
                 for (int n = 0; n < row.n; ++n) {
                     cols[ps+n] = static_cast<HypreNodeLap::Int>
                         (gid(row.node[n].x, row.node[n].y, row.node[n].z));
                     mat[ps+n] = row.val[n];
                 }
             },
             Scan::Type::exclusive);
    } else
#endif
    {
        // The nodes are visited in the same order in which local ids were
        // assigned, so the rows come out sorted by local id.
        int nelems = 0;
        amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
        {
            if (lid(i,j,k) >= 0) {
                auto const& row = mlebndfdlap_ijmat_row(i, j, k, gid, has_eb, has_sig,
                                                        bcoef, sig, levset, ebp,
                                                        ndlo, ndhi, reflect_lo, reflect_hi);
                ncols[lid(i,j,k)] = row.n;
                for (int n = 0; n < row.n; ++n) {
                    cols[nelems] = static_cast<HypreNodeLap::Int>
                        (gid(row.node[n].x, row.node[n].y, row.node[n].z));
                    mat[nelems] = row.val[n];
                    ++nelems;
                }
            }
        });
    }
}

void
MLEBNodeFDLaplacian::fillRHS (MFIter const& mfi, Array4<int const> const& lid,
                              Real* rhs, Array4<Real const> const& bfab) const
{
    // Unlike MLNodeLaplacian, this is a finite-difference operator, so nodes on
    // a Neumann boundary need no volume factor here.  fillIJMatrix folds the
    // ghost node onto its mirror image instead.
    const Box& bx = mfi.validbox();
    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
    {
        if (lid(i,j,k) >= 0) {
            rhs[lid(i,j,k)] = bfab(i,j,k);
        }
    });
}
#endif

void
MLEBNodeFDLaplacian::postSolve (Vector<MultiFab*> const& sol) const
{
#ifdef AMREX_USE_EB
    if (this->m_precond_mode) { return; }
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        const auto phieb = m_s_phi_eb;
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        if (!factory || factory->isAllRegular()) { return; }
        auto const& levset_mf = factory->getLevelSet();
        auto const& levset_ar = levset_mf.const_arrays();
        MultiFab& mf = *sol[amrlev];
        auto const& sol_ar = mf.arrays();
        if (phieb == std::numeric_limits<Real>::lowest()) {
            auto const& phieb_ar = m_phi_eb[amrlev].const_arrays();
            amrex::ParallelFor(mf, IntVect(1),
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
            {
                if (levset_ar[bi](i,j,k) >= Real(0.0)) {
                    sol_ar[bi](i,j,k) = phieb_ar[bi](i,j,k);
                }
            });
        } else {
            amrex::ParallelFor(mf, IntVect(1),
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
            {
                if (levset_ar[bi](i,j,k) >= Real(0.0)) {
                    sol_ar[bi](i,j,k) = phieb;
                }
            });
        }
    }
#else
    amrex::ignore_unused(sol);
#endif
}

void
MLEBNodeFDLaplacian::update ()
{
    if (MLNodeLinOp::needsUpdate()) {
        MLNodeLinOp::update();
    }

    if (m_needs_update && m_has_sigma_mf) {
        update_sigma();
    }
    m_needs_update = false;
}

void
MLEBNodeFDLaplacian::update_sigma ()
{
    BL_PROFILE("MLEBNodeFDLaplacian::update_sigma()");

    AMREX_D_TERM(m_sigma[0] = Real(1.0);,
                 m_sigma[1] = Real(1.0);,
                 m_sigma[2] = Real(1.0));
    AMREX_ALWAYS_ASSERT(this->m_num_amr_levels == 1);
    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 0; mglev < this->m_num_mg_levels[amrlev]; ++mglev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto& se = m_sigma_edge[amrlev][mglev][idim];
                if (se.empty()) {
                    se.define(amrex::convert(this->m_grids[amrlev][mglev],
                                             IntVect::TheEdgeVector(idim)),
                              this->m_dmap[amrlev][mglev], 1, 1);
                }
            }
        }

        // Level 0: cell-centered to edge-centered
        {
            auto const& geom = this->m_geom[amrlev][0];
            auto& sigma = *m_sigma_mf[amrlev];
            sigma.FillBoundary(geom.periodicity());

            const Box& domain = geom.Domain();
            const auto lobc = LoBC();
            const auto hibc = HiBC();

            MFItInfo mfi_info;
            if (Gpu::notInLaunchRegion()) { mfi_info.SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(sigma, mfi_info); mfi.isValid(); ++mfi)
            {
                Array4<Real> const& sfab = sigma.array(mfi);
                mlndlap_fillbc_cc<Real>(mfi.validbox(),sfab,domain,lobc,hibc);
            }

            bool const rz = m_rz;
            Real const dr = geom.CellSize(0);
            Real const rlo = geom.ProbLo(0);
            MultiFab const* vfrac = nullptr;
#ifdef AMREX_USE_EB
            auto const* factory = dynamic_cast<EBFArrayBoxFactory const*>
                (m_factory[amrlev][0].get());
            if (factory && !m_levset[amrlev].empty()) {
                vfrac = &(factory->getVolFrac());
            }
#endif
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto& se = m_sigma_edge[amrlev][0][idim];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(se, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box const& bx = mfi.tilebox();
                    Array4<Real> const& sea = se.array(mfi);
                    Array4<Real const> const& sca = sigma.const_array(mfi);
                    Array4<Real const> const vfa = vfrac ? vfrac->const_array(mfi)
                                                         : Array4<Real const>{};
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlebndfdlap_cc_to_edge_sigma(i,j,k,sea,sca,vfa,idim,rz,dr,rlo);
                    });
                }
                fill_domain_ghost(se, geom, -1);
            }
        }

        // Coarse levels
        for (int mglev = 1; mglev < this->m_num_mg_levels[amrlev]; ++mglev)
        {
            IntVect const ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[mglev-1];
            Dim3 const rr = ratio.dim3(1);
            auto const& fse = m_sigma_edge[amrlev][mglev-1];
            auto& cse = m_sigma_edge[amrlev][mglev];

            bool const need_parallel_copy = !amrex::isMFIterSafe(cse[0], fse[0]);
            Array<MultiFab,AMREX_SPACEDIM> cse_tmp;
            Array<MultiFab*,AMREX_SPACEDIM> pcse = GetArrOfPtrs(cse);
            if (need_parallel_copy) {
                BoxArray const& cba = amrex::coarsen(m_grids[amrlev][mglev-1], ratio);
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    cse_tmp[idim].define(amrex::convert(cba,IntVect::TheEdgeVector(idim)),
                                         fse[idim].DistributionMap(), 1, 0,
                                         MFInfo().SetArena(The_Async_Arena()));
                    pcse[idim] = &cse_tmp[idim];
                }
            }

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto const off = IntVect::TheDimensionVector(idim).dim3();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*pcse[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box const& bx = mfi.tilebox();
                    Array4<Real> const& csa = pcse[idim]->array(mfi);
                    Array4<Real const> const& fsa = fse[idim].const_array(mfi);
                    Array4<Real const> febpa;
                    Array4<Real const> flsa;
#ifdef AMREX_USE_EB
                    if (!m_levset[amrlev].empty()) {
                        febpa = m_eb_pos[amrlev][mglev-1][idim].const_array(mfi);
                        flsa = m_levset[amrlev][mglev-1].const_array(mfi);
                    }
#endif
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlebndfdlap_coarsen_edge_sigma(i,j,k,csa,fsa,febpa,flsa,off,rr);
                    });
                }
                if (need_parallel_copy) {
                    cse[idim].ParallelCopy(cse_tmp[idim]);
                }
                fill_domain_ghost(cse[idim], m_geom[amrlev][mglev], -1);
            }
        }
    }
}

namespace {
    struct LPBase
    {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real xdoty (IntVect const& iv, int, Real vx, Real vy) const
        {
            return dotmsk(iv)*vx*vy;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void normalize (IntVect const&, int, Real&) const {}

        // The neighbor index returned by lowerNeighbor/upperNeighbor is for
        // the solution data only, which lives in a Box without ghost cells.
        // EB data (level set and EB positions) does have ghost cells.
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        int lowerNeighbor (int i, int idim) const
        {
            // There is only a single Box. Thus a box boundary node is
            // either Dirichlet (including domain Dirichlet boundary or
            // coarse/fine Dirichlet boundary) or on the periodic or Neumann
            // domain boundary. Because this is called on non-Dirichlet
            // nodes only, if a boundary node (e.g., i == dlo[idim]) gets
            // here, it's either periodic or Neumann.
            return (i == dlo[idim])
                ? (is_periodic[idim] ? dhi[idim]-1 : i+1)
                : i-1;
        }

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        int upperNeighbor (int i, int idim) const
        {
            return (i == dhi[idim])
                ? (is_periodic[idim] ? dlo[idim]+1 : i-1)
                : i+1;
        }

        Array4<Real const> dotmsk;
        Array4<int const> dirmsk;
        GpuArray<Real,AMREX_SPACEDIM> beta;
        GpuArray<int,AMREX_SPACEDIM> dlo, dhi;
        GpuArray<bool,AMREX_SPACEDIM> is_periodic;
    };

    template <typename S>
    struct LP
        : public LPBase
    {
        LP (LPBase const& a_lpbase, S const& a_sigma)
            : LPBase(a_lpbase), sigma(a_sigma)
            {}

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real apply (IntVect const& iv, int, Array4<Real> const& xa, int n) const
        {
            int const i = iv[0];
            int const j = iv[1];
#if (AMREX_SPACEDIM == 3)
            int const k = iv[2];
#else
            int const k = 0;
#endif
            if (dirmsk(i,j,k)) {
                return Real(0.0);
            }

            Real const xc = xa(i,j,k,n);
            int const im = lowerNeighbor(i,0);
            int const ip = upperNeighbor(i,0);
            int const jm = lowerNeighbor(j,1);
            int const jp = upperNeighbor(j,1);
            Real y = beta[0] * (sigma.x(i-1,j,k)*(xa(im,j,k,n)-xc) + sigma.x(i,j,k)*(xa(ip,j,k,n)-xc))
                +    beta[1] * (sigma.y(i,j-1,k)*(xa(i,jm,k,n)-xc) + sigma.y(i,j,k)*(xa(i,jp,k,n)-xc));
#if (AMREX_SPACEDIM == 3)
            int const km = lowerNeighbor(k,2);
            int const kp = upperNeighbor(k,2);
            y += beta[2] * (sigma.z(i,j,k-1)*(xa(i,j,km,n)-xc) + sigma.z(i,j,k)*(xa(i,j,kp,n)-xc));
#endif
            return y;
        }

        S sigma;
    };

#if (AMREX_SPACEDIM == 2)
    template <bool UseEB>
    struct RZEBData {};

#ifdef AMREX_USE_EB
    template <>
    struct RZEBData<true>
    {
        Array4<Real const> levset;
        GpuArray<Array4<Real const>,AMREX_SPACEDIM> ebp;
    };
#endif

    template <bool UseEB, typename S>
    struct LPRZ
        : public LPBase, public RZEBData<UseEB>
    {
        LPRZ (LPBase const& a_lpbase, S const& a_sigma,
              RZEBData<UseEB> const& a_eb_data, Real a_dr, Real a_dz,
              Real a_rlo, Real a_alpha)
            : LPBase(a_lpbase), RZEBData<UseEB>(a_eb_data), sigma(a_sigma),
              dr(a_dr), dz(a_dz), rlo(a_rlo), alpha(a_alpha)
            {}

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real apply (IntVect const& iv, int, Array4<Real> const& xa, int n) const
        {
            int const i = iv[0];
            int const j = iv[1];
            int const k = 0;
            Real const r = rlo + Real(i)*dr;
            if (dirmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
                return Real(0.0);
            }

            int const im = lowerNeighbor(i,0);
            int const ip = upperNeighbor(i,0);
            int const jm = lowerNeighbor(j,1);
            int const jp = upperNeighbor(j,1);
            Real const xc = xa(i,j,k,n);

            Real hpz = Real(1.0);
            Real hmz = Real(1.0);
            Real hp = Real(1.0);
            Real hm = Real(1.0);
            bool open_pz = true, open_mz = true, open_p = true, open_m = true;
            if constexpr (UseEB) {
                hpz = this->ebp[1](i,j  ,k);
                hmz = mlebndfdlap_hm(this->ebp[1](i,j-1,k));
                hp  = this->ebp[0](i  ,j,k);
                hm  = mlebndfdlap_hm(this->ebp[0](i-1,j,k));
                open_pz = this->levset(i,j+1,k) < Real(0.0) && hpz == Real(1.0);
                open_mz = this->levset(i,j-1,k) < Real(0.0) && hmz == Real(1.0);
                open_p  = this->levset(i+1,j,k) < Real(0.0) && hp  == Real(1.0);
                open_m  = this->levset(i-1,j,k) < Real(0.0) && hm  == Real(1.0);
            }

            Real out;
            Real scale;

            if (r == Real(0.0)) {
                Real const sigp = sigma.x(i,j,k);
                if (open_p) {
                    scale = amrex::min(Real(1.0),hmz,hpz);
                    out = scale*Real(4.0)*sigp*(xa(ip,j,k,n)-xc)/(dr*dr);
                } else {
                    scale = amrex::min(hp,hmz,hpz);
                    out = -Real(4.0)*sigp*mlebndfdlap_scaled_h2inv(scale,hp)*xc/(dr*dr);
                }
            } else {
                scale = amrex::min(hm,hp,hmz,hpz);
                Real const sigp = sigma.x(i  ,j,k);
                Real const sigm = sigma.x(i-1,j,k);
                Real tmp = open_p
                    ? scale*sigp*(xa(ip,j,k,n)-xc)*(r+Real(0.5)*dr)
                    : -sigp*mlebndfdlap_scaled_hinv(scale,hp)*xc*(r+Real(0.5)*hp*dr);
                tmp += open_m
                    ? scale*sigm*(xa(im,j,k,n)-xc)*(r-Real(0.5)*dr)
                    : -sigm*mlebndfdlap_scaled_hinv(scale,hm)*xc*(r-Real(0.5)*hm*dr);
                out = tmp*Real(2.0)/((hp+hm)*r*dr*dr);
            }

            Real const sigp = sigma.y(i,j  ,k);
            Real const sigm = sigma.y(i,j-1,k);
            Real tmp = open_pz
                ? scale*sigp*(xa(i,jp,k,n)-xc)
                : -sigp*mlebndfdlap_scaled_hinv(scale,hpz)*xc;
            tmp += open_mz
                ? scale*sigm*(xa(i,jm,k,n)-xc)
                : -sigm*mlebndfdlap_scaled_hinv(scale,hmz)*xc;
            out += tmp*Real(2.0)/((hpz+hmz)*dz*dz);

            if (r != Real(0.0)) {
                out -= scale*alpha*xc/(r*r);
            }
            return out;
        }

        S sigma;
        Real dr;
        Real dz;
        Real rlo;
        Real alpha;
    };
#endif

#ifdef AMREX_USE_EB
    template <typename S>
    struct LPEB
        : public LPBase
    {
        LPEB (LPBase const& a_lpbase, Array4<Real const> const& a_levset,
              GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& a_el,
              S const& a_sigma)
            : LPBase(a_lpbase), sigma(a_sigma), levset(a_levset), ebp(a_el)
            {}

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real apply (IntVect const& iv, int, Array4<Real> const& xa, int n) const
        {
            int const i = iv[0];
            int const j = iv[1];
#if (AMREX_SPACEDIM == 3)
            int const k = iv[2];
#else
            int const k = 0;
#endif
            if (dirmsk(i,j,k)) {
                return Real(0.0);
            }

            Real const xc = xa(i,j,k,n);
            Real const hpx = ebp[0](i  ,j  ,k  );
            Real const hmx = mlebndfdlap_hm(ebp[0](i-1,j  ,k  ));
            Real const hpy = ebp[1](i  ,j  ,k  );
            Real const hmy = mlebndfdlap_hm(ebp[1](i  ,j-1,k  ));
#if (AMREX_SPACEDIM == 3)
            Real const hpz = ebp[2](i  ,j  ,k  );
            Real const hmz = mlebndfdlap_hm(ebp[2](i  ,j  ,k-1));
            Real const scale = amrex::min(hmx,hpx,hmy,hpy,hmz,hpz);
#else
            Real const scale = amrex::min(hmx,hpx,hmy,hpy);
#endif

            int const im = lowerNeighbor(i,0);
            int const ip = upperNeighbor(i,0);
            Real const sigxp = sigma.x(i  ,j,k);
            Real const sigxm = sigma.x(i-1,j,k);
            Real tmp = (levset(i+1,j,k) < Real(0.0) && hpx == Real(1.0))
                ? sigxp*scale*(xa(ip,j,k,n)-xc)
                : -sigxp*mlebndfdlap_scaled_hinv(scale,hpx)*xc;
            tmp += (levset(i-1,j,k) < Real(0.0) && hmx == Real(1.0))
                ? sigxm*scale*(xa(im,j,k,n)-xc)
                : -sigxm*mlebndfdlap_scaled_hinv(scale,hmx)*xc;
            Real y = beta[0]*tmp*Real(2.0)/(hpx+hmx);

            int const jm = lowerNeighbor(j,1);
            int const jp = upperNeighbor(j,1);
            Real const sigyp = sigma.y(i,j  ,k);
            Real const sigym = sigma.y(i,j-1,k);
            tmp = (levset(i,j+1,k) < Real(0.0) && hpy == Real(1.0))
                ? sigyp*scale*(xa(i,jp,k,n)-xc)
                : -sigyp*mlebndfdlap_scaled_hinv(scale,hpy)*xc;
            tmp += (levset(i,j-1,k) < Real(0.0) && hmy == Real(1.0))
                ? sigym*scale*(xa(i,jm,k,n)-xc)
                : -sigym*mlebndfdlap_scaled_hinv(scale,hmy)*xc;
            y += beta[1]*tmp*Real(2.0)/(hpy+hmy);

#if (AMREX_SPACEDIM == 3)
            int const km = lowerNeighbor(k,2);
            int const kp = upperNeighbor(k,2);
            Real const sigzp = sigma.z(i,j,k  );
            Real const sigzm = sigma.z(i,j,k-1);
            tmp = (levset(i,j,k+1) < Real(0.0) && hpz == Real(1.0))
                ? sigzp*scale*(xa(i,j,kp,n)-xc)
                : -sigzp*mlebndfdlap_scaled_hinv(scale,hpz)*xc;
            tmp += (levset(i,j,k-1) < Real(0.0) && hmz == Real(1.0))
                ? sigzm*scale*(xa(i,j,km,n)-xc)
                : -sigzm*mlebndfdlap_scaled_hinv(scale,hmz)*xc;
            y += beta[2]*tmp*Real(2.0)/(hpz+hmz);
#endif

            return y;
        }

        S sigma;
        Array4<Real const> levset;
        GpuArray<Array4<Real const>,AMREX_SPACEDIM> ebp;
    };
#endif
}

void
MLEBNodeFDLaplacian::customBottomSolve (MLMGT<MultiFab>* mlmg, MultiFab& x, const MultiFab& b,
                                        Real eps_rel, Real eps_abs, int maxiter)
{
    amrex::ignore_unused(maxiter, eps_rel, eps_abs);

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    bool use_custom_solver = (x.size() == 1);
    if (use_custom_solver)
    {
        int const amrlev = 0;
        int const mglev = NMGLevels(0) - 1;
        int const bottom_verbose = mlmg->getBottomVerbose();
        int niters = 0;

        int ret = 0;
        if (ParallelDescriptor::MyProc() == x.DistributionMap()[0])
        {
#ifdef AMREX_USE_EB
            bool const use_eb = !m_levset[amrlev].empty() && m_has_eb[amrlev][mglev][0];
#endif
            auto const& geom = m_geom[amrlev][mglev];
            const auto dxinv = geom.InvCellSizeArray();
#if (AMREX_SPACEDIM == 2)
            const auto sig0 = m_sigma[0];
            const auto dx0 = geom.CellSize(0);
            const auto dx1 = geom.CellSize(1)/std::sqrt(m_sigma[1]);
            const auto xlo = geom.ProbLo(0);
            const auto alpha = m_rz_alpha;
#endif
            AMREX_D_TERM(const Real bx = m_sigma[0]*dxinv[0]*dxinv[0];,
                         const Real by = m_sigma[1]*dxinv[1]*dxinv[1];,
                         const Real bz = m_sigma[2]*dxinv[2]*dxinv[2];)

            auto const& dotmsk = m_bottom_dot_mask[0].const_array();
            auto const& dirmsk = (*m_dirichlet_mask[amrlev][mglev])[0].const_array();
            Box box = x.boxArray()[0];
            Box dbox = amrex::convert(geom.Domain(),IntVect(1));
            LPBase lpbase{dotmsk,dirmsk,
                          GpuArray<Real,AMREX_SPACEDIM>{AMREX_D_DECL(bx,by,bz)},
                          GpuArray<int,AMREX_SPACEDIM>
                              {AMREX_D_DECL(dbox.smallEnd(0),
                                            dbox.smallEnd(1),
                                            dbox.smallEnd(2))},
                          GpuArray<int,AMREX_SPACEDIM>
                              {AMREX_D_DECL(dbox.bigEnd(0),
                                            dbox.bigEnd(1),
                                            dbox.bigEnd(2))},
                          GpuArray<bool,AMREX_SPACEDIM>{AMREX_D_DECL(geom.isPeriodic(0),
                                                                     geom.isPeriodic(1),
                                                                     geom.isPeriodic(2))}};

            EBNodeFDEdgeSigma esig;
            if (m_has_sigma_mf) {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    esig.s[idim] = m_sigma_edge[amrlev][mglev][idim][0].const_array();
                }
            }
            EBNodeFDConstSigma const csig{};
#if (AMREX_SPACEDIM == 2)
            EBNodeFDRZConstSigma const rzsig{sig0};
#endif

#ifdef AMREX_USE_EB
            Array4<Real const> levset;
            GpuArray<Array4<Real const>,AMREX_SPACEDIM> ebp;
            if (use_eb) {
                levset = m_levset[amrlev][mglev][0].const_array();
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    ebp[idim] = m_eb_pos[amrlev][mglev][idim][0].const_array();
                }
            }
#endif

#if (AMREX_SPACEDIM == 2)
            if (m_rz) {
#ifdef AMREX_USE_EB
                if (use_eb) {
                    if (m_has_sigma_mf) {
                        LPRZ<true,EBNodeFDEdgeSigma> lp
                            (lpbase, esig, RZEBData<true>{levset,ebp}, dx0, dx1, xlo, alpha);
                        ret = bicgstab_solve(box, x[0], b[0], lp,
                                             eps_rel, eps_abs, maxiter, bottom_verbose, niters);
                    } else {
                        LPRZ<true,EBNodeFDRZConstSigma> lp
                            (lpbase, rzsig, RZEBData<true>{levset,ebp}, dx0, dx1, xlo, alpha);
                        ret = bicgstab_solve(box, x[0], b[0], lp,
                                             eps_rel, eps_abs, maxiter, bottom_verbose, niters);
                    }
                } else
#endif
                {
                    if (m_has_sigma_mf) {
                        LPRZ<false,EBNodeFDEdgeSigma> lp
                            (lpbase, esig, RZEBData<false>{}, dx0, dx1, xlo, alpha);
                        ret = bicgstab_solve(box, x[0], b[0], lp,
                                             eps_rel, eps_abs, maxiter, bottom_verbose, niters);
                    } else {
                        LPRZ<false,EBNodeFDRZConstSigma> lp
                            (lpbase, rzsig, RZEBData<false>{}, dx0, dx1, xlo, alpha);
                        ret = bicgstab_solve(box, x[0], b[0], lp,
                                             eps_rel, eps_abs, maxiter, bottom_verbose, niters);
                    }
                }
            } else
#endif
#ifdef AMREX_USE_EB
            if (use_eb) {
                if (m_has_sigma_mf) {
                    LPEB<EBNodeFDEdgeSigma> lp(lpbase, levset, ebp, esig);
                    ret = bicgstab_solve(box, x[0], b[0], lp, eps_rel, eps_abs, maxiter, bottom_verbose, niters);
                } else {
                    LPEB<EBNodeFDConstSigma> lp(lpbase, levset, ebp, csig);
                    ret = bicgstab_solve(box, x[0], b[0], lp, eps_rel, eps_abs, maxiter, bottom_verbose, niters);
                }
            } else
#endif
            if (m_has_sigma_mf) {
                LP<EBNodeFDEdgeSigma> lp(lpbase, esig);
                ret = bicgstab_solve(box, x[0], b[0], lp, eps_rel, eps_abs, maxiter, bottom_verbose, niters);
            } else {
                LP<EBNodeFDConstSigma> lp(lpbase, csig);
                ret = bicgstab_solve(box, x[0], b[0], lp, eps_rel, eps_abs, maxiter, bottom_verbose, niters);
            }
        }

        if (ParallelContext::NProcsSub() > 1) {
            int root = ParallelContext::global_to_local_rank(x.DistributionMap()[0]);
            int buf[2] = {ret, niters};
            ParallelDescriptor::Bcast(buf, 2, root, ParallelContext::CommunicatorSub());
            ret = buf[0];
            niters = buf[1];
        }

        if (ret != 0 && mlmg->getVerbose() > 1) {
            amrex::Print() << "MLMG: Bottom solve failed.\n";
        }

        mlmg->postCG(ret, niters);
    } else
#endif
    {
        int ret = mlmg->bottomSolveWithCG(x, b, MLCGSolverT<MultiFab>::Type::BiCGStab);
        mlmg->postCG(ret);
    }
}

}
