#include <AMReX_OpenBC.H>
#include <AMReX_OpenBC_K.H>
#include <AMReX_Algorithm.H>

namespace amrex
{

OpenBCSolver::OpenBCSolver (const Vector<Geometry>& a_geom,
                            const Vector<BoxArray>& a_grids,
                            const Vector<DistributionMapping>& a_dmap,
                            const LPInfo& a_info)
{
    define(a_geom, a_grids, a_dmap, a_info);
}

void OpenBCSolver::define (const Vector<Geometry>& a_geom,
                           const Vector<BoxArray>& a_grids,
                           const Vector<DistributionMapping>& a_dmap,
                           const LPInfo& a_info)
{
    BL_PROFILE("OpenBCSoler::define()");

    m_geom = a_geom;
    m_grids = a_grids;
    m_dmap = a_dmap;
    m_info = a_info;
    for (auto& grids : m_grids) {
        grids.enclosedCells();
    }

    int nlevels = static_cast<int>(m_geom.size());
    m_ba_all.resize(nlevels);
    m_dm_all.resize(nlevels);
    m_geom_all.resize(nlevels);

    Box const domain0 = m_geom[0].Domain();
    m_coarsen_ratio = 8;
    AMREX_ALWAYS_ASSERT(domain0.coarsenable(m_coarsen_ratio));
    int N1d = static_cast<int>(std::round(std::pow(domain0.d_numPts(),1./3.)));
    while (domain0.coarsenable(m_coarsen_ratio*2)
           && 4*m_coarsen_ratio*m_coarsen_ratio <= N1d) {
        m_coarsen_ratio *= 2;
    }

    int ntags = 0;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Box lo = amrex::coarsen(amrex::bdryLo(domain0, idim), m_coarsen_ratio);
        Box hi = amrex::coarsen(amrex::bdryHi(domain0, idim), m_coarsen_ratio);
        BoxList bl({lo,hi});
        IntVect chunk = lo.length();
        while (bl.size() < ParallelContext::NProcsSub()) {
            IntVect chunk_prev = chunk;
            for (int jdim = AMREX_SPACEDIM-1; jdim >= 0; --jdim) {
                if (jdim != idim) {
                    int new_chunk_size = chunk[jdim] / 2;
                    if (bl.size() < ParallelContext::NProcsSub()
                        && new_chunk_size > 0) {
                        chunk[jdim] = new_chunk_size;
                        bl.maxSize(chunk);
                    }
                }
            }
            if (chunk == chunk_prev) {
                break;
            }
        }
        int mgs = std::max(1, 256/m_coarsen_ratio);
        bl.maxSize(mgs);
        bl.refine(m_coarsen_ratio);
        BoxArray ba2d(std::move(bl));
        DistributionMapping dm2d{ba2d};
        m_dpdn[idim].define(ba2d, dm2d, 1, 0);
        ntags += m_dpdn[idim].local_size();
    }

    m_momtags_h.reserve(ntags);
    int nblocks = 0;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (MFIter mfi(m_dpdn[idim]); mfi.isValid(); ++mfi) {
            Box const& b2d = mfi.validbox();
            Orientation::Side side = (b2d.smallEnd(idim) == domain0.smallEnd(idim))
                ? Orientation::low : Orientation::high;
            Orientation face(idim, side);
            m_momtags_h.push_back({m_dpdn[idim].const_array(mfi), b2d, face,
                                   nblocks});
            nblocks += static_cast<int>(b2d.numPts())
                / (m_coarsen_ratio*m_coarsen_ratio);
        }
    }
    m_nblocks_local = nblocks;

#ifdef AMREX_USE_GPU
    if (ntags > 0) {
        m_momtags_d.resize(ntags);
        Gpu::copyAsync(Gpu::hostToDevice, m_momtags_h.begin(), m_momtags_h.end(), m_momtags_d.begin());

        m_nthreads_momtag = (m_coarsen_ratio == 8) ? 64 : 128;
        int ntotgpublocks = 0;
        m_ngpublocks_h.reserve(ntags+1);
        for (auto const& tag : m_momtags_h) {
            m_ngpublocks_h.push_back(ntotgpublocks);
            Box cb2d = amrex::coarsen(tag.b2d, m_coarsen_ratio);
            ntotgpublocks += static_cast<int>(cb2d.numPts());
        }
        m_ngpublocks_h.push_back(ntotgpublocks);
        m_ngpublocks_d.resize(m_ngpublocks_h.size());
        Gpu::copyAsync(Gpu::hostToDevice, m_ngpublocks_h.begin(), m_ngpublocks_h.end(),
                       m_ngpublocks_d.begin());
    }
#endif

    const auto *const dx = m_geom[0].CellSize();
    Real dmax = amrex::max(std::sqrt(dx[0]*dx[0]+dx[1]*dx[1]),
                           std::sqrt(dx[0]*dx[0]+dx[2]*dx[2]),
                           std::sqrt(dx[1]*dx[1]+dx[2]*dx[2]));
    m_ngrowdomain[0] = static_cast<int>(std::ceil(dmax/dx[0])) * m_coarsen_ratio;
    m_ngrowdomain[1] = static_cast<int>(std::ceil(dmax/dx[1])) * m_coarsen_ratio;
    m_ngrowdomain[2] = static_cast<int>(std::ceil(dmax/dx[2])) * m_coarsen_ratio;
    // This is the minimal size we need to embiggen the domain.

    Box const domain1 = amrex::grow(domain0, m_ngrowdomain);
    BoxList bl_crse_grown_faces(IndexType::TheNodeType());
    for (OrientationIter oit; oit.isValid(); ++oit) {
        Orientation face = oit();
        Box face_box = amrex::surroundingNodes(amrex::bdryNode(domain1,face));
        face_box.coarsen(m_coarsen_ratio);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (idim != face.coordDir()) {
                face_box.grow(idim,openbc::P);
            }
        }
        bl_crse_grown_faces.push_back(face_box);
    }

    bl_crse_grown_faces.maxSize(16); // xxxxx make this a parameter?
    BoxArray ba_crse_grown_faces(std::move(bl_crse_grown_faces));
    DistributionMapping dm_crse_grown_faces(ba_crse_grown_faces);
    m_crse_grown_faces_phi.define(ba_crse_grown_faces, dm_crse_grown_faces, 1, 0);

    BoxList blg = amrex::boxDiff(domain1, domain0);
    blg.maxSize(std::max(64,m_coarsen_ratio)); // xxxxx make this a parameter?
    m_bag = BoxArray(std::move(blg));
    DistributionMapping dmg(m_bag);
    m_phind.define(amrex::coarsen(amrex::convert(m_bag,IntVect(1)),m_coarsen_ratio),
                   dmg, 1, openbc::P);

    BoxList bl0 = m_grids[0].boxList();
    BoxList bl1 = m_bag.boxList();
    Vector<int> p0 = m_dmap[0].ProcessorMap();
    Vector<int> p1 = dmg.ProcessorMap();
    bl0.join(bl1);
    p0.insert(p0.end(), p1.begin(), p1.end());
    IntVect offset = -domain1.smallEnd();
    for (auto& b : bl0) {
        b.shift(offset);
    }
    m_ba_all[0] = BoxArray(std::move(bl0));
    m_dm_all[0] = DistributionMapping(std::move(p0));
    m_box_offset.push_back(offset);

    const auto *const problo = m_geom[0].ProbLo();
    const auto *const probhi = m_geom[0].ProbHi();
    std::array<Real,AMREX_SPACEDIM> problo_all, probhi_all;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        problo_all[idim] = problo[idim] - static_cast<Real>(m_ngrowdomain[idim])*dx[idim];
        probhi_all[idim] = probhi[idim] + static_cast<Real>(m_ngrowdomain[idim])*dx[idim];
    }
    m_geom_all[0] = Geometry(amrex::shift(domain1,offset),
                             RealBox(problo_all,probhi_all),
                             m_geom[0].Coord(), m_geom[0].isPeriodic());

    for (int ilev = 1; ilev < nlevels; ++ilev) {
        IntVect rr = m_geom[ilev].Domain().length()
                   / m_geom[ilev-1].Domain().length();
        offset *= rr;
        m_box_offset.push_back(offset);
        m_geom_all[ilev] = amrex::refine(m_geom_all[ilev-1], rr);
        m_ba_all[ilev] = a_grids[ilev];
        m_ba_all[ilev].shift(offset);
        m_dm_all[ilev] = a_dmap[ilev];
    }
}

void OpenBCSolver::setVerbose (int v) noexcept
{
    m_verbose = v;
}

void OpenBCSolver::setBottomVerbose (int v) noexcept
{
    m_bottom_verbose = v;
}

void OpenBCSolver::useHypre (bool use_hypre) noexcept
{
    if (use_hypre) {
        m_bottom_solver_type = BottomSolver::hypre;
        m_info.setMaxCoarseningLevel(0);
#ifndef AMREX_USE_HYPRE
        amrex::Abort("OpenBCSolver: Must enable Hypre support to use it.");
#endif
    }
}

Real OpenBCSolver::solve (const Vector<MultiFab*>& a_sol,
                          const Vector<MultiFab const*>& a_rhs,
                          Real a_tol_rel, Real a_tol_abs)
{
    BL_PROFILE("OpenBCSolver::solve()");

    auto solve_start_time = amrex::second();

    int nlevels = static_cast<int>(m_geom.size());

    BL_PROFILE_VAR("OpenBCSolver::MG1", blp_mg1);

    if (m_poisson_1 == nullptr) {
        m_poisson_1 = std::make_unique<MLPoisson>(m_geom, m_grids, m_dmap, m_info);
        m_poisson_1->setVerbose(m_verbose);
        m_poisson_1->setMaxOrder(4);
        m_poisson_1->setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                               LinOpBCType::Dirichlet,
                                               LinOpBCType::Dirichlet)},
                                 {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                               LinOpBCType::Dirichlet,
                                               LinOpBCType::Dirichlet)});
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            m_poisson_1->setLevelBC(ilev, nullptr);
        }

        m_mlmg_1 = std::make_unique<MLMG>(*m_poisson_1);
        m_mlmg_1->setVerbose(m_verbose);
        m_mlmg_1->setBottomVerbose(m_bottom_verbose);
        m_mlmg_1->setBottomSolver(m_bottom_solver_type);
#ifdef AMREX_USE_HYPRE
        if (m_bottom_solver_type == BottomSolver::hypre) {
            m_mlmg_1->setHypreInterface(Hypre::Interface::structed);
        }
#endif
    }
    m_mlmg_1->solve(a_sol, a_rhs, a_tol_rel, a_tol_abs);

    BL_PROFILE_VAR_STOP(blp_mg1);

    Array<MultiFab,AMREX_SPACEDIM> dpdn_tmp;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        dpdn_tmp[idim].define(amrex::convert(m_grids[0],
                                             IntVect::TheDimensionVector(idim)),
                          m_dmap[0], 1, 0);
    }
    m_poisson_1->get_dpdn_on_domain_faces(GetArrOfPtrs(dpdn_tmp), *a_sol[0]);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_dpdn[idim].ParallelCopy(dpdn_tmp[idim]);
    }

    {
        Gpu::DeviceVector<openbc::Moments> moments(m_nblocks_local);
        compute_moments(moments);
        compute_potential(moments);
    }

    MultiFab rhsg(m_bag, m_phind.DistributionMap(), 1, a_rhs[0]->nGrowVect());
    rhsg.setVal(0._rt);

    MultiFab solg(m_bag, m_phind.DistributionMap(), 1, 1);
    solg.setVal(0._rt);
    interpolate_potential(solg);

    const int nboxes0 = static_cast<int>(m_grids[0].size());
    Vector<MultiFab> sol_all(nlevels);
    Vector<MultiFab> rhs_all(nlevels);

    sol_all[0].define(m_ba_all[0], m_dm_all[0], 1, solg.nGrowVect(),
                      MFInfo().SetAlloc(false));
    rhs_all[0].define(m_ba_all[0], m_dm_all[0], 1, rhsg.nGrowVect(),
                      MFInfo().SetAlloc(false));

    for (MFIter mfi(sol_all[0]); mfi.isValid(); ++mfi) {
        const int index = mfi.index();
        FArrayBox solfab, rhsfab;
        if (index < nboxes0) {
            FArrayBox& sfab0 = (*a_sol[0])[index];
            if (sol_all[0].nGrowVect() == a_sol[0]->nGrowVect()) {
                solfab = FArrayBox(sfab0, amrex::make_alias, 0, 1);
            } else {
                Box b = sfab0.box();
                b.grow(sol_all[0].nGrowVect()-a_sol[0]->nGrowVect());
                solfab.resize(b,1);
                solfab.template setVal<RunOn::Device>(0._rt);
            }
            rhsfab = FArrayBox((*a_rhs[0])[index], amrex::make_alias, 0, 1);
        } else {
            solfab = FArrayBox(solg[index-nboxes0], amrex::make_alias, 0, 1);
            rhsfab = FArrayBox(rhsg[index-nboxes0], amrex::make_alias, 0, 1);
        }
        solfab.shift(m_box_offset[0]);
        rhsfab.shift(m_box_offset[0]);
        sol_all[0].setFab(mfi, std::move(solfab));
        rhs_all[0].setFab(mfi, std::move(rhsfab));
    }

    for (int ilev = 1; ilev < nlevels; ++ilev) {
        sol_all[ilev].define(m_ba_all[ilev], m_dm_all[ilev], 1,
                             a_sol[ilev]->nGrowVect(), MFInfo().SetAlloc(false));
        rhs_all[ilev].define(m_ba_all[ilev], m_dm_all[ilev], 1,
                             a_rhs[ilev]->nGrowVect(), MFInfo().SetAlloc(false));
        for (MFIter mfi(sol_all[ilev]); mfi.isValid(); ++mfi) {
            const int index = mfi.index();
            auto const& a_sol_fab = (*a_sol[ilev])[index];
            auto const& a_rhs_fab = (*a_rhs[ilev])[index];
            FArrayBox solfab(a_sol_fab, amrex::make_alias, 0, 1);
            FArrayBox rhsfab(a_rhs_fab, amrex::make_alias, 0, 1);
            solfab.shift(m_box_offset[ilev]);
            rhsfab.shift(m_box_offset[ilev]);
            sol_all[ilev].setFab(mfi, std::move(solfab));
            rhs_all[ilev].setFab(mfi, std::move(rhsfab));
        }
    }

    BL_PROFILE_VAR("OpenBCSolver::MG2", blp_mg2);

    if (m_poisson_2 == nullptr) {
        m_poisson_2 = std::make_unique<MLPoisson>(m_geom_all, m_ba_all,
                                                  m_dm_all, m_info);
        m_poisson_2->setVerbose(m_verbose);
        m_poisson_2->setMaxOrder(4);
        m_poisson_2->setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                               LinOpBCType::Dirichlet,
                                               LinOpBCType::Dirichlet)},
                                 {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                               LinOpBCType::Dirichlet,
                                               LinOpBCType::Dirichlet)});
        m_poisson_2->setLevelBC(0, &sol_all[0]);
        for (int ilev = 1; ilev < nlevels; ++ilev) {
            m_poisson_2->setLevelBC(ilev, nullptr);
        }

        m_mlmg_2 = std::make_unique<MLMG>(*m_poisson_2);
        m_mlmg_2->setVerbose(m_verbose);
        m_mlmg_2->setBottomVerbose(m_bottom_verbose);
        m_mlmg_2->setBottomSolver(m_bottom_solver_type);
#ifdef AMREX_USE_HYPRE
        if (m_bottom_solver_type == BottomSolver::hypre) {
            m_mlmg_2->setHypreInterface(Hypre::Interface::structed);
        }
#endif
    } else {
        m_poisson_2->setLevelBC(0, &sol_all[0]);
    }

    Real err = m_mlmg_2->solve(GetVecOfPtrs(sol_all), GetVecOfConstPtrs(rhs_all),
                               a_tol_rel, a_tol_abs);

    BL_PROFILE_VAR_STOP(blp_mg2);

    if (sol_all[0].nGrowVect() != a_sol[0]->nGrowVect()) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*a_sol[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real const> const& sall = sol_all[0].const_array(mfi.index());
            Array4<Real> const& s = a_sol[0]->array(mfi);
            auto const offset = m_box_offset[0].dim3();
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                s(i,j,k) = sall(i+offset.x,j+offset.y,k+offset.z);
            });
        }
    }

    auto solve_stop_time = amrex::second();
    if (m_verbose >= 1) {
        amrex::Print() << "OpenBCSolver time = "
                       << solve_stop_time - solve_start_time << "\n";
    }

    return err;
}

void OpenBCSolver::compute_moments (Gpu::DeviceVector<openbc::Moments>& moments)
{
    BL_PROFILE("OpenBCSolver::comp_mom()");

    auto const problo = m_geom[0].ProbLoArray();
    auto const probhi = m_geom[0].ProbHiArray();
    auto const dx     = m_geom[0].CellSizeArray();

#ifdef AMREX_USE_GPU
    if (m_momtags_h.size() > 0)
    {
        int crse_ratio = m_coarsen_ratio;
        int ntags = m_momtags_h.size();
        openbc::Moments* pm = moments.data();
        openbc::MomTag const* ptag = m_momtags_d.data();
        int const* pnblks = m_ngpublocks_d.data();
        std::size_t shared_mem_bytes = m_nthreads_momtag * sizeof(openbc::Moments::array_type);

#ifdef AMREX_USE_SYCL
        amrex::ignore_unused(problo,probhi,dx,crse_ratio,ntags,pm,ptag,pnblks,
                             shared_mem_bytes);
        amrex::Abort("xxxx SYCL todo: openbc compute_moments");
#else
        amrex::launch(m_ngpublocks_h.back(), m_nthreads_momtag, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept
        {
            Gpu::SharedMemory<openbc::Moments::array_type> gsm;
            openbc::Moments::array_type* const shared = gsm.dataPtr();
            openbc::Moments::array_type& tmom = shared[threadIdx.x];
            for (int i = 0; i < (openbc::M+1)*(openbc::M+2)/2; ++i) {
                tmom[i] = Real(0.);
            }

            int tag_id = amrex::bisect(pnblks, 0, ntags, static_cast<int>(blockIdx.x));
            int iblock = blockIdx.x - pnblks[tag_id]; // iblock'th gpublock on this box.
            auto const& tag = ptag[tag_id];
            openbc::Moments& mom = pm[tag.offset+iblock];
            if (tag.face.coordDir() == 0) {
                int const nby = tag.b2d.length(1) / crse_ratio;
                int const kb = iblock / nby;
                int const jb = iblock - kb*nby;
                int const i = tag.b2d.smallEnd(0);
                int const jlo = tag.b2d.smallEnd(1) + jb*crse_ratio;
                int const klo = tag.b2d.smallEnd(2) + kb*crse_ratio;
                Real const fac = dx[1]*dx[2];
                Real const xc = tag.face.isLow() ? problo[0] : probhi[0];
                for (int icell = threadIdx.x; icell < crse_ratio*crse_ratio; icell += blockDim.x) {
                    int k = icell/crse_ratio;
                    int j = icell - k*crse_ratio;
                    Real const yy = (j-crse_ratio/2+Real(0.5))*dx[1];
                    Real const zz = (k-crse_ratio/2+Real(0.5))*dx[2];
                    j += jlo;
                    k += klo;
                    Real const charge = tag.gp(i,j,k) * fac;
                    Real zpow = Real(1.);
                    int m = 0;
                    for (int q = 0; q <= openbc::M; ++q) {
                        Real ypow = Real(1.);
                        for (int p = 0; p <= openbc::M-q; ++p) {
                            tmom[m++] += charge*ypow*zpow;
                            ypow *= yy;
                        }
                        zpow *= zz;
                    }
                }
                if (threadIdx.x == 0) {
                    mom.x = xc;
                    mom.y = problo[1] + dx[1]*(jlo + crse_ratio/2);
                    mom.z = problo[2] + dx[2]*(klo + crse_ratio/2);
                    mom.face = tag.face;
                }
            } else if (tag.face.coordDir() == 1) {
                int const nbx = tag.b2d.length(0) / crse_ratio;
                int const kb = iblock / nbx;
                int const ib = iblock - kb*nbx;
                int const j = tag.b2d.smallEnd(1);
                int const ilo = tag.b2d.smallEnd(0) + ib*crse_ratio;
                int const klo = tag.b2d.smallEnd(2) + kb*crse_ratio;
                Real const fac = dx[0]*dx[2];
                Real const yc = tag.face.isLow() ? problo[1] : probhi[1];
                for (int icell = threadIdx.x; icell < crse_ratio*crse_ratio; icell += blockDim.x) {
                    int k = icell/crse_ratio;
                    int i = icell - k*crse_ratio;
                    Real const xx = (i-crse_ratio/2+Real(0.5))*dx[0];
                    Real const zz = (k-crse_ratio/2+Real(0.5))*dx[2];
                    i += ilo;
                    k += klo;
                    Real const charge = tag.gp(i,j,k) * fac;
                    Real zpow = Real(1.);
                    int m = 0;
                    for (int q = 0; q <= openbc::M; ++q) {
                        Real xpow = Real(1.);
                        for (int p = 0; p <= openbc::M-q; ++p) {
                            tmom[m++] += charge*xpow*zpow;
                            xpow *= xx;
                        }
                        zpow *= zz;
                    }
                }
                if (threadIdx.x == 0) {
                    mom.x = problo[0] + dx[0]*(ilo + crse_ratio/2);
                    mom.y = yc;
                    mom.z = problo[2] + dx[2]*(klo + crse_ratio/2);
                    mom.face = tag.face;
                }
            } else {
                int const nbx = tag.b2d.length(0) / crse_ratio;
                int const jb = iblock / nbx;
                int const ib = iblock - jb*nbx;
                int const k = tag.b2d.smallEnd(2);
                int const ilo = tag.b2d.smallEnd(0) + ib*crse_ratio;
                int const jlo = tag.b2d.smallEnd(1) + jb*crse_ratio;
                Real const fac = dx[0]*dx[1];
                Real const zc = tag.face.isLow() ? problo[2] : probhi[2];
                for (int icell = threadIdx.x; icell < crse_ratio*crse_ratio; icell += blockDim.x) {
                    int j = icell/crse_ratio;
                    int i = icell - j*crse_ratio;
                    Real const xx = (i-crse_ratio/2+Real(0.5))*dx[0];
                    Real const yy = (j-crse_ratio/2+Real(0.5))*dx[1];
                    i += ilo;
                    j += jlo;
                    Real const charge = tag.gp(i,j,k) * fac;
                    Real ypow = Real(1.);
                    int m = 0;
                    for (int q=0; q <= openbc::M; ++q) {
                        Real xpow = Real(1.);
                        for (int p = 0; p <= openbc::M-q; ++p) {
                            tmom[m++] += charge*xpow*ypow;
                            xpow *= xx;
                        }
                        ypow *= yy;
                    }
                }
                if (threadIdx.x == 0) {
                    mom.x = problo[0] + dx[0]*(ilo + crse_ratio/2);
                    mom.y = problo[1] + dx[1]*(jlo + crse_ratio/2);
                    mom.z = zc;
                    mom.face = tag.face;
                }
            }
            openbc::scale_moments(tmom);

            __syncthreads();

            if (threadIdx.x < (openbc::M+1)*(openbc::M+2)/2) {
                mom.mom[threadIdx.x] = Real(0.);
                for (unsigned int i = 0; i < blockDim.x; ++i) {
                    mom.mom[threadIdx.x] += shared[i][threadIdx.x];
                }
            }
        });
#endif
    }
#else
    for (auto const& tag : m_momtags_h) {
        if (tag.face.coordDir() == 0) {
            int nby = tag.b2d.length(1) / m_coarsen_ratio;
            int nbz = tag.b2d.length(2) / m_coarsen_ratio;
            int i = tag.b2d.smallEnd(0);
            int jlo = tag.b2d.smallEnd(1);
            int klo = tag.b2d.smallEnd(2);
            Real fac = dx[1]*dx[2];
            Real xc = tag.face.isLow() ? problo[0] : probhi[0];
            for (int kb = 0; kb < nbz; ++kb) {
            for (int jb = 0; jb < nby; ++jb) {
                openbc::Moments& mom = moments[tag.offset+jb+kb*nby];
                for (auto& m : mom.mom) {
                    m = 0._rt;
                }
                for (int kk = 0; kk < m_coarsen_ratio; ++kk) {
                for (int jj = 0; jj < m_coarsen_ratio; ++jj) {
                    Real charge = tag.gp(i, jlo+jb*m_coarsen_ratio+jj,
                                         klo+kb*m_coarsen_ratio+kk) * fac;
                    Real yy = (jj-m_coarsen_ratio/2+0.5_rt)*dx[1]; // NOLINT
                    Real zz = (kk-m_coarsen_ratio/2+0.5_rt)*dx[2]; // NOLINT
                    Real zpow = 1._rt;
                    int m = 0;
                    for (int q = 0; q <= openbc::M; ++q) {
                        Real ypow = 1._rt;
                        for (int p = 0; p <= openbc::M-q; ++p) {
                            mom.mom[m++] += charge*ypow*zpow;
                            ypow *= yy;
                        }
                        zpow *= zz;
                    }
                }}
                openbc::scale_moments(mom.mom);
                // center of the block
                mom.x = xc;
                mom.y = problo[1] + dx[1]*static_cast<Real>(tag.b2d.smallEnd(1)
                                           + jb*m_coarsen_ratio
                                           + m_coarsen_ratio/2); // NOLINT
                mom.z = problo[2] + dx[2]*static_cast<Real>(tag.b2d.smallEnd(2)
                                           + kb*m_coarsen_ratio
                                           + m_coarsen_ratio/2); // NOLINT
                mom.face = tag.face;
            }}
        } else if (tag.face.coordDir() == 1) {
            int nbx = tag.b2d.length(0) / m_coarsen_ratio;
            int nbz = tag.b2d.length(2) / m_coarsen_ratio;
            int j = tag.b2d.smallEnd(1);
            int ilo = tag.b2d.smallEnd(0);
            int klo = tag.b2d.smallEnd(2);
            Real fac = dx[0]*dx[2];
            Real yc = tag.face.isLow() ? problo[1] : probhi[1];
            for (int kb = 0; kb < nbz; ++kb) {
            for (int ib = 0; ib < nbx; ++ib) {
                openbc::Moments& mom = moments[tag.offset+ib+kb*nbx];
                for (auto& m : mom.mom) {
                    m = 0._rt;
                }
                for (int kk = 0; kk < m_coarsen_ratio; ++kk) {
                for (int ii = 0; ii < m_coarsen_ratio; ++ii) {
                    Real charge = tag.gp(ilo+ib*m_coarsen_ratio+ii, j,
                                         klo+kb*m_coarsen_ratio+kk) * fac;
                    Real xx = (ii-m_coarsen_ratio/2+0.5_rt)*dx[0]; // NOLINT
                    Real zz = (kk-m_coarsen_ratio/2+0.5_rt)*dx[2]; // NOLINT
                    Real zpow = 1._rt;
                    int m = 0;
                    for (int q = 0; q <= openbc::M; ++q) {
                        Real xpow = 1._rt;
                        for (int p = 0; p <= openbc::M-q; ++p) {
                            mom.mom[m++] += charge*xpow*zpow;
                            xpow *= xx;
                        }
                        zpow *= zz;
                    }
                }}
                openbc::scale_moments(mom.mom);
                mom.x = problo[0] + dx[0]*static_cast<Real>(tag.b2d.smallEnd(0)
                                           + ib*m_coarsen_ratio
                                           + m_coarsen_ratio/2); // NOLINT
                mom.y = yc;
                mom.z = problo[2] + dx[2]*static_cast<Real>(tag.b2d.smallEnd(2)
                                           + kb*m_coarsen_ratio
                                           + m_coarsen_ratio/2); // NOLINT
                mom.face = tag.face;
            }}
        } else {
            int nbx = tag.b2d.length(0) / m_coarsen_ratio;
            int nby = tag.b2d.length(1) / m_coarsen_ratio;
            int k = tag.b2d.smallEnd(2);
            int ilo = tag.b2d.smallEnd(0);
            int jlo = tag.b2d.smallEnd(1);
            Real fac = dx[0]*dx[1];
            Real zc = tag.face.isLow() ? problo[2] : probhi[2];
            for (int jb = 0; jb < nby; ++jb) {
            for (int ib = 0; ib < nbx; ++ib) {
                openbc::Moments& mom = moments[tag.offset+ib+jb*nbx];
                for (auto& m : mom.mom) {
                    m = 0._rt;
                }
                for (int jj = 0; jj < m_coarsen_ratio; ++jj) {
                for (int ii = 0; ii < m_coarsen_ratio; ++ii) {
                    Real charge = tag.gp(ilo+ib*m_coarsen_ratio+ii,
                                         jlo+jb*m_coarsen_ratio+jj, k) * fac;
                    Real xx = (ii-m_coarsen_ratio/2+0.5_rt)*dx[0]; // NOLINT
                    Real yy = (jj-m_coarsen_ratio/2+0.5_rt)*dx[1]; // NOLINT
                    Real ypow = 1._rt;
                    int m = 0;
                    for (int q = 0; q <= openbc::M; ++q) {
                        Real xpow = 1._rt;
                        for (int p = 0; p <= openbc::M-q; ++p) {
                            mom.mom[m++] += charge*xpow*ypow;
                            xpow *= xx;
                        }
                        ypow *= yy;
                    }
                }}
                openbc::scale_moments(mom.mom);
                mom.x = problo[0] + dx[0]*static_cast<Real>(tag.b2d.smallEnd(0)
                                           + ib*m_coarsen_ratio
                                           + m_coarsen_ratio/2); // NOLINT
                mom.y = problo[1] + dx[1]*static_cast<Real>(tag.b2d.smallEnd(1)
                                           + jb*m_coarsen_ratio
                                           + m_coarsen_ratio/2); // NOLINT
                mom.z = zc;
                mom.face = tag.face;
            }}
        }
    }
#endif

#ifdef AMREX_USE_MPI
    bcast_moments(moments);
#endif
    m_nblocks = static_cast<int>(moments.size());
}

#ifdef AMREX_USE_MPI
void OpenBCSolver::bcast_moments (Gpu::DeviceVector<openbc::Moments>& moments)
{
    if (ParallelContext::NProcsSub() > 1)
    {
        MPI_Comm comm = ParallelContext::CommunicatorSub();
        if (m_nblocks == 0) {
            int count = static_cast<int>(moments.size() * sizeof(openbc::Moments));
            m_countvec.resize(ParallelContext::NProcsSub());
            MPI_Allgather(&count, 1, MPI_INT, m_countvec.data(), 1, MPI_INT, comm);

            m_offset.resize(m_countvec.size(), 0);
            Long count_tot = m_countvec[0];
            for (int i = 1, N = static_cast<int>(m_offset.size()); i < N; ++i) {
                m_offset[i] = m_offset[i-1] + m_countvec[i-1];
                count_tot += m_countvec[i];
            }

            if (count_tot > static_cast<Long>(std::numeric_limits<int>::max())) {
                amrex::Abort("OpenBC: integer overflow. Let us know and we will fix this.");
            }

            m_nblocks = static_cast<int>(count_tot/sizeof(openbc::Moments));
        }

        Gpu::DeviceVector<openbc::Moments> moments_all(m_nblocks);

#ifdef AMREX_USE_GPU
        Gpu::PinnedVector<openbc::Moments> h_moments(moments.size());
        Gpu::PinnedVector<openbc::Moments> h_moments_all(moments_all.size());
        Gpu::copyAsync(Gpu::deviceToHost, moments.begin(), moments.end(),
                       h_moments.begin());
        Gpu::streamSynchronize();
#else
        auto const& h_moments = moments;
        auto& h_moments_all = moments_all;
#endif

        int count = m_nblocks_local*static_cast<int>(sizeof(openbc::Moments));
        MPI_Allgatherv(h_moments.data(), count, MPI_CHAR, h_moments_all.data(),
                       m_countvec.data(), m_offset.data(), MPI_CHAR, comm);

#ifdef AMREX_USE_GPU
        Gpu::copyAsync(Gpu::hostToDevice, h_moments_all.begin(), h_moments_all.end(),
                       moments_all.begin());
        Gpu::streamSynchronize();
#endif

        std::swap(moments, moments_all);
    }
}
#endif

void OpenBCSolver::compute_potential (Gpu::DeviceVector<openbc::Moments> const& moments)
{
    BL_PROFILE("OpenBCSolver::comp_phi()");

    auto const problo = m_geom[0].ProbLoArray();
    auto const dx     = m_geom[0].CellSizeArray();

    int crse_ratio = m_coarsen_ratio;
    int nblocks = m_nblocks;
    openbc::Moments const* pmom = moments.data();
    for (MFIter mfi(m_crse_grown_faces_phi); mfi.isValid(); ++mfi) {
        Box const& b = mfi.validbox();
        Array4<Real> const& phi_arr = m_crse_grown_faces_phi.array(mfi);
#if defined(AMREX_USE_GPU)
        const auto lo  = amrex::lbound(b);
        const auto len = amrex::length(b);
        const auto lenxy = len.x*len.y;
        const auto lenx = len.x;
#ifdef AMREX_USE_SYCL
        amrex::ignore_unused(problo,dx,crse_ratio,nblocks,pmom,b,phi_arr,lo,
                             lenxy,lenx);
        amrex::Abort("xxxxx SYCL todo: openbc compute_potential");
#else
        amrex::launch(b.numPts(), AMREX_GPU_MAX_THREADS, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept
        {
            int icell = blockIdx.x;
            int k =  icell /   lenxy;
            int j = (icell - k*lenxy) /   lenx;
            int i = (icell - k*lenxy) - j*lenx;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            Real xb = problo[0] + i*crse_ratio*dx[0];
            Real yb = problo[1] + j*crse_ratio*dx[1];
            Real zb = problo[2] + k*crse_ratio*dx[2];
            Real phi = Real(0.);
            for (int iblock = threadIdx.x; iblock < nblocks; iblock += blockDim.x) {
                phi += openbc::block_potential(pmom[iblock], xb, yb, zb);
            }
            Real phitot = Gpu::blockReduceSum<AMREX_GPU_MAX_THREADS>(phi);
            if (threadIdx.x == 0) {
                phi_arr(i,j,k) = phitot;
            }
        });
#endif
#else
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) noexcept
        {
            Real xb = problo[0] + static_cast<Real>(i*crse_ratio)*dx[0];
            Real yb = problo[1] + static_cast<Real>(j*crse_ratio)*dx[1];
            Real zb = problo[2] + static_cast<Real>(k*crse_ratio)*dx[2];
            Real phi = 0._rt;
            for (int iblock = 0; iblock < nblocks; ++iblock) {
                phi += openbc::block_potential(pmom[iblock], xb, yb, zb);
            }
            phi_arr(i,j,k) = phi;
        });
#endif
    }

    m_phind.ParallelCopy(m_crse_grown_faces_phi, 0, 0, 1, IntVect(0),
                         m_phind.nGrowVect());
}

void OpenBCSolver::interpolate_potential (MultiFab& solg)
{
    BL_PROFILE("OpenBCSolver::interp_phi");

    Box const domain1 = amrex::grow(m_geom[0].Domain(), m_ngrowdomain);
    int crse_ratio = m_coarsen_ratio;

    for (MFIter mfi(solg); mfi.isValid(); ++mfi) {
        Box const& vbx = mfi.validbox();
        for (OrientationIter oit; oit.isValid(); ++oit) {
            Orientation face = oit();
            if (vbx[face] == domain1[face]) {
                Array4<Real> const& solg_arr = solg.array(mfi);
                Array4<Real const> const& phi_arr = m_phind.const_array(mfi);
                Box const& b2d = amrex::bdryNode(vbx, face);
                int offset = face.isLow() ? -1 : 0;
                if (face.coordDir() == 0) {
                    Box b = amrex::coarsen(b2d,IntVect(crse_ratio,crse_ratio,1));
                    b.grow(1,openbc::P).surroundingNodes(1);
                    FArrayBox tmpfab(b,1,The_Async_Arena());
                    Array4<Real> const& tmp = tmpfab.array();
                    Array4<Real const> const& ctmp = tmpfab.const_array();
                    amrex::ParallelFor(b,
                    [=] AMREX_GPU_DEVICE (int ic, int jc, int k) noexcept
                    {
                        tmp(ic,jc,k) = openbc::interpccz(ic,jc,k,phi_arr,crse_ratio);
                    });
                    b = amrex::coarsen(b2d,IntVect(crse_ratio,1,1));
                    amrex::ParallelFor(b,
                    [=] AMREX_GPU_DEVICE (int ic, int j, int k) noexcept
                    {
                        int i = ic*crse_ratio+offset;
                        solg_arr(i,j,k) = openbc::interpccy(ic,j,k,ctmp,crse_ratio);
                    });
                } else if (face.coordDir() == 1) {
                    Box b = amrex::coarsen(b2d,IntVect(crse_ratio,crse_ratio,1));
                    b.grow(0,openbc::P).surroundingNodes(0);
                    FArrayBox tmpfab(b,1,The_Async_Arena());
                    Array4<Real> const& tmp = tmpfab.array();
                    Array4<Real const> const& ctmp = tmpfab.const_array();
                    amrex::ParallelFor(b,
                    [=] AMREX_GPU_DEVICE (int ic, int jc, int k) noexcept
                    {
                        tmp(ic,jc,k) = openbc::interpccz(ic,jc,k,phi_arr,crse_ratio);
                    });
                    b = amrex::coarsen(b2d,IntVect(1,crse_ratio,1));
                    amrex::ParallelFor(b,
                    [=] AMREX_GPU_DEVICE (int i, int jc, int k) noexcept
                    {
                        int j = jc*crse_ratio+offset;
                        solg_arr(i,j,k) = openbc::interpccx(i,jc,k,ctmp,crse_ratio);
                    });
                } else {
                    Box b = amrex::coarsen(b2d,IntVect(crse_ratio,1,crse_ratio));
                    b.grow(0,openbc::P).surroundingNodes(0);
                    FArrayBox tmpfab(b,1,The_Async_Arena());
                    Array4<Real> const& tmp = tmpfab.array();
                    Array4<Real const> const& ctmp = tmpfab.const_array();
                    amrex::ParallelFor(b,
                    [=] AMREX_GPU_DEVICE (int ic, int j, int kc) noexcept
                    {
                        tmp(ic,j,kc) = openbc::interpccy(ic,j,kc,phi_arr,crse_ratio);
                    });
                    b = amrex::coarsen(b2d,IntVect(1,1,crse_ratio));
                    amrex::ParallelFor(b,
                    [=] AMREX_GPU_DEVICE (int i, int j, int kc) noexcept
                    {
                        int k = kc*crse_ratio+offset;
                        solg_arr(i,j,k) = openbc::interpccx(i,j,kc,ctmp,crse_ratio);
                    });
                }
            }
        }
    }
}

namespace openbc {
std::ostream& operator<< (std::ostream& os, Moments const& mom)
{
    os << "Face " << mom.face << ", x = " << mom.x << ", y = " << mom.y
       << ", z = " << mom.z << "\n"
       << "  " << mom.mom[0] << "\n"
       << "  " << mom.mom[1] << ", " << mom.mom[8] << "\n"
       << "  " << mom.mom[2] << ", " << mom.mom[9] << ", " << mom.mom[15] << "\n"
       << "  " << mom.mom[3] << ", " << mom.mom[10] << ", " << mom.mom[16]
       << ", " << mom.mom[21] << "\n"
       << "  " << mom.mom[4] << ", " << mom.mom[11] << ", " << mom.mom[17]
       << ", " << mom.mom[22] << ", " << mom.mom[26] << "\n"
       << "  " << mom.mom[5] << ", " << mom.mom[12] << ", " << mom.mom[18]
       << ", " << mom.mom[23] << ", " << mom.mom[27] << ", " << mom.mom[30] << "\n"
       << "  " << mom.mom[6] << ", " << mom.mom[13] << ", " << mom.mom[19]
       << ", " << mom.mom[24] << ", " << mom.mom[28] << ", " << mom.mom[31]
       << ", " << mom.mom[33] << "\n"
       << "  " << mom.mom[7] << ", " << mom.mom[14] << ", " << mom.mom[20]
       << ", " << mom.mom[25] << ", " << mom.mom[29] << ", " << mom.mom[32]
       << ", " << mom.mom[34] << ", " << mom.mom[35] << "\n";
    return os;
}
}

}
