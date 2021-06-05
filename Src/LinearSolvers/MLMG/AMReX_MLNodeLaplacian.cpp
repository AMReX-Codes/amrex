#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_algoim.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex {

MLNodeLaplacian::MLNodeLaplacian (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap,
                                  const LPInfo& a_info,
                                  const Vector<FabFactory<FArrayBox> const*>& a_factory,
                                  Real  a_const_sigma)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory, a_const_sigma);
}

#ifdef AMREX_USE_EB
MLNodeLaplacian::MLNodeLaplacian (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap,
                                  const LPInfo& a_info,
                                  const Vector<EBFArrayBoxFactory const*>& a_factory,
                                  Real  a_const_sigma)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory, a_const_sigma);
}
#endif

MLNodeLaplacian::~MLNodeLaplacian ()
{}

void
MLNodeLaplacian::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<FabFactory<FArrayBox> const*>& a_factory,
                         Real  a_const_sigma)
{
    BL_PROFILE("MLNodeLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, a_factory);

    m_const_sigma = a_const_sigma;
    m_sigma.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_sigma[amrlev].resize(m_num_mg_levels[amrlev]);
        const int mglev = 0;
        const int idim = 0;
#ifdef AMREX_USE_EB
        bool allocate_sigma_mfs = true;
#else
        bool allocate_sigma_mfs = m_const_sigma == Real(0.0);
#endif
        if (allocate_sigma_mfs) {
            m_sigma[amrlev][mglev][idim] = std::make_unique<MultiFab>(m_grids[amrlev][mglev],
                                                                      m_dmap[amrlev][mglev],
                                                                      1, 1, MFInfo(),
                                                                      *m_factory[amrlev][0]);
            m_sigma[amrlev][mglev][idim]->setVal(m_const_sigma);
#ifdef AMREX_USE_EB
            m_sigma[amrlev][mglev][idim]->setDomainBndry(0.0, m_geom[amrlev][mglev]);
#endif
        }
    }

#ifdef AMREX_USE_EB
#if (AMREX_SPACEDIM == 2)
    const int ncomp_i = 5;
#else
    const int ncomp_i = algoim::numIntgs;
#endif
    m_integral.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
#ifdef AMREX_USE_EB
        m_integral[amrlev] = std::make_unique<MultiFab>(m_grids[amrlev][0],
                                                        m_dmap[amrlev][0],
                                                        ncomp_i, 1, MFInfo(),
                                                        *m_factory[amrlev][0]);
#else
        m_integral[amrlev] = std::make_unique<MultiFab>(m_grids[amrlev][0],
                                                        m_dmap[amrlev][0], ncomp_i, 1));
#endif
    }
#endif

#if (AMREX_SPACEDIM == 2)
    m_is_rz = m_geom[0][0].IsRZ();
#endif
}

#ifdef AMREX_USE_EB
void
MLNodeLaplacian::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<EBFArrayBoxFactory const*>& a_factory,
                         Real  a_const_sigma)
{
    Vector<FabFactory<FArrayBox> const*> _factory;
    for (auto x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }
    define(a_geom, a_grids, a_dmap, a_info, _factory, a_const_sigma);
}
#endif

void
MLNodeLaplacian::unimposeNeumannBC (int amrlev, MultiFab& rhs) const
{
    if (m_coarsening_strategy == CoarseningStrategy::RAP) {
        const Box& nddom = amrex::surroundingNodes(Geom(amrlev).Domain());
        const auto lobc = LoBC();
        const auto hibc = HiBC();

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs,mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& rhsarr = rhs.array(mfi);
            mlndlap_unimpose_neumann_bc(bx, rhsarr, nddom, lobc, hibc);
        }
    }
}

void
MLNodeLaplacian::setSigma (int amrlev, const MultiFab& a_sigma)
{
    AMREX_ALWAYS_ASSERT(m_sigma[amrlev][0][0]);
    MultiFab::Copy(*m_sigma[amrlev][0][0], a_sigma, 0, 0, 1, 0);
}

void
MLNodeLaplacian::FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom)
{
    BL_PROFILE("MLNodeLaplacian::FillBoundaryCoeff()");

    sigma.FillBoundary(geom.periodicity());

    if (m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
        const Box& domain = geom.Domain();
        const auto lobc = LoBC();
        const auto hibc = HiBC();

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sigma, mfi_info); mfi.isValid(); ++mfi)
        {
            Array4<Real> const& sfab = sigma.array(mfi);
            mlndlap_fillbc_cc<Real>(mfi.validbox(),sfab,domain,lobc,hibc);
        }
    }
}

void
MLNodeLaplacian::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
    if (!m_masks_built) buildMasks();

    const iMultiFab& cfmask = *m_nd_fine_mask[amrlev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(resmsk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<int> const& rmsk = resmsk.array(mfi);
        Array4<int const> const& fmsk = cfmask.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
        {
            if (fmsk(i,j,k) == crse_fine_node) rmsk(i,j,k) = 1;
        });
    }
}

void
MLNodeLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();

    averageDownCoeffs();

#ifdef AMREX_USE_EB
    buildIntegral();
#endif

    buildStencil();
}

void
MLNodeLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLNodeLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][cmglev-1];

    const auto& stencil = m_stencil[amrlev][cmglev-1];

    bool regular_coarsening = true; int idir = 2;
    if (cmglev > 0) {
        regular_coarsening = mg_coarsen_ratio_vec[cmglev-1] == mg_coarsen_ratio;
        IntVect ratio = mg_coarsen_ratio_vec[cmglev-1];
        if (ratio[1] == 1) {
            idir = 1;
        } else if (ratio[0] == 1) {
            idir = 0;
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> cfab = pcrse->array(mfi);
        Array4<Real const> const& ffab = fine.const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::Sigma)
        {
            if (regular_coarsening)
            {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_restriction(i,j,k,cfab,ffab,mfab);
                });
            }
            else
            {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_semi_restriction(i,j,k,cfab,ffab,mfab,idir);
                });
            }
        }
        else
        {
            Array4<Real const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
            });
        }
    }

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }
}

void
MLNodeLaplacian::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    BL_PROFILE("MLNodeLaplacian::interpolation()");

    const auto& sigma = m_sigma[amrlev][fmglev];
    const auto& stencil = m_stencil[amrlev][fmglev];

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    const MultiFab* cmf = &crse;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
        cfine.ParallelCopy(crse);
        cmf = &cfine;
    }

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][fmglev];

    bool regular_coarsening = true; int idir = 2;
    if (fmglev > 0) {
        regular_coarsening = mg_coarsen_ratio_vec[fmglev] == mg_coarsen_ratio;
        IntVect ratio = mg_coarsen_ratio_vec[fmglev];
        if (ratio[1] == 1) {
            idir = 1;
        } else if (ratio[0] == 1) {
            idir = 0;
        }
    }
    if (sigma[0] == nullptr) {
        AMREX_ALWAYS_ASSERT(regular_coarsening);
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& ffab = fine.array(mfi);
        Array4<Real const> const& cfab = cmf->const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            Array4<Real const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_interpadd_rap(i,j,k,ffab,cfab,stfab,mfab);
            });
        }
        else if (sigma[0] == nullptr)
        {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_interpadd_c(i,j,k,ffab,cfab,mfab);
            });
        }
        else if (m_use_harmonic_average && fmglev > 0)
        {
            AMREX_D_TERM(Array4<Real const> const& sxfab = sigma[0]->const_array(mfi);,
                         Array4<Real const> const& syfab = sigma[1]->const_array(mfi);,
                         Array4<Real const> const& szfab = sigma[2]->const_array(mfi););
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_interpadd_ha(i,j,k,ffab,cfab,AMREX_D_DECL(sxfab,syfab,szfab),mfab);
            });
        }
        else
        {
            Array4<Real const> const& sfab = sigma[0]->const_array(mfi);
            if (regular_coarsening)
            {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_interpadd_aa(i,j,k,ffab,cfab,sfab,mfab);
                });
            }
            else
            {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_semi_interpadd_aa(i,j,k,ffab,cfab,sfab,mfab,idir);
                });
            }
        }
    }
}

void
MLNodeLaplacian::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
    const auto& amrrr = AMRRefRatio(camrlev);
    amrex::average_down(fine_sol, crse_sol, 0, 1, amrrr);

    if (isSingular(0))
    {
        MultiFab frhs(fine_rhs.boxArray(), fine_rhs.DistributionMap(), 1, amrrr-1);
        MultiFab::Copy(frhs, fine_rhs, 0, 0, 1, 0);
        restrictInteriorNodes(camrlev, crse_rhs, frhs);
    }
}

void
MLNodeLaplacian::restrictInteriorNodes (int camrlev, MultiFab& crhs, MultiFab& a_frhs) const
{
    const BoxArray& fba = a_frhs.boxArray();
    const DistributionMapping& fdm = a_frhs.DistributionMap();
    const int amrrr = AMRRefRatio(camrlev);

    MultiFab* frhs = nullptr;
    std::unique_ptr<MultiFab> mf;
    if (a_frhs.nGrowVect().allGE(IntVect(amrrr-1)))
    {
        frhs = &a_frhs;
    }
    else
    {
        mf = std::make_unique<MultiFab>(fba, fdm, 1, amrrr-1);
        frhs = mf.get();
        MultiFab::Copy(*frhs, a_frhs, 0, 0, 1, 0);
    }

    const Geometry& cgeom = m_geom[camrlev  ][0];
    const Geometry& fgeom = m_geom[camrlev+1][0];

    const Box& f_nd_domain = amrex::surroundingNodes(fgeom.Domain());

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    const iMultiFab& fdmsk = *m_dirichlet_mask[camrlev+1][0];
    const auto& stencil    =  m_stencil[camrlev+1][0];

    MultiFab cfine(amrex::coarsen(fba, amrrr), fdm, 1, 0);

    frhs->setBndry(0.0);
    frhs->FillBoundary(fgeom.periodicity());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cfine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& cfab = cfine.array(mfi);
        Array4<Real const> const& ffab = frhs->const_array(mfi);
        Array4<int const> const& mfab = fdmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::Sigma) {
            if (amrrr == 2) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_restriction<2>(i,j,k,cfab,ffab,mfab,f_nd_domain,lobc,hibc);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_restriction<4>(i,j,k,cfab,ffab,mfab,f_nd_domain,lobc,hibc);
                });
            }
        } else {
            Array4<Real const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
            });
        }
    }

    MultiFab tmp_crhs(crhs.boxArray(), crhs.DistributionMap(), 1, 0);
    tmp_crhs.setVal(0.0);
    tmp_crhs.ParallelCopy(cfine, cgeom.periodicity());

    const iMultiFab& c_nd_mask = *m_nd_fine_mask[camrlev];
    const auto& has_fine_bndry = *m_has_fine_bndry[camrlev];

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(crhs, mfi_info); mfi.isValid(); ++mfi)
    {
        if (has_fine_bndry[mfi])
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& dfab = crhs.array(mfi);
            Array4<Real const> const& sfab = tmp_crhs.const_array(mfi);
            Array4<int const> const& mfab = c_nd_mask.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                if (mfab(i,j,k) == fine_node) dfab(i,j,k) = sfab(i,j,k);
            });
        }
    }
}

void
MLNodeLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    BL_PROFILE("MLNodeLaplacian::normalize()");

    if (m_sigma[0][0][0] == nullptr) return;

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& stencil = m_stencil[amrlev][mglev];
    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];
    const Real s0_norm0 = m_s0_norm0[amrlev][mglev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& arr = mf.array(mfi);
        Array4<int const> const& dmskarr = dmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            Array4<Real const> const& stenarr = stencil->const_array(mfi);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndlap_normalize_sten(tbx,arr,stenarr,dmskarr,s0_norm0);
            });
        }
        else if (m_use_harmonic_average && mglev > 0)
        {
            AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                         Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                         Array4<Real const> const& szarr = sigma[2]->const_array(mfi););

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndlap_normalize_ha(tbx,arr,AMREX_D_DECL(sxarr,syarr,szarr),dmskarr,dxinv);
            });
        }
        else
        {
            Array4<Real const> const& sarr = sigma[0]->const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndlap_normalize_aa(tbx,arr,sarr,dmskarr,dxinv);
            });
        }
    }
}

void
MLNodeLaplacian::checkPoint (std::string const& file_name) const
{
    if (ParallelContext::IOProcessorSub())
    {
        UtilCreateCleanDirectory(file_name, false);
        {
            std::string HeaderFileName(file_name+"/Header");
            std::ofstream HeaderFile;
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) {
                FileOpenFailed(HeaderFileName);
            }

            HeaderFile.precision(17);

            // MLLinop stuff
            HeaderFile << "verbose = " << verbose << "\n"
                       << "nlevs = " << NAMRLevels() << "\n"
                       << "do_agglomeration = " << info.do_agglomeration << "\n"
                       << "do_consolidation = " << info.do_consolidation << "\n"
                       << "agg_grid_size = " << info.agg_grid_size << "\n"
                       << "con_grid_size = " << info.con_grid_size << "\n"
                       << "has_metric_term = " << info.has_metric_term << "\n"
                       << "max_coarsening_level = " << info.max_coarsening_level << "\n";
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0])
                       << " "       << static_cast<int>(m_lobc[0][1]) << "\n";
#else
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0])
                       << " "       << static_cast<int>(m_lobc[0][1])
                       << " "       << static_cast<int>(m_lobc[0][2]) << "\n";
#endif
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0])
                       << " "       << static_cast<int>(m_hibc[0][1]) << "\n";
#else
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0])
                       << " "       << static_cast<int>(m_hibc[0][1])
                       << " "       << static_cast<int>(m_hibc[0][2]) << "\n";
#endif
            // m_coarse_data_for_bc: not used
            HeaderFile << "maxorder = " << getMaxOrder() << "\n";

            // MLNodeLaplacian stuff
            HeaderFile << "is_rz = " << m_is_rz << "\n";
            HeaderFile << "m_const_sigma = " << m_const_sigma << "\n";
            HeaderFile << "use_gauss_seidel = " << m_use_gauss_seidel << "\n";
            HeaderFile << "use_harmonic_average = " << m_use_harmonic_average << "\n";
            HeaderFile << "coarsen_strategy = " << static_cast<int>(m_coarsening_strategy) << "\n";
            // No level bc multifab
        }

        for (int ilev = 0; ilev < NAMRLevels(); ++ilev)
        {
            UtilCreateCleanDirectory(file_name+"/Level_"+std::to_string(ilev), false);
            std::string HeaderFileName(file_name+"/Level_"+std::to_string(ilev)+"/Header");
            std::ofstream HeaderFile;
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) {
                FileOpenFailed(HeaderFileName);
            }

            HeaderFile.precision(17);

            HeaderFile << Geom(ilev) << "\n";
            m_grids[ilev][0].writeOn(HeaderFile);  HeaderFile << "\n";
        }
    }

    ParallelContext::BarrierSub();

    for (int ilev = 0; ilev < NAMRLevels(); ++ilev)
    {
        if (m_sigma[ilev][0][0]) {
            VisMF::Write(*m_sigma[ilev][0][0], file_name+"/Level_"+std::to_string(ilev)+"/sigma");
        }
    }
}

}
