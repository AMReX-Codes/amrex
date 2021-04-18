
#include <AMReX_MLCellABecLap.H>
#include <AMReX_MLLinOp_K.H>
#include <AMReX_MLCellABecLap_K.H>

#ifdef AMREX_USE_PETSC
#include <petscksp.h>
#include <AMReX_PETSc.H>
#endif

namespace amrex {

MLCellABecLap::MLCellABecLap ()
{
}

MLCellABecLap::~MLCellABecLap () {}

void
MLCellABecLap::define (const Vector<Geometry>& a_geom,
                       const Vector<BoxArray>& a_grids,
                       const Vector<DistributionMapping>& a_dmap,
                       const LPInfo& a_info,
                       const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    m_overset_mask.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        m_overset_mask[amrlev].resize(m_num_mg_levels[amrlev]);
    }
}

void
MLCellABecLap::define (const Vector<Geometry>& a_geom,
                       const Vector<BoxArray>& a_grids,
                       const Vector<DistributionMapping>& a_dmap,
                       const Vector<iMultiFab const*>& a_overset_mask,
                       const LPInfo& a_info,
                       const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLCellABecLap::define(overset)");

    AMREX_ALWAYS_ASSERT(!hasHiddenDimension());

    m_lpinfo_arg = a_info;

    int namrlevs = a_geom.size();
    m_overset_mask.resize(namrlevs);
    for (int amrlev = 0; amrlev < namrlevs; ++amrlev)
    {
        m_overset_mask[amrlev].push_back(std::make_unique<iMultiFab>(a_grids[amrlev],
                                                                     a_dmap[amrlev], 1, 1));
        iMultiFab::Copy(*m_overset_mask[amrlev][0], *a_overset_mask[amrlev], 0, 0, 1, 0);
        if (amrlev > 1) {
            AMREX_ALWAYS_ASSERT(amrex::refine(a_geom[amrlev-1].Domain(),2)
                                == a_geom[amrlev].Domain());
        }
    }

    int amrlev = 0;
    Box dom = a_geom[0].Domain();
    for (int mglev = 1; mglev <= a_info.max_coarsening_level; ++mglev)
    {
        AMREX_ALWAYS_ASSERT(mg_coarsen_ratio == 2);
        iMultiFab const& fine = *m_overset_mask[amrlev][mglev-1];
        if (dom.coarsenable(2) && fine.boxArray().coarsenable(2)) {
            dom.coarsen(2);
            auto crse = std::make_unique<iMultiFab>(amrex::coarsen(fine.boxArray(),2),
                                                    fine.DistributionMap(), 1, 1);
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*crse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<int const> const& fmsk = fine.const_array(mfi);
                Array4<int> const& cmsk = crse->array(mfi);
                reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_HOST_DEVICE (Box const& b) -> ReduceTuple
                {
                    return { coarsen_overset_mask(b, cmsk, fmsk) };
                });
            }
            ReduceTuple hv = reduce_data.value();
            if (amrex::get<0>(hv) == 0) {
                m_overset_mask[amrlev].push_back(std::move(crse));
            } else {
                break;
            }
        } else {
            break;
        }
    }
    int max_overset_mask_coarsening_level = m_overset_mask[amrlev].size()-1;
    ParallelAllReduce::Min(max_overset_mask_coarsening_level, ParallelContext::CommunicatorSub());
    m_overset_mask[amrlev].resize(max_overset_mask_coarsening_level+1);

    LPInfo linfo = a_info;
    linfo.max_coarsening_level = std::min(a_info.max_coarsening_level,
                                          max_overset_mask_coarsening_level);

    MLCellLinOp::define(a_geom, a_grids, a_dmap, linfo, a_factory);

    amrlev = 0;
    for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev) {
        MultiFab foo(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 0, MFInfo().SetAlloc(false));
        if (! isMFIterSafe(*m_overset_mask[amrlev][mglev], foo)) {
            auto osm = std::make_unique<iMultiFab>(m_grids[amrlev][mglev],
                                                   m_dmap[amrlev][mglev], 1, 1);
            osm->ParallelCopy(*m_overset_mask[amrlev][mglev]);
            std::swap(osm, m_overset_mask[amrlev][mglev]);
        }
    }

    for (amrlev = 1; amrlev < m_num_amr_levels; ++amrlev) {
        for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev) { // for ref_ratio 4
            m_overset_mask[amrlev].push_back(std::make_unique<iMultiFab>(m_grids[amrlev][mglev],
                                                                         m_dmap[amrlev][mglev],
                                                                         1, 1));
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*m_overset_mask[amrlev][mglev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<int> const& cmsk = m_overset_mask[amrlev][mglev]->array(mfi);
                Array4<int const> const fmsk = m_overset_mask[amrlev][mglev-1]->const_array(mfi);
                AMREX_LAUNCH_HOST_DEVICE_FUSIBLE_LAMBDA(bx, tbx,
                {
                    coarsen_overset_mask(tbx, cmsk, fmsk);
                });
            }
        }
    }

    for (amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            m_overset_mask[amrlev][mglev]->setBndry(1);
            m_overset_mask[amrlev][mglev]->FillBoundary(m_geom[amrlev][mglev].periodicity());
        }
    }
}

void
MLCellABecLap::prepareForSolve ()
{
    MLCellLinOp::prepareForSolve();
}

void
MLCellABecLap::update ()
{
    if (MLCellLinOp::needsUpdate()) MLCellLinOp::update();
}

void
MLCellABecLap::getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_flux,
                          const Vector<MultiFab*>& a_sol,
                          Location a_loc) const
{
    BL_PROFILE("MLMG::getFluxes()");

    const Real betainv = Real(1.0) / getBScalar();
    const int nlevs = NAMRLevels();
    for (int alev = 0; alev < nlevs; ++alev) {
        compFlux(alev, a_flux[alev], *a_sol[alev], a_loc);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            unapplyMetricTerm(alev, 0, *a_flux[alev][idim]);
            if (betainv != Real(1.0)) {
                a_flux[alev][idim]->mult(betainv);
            }
        }
    }
}

void
MLCellABecLap::applyInhomogNeumannTerm (int amrlev, MultiFab& rhs) const
{
    bool has_inhomog_neumann = hasInhomogNeumannBC();
    bool has_robin = hasRobinBC();

    if (!has_inhomog_neumann && !has_robin) return;

    int ncomp = getNComp();
    const int mglev = 0;

    const auto problo = m_geom[amrlev][mglev].ProbLoArray();
    const auto probhi = m_geom[amrlev][mglev].ProbHiArray();
    amrex::ignore_unused(probhi);
    const Real dxi = m_geom[amrlev][mglev].InvCellSize(0);
    const Real dyi = (AMREX_SPACEDIM >= 2) ? m_geom[amrlev][mglev].InvCellSize(1) : Real(1.0);
    const Real dzi = (AMREX_SPACEDIM == 3) ? m_geom[amrlev][mglev].InvCellSize(2) : Real(1.0);
    const Real xlo = problo[0];
    const Real dx = m_geom[amrlev][mglev].CellSize(0);
    const Box& domain = m_geom[amrlev][mglev].Domain();

    const Real beta = getBScalar();
    Array<MultiFab const*, AMREX_SPACEDIM> const& bcoef = getBCoeffs(amrlev,mglev);
    FArrayBox foo(Box(IntVect(0),IntVect(1)));
    bool has_bcoef = (bcoef[0] != nullptr);

    const auto& maskvals = m_maskvals[amrlev][mglev];
    const auto& bcondloc = *m_bcondloc[amrlev][mglev];
    const auto& bndry = *m_bndry_sol[amrlev];

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        Array4<Real> const& rhsfab = rhs.array(mfi);

        const auto & bdlv = bcondloc.bndryLocs(mfi);
        const auto & bdcv = bcondloc.bndryConds(mfi);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            Array4<Real const> const bfab = (has_bcoef)
                ? bcoef[idim]->const_array(mfi) : foo.const_array();
            const Orientation olo(idim,Orientation::low);
            const Orientation ohi(idim,Orientation::high);
            const Box blo = amrex::adjCellLo(vbx, idim);
            const Box bhi = amrex::adjCellHi(vbx, idim);
            const auto& mlo = maskvals[olo].array(mfi);
            const auto& mhi = maskvals[ohi].array(mfi);
            const auto& bvlo = bndry.bndryValues(olo).array(mfi);
            const auto& bvhi = bndry.bndryValues(ohi).array(mfi);
            bool outside_domain_lo = !(domain.contains(blo));
            bool outside_domain_hi = !(domain.contains(bhi));
            if ((!outside_domain_lo) && (!outside_domain_hi)) continue;
            for (int icomp = 0; icomp < ncomp; ++icomp) {
                const BoundCond bctlo = bdcv[icomp][olo];
                const BoundCond bcthi = bdcv[icomp][ohi];
                const Real bcllo = bdlv[icomp][olo];
                const Real bclhi = bdlv[icomp][ohi];
                if (m_lobc_orig[icomp][idim] == LinOpBCType::inhomogNeumann && outside_domain_lo)
                {
                    if (idim == 0) {
                        Real fac = beta*dxi;
                        if (m_has_metric_term && !has_bcoef) {
#if (AMREX_SPACEDIM == 1)
                            fac *= problo[0]*problo[0];
#elif (AMREX_SPACEDIM == 2)
                            fac *= problo[0];
#endif
                        }
                        AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                        {
                            mllinop_apply_innu_xlo(i,j,k, rhsfab, mlo, bfab,
                                                   bctlo, bcllo, bvlo,
                                                   fac, has_bcoef, icomp);
                        });
                    } else if (idim == 1) {
                        Real fac = beta*dyi;
                        if (m_has_metric_term && !has_bcoef) {
                            AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                            {
                                mllinop_apply_innu_ylo_m(i,j,k, rhsfab, mlo,
                                                         bctlo, bcllo, bvlo,
                                                         fac, xlo, dx, icomp);
                            });
                        }
                        else {
                            AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                            {
                                mllinop_apply_innu_ylo(i,j,k, rhsfab, mlo, bfab,
                                                       bctlo, bcllo, bvlo,
                                                       fac, has_bcoef, icomp);
                            });
                        }
                    } else {
                        Real fac = beta*dzi;
                        AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                        {
                            mllinop_apply_innu_zlo(i,j,k, rhsfab, mlo, bfab,
                                                   bctlo, bcllo, bvlo,
                                                   fac, has_bcoef, icomp);
                        });
                    }
                }
                if (m_hibc_orig[icomp][idim] == LinOpBCType::inhomogNeumann && outside_domain_hi)
                {
                    if (idim == 0) {
                        Real fac = beta*dxi;
                        if (m_has_metric_term && !has_bcoef) {
#if (AMREX_SPACEDIM == 1)
                            fac *= probhi[0]*probhi[0];
#elif (AMREX_SPACEDIM == 2)
                            fac *= probhi[0];
#endif
                        }
                        AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                        {
                            mllinop_apply_innu_xhi(i,j,k, rhsfab, mhi, bfab,
                                                   bcthi, bclhi, bvhi,
                                                   fac, has_bcoef, icomp);
                        });
                    } else if (idim == 1) {
                        Real fac = beta*dyi;
                        if (m_has_metric_term && !has_bcoef) {
                            AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                            {
                                mllinop_apply_innu_yhi_m(i,j,k, rhsfab, mhi,
                                                         bcthi, bclhi, bvhi,
                                                         fac, xlo, dx, icomp);
                            });
                        } else {
                            AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                            {
                                mllinop_apply_innu_yhi(i,j,k, rhsfab, mhi, bfab,
                                                       bcthi, bclhi, bvhi,
                                                       fac, has_bcoef, icomp);
                            });
                        }
                    } else {
                        Real fac = beta*dzi;
                        AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                        {
                            mllinop_apply_innu_zhi(i,j,k, rhsfab, mhi, bfab,
                                                   bcthi, bclhi, bvhi,
                                                   fac, has_bcoef, icomp);
                        });
                    }
                }

                if (has_robin) {
                    // For Robin BC, see comments in AMReX_MLABecLaplacian.cpp above
                    // function applyRobinBCTermsCoeffs.
                    Array4<Real const> const& rbc = (*m_robin_bcval[amrlev])[mfi].const_array(icomp*3);
                    if (m_lobc_orig[icomp][idim] == LinOpBCType::Robin && outside_domain_lo)
                    {
                        if (idim == 0) {
                            Real fac = beta*dxi*dxi;
                            AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                            {
                                Real A = rbc(i,j,k,2)
                                    / (rbc(i,j,k,1)*dxi + rbc(i,j,k,0)*Real(0.5));
                                rhsfab(i+1,j,k,icomp) += fac*bfab(i+1,j,k,icomp)*A;
                            });
                        } else if (idim == 1) {
                            Real fac = beta*dyi*dyi;
                            AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                            {
                                Real A = rbc(i,j,k,2)
                                    / (rbc(i,j,k,1)*dyi + rbc(i,j,k,0)*Real(0.5));
                                rhsfab(i,j+1,k,icomp) += fac*bfab(i,j+1,k,icomp)*A;
                            });
                        } else {
                            Real fac = beta*dzi*dzi;
                            AMREX_HOST_DEVICE_FOR_3D(blo, i, j, k,
                            {
                                Real A = rbc(i,j,k,2)
                                    / (rbc(i,j,k,1)*dzi + rbc(i,j,k,0)*Real(0.5));
                                rhsfab(i,j,k+1,icomp) += fac*bfab(i,j,k+1,icomp)*A;
                            });
                        }
                    }
                    if (m_hibc_orig[icomp][idim] == LinOpBCType::Robin && outside_domain_hi)
                    {
                        if (idim == 0) {
                            Real fac = beta*dxi*dxi;
                            AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                            {
                                Real A = rbc(i,j,k,2)
                                    / (rbc(i,j,k,1)*dxi + rbc(i,j,k,0)*Real(0.5));
                                rhsfab(i-1,j,k,icomp) += fac*bfab(i,j,k,icomp)*A;
                            });
                        } else if (idim == 1) {
                            Real fac = beta*dyi*dyi;
                            AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                            {
                                Real A = rbc(i,j,k,2)
                                    / (rbc(i,j,k,1)*dyi + rbc(i,j,k,0)*Real(0.5));
                                rhsfab(i,j-1,k,icomp) += fac*bfab(i,j,k,icomp)*A;
                            });
                        } else {
                            Real fac = beta*dzi*dzi;
                            AMREX_HOST_DEVICE_FOR_3D(bhi, i, j, k,
                            {
                                Real A = rbc(i,j,k,2)
                                    / (rbc(i,j,k,1)*dzi + rbc(i,j,k,0)*Real(0.5));
                                rhsfab(i,j,k-1,icomp) += fac*bfab(i,j,k,icomp)*A;
                            });
                        }
                    }
                }
            }
        }

    }
}

void
MLCellABecLap::applyOverset (int amrlev, MultiFab& rhs) const
{
    if (m_overset_mask[amrlev][0]) {
        const int ncomp = getNComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*m_overset_mask[amrlev][0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& rfab = rhs.array(mfi);
            Array4<int const> const& osm = m_overset_mask[amrlev][0]->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
            {
                if (osm(i,j,k) == 0) rfab(i,j,k,n) = 0.0;
            });
        }
    }
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
std::unique_ptr<Hypre>
MLCellABecLap::makeHypre (Hypre::Interface hypre_interface) const
{
    const BoxArray& ba = m_grids[0].back();
    const DistributionMapping& dm = m_dmap[0].back();
    const Geometry& geom = m_geom[0].back();
    const auto& factory = *(m_factory[0].back());
    MPI_Comm comm = BottomCommunicator();

    const int mglev = NMGLevels(0)-1;

    auto om = getOversetMask(0, mglev);

    auto hypre_solver = amrex::makeHypre(ba, dm, geom, comm, hypre_interface, om);

    hypre_solver->setScalars(getAScalar(), getBScalar());

    auto ac = getACoeffs(0, mglev);
    if (ac)
    {
        hypre_solver->setACoeffs(*ac);
    }
    else
    {
        MultiFab alpha(ba,dm,1,0,MFInfo(),factory);
        alpha.setVal(0.0);
        hypre_solver->setACoeffs(alpha);
    }

    auto bc = getBCoeffs(0, mglev);
    if (bc[0])
    {
        hypre_solver->setBCoeffs(bc);
    }
    else
    {
        Array<MultiFab,AMREX_SPACEDIM> beta;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            beta[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                              dm, 1, 0, MFInfo(), factory);
            beta[idim].setVal(1.0);
        }
        hypre_solver->setBCoeffs(amrex::GetArrOfConstPtrs(beta));
    }
    hypre_solver->setIsMatrixSingular(this->isBottomSingular());

    return hypre_solver;
}
#endif

#ifdef AMREX_USE_PETSC
std::unique_ptr<PETScABecLap>
MLCellABecLap::makePETSc () const
{
    const BoxArray& ba = m_grids[0].back();
    const DistributionMapping& dm = m_dmap[0].back();
    const Geometry& geom = m_geom[0].back();
    const auto& factory = *(m_factory[0].back());
    MPI_Comm comm = BottomCommunicator();

    auto petsc_solver = makePetsc(ba, dm, geom, comm);

    petsc_solver->setScalars(getAScalar(), getBScalar());

    const int mglev = NMGLevels(0)-1;
    auto ac = getACoeffs(0, mglev);
    if (ac)
    {
        petsc_solver->setACoeffs(*ac);
    }
    else
    {
        MultiFab alpha(ba,dm,1,0,MFInfo(),factory);
        alpha.setVal(0.0);
        petsc_solver->setACoeffs(alpha);
    }

    auto bc = getBCoeffs(0, mglev);
    if (bc[0])
    {
        petsc_solver->setBCoeffs(bc);
    }
    else
    {
        Array<MultiFab,AMREX_SPACEDIM> beta;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            beta[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                              dm, 1, 0, MFInfo(), factory);
            beta[idim].setVal(1.0);
        }
        petsc_solver->setBCoeffs(amrex::GetArrOfConstPtrs(beta));
    }
    return petsc_solver;
}
#endif

}
