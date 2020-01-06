
#include <AMReX_MLCellABecLap.H>
#include <AMReX_MLLinOp_K.H>

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

    const Real betainv = 1.0 / getBScalar();
    const int nlevs = NAMRLevels();
    for (int alev = 0; alev < nlevs; ++alev) {
        compFlux(alev, a_flux[alev], *a_sol[alev], a_loc);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            unapplyMetricTerm(alev, 0, *a_flux[alev][idim]);
            if (betainv != 1.0) {
                a_flux[alev][idim]->mult(betainv);
            }
        }
    }
}

void
MLCellABecLap::applyInhomogNeumannTerm (int amrlev, MultiFab& rhs) const
{
    int ncomp = rhs.nComp();
    bool has_inhomog_neumann = false;
    for (int n = 0; n < ncomp; ++n)
    {
        auto itlo = std::find(m_lo_inhomog_neumann[n].begin(),
                              m_lo_inhomog_neumann[n].end(),   1);
        auto ithi = std::find(m_hi_inhomog_neumann[n].begin(),
                              m_hi_inhomog_neumann[n].end(),   1);
        if (itlo != m_lo_inhomog_neumann[n].end() or
            ithi != m_hi_inhomog_neumann[n].end())
        {
            has_inhomog_neumann = true;
        }
    }

    if (!has_inhomog_neumann) return;

    const int mglev = 0;

    const auto problo = m_geom[amrlev][mglev].ProbLoArray();
    const auto probhi = m_geom[amrlev][mglev].ProbHiArray();
    amrex::ignore_unused(probhi);
    const Real dxi = m_geom[amrlev][mglev].InvCellSize(0);
    const Real dyi = (AMREX_SPACEDIM >= 2) ? m_geom[amrlev][mglev].InvCellSize(1) : 1.0;
    const Real dzi = (AMREX_SPACEDIM == 3) ? m_geom[amrlev][mglev].InvCellSize(2) : 1.0;
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

#ifdef _OPENMP
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
            Array4<Real const> const bfab = (has_bcoef) ? bcoef[idim]->array(mfi) : foo.array();
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
            if ((!outside_domain_lo) and (!outside_domain_hi)) continue;
            for (int icomp = 0; icomp < ncomp; ++icomp) {
                const BoundCond bctlo = bdcv[icomp][olo];
                const BoundCond bcthi = bdcv[icomp][ohi];
                const Real bcllo = bdlv[icomp][olo];
                const Real bclhi = bdlv[icomp][ohi];
                if (m_lo_inhomog_neumann[icomp][idim] and outside_domain_lo)
                {
                    if (idim == 0) {
                        Real fac = beta*dxi;
                        if (m_has_metric_term and !has_bcoef) {
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
                        if (m_has_metric_term and !has_bcoef) {
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
                if (m_hi_inhomog_neumann[icomp][idim] and outside_domain_hi)
                {
                    if (idim == 0) {
                        Real fac = beta*dxi;
                        if (m_has_metric_term and !has_bcoef) {
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
                        if (m_has_metric_term and !has_bcoef) {
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

            }
        }

    }
}

#ifdef AMREX_USE_HYPRE
std::unique_ptr<Hypre>
MLCellABecLap::makeHypre (Hypre::Interface hypre_interface) const
{
    const BoxArray& ba = m_grids[0].back();
    const DistributionMapping& dm = m_dmap[0].back();
    const Geometry& geom = m_geom[0].back();
    const auto& factory = *(m_factory[0].back());
    MPI_Comm comm = BottomCommunicator();

    auto hypre_solver = amrex::makeHypre(ba, dm, geom, comm, hypre_interface);

    hypre_solver->setScalars(getAScalar(), getBScalar());

    const int mglev = NMGLevels(0)-1;
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
