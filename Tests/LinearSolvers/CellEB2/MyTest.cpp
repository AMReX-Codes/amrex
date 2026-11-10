#include "MyTest.H"
#include "MyTest_K.H"

#include <AMReX_MLEBABecLap.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#ifdef AMREX_USE_HYPRE
#include <AMReX_HypreMLABecLap.H>
#endif

#include <cmath>
#include <limits>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEB();

    addFineGrids();

    initData();
}

void
MyTest::solve ()
{
    if (verbose > 0) {
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            const MultiFab& vfrc = factory[ilev]->getVolFrac();
            MultiFab v(vfrc.boxArray(), vfrc.DistributionMap(), 1, 0,
                       MFInfo(), *factory[ilev]);
            MultiFab::Copy(v, vfrc, 0, 0, 1, 0);
            amrex::EB_set_covered(v, 1.0);
            amrex::Print() << "Level " << ilev << ": vfrc min = " << v.min(0) << '\n';
        }
    }

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom[0].isPeriodic(idim)) {
            mlmg_lobc[idim] = LinOpBCType::Periodic;
            mlmg_hibc[idim] = LinOpBCType::Periodic;
        } else if (domain_bc == "neumann") {
            mlmg_lobc[idim] = LinOpBCType::Neumann;
            mlmg_hibc[idim] = LinOpBCType::Neumann;
        } else {
            mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            mlmg_hibc[idim] = LinOpBCType::Dirichlet;
        }
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setAgglomerationGridSize(agg_grid_size);
    info.setConsolidationGridSize(con_grid_size);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        phi[ilev].setVal(0.0,0,1,IntVect(0));
    }

    static int ipass = 0;
    ++ipass;

    if (composite_solve)
    {
        BL_PROFILE_REGION("COMPOSITE_SOLVE-pass"+std::to_string(ipass));

        MLEBABecLap mleb (geom, grids, dmap, info, amrex::GetVecOfConstPtrs(factory));
        mleb.setMaxOrder(linop_maxorder);

        mleb.setDomainBC(mlmg_lobc, mlmg_hibc);

        for (int ilev = 0; ilev <= max_level; ++ilev) {
            mleb.setLevelBC(ilev, &phi[ilev]);
        }

        mleb.setScalars(scalars[0], scalars[1]);

        for (int ilev = 0; ilev <= max_level; ++ilev) {
            mleb.setACoeffs(ilev, acoef[ilev]);
            mleb.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcoef[ilev]));
        }

        if (!eb_neumann) {
            for (int ilev = 0; ilev <= max_level; ++ilev) {
                mleb.setEBDirichlet(ilev, phieb[ilev], bcoef_eb[ilev]);
            }
        }
        if (use_mlhypre) { mleb.setInterpBndryHalfWidth(1); }

        MLMG mlmg(mleb);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setBottomMaxIter(max_bottom_iter);
        mlmg.setBottomTolerance(bottom_reltol);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
        if (use_hypre) {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
        } else if (use_petsc) {
            mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
        }

        // The boundary values, which the solve overwrites
        Vector<MultiFab> bcdata(max_level+1);
        for (int ilev = 0; ilev <= max_level && use_mlhypre; ++ilev) {
            bcdata[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
            MultiFab::Copy(bcdata[ilev], phi[ilev], 0, 0, 1, 1);
        }

        const Real tol_rel = reltol;
        const Real tol_abs = 0.0;
        mlmg.solve(amrex::GetVecOfPtrs(phi), amrex::GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

        // Test getEBFluxes
        mleb.getEBFluxes(amrex::GetVecOfPtrs(fluxeb_sol), amrex::GetVecOfPtrs(phi));

#ifdef AMREX_USE_HYPRE
        if (use_mlhypre) { compareWithMLHypre(mleb, mlmg, bcdata, mlmg_lobc, mlmg_hibc); }
#endif
    }
    else
    {
        for (int ilev = 0; ilev <= max_level; ++ilev)
        {
            BL_PROFILE_REGION("LEVEL-SOLVE-lev"+std::to_string(ilev)+"-pass"+std::to_string(ipass));

            MLEBABecLap mleb({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info, {factory[ilev].get()});
            mleb.setMaxOrder(linop_maxorder);

            mleb.setDomainBC(mlmg_lobc, mlmg_hibc);

            if (ilev > 0) {
                mleb.setCoarseFineBC(&phi[ilev-1], ref_ratio);
            }
            mleb.setLevelBC(0, &phi[ilev]);

            mleb.setScalars(scalars[0], scalars[1]);

            mleb.setACoeffs(0, acoef[ilev]);
            mleb.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoef[ilev]));

            if (!eb_neumann) {
                mleb.setEBDirichlet(0, phieb[ilev], bcoef_eb[ilev]);
            }

            MLMG mlmg(mleb);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setBottomMaxIter(max_bottom_iter);
            mlmg.setBottomTolerance(bottom_reltol);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);
            if (use_hypre) {
                mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
            } else if (use_petsc) {
                mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
            }

            const Real tol_rel = reltol;
            const Real tol_abs = 0.0;
            mlmg.solve({&phi[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);

            // Test getEBFluxes
            mleb.getEBFluxes({&fluxeb_sol[ilev]}, {&phi[ilev]});
        }
    }

    if (verbose > 0) {
        for (int ilev = 0; ilev <= max_level; ++ilev)
        {
            MultiFab mf(phi[ilev].boxArray(),phi[ilev].DistributionMap(), 1, 0);
            MultiFab::Copy(mf,phi[ilev],0,0,1,0);
            MultiFab::Subtract(mf, phiexact[ilev], 0, 0, 1, 0);

            const MultiFab& vfrc = factory[ilev]->getVolFrac();

            MultiFab::Multiply(mf, vfrc, 0, 0, 1, 0);

            Real norminf = mf.norm0();
            Real norm1 = mf.norm1()*AMREX_D_TERM((1.0/n_cell), *(1.0/n_cell), *(1.0/n_cell));
            amrex::Print() << "Level " << ilev << ": weighted max and 1 norms " << norminf << ", " << norm1 << '\n';
        }
        for (int ilev = 0; ilev <= max_level; ++ilev)
        {
            const MultiFab& barea = factory[ilev]->getBndryArea().ToMultiFab(0.0, 0.0);
            const Real nxyzi = AMREX_D_TERM((1.0/n_cell), *(1.0/n_cell), *(1.0/n_cell));

            MultiFab mf_feb_err(fluxeb_sol[ilev].boxArray(), fluxeb_sol[ilev].DistributionMap(), 1, 0);
            MultiFab::Copy(mf_feb_err, fluxeb_sol[ilev], 0, 0, 1, 0);
            MultiFab::Subtract(mf_feb_err, fluxeb_exact[ilev], 0, 0, 1, 0);

            amrex::Print() << "Level " << ilev << ": EB flux error norms "
                << mf_feb_err.norm0() << ", " << mf_feb_err.norm1()*nxyzi << '\n';

            MultiFab::Multiply(mf_feb_err, barea, 0, 0, 1, 0);

            amrex::Print() << "Level " << ilev << ": EB flux error weighted norms "
                << mf_feb_err.norm0() << ", " << mf_feb_err.norm1()*nxyzi << '\n';
        }
    }
}

#ifdef AMREX_USE_HYPRE
void
MyTest::compareWithMLHypre (MLEBABecLap& mleb, MLMG& mlmg, Vector<MultiFab> const& bcdata,
                            Array<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                            Array<LinOpBCType,AMREX_SPACEDIM> const& hibc)
{
    int const nlevels = max_level+1;
    bool const singular = mleb.isSingular(0);

    // A singular problem gets MLMG's consistent rhs
    Vector<MultiFab> hrhs(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        hrhs[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        MultiFab::Copy(hrhs[ilev], rhs[ilev], 0, 0, 1, 0);
    }
    if (singular) {
        for (int ilev = max_level; ilev > 0; --ilev) {
            amrex::EB_average_down(hrhs[ilev], hrhs[ilev-1], 0, 1, ref_ratio);
        }
        auto const offset = mleb.getSolvabilityOffset(0, 0, hrhs[0]);
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            mleb.fixSolvabilityByOffset(ilev, 0, hrhs[ilev], offset);
        }
    }

    Vector<FabFactory<FArrayBox> const*> fact;
    for (auto const& f : factory) { fact.push_back(f.get()); }
    HypreMLABecLap hmleb(geom, grids, dmap, HypreSolverID::BoomerAMG, fact);
    hmleb.setVerbose(verbose);
    hmleb.setIsSingular(singular);
    hmleb.setReuseSolverSetup(mlhypre_num_setups > 1);
    if (!eb_neumann) {
        hmleb.setEBDirichlet(amrex::GetVecOfConstPtrs(phieb), amrex::GetVecOfConstPtrs(bcoef_eb));
    }
    Vector<Array<MultiFab const*,AMREX_SPACEDIM>> b(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        if (null_fine_bcoef == 1 && ilev > 0) {
            b[ilev].fill(nullptr);
        } else {
            b[ilev] = amrex::GetArrOfConstPtrs(bcoef[ilev]);
        }
    }

    Vector<MultiFab> hphi(nlevels), srhs(nlevels), res(nlevels), diff(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        hphi[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        MultiFab::Copy(hphi[ilev], bcdata[ilev], 0, 0, 1, 1);
        srhs[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        res[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        diff[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
    }
    for (int isetup = 0; isetup < mlhypre_num_setups; ++isetup) {
        // Scaled operator and rhs, same solution
        Real const s = isetup+1;
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            MultiFab::Copy(srhs[ilev], hrhs[ilev], 0, 0, 1, 0);
            srhs[ilev].mult(s);
            hphi[ilev].setVal(0.0, 0, 1, 0);
        }
        if (nan_in_covered) {
            // The cells the system does not use (EB-covered, covered by a
            // finer level) must not matter, whatever they hold
            Real const nan = std::numeric_limits<Real>::quiet_NaN();
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                amrex::EB_set_covered(srhs[ilev], 0, 1, 0, nan);
                amrex::EB_set_covered(hphi[ilev], 0, 1, 0, nan);
                if (ilev == max_level) { continue; }
                auto const fmask = amrex::makeFineMask(grids[ilev], dmap[ilev], grids[ilev+1],
                                                       IntVect(ref_ratio), 0, 1);
                for (MFIter mfi(srhs[ilev]); mfi.isValid(); ++mfi) {
                    auto const& m = fmask.const_array(mfi);
                    auto const& r = srhs[ilev].array(mfi);
                    auto const& x = hphi[ilev].array(mfi);
                    amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (m(i,j,k)) { r(i,j,k) = nan; x(i,j,k) = nan; }
                    });
                }
            }
        }
        Vector<MultiFab const*> acoef_ptrs = amrex::GetVecOfConstPtrs(acoef);
        if (null_acoef) { for (auto& p : acoef_ptrs) { p = nullptr; } }
        hmleb.setup(s*scalars[0], s*scalars[1], acoef_ptrs, b,
                    lobc, hibc, amrex::GetVecOfConstPtrs(bcdata));
        hmleb.solve(amrex::GetVecOfPtrs(hphi), amrex::GetVecOfConstPtrs(srhs), reltol, 0.0);

        // Weighted by the volume fraction, as the rows of HypreMLABecLap's
        // system (res), and by its square root (ress), as the rows of the
        // scaled system the solver converges: an error in a small cell's row
        // appears there as sqrt(vfrac) times the right-hand side, not vfrac
        // times it
        mlmg.compResidual(amrex::GetVecOfPtrs(res), amrex::GetVecOfPtrs(hphi),
                          amrex::GetVecOfConstPtrs(hrhs));
        Vector<MultiFab> ress(nlevels);
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            auto const& vfrac = factory[ilev]->getVolFrac();
            ress[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
            MultiFab::Copy(ress[ilev], vfrac, 0, 0, 1, 0);
            for (MFIter mfi(ress[ilev]); mfi.isValid(); ++mfi) {
                auto const& a = ress[ilev].array(mfi);
                amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    a(i,j,k) = std::sqrt(a(i,j,k));
                });
            }
            MultiFab::Multiply(ress[ilev], res[ilev], 0, 0, 1, 0);
            amrex::EB_set_covered(ress[ilev], 0.0);
            MultiFab::Multiply(res[ilev], vfrac, 0, 0, 1, 0);
        }
        // Zero the cells covered by the finer level and their face neighbors:
        // there the composite operators differ by design (MLMG's reflux is
        // not volume weighted)
        auto zero_near_cf = [&] (MultiFab& mf, int ilev)
        {
            auto const fmask = amrex::makeFineMask(grids[ilev], dmap[ilev], IntVect(1), grids[ilev+1],
                                                   IntVect(ref_ratio), geom[ilev].periodicity(), 0, 1);
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                auto const& m = fmask.const_array(mfi);
                auto const& r = mf.array(mfi);
                amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    IntVect const iv(AMREX_D_DECL(i,j,k));
                    bool near = m(iv);
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        IntVect const ed = IntVect::TheDimensionVector(d);
                        near = near || m(iv-ed) || m(iv+ed);
                    }
                    if (near) { r(i,j,k) = 0.0; }
                });
            }
        };
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            MultiFab::LinComb(diff[ilev], 1.0, hphi[ilev], 0, -1.0, phi[ilev], 0, 0, 1, 0);
        }
        if (singular) { // up to a constant
            auto const& vfrac = factory[max_level]->getVolFrac();
            Real const mean = MultiFab::Dot(diff[max_level], 0, vfrac, 0, 1, 0)
                / vfrac.sum(0);
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                diff[ilev].plus(-mean, 0, 1);
                amrex::EB_set_covered(diff[ilev], 0.0);
            }
        }
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            Real const dmax = diff[ilev].norm0();
            Real const rmax = res[ilev].norm0();
            amrex::Print() << "Setup " << isetup << ", level " << ilev
                           << ": max |phi_hypre - phi_mlmg| = " << dmax
                           << ", max |MLMG residual of phi_hypre| * vfrac = " << rmax << '\n';
            // On a coarser level the two operators differ only in the cells
            // next to the coarse/fine boundary: elsewhere the residual is
            // checked like the finest level's, which cannot see an error in
            // the coarse rows
            Real rmaxs = ress[ilev].norm0();
            if (ilev < max_level) {
                MultiFab away(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
                MultiFab::Copy(away, res[ilev], 0, 0, 1, 0);
                zero_near_cf(away, ilev);
                Real const amax = away.norm0();
                MultiFab::Copy(away, ress[ilev], 0, 0, 1, 0);
                zero_near_cf(away, ilev);
                rmaxs = away.norm0();
                amrex::Print() << "    away from the coarse/fine boundary: max |MLMG residual of phi_hypre| * vfrac = "
                               << amax << '\n';
                if (mlhypre_max_fine_res >= 0.0 && !(amax <= mlhypre_max_fine_res)) {
                    amrex::Abort("mlhypre: a coarser level's residual away from the coarse/fine boundary exceeds mlhypre_max_fine_res");
                }
            }
            amrex::Print() << "    scaled rows" << ((ilev < max_level) ? " (away from the c/f boundary)" : "")
                           << ": max |MLMG residual of phi_hypre| * sqrt(vfrac) = " << rmaxs << '\n';
            if (mlhypre_max_res_scaled >= 0.0 && !(rmaxs <= mlhypre_max_res_scaled)) {
                amrex::Abort("mlhypre: the residual of the scaled rows exceeds mlhypre_max_res_scaled");
            }
            // Against the exact solution: an error in the rows (which the
            // residual checks above may not see) makes HypreMLABecLap's error
            // orders of magnitude larger than MLMG's, while the operators'
            // difference at the coarse/fine boundary is below the
            // discretization error; not for singular problems (constant)
            if (!singular) {
                MultiFab err(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
                MultiFab::LinComb(err, 1.0, hphi[ilev], 0, -1.0, phiexact[ilev], 0, 0, 1, 0);
                MultiFab::Multiply(err, factory[ilev]->getVolFrac(), 0, 0, 1, 0);
                Real const emax_h = err.norm0();
                MultiFab::LinComb(err, 1.0, phi[ilev], 0, -1.0, phiexact[ilev], 0, 0, 1, 0);
                MultiFab::Multiply(err, factory[ilev]->getVolFrac(), 0, 0, 1, 0);
                Real const emax_m = err.norm0();
                amrex::Print() << "    max |phi - exact| * vfrac: HypreMLABecLap " << emax_h
                               << ", MLMG " << emax_m << '\n';
                if (mlhypre_max_exact_ratio >= 0.0 && !(emax_h <= mlhypre_max_exact_ratio*emax_m)) {
                    amrex::Abort("mlhypre: the error against the exact solution exceeds mlhypre_max_exact_ratio times MLMG's");
                }
                if (mlhypre_max_exact_err >= 0.0 && !(emax_h <= mlhypre_max_exact_err)) {
                    amrex::Abort("mlhypre: the error against the exact solution exceeds mlhypre_max_exact_err");
                }
            }
            if (verbose > 0) {
                // Where: the cell of the largest residual, its volume fraction
                // and whether it is at the domain boundary or at the
                // coarse/fine boundary
                MultiFab absres(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
                MultiFab::Copy(absres, res[ilev], 0, 0, 1, 0);
                absres.abs(0, 1);
                IntVect const iv = absres.maxIndex(0);
                Real vf = -1.0;
                auto const& vfrac = factory[ilev]->getVolFrac();
                for (MFIter mfi(vfrac); mfi.isValid(); ++mfi) {
                    if (mfi.validbox().contains(iv)) {
                        Gpu::dtoh_memcpy(&vf, vfrac[mfi].dataPtr() + vfrac[mfi].box().index(iv), sizeof(Real));
                    }
                }
                ParallelDescriptor::ReduceRealMax(vf);
                Box const& domain = geom[ilev].Domain();
                bool const at_wall = !Box(domain).grow(-1).contains(iv);
                bool at_cf = false;
                if (ilev > 0) {
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        for (int side = -1; side <= 1; side += 2) {
                            IntVect nb = iv; nb[idim] += side;
                            at_cf = at_cf || (domain.contains(nb) && !grids[ilev].contains(nb));
                        }
                    }
                }
                amrex::Print() << "    largest residual at cell " << iv << ", vfrac " << vf
                               << (at_wall ? ", at the domain boundary" : "")
                               << (at_cf ? ", at the coarse/fine boundary" : "") << '\n';
            }
            if (mlhypre_max_diff >= 0.0 && !(dmax <= mlhypre_max_diff)) {
                amrex::Abort("mlhypre: the solutions differ by more than mlhypre_max_diff");
            }
            if (mlhypre_max_fine_res >= 0.0 && ilev == max_level && !(rmax <= mlhypre_max_fine_res)) {
                amrex::Abort("mlhypre: the finest level's residual exceeds mlhypre_max_fine_res");
            }
        }
        // The coarse cells under a finer level hold the EB average of the
        // fine solution (AMReX's), not just the average-down rows' solve
        for (int ilev = 0; ilev < max_level; ++ilev) {
            MultiFab avg(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
            MultiFab::Copy(avg, hphi[ilev], 0, 0, 1, 0);
            amrex::EB_average_down(hphi[ilev+1], avg, 0, 1, ref_ratio);
            MultiFab::Subtract(avg, hphi[ilev], 0, 0, 1, 0);
            Real const amax = avg.norm0();
            amrex::Print() << "Setup " << isetup << ", level " << ilev
                           << ": max |phi_hypre - EB average of the finer phi_hypre| = " << amax << '\n';
            if (mlhypre_max_avg_diff >= 0.0 && !(amax <= mlhypre_max_avg_diff)) {
                amrex::Abort("mlhypre: the coarse solution is not the average of the fine one");
            }
        }
    }
}
#endif

void
MyTest::writePlotfile ()
{
    ParmParse pp;
    bool gpu_regtest = false;
#ifdef AMREX_USE_GPU
    pp.query("gpu_regtest", gpu_regtest);
#endif
    Vector<MultiFab> plotmf(max_level+1);
    if (gpu_regtest) {
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            const MultiFab& vfrc = factory[ilev]->getVolFrac();
            const MultiFab& barea = factory[ilev]->getBndryArea().ToMultiFab(0.0, 0.0);
            plotmf[ilev].define(grids[ilev],dmap[ilev], 6, 0);
            MultiFab::Copy(plotmf[ilev], phi[ilev], 0, 0, 1, 0);
            MultiFab::Copy(plotmf[ilev], phiexact[ilev], 0, 1, 1, 0);
            MultiFab::Copy(plotmf[ilev], vfrc, 0, 2, 1, 0);
            MultiFab::Copy(plotmf[ilev], fluxeb_sol[ilev], 0, 3, 1, 0);
            MultiFab::Copy(plotmf[ilev], fluxeb_exact[ilev], 0, 4, 1, 0);
            MultiFab::Copy(plotmf[ilev], barea, 0, 5, 1, 0);
        }
        WriteMultiLevelPlotfile(plot_file_name, max_level+1,
                                amrex::GetVecOfConstPtrs(plotmf),
                                {"phi","exact","vfrac","eb flux","eb flux exact","barea"},
                                geom, 0.0, Vector<int>(max_level+1,0),
                                Vector<IntVect>(max_level,IntVect{ref_ratio}));
    } else {
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            const MultiFab& vfrc = factory[ilev]->getVolFrac();
            const MultiFab& barea = factory[ilev]->getBndryArea().ToMultiFab(0.0, 0.0);
            plotmf[ilev].define(grids[ilev],dmap[ilev], 10, 0);

            MultiFab::Copy(plotmf[ilev], phi[ilev], 0, 0, 1, 0);

            MultiFab::Copy(plotmf[ilev], phiexact[ilev], 0, 1, 1, 0);

            MultiFab::Copy(plotmf[ilev], phi[ilev], 0, 2, 1, 0);
            MultiFab::Subtract(plotmf[ilev], phiexact[ilev], 0, 2, 1, 0);

            MultiFab::Copy(plotmf[ilev], phi[ilev], 0, 3, 1, 0);
            MultiFab::Subtract(plotmf[ilev], phiexact[ilev], 0, 3, 1, 0);
            MultiFab::Multiply(plotmf[ilev], vfrc, 0, 3, 1, 0);

            MultiFab::Copy(plotmf[ilev], vfrc, 0, 4, 1, 0);

            // EB fluxes
            MultiFab::Copy(plotmf[ilev], fluxeb_sol[ilev], 0, 5, 1, 0);

            MultiFab::Copy(plotmf[ilev], fluxeb_exact[ilev], 0, 6, 1, 0);

            MultiFab::Copy(    plotmf[ilev], fluxeb_sol[ilev], 0, 7, 1, 0);
            MultiFab::Subtract(plotmf[ilev], fluxeb_exact[ilev], 0, 7, 1, 0);

            MultiFab::Copy(    plotmf[ilev], fluxeb_sol[ilev], 0, 8, 1, 0);
            MultiFab::Subtract(plotmf[ilev], fluxeb_exact[ilev], 0, 8, 1, 0);
            MultiFab::Multiply(plotmf[ilev], barea, 0, 8, 1, 0);

            MultiFab::Copy(plotmf[ilev], barea, 0, 9, 1, 0);
        }
        WriteMultiLevelPlotfile(plot_file_name, max_level+1,
                                amrex::GetVecOfConstPtrs(plotmf),
                                {
                                    "phi", "exact", "error", "error*vfrac", "vfrac",
                                    "eb flux", "eb flux exact", "flux error", "(flux error)*barea", "barea"
                                },
                                geom, 0.0, Vector<int>(max_level+1,0),
                                Vector<IntVect>(max_level,IntVect{ref_ratio}));
    }

}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);
    pp.query("prob_type", prob_type);
    pp.query("ref_ratio", ref_ratio);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ref_ratio == 2 || ref_ratio == 4,
                                     "ref_ratio must be either 2 or 4");

    pp.query("plot_file", plot_file_name);

    scalars.resize(2);
    scalars[0] = 1.0;
    scalars[1] = 1.0;

    pp.queryarr("scalars", scalars);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("max_bottom_iter", max_bottom_iter);
    pp.query("bottom_reltol", bottom_reltol);
    pp.query("reltol", reltol);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("max_coarsening_level", max_coarsening_level);
#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
#endif
#ifdef AMREX_USE_PETSC
    pp.query("use_petsc", use_petsc);
#endif
    pp.query("agg_grid_size", agg_grid_size);
    pp.query("con_grid_size", con_grid_size);

    pp.query("composite_solve", composite_solve);
    pp.query("eb_neumann", eb_neumann);
    pp.query("acoef", acoef_value);
    pp.query("ghost_value", ghost_value);
    // Exactly AMREX_SPACEDIM values: a shorter list is an error
    pp.queryarr("is_periodic", is_periodic, 0, AMREX_SPACEDIM);
    pp.queryarr("fine_lo", fine_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("fine_hi", fine_hi, 0, AMREX_SPACEDIM);
    pp.query("fine_shrink", fine_shrink);
    if (pp.contains("fine_notch_lo") || pp.contains("fine_notch_hi")) {
        pp.getarr("fine_notch_lo", fine_notch_lo, 0, AMREX_SPACEDIM);
        pp.getarr("fine_notch_hi", fine_notch_hi, 0, AMREX_SPACEDIM);
    }
    pp.query("null_fine_bcoef", null_fine_bcoef);
    pp.query("null_acoef", null_acoef);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!null_acoef || acoef_value == 0.0, "null_acoef needs acoef = 0");
    pp.query("nan_in_covered", nan_in_covered);
    pp.query("domain_bc", domain_bc);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(domain_bc == "dirichlet" || domain_bc == "neumann",
                                     "domain_bc must be dirichlet or neumann");
    pp.queryarr("prob_hi", prob_hi, 0, AMREX_SPACEDIM);
    pp.query("use_mlhypre", use_mlhypre);
#ifndef AMREX_USE_HYPRE
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!use_mlhypre, "use_mlhypre needs USE_HYPRE=TRUE");
#endif
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!use_mlhypre || composite_solve,
                                     "use_mlhypre needs composite_solve");
    pp.query("mlhypre_num_setups", mlhypre_num_setups);
    pp.query("mlhypre_max_diff", mlhypre_max_diff);
    pp.query("mlhypre_max_fine_res", mlhypre_max_fine_res);
    pp.query("mlhypre_max_avg_diff", mlhypre_max_avg_diff);
    pp.query("mlhypre_max_exact_ratio", mlhypre_max_exact_ratio);
    pp.query("mlhypre_max_res_scaled", mlhypre_max_res_scaled);
    pp.query("mlhypre_max_exact_err", mlhypre_max_exact_err);
}

void
MyTest::initGrids ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    // Fine levels will be added later
    grids[0].define(domain0);
    grids[0].maxSize(max_grid_size);
}

void
MyTest::addFineGrids ()
{
    for (int ilev = 1; ilev <= max_level; ++ilev)
    {
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        BoxList bl = eb_level.boxArray().boxList();
        if (bl.isEmpty()) { bl.push_back(geom[ilev].Domain()); } // all regular
        Box domain = geom[ilev].Domain();
        auto const dx = geom[ilev].CellSizeArray();
        Real const shrink = (ilev-1)*fine_shrink;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            domain.setSmall(idim, int(std::round((fine_lo[idim]+shrink)/dx[idim])));
            domain.setBig(idim, int(std::round((fine_hi[idim]-shrink)/dx[idim]))-1);
        }
        for (Box& b : bl) {
            b &= domain;
        }
        if (!fine_notch_lo.empty()) {
            Box notch;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                notch.setSmall(idim, int(std::round(fine_notch_lo[idim]/dx[idim])));
                notch.setBig(idim, int(std::round(fine_notch_hi[idim]/dx[idim]))-1);
            }
            BoxList notched;
            for (Box const& b : bl) {
                notched.join(amrex::boxDiff(b, notch));
            }
            bl = notched;
        }
        bl.simplify();
        grids[ilev].define(bl);
        grids[ilev].maxSize(max_grid_size);
    }
}

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    dmap.resize(nlevels);
    factory.resize(nlevels);
    phi.resize(nlevels);
    phiexact.resize(nlevels);
    phieb.resize(nlevels);
    rhs.resize(nlevels);
    acoef.resize(nlevels);
    bcoef.resize(nlevels);
    bcoef_eb.resize(nlevels);
    fluxeb_sol.resize(nlevels);
    fluxeb_exact.resize(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        factory[ilev] = std::make_unique<EBFArrayBoxFactory>
            (eb_level, geom[ilev], grids[ilev], dmap[ilev], Vector<int>{2,2,2}, EBSupport::full);

        phi[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        phiexact[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        phieb[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].define(amrex::convert(grids[ilev],IntVect::TheDimensionVector(idim)),
                                     dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        }
        bcoef_eb[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        bcoef_eb[ilev].setVal(1.0);
        fluxeb_sol[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        fluxeb_exact[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);

        phi[ilev].setVal(ghost_value); // the valid and domain boundary values are set below
        rhs[ilev].setVal(0.0);
        acoef[ilev].setVal(acoef_value);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].setVal(1.0);
        }
        fluxeb_sol[ilev].setVal(0.0);

        const auto dx = geom[ilev].CellSizeArray();
        const Box& domainbox = geom[ilev].Domain();
        const auto lprob_type = this->prob_type;

        const FabArray<EBCellFlagFab>& flags = factory[ilev]->getMultiEBCellFlagFab();
        const MultiCutFab& bcent = factory[ilev]->getBndryCent();
        const MultiCutFab& cent = factory[ilev]->getCentroid();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(phiexact[ilev]); mfi.isValid(); ++mfi)
        {
            const amrex::EBData& ebdata = factory[ilev]->getEBData(mfi);
            const Box& bx = mfi.validbox();
            const Box& nbx = amrex::surroundingNodes(bx);
            Array4<Real> const& phi_arr = phi[ilev].array(mfi);
            Array4<Real> const& phi_ex_arr = phiexact[ilev].array(mfi);
            Array4<Real> const& phi_eb_arr = phieb[ilev].array(mfi);
            Array4<Real> const& rhs_arr = rhs[ilev].array(mfi);
            AMREX_D_TERM(Array4<Real> const& bx_arr = bcoef[ilev][0].array(mfi);,
                         Array4<Real> const& by_arr = bcoef[ilev][1].array(mfi);,
                         Array4<Real> const& bz_arr = bcoef[ilev][2].array(mfi););
            Array4<Real> const& feb_ex_arr = fluxeb_exact[ilev].array(mfi);

            auto fabtyp = flags[mfi].getType(bx);
            if (FabType::covered == fabtyp) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    phi_ex_arr(i,j,k) = 0.0;
                    phi_eb_arr(i,j,k) = 0.0;
                    feb_ex_arr(i,j,k) = 0.0;
                });
            } else if (FabType::regular == fabtyp) {
                amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    mytest_set_phi_reg(i,j,k,phi_ex_arr,rhs_arr,
                                       AMREX_D_DECL(bx_arr,by_arr,bz_arr),
                                       dx, lprob_type, bx);
                    if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                        feb_ex_arr(i,j,k) = 0.0;
                    }
                });
            } else {
                Array4<Real> const& beb_arr = bcoef_eb[ilev].array(mfi);
                Array4<EBCellFlag const> const& flag_arr = flags.const_array(mfi);
                Array4<Real const> const& cent_arr = cent.const_array(mfi);
                Array4<Real const> const& bcent_arr = bcent.const_array(mfi);
                amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    mytest_set_phi_eb(i,j,k,phi_ex_arr,phi_eb_arr,rhs_arr,
                                      AMREX_D_DECL(bx_arr,by_arr,bz_arr),
                                      beb_arr,flag_arr,cent_arr,bcent_arr,
                                      dx, lprob_type, bx);
                    mytest_set_fluxeb_exact(i,j,k,feb_ex_arr,flag_arr,ebdata,dx,lprob_type,bx);
                });
            }

            const Box& gbx = mfi.growntilebox(1);
            if (!domainbox.contains(gbx)) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    mytest_set_phi_boundary(i,j,k,phi_arr,dx,domainbox,lprob_type);
                });
            }
        }
        // b = 1 on the finer levels: given to HypreMLABecLap as null (1) or
        // explicitly (2), and to MLMG
        if (null_fine_bcoef && ilev > 0) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                bcoef[ilev][idim].setVal(1.0);
            }
        }
    }
}
