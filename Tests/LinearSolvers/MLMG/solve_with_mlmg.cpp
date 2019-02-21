#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <prob_par.H>

using namespace amrex;

namespace {
static bool composite_solve = true;
static bool fine_leve_solve_only = false;
static int max_iter = 100;
static int max_fmg_iter = 20;
static int max_coarsening_level = 30;
static int verbose  = 2;
static int cg_verbose = 0;
static int linop_maxorder = 2;
static bool agglomeration = false;
static bool consolidation = false;
static int  use_hypre = 0;
}

void solve_with_mlmg(const Vector<Geometry>& geom, int ref_ratio,
                      Vector<MultiFab>& soln,
                      const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                      Vector<MultiFab>& rhs, const Vector<MultiFab>& exact) {
  BL_PROFILE("solve_with_mlmg");

  Real tol_rel = 1.e-10;
  Real tol_abs = 0.0;

  {
    ParmParse pp;
    pp.query("composite_solve", composite_solve);
    pp.query("fine_leve_solve_only", fine_leve_solve_only);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("max_coarsening_level", max_coarsening_level);
    pp.query("verbose", verbose);
    pp.query("cg_verbose", cg_verbose);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("use_hypre", use_hypre);
    pp.query("tol_rel", tol_rel);
    pp.query("tol_abs", tol_abs);
  }

  LPInfo info;
  info.setAgglomeration(agglomeration);
  info.setConsolidation(consolidation);
  info.setMaxCoarseningLevel(max_coarsening_level);

  const int nlevels = geom.size();

  if (composite_solve) {
    Vector<BoxArray> grids;
    Vector<DistributionMapping> dmap;
    Vector<MultiFab*> psoln;
    Vector<MultiFab const*> prhs;
    for (int ilev = 0; ilev < nlevels; ++ilev) {
      grids.push_back(soln[ilev].boxArray());
      dmap.push_back(soln[ilev].DistributionMap());
      psoln.push_back(&(soln[ilev]));
      prhs.push_back(&(rhs[ilev]));
    }

    MLABecLaplacian mlabec(geom, grids, dmap, info);
    mlabec.setMaxOrder(linop_maxorder);
    // BC
    mlabec.setDomainBC({prob::bc_type, prob::bc_type, prob::bc_type},
                       {prob::bc_type, prob::bc_type, prob::bc_type});
    for (int ilev = 0; ilev < nlevels; ++ilev) {
      mlabec.setLevelBC(ilev, psoln[ilev]);
    }
    mlabec.setScalars(prob::a, prob::b);
    for (int ilev = 0; ilev < nlevels; ++ilev) {
      mlabec.setACoeffs(ilev, alpha[ilev]);
      std::array<MultiFab, AMREX_SPACEDIM> bcoefs;
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        const BoxArray& ba = amrex::convert(beta[ilev].boxArray(),
                                            IntVect::TheDimensionVector(idim));
        bcoefs[idim].define(ba, beta[ilev].DistributionMap(), 1, 0);
      }
      amrex::average_cellcenter_to_face(amrex::GetArrOfPtrs(bcoefs),
                                        beta[ilev], geom[ilev]);
      mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcoefs));
    }

    MLMG mlmg(mlabec);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    if (use_hypre) mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(cg_verbose);

    mlmg.solve(psoln, prhs, tol_rel, tol_abs);
  } else {
    const int levbegin = (fine_leve_solve_only) ? nlevels-1 : 0;
    for (int ilev = 0; ilev < levbegin; ++ilev) {
      MultiFab::Copy(soln[ilev], exact[ilev], 0, 0, 1, 0);
    }

    for (int ilev = levbegin; ilev < nlevels; ++ilev) {
      MLABecLaplacian mlabec({geom[ilev]},
                             {soln[ilev].boxArray()},
                             {soln[ilev].DistributionMap()},
                             info);

      mlabec.setMaxOrder(linop_maxorder);

      mlabec.setDomainBC({prob::bc_type, prob::bc_type, prob::bc_type},
                         {prob::bc_type, prob::bc_type, prob::bc_type});
      const int solver_level = 0;  // 0 even though ilev may be > 0
      if (ilev > 0) {
        mlabec.setCoarseFineBC(&soln[ilev-1], ref_ratio);
      }
      mlabec.setLevelBC(solver_level, &soln[ilev]);

      mlabec.setScalars(prob::a, prob::b);
      mlabec.setACoeffs(solver_level, alpha[ilev]);

      std::array<MultiFab, AMREX_SPACEDIM> bcoefs;
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        const BoxArray& ba = amrex::convert(beta[ilev].boxArray(),
                                            IntVect::TheDimensionVector(idim));
        bcoefs[idim].define(ba, beta[ilev].DistributionMap(), 1, 0);
      }
      amrex::average_cellcenter_to_face(amrex::GetArrOfPtrs(bcoefs),
                                        beta[ilev], geom[ilev]);
      mlabec.setBCoeffs(solver_level, amrex::GetArrOfConstPtrs(bcoefs));

      MLMG mlmg(mlabec);
      mlmg.setMaxIter(max_iter);
      mlmg.setMaxFmgIter(max_fmg_iter);
      mlmg.setVerbose(verbose);
      mlmg.setBottomVerbose(cg_verbose);

      mlmg.solve({&soln[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
    }
  }
}

