
 #include <AMReX_MLMG.H>

 using namespace amrex;

 extern "C" {

     void amrex_fi_new_multigrid (MLMG*& mlmg, MLLinOp* lp)
     {
         mlmg = new MLMG(*lp);
     }

     void amrex_fi_delete_multigrid (MLMG* mlmg)
     {
         delete mlmg;
     }

     Real amrex_fi_multigrid_solve (MLMG* mlmg, MultiFab* a_sol[], MultiFab* a_rhs[],
                                    Real a_tol_rel, Real a_tol_abs)
     {
         const int n = mlmg->numAMRLevels();
         return mlmg->solve(Vector<MultiFab*>{a_sol, a_sol+n},
                            Vector<const MultiFab*>{a_rhs, a_rhs+n},
                            a_tol_rel, a_tol_abs);
     }

     void amrex_fi_multigrid_get_grad_solution (MLMG* mlmg, MultiFab* a_grad_sol[])
     {
         MultiFab** p = a_grad_sol;
         const int n = mlmg->numAMRLevels();
         Vector<std::array<MultiFab*,AMREX_SPACEDIM> > grad(n);
         for (auto& v : grad) {
             for (auto& a : v) {
                 a = *p++;
             }
         }
         mlmg->getGradSolution(grad);
     }

     void amrex_fi_multigrid_get_fluxes (MLMG* mlmg, MultiFab* a_fluxes[])
     {
         MultiFab** p = a_fluxes;
         const int n = mlmg->numAMRLevels();
         Vector<std::array<MultiFab*,AMREX_SPACEDIM> > fluxes(n);
         for (auto& v : fluxes) {
             for (auto& a : v) {
                 a = *p++;
             }
         }
         mlmg->getFluxes(fluxes);
     }

     void amrex_fi_multigrid_comp_residual (MLMG* mlmg, MultiFab* a_res[], MultiFab* a_sol[],
                                            MultiFab* a_rhs[])
     {
         const int n = mlmg->numAMRLevels();
         mlmg->compResidual(Vector<MultiFab*>{a_res, a_res+n},
                            Vector<MultiFab*>{a_sol, a_sol+n},
                            Vector<const MultiFab*>{a_rhs, a_rhs+n});
     }

     void amrex_fi_multigrid_set_verbose (MLMG* mlmg, int v)
     {
         mlmg->setVerbose(v);
     }

     void amrex_fi_multigrid_set_max_iter (MLMG* mlmg, int n)
     {
         mlmg->setMaxIter(n);
     }

     void amrex_fi_multigrid_set_max_fmg_iter (MLMG* mlmg, int n)
     {
         mlmg->setMaxFmgIter(n);
     }

     void amrex_fi_multigrid_set_bottom_solver (MLMG* mlmg, int s)
     {
         if (s == 0) {
             mlmg->setBottomSolver(MLMG::BottomSolver::smoother);
         } else if (s == 1) {
             mlmg->setBottomSolver(MLMG::BottomSolver::bicgstab);
         } else if (s == 2) {
             mlmg->setBottomSolver(MLMG::BottomSolver::cg);
         } else if (s == 3) {
             mlmg->setBottomSolver(MLMG::BottomSolver::hypre);
         } else if (s == 4) {
             mlmg->setBottomSolver(MLMG::BottomSolver::petsc);
         } else {
             amrex::Abort("amrex_fi_multigrid_set_bottom_solver: unknown bottom solver");
         }
     }

     void amrex_fi_multigrid_set_cg_verbose (MLMG* mlmg, int n)
     {
         mlmg->setCGVerbose(n);
     }
     
     void amrex_fi_multigrid_set_always_use_bnorm (MLMG* mlmg, int f)
     {
         mlmg->setAlwaysUseBNorm(f);
     }
     
     void amrex_fi_multigrid_set_final_fill_bc (MLMG* mlmg, int f)
     {
         mlmg->setFinalFillBC(f);
     }
}
