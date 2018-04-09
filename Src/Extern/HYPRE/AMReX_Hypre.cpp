
#include <AMReX_Hypre.H>

namespace amrex
{

Hypre::Hypre (const BoxArray& grids,
              const DistributionMapping& dmap,
              const Geometry& geom,
              MPI_Comm comm_)
{
    ParmParse pp("hypre");
    
    // solver_flag = 0 for SMG
    // solver_flag = 1 for PFMG
    // solver_falg = 2 for BoomerAMG
    int solver_flag = 1;
    {
        std::string solver_flag_s {"null"};
        pp.query("solver_flag", solver_flag_s);
        std::transform(solver_flag_s.begin(), solver_flag_s.end(), solver_flag_s.begin(), ::tolower);
        if (solver_flag_s == "smg") {
            solver_flag = 0;
        } else if (solver_flag_s == "pfmg") {
            solver_flag = 1;
        } else if (solver_flag_s == "boomeramg") {
            solver_flag = 2;
        } else if (solver_flag_s == "boomeramgij") {
            solver_flag = 3;
        } else if (solver_flag_s == "none") {
            pp.query("solver_flag", solver_flag);
        } else {
            amrex::Abort("Hypre: unknown solver flag");
        }
    }

    if (solver_flag == 2) {
        semi_struct_solver.reset(new HypreABecLap2(grids, dmap, geom, comm_));
    } else if (solver_flag == 3) {
        IJ_solver.reset(new HypreABecLap3(grids, dmap, geom, comm_));
    }else {
        struct_solver.reset(new HypreABecLap(grids, dmap, geom, comm_));
    }
}

Hypre::~Hypre ()
{

}

void
Hypre::setScalars (Real sa, Real sb)
{
    if (struct_solver) {
        struct_solver->setScalars(sa,sb);
    } else if (semi_struct_solver) {
        semi_struct_solver->setScalars(sa,sb);
    }else if (IJ_solver) {
      IJ_solver->setScalars(sa,sb);
    }else {
        amrex::Abort("Hypre::setScalars: How did this happen?");
    }
}

void
Hypre::setACoeffs (const MultiFab& alpha)
{
    if (struct_solver) {
        struct_solver->setACoeffs(alpha);
    } else if (semi_struct_solver) {
        semi_struct_solver->setACoeffs(alpha);
    }else if (IJ_solver) {
        IJ_solver->setACoeffs(alpha);
    } else {
        amrex::Abort("Hypre::setACoeffs: How did this happen?");
    }
}

void
Hypre::setBCoeffs (const std::array<const MultiFab*,BL_SPACEDIM>& beta)
{
    if (struct_solver) {
        struct_solver->setBCoeffs(beta);
    } else if (semi_struct_solver) {
        semi_struct_solver->setBCoeffs(beta);
    }else if (IJ_solver) {
        IJ_solver->setBCoeffs(beta);
    } else {
        amrex::Abort("Hypre::setBCoeffs: How did this happen?");        
    }
}

void
Hypre::setVerbose (int _verbose)
{
    if (struct_solver) {
        struct_solver->setVerbose(_verbose);
    } else if (semi_struct_solver) {
        semi_struct_solver->setVerbose(_verbose);
    }else if (IJ_solver) {
        IJ_solver->setVerbose(_verbose);
    } else {
        amrex::Abort("Hypre::setVerbose: How did this happen?");
    }
}

void
Hypre::solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
              int max_iter, LinOpBCType bc_type, Real bc_value)
{
    if (struct_solver) {
        struct_solver->solve(soln, rhs, rel_tol, abs_tol, max_iter, bc_type, bc_value);
    } else if (semi_struct_solver) {
        semi_struct_solver->solve(soln, rhs, rel_tol, abs_tol, max_iter, bc_type, bc_value);
    } else if (IJ_solver) {
      IJ_solver->solve(soln, rhs, rel_tol, abs_tol, max_iter, bc_type, bc_value);
    } else {
        amrex::Abort("Hypre::solve: How did this happen?");
    }
}

}
