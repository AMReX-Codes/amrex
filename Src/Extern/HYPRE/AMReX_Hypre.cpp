
#include <AMReX_Hypre.H>
#include <HYPRE_config.h>
#include <type_traits>

namespace amrex {

Hypre::Hypre (const BoxArray& grids,
              const DistributionMapping& dmap,
              const Geometry& geom,
              MPI_Comm comm_, Interface interface)
{
    static_assert(std::is_same<Real, HYPRE_Real>::value, "amrex::Real != HYPRE_Real");
#ifdef HYPRE_BIGINT
    static_assert(std::is_same<long long int, HYPRE_Int>::value, "long long int != HYPRE_Int");
#else
    static_assert(std::is_same<int, HYPRE_Int>::value, "int != HYPRE_Int");
#endif
    
    if (interface == Interface::structed) {
        struct_solver.reset(new HypreABecLap(grids, dmap, geom, comm_));
    } else if (interface == Interface::semi_structed) {
        semi_struct_solver.reset(new HypreABecLap2(grids, dmap, geom, comm_));
    } else if (interface == Interface::ij) {
        IJ_solver.reset(new HypreABecLap3(grids, dmap, geom, comm_));
    } else {
        amrex::Abort("Hypre: unknown interface");
    }
}

Hypre::~Hypre () {
}

void
Hypre::setScalars (Real sa, Real sb)
{
    if (struct_solver) {
        struct_solver->setScalars(sa, sb);
    } else if (semi_struct_solver) {
        semi_struct_solver->setScalars(sa, sb);
    } else if (IJ_solver) {
        IJ_solver->setScalars(sa, sb);
    } else {
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
    } else if (IJ_solver) {
        IJ_solver->setACoeffs(alpha);
    } else {
        amrex::Abort("Hypre::setACoeffs: How did this happen?");
    }
}

void
Hypre::setBCoeffs (const Array<const MultiFab*, BL_SPACEDIM>& beta)
{
    if (struct_solver) {
        struct_solver->setBCoeffs(beta);
    } else if (semi_struct_solver) {
        semi_struct_solver->setBCoeffs(beta);
    } else if (IJ_solver) {
        IJ_solver->setBCoeffs(beta);
    } else {
        amrex::Abort("Hypre::setBCoeffs: How did this happen?");        
    }
}

void
Hypre::setVerbose (int _verbose) {
    if (struct_solver) {
        struct_solver->setVerbose(_verbose);
    } else if (semi_struct_solver) {
        semi_struct_solver->setVerbose(_verbose);
    } else if (IJ_solver) {
        IJ_solver->setVerbose(_verbose);
    } else {
        amrex::Abort("Hypre::setVerbose: How did this happen?");
    }
}

void
Hypre::solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
              int max_iter, const BndryData& bndry, int max_bndry_order)
{
    if (struct_solver) {
        struct_solver->solve(soln, rhs, rel_tol, abs_tol, max_iter, bndry, max_bndry_order);
    } else if (semi_struct_solver) {
        semi_struct_solver->solve(soln, rhs, rel_tol, abs_tol, max_iter, bndry, max_bndry_order);
    } else if (IJ_solver) {
        IJ_solver->solve(soln, rhs, rel_tol, abs_tol, max_iter, bndry, max_bndry_order);
    } else {
        amrex::Abort("Hypre::solve: How did this happen?");
    }
}

}  // namespace amrex
