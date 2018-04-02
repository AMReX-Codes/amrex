#ifndef AMREX_HYPRE_H_
#define AMREX_HYPRE_H_

#include <AMReX_HypreABecLap.H>
#include <AMReX_HypreABecLap2.H>
#include <AMReX_HypreABecLap3.H>

#include <memory>

namespace amrex
{

class Hypre
{
public:
    
    Hypre (const BoxArray& grids,
           const DistributionMapping& dmap,
           const Geometry& geom,
           MPI_Comm comm_);

    ~Hypre ();

    void setScalars(Real sa, Real sb);
    void setACoeffs(const MultiFab& alpha);
    void setBCoeffs(const std::array<const MultiFab*,BL_SPACEDIM>& beta);
    void setVerbose(int _verbose);
    void solve(MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
               int max_iter, LinOpBCType bc_type, Real bc_value);

private:

    std::unique_ptr<HypreABecLap> struct_solver;
    std::unique_ptr<HypreABecLap2> semi_struct_solver;
    std::unique_ptr<HypreABecLap3> IJ_solver;
};

}

#endif
