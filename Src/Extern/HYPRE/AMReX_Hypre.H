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

    enum class Interface : int { structed, semi_structed, ij };
    
    Hypre (const BoxArray& grids,
           const DistributionMapping& dmap,
           const Geometry& geom,
           MPI_Comm comm_,
           Interface interface);

    ~Hypre ();

    void setScalars (Real sa, Real sb);
    void setACoeffs (const MultiFab& alpha);
    void setBCoeffs (const Array<const MultiFab*,BL_SPACEDIM>& beta);
    void setVerbose (int _verbose);
    void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                int max_iter, const BndryData& bndry, int max_bndry_order);

    static HYPRE_Int ispow2 (HYPRE_Int i)
    {
        return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
    }

    static Array<HYPRE_Int,AMREX_SPACEDIM> loV (const Box& b) {
        const auto& v = b.loVect();
        return {AMREX_D_DECL(static_cast<HYPRE_Int>(v[0]),
                             static_cast<HYPRE_Int>(v[1]),
                             static_cast<HYPRE_Int>(v[2]))};
    }

    static Array<HYPRE_Int,AMREX_SPACEDIM> hiV (const Box& b) {
        const auto& v = b.hiVect();
        return {AMREX_D_DECL(static_cast<HYPRE_Int>(v[0]),
                             static_cast<HYPRE_Int>(v[1]),
                             static_cast<HYPRE_Int>(v[2]))};
    }
    
private:

    std::unique_ptr<HypreABecLap> struct_solver;
    std::unique_ptr<HypreABecLap2> semi_struct_solver;
    std::unique_ptr<HypreABecLap3> IJ_solver;
};

}

#endif
