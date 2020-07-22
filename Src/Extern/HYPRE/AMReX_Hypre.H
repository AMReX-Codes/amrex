#ifndef AMREX_HYPRE_H_
#define AMREX_HYPRE_H_

#include <memory>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BndryData.H>

#include "HYPRE.h"
#include "_hypre_utilities.h"

namespace amrex
{

class Hypre
{
public:

    enum class Interface : int { structed, semi_structed, ij };
    
    Hypre (const BoxArray& grids,
           const DistributionMapping& dmap,
           const Geometry& geom,
           MPI_Comm comm_);

    virtual ~Hypre ();

    void setScalars (Real sa, Real sb);
    void setACoeffs (const MultiFab& alpha);
    void setBCoeffs (const Array<const MultiFab*,BL_SPACEDIM>& beta);
    void setVerbose (int _verbose);
    virtual void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                        int max_iter, const BndryData& bndry, int max_bndry_order) = 0;

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
    
    void setHypreOldDefault (bool l) noexcept {old_default = l;}
    void setHypreRelaxType (int n) noexcept {relax_type = n;}
    void setHypreRelaxOrder (int n) noexcept {relax_order = n;}
    void setHypreNumSweeps (int n) noexcept {num_sweeps = n;}
    void setHypreStrongThreshold (Real t) noexcept {strong_threshold = t;}

protected:

    static constexpr HYPRE_Int regular_stencil_size = 2*AMREX_SPACEDIM + 1;
    static constexpr HYPRE_Int eb_stencil_size = AMREX_D_TERM(3, *3, *3);
    
    MPI_Comm comm = MPI_COMM_NULL;
    Geometry geom;

    int verbose = 0;
    bool old_default = true; // Falgout coarsening with modified classical interpolation
    int relax_type = 6;  // G-S/Jacobi hybrid relaxation
    int relax_order = 1; // uses C/F relaxation
    int num_sweeps = 2;  // Sweeeps on each level
    Real strong_threshold = 0.25; // Hypre default is 0.25

    MultiFab acoefs;
    Array<MultiFab,AMREX_SPACEDIM> bcoefs;
    Real scalar_a, scalar_b;

    MultiFab diaginv;
    
    FabFactory<FArrayBox> const* m_factory = nullptr;
    BndryData const* m_bndry = nullptr;
    int m_maxorder = -1;
};

std::unique_ptr<Hypre> makeHypre (const BoxArray& grids, const DistributionMapping& damp,
                                  const Geometry& geom, MPI_Comm comm_, Hypre::Interface interface);

}

#endif
