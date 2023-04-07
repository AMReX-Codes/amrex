#ifndef AMREX_MLEBNODEFDLAPLACIAN_H_
#define AMREX_MLEBNODEFDLAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif
#include <AMReX_MLNodeLinOp.H>

#include <limits>

namespace amrex {

// Although the class has EB in the name, it works for non-EB build too.
//
// del dot (sigma grad phi) = rhs, for non-RZ
// where phi and rhs are nodal multifab, and sigma is a tensor constant
// with only diagonal components.  The EB is assumed to be Dirichlet.
//
// del dot (simga grad phi) - alpha/r^2 phi = rhs, for RZ where alpha is a
// scalar constant that is zero by default.  For now the `alpha` term has
// not been implemented yet

class MLEBNodeFDLaplacian
    : public MLNodeLinOp
{
public:

    MLEBNodeFDLaplacian () = default;

#ifdef AMREX_USE_EB
    MLEBNodeFDLaplacian (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<EBFArrayBoxFactory const*>& a_factory);
#endif

    MLEBNodeFDLaplacian (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info);

    virtual ~MLEBNodeFDLaplacian () = default;

    MLEBNodeFDLaplacian (const MLEBNodeFDLaplacian&) = delete;
    MLEBNodeFDLaplacian (MLEBNodeFDLaplacian&&) = delete;
    MLEBNodeFDLaplacian& operator= (const MLEBNodeFDLaplacian&) = delete;
    MLEBNodeFDLaplacian& operator= (MLEBNodeFDLaplacian&&) = delete;

    void setSigma (Array<Real,AMREX_SPACEDIM> const& a_sigma) noexcept;

    void setRZ (bool flag);

    void setAlpha (Real a_alpha);

#ifdef AMREX_USE_EB

    // Phi on EB
    void setEBDirichlet (Real a_phi_eb);
    //
    template <typename F>
    std::enable_if_t<IsCallableR<Real,F,AMREX_D_DECL(Real,Real,Real)>::value>
    setEBDirichlet (F&& f);

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

    virtual std::unique_ptr<FabFactory<FArrayBox> > makeFactory (int amrlev, int mglev) const final override;

    virtual void scaleRHS (int amrlev, MultiFab& rhs) const final;

#endif

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info);

    virtual std::string name () const override { return std::string("MLEBNodeFDLaplacian"); }

    virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final override;
    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final override;
    virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final override;

    virtual void reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                         MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final override;

    virtual void prepareForSolve () final override;
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final override;
    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const final override;

    virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final override;

    virtual bool isSingular (int) const final override { return false; }
    virtual bool isBottomSingular () const final override { return false; }

    virtual void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                           MultiFab& sol, Location /*loc*/) const override;

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    virtual void fillIJMatrix (MFIter const& mfi,
                               Array4<HypreNodeLap::AtomicInt const> const& gid,
                               Array4<int const> const& lid,
                               HypreNodeLap::Int* const ncols,
                               HypreNodeLap::Int* const cols,
                               Real* const mat) const override;

    virtual void fillRHS (MFIter const& mfi,
                          Array4<int const> const& lid,
                          Real* const rhs,
                          Array4<Real const> const& bfab) const override;
#endif

    virtual void postSolve (Vector<MultiFab>& sol) const override;

private:
    GpuArray<Real,AMREX_SPACEDIM> m_sigma{{AMREX_D_DECL(1_rt,1_rt,1_rt)}};
    Real m_s_phi_eb = std::numeric_limits<Real>::lowest();
    Vector<MultiFab> m_phi_eb;
    int m_rz = false;
    Real m_rz_alpha = 0._rt;
};

#ifdef AMREX_USE_EB

template <typename F>
std::enable_if_t<IsCallableR<Real,F,AMREX_D_DECL(Real,Real,Real)>::value>
MLEBNodeFDLaplacian::setEBDirichlet (F&& f)
{
    m_phi_eb.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        if (factory) {
            Geometry const& geom = m_geom[amrlev][0];
            auto const problo = geom.ProbLoArray();
            auto const cellsize = geom.CellSizeArray();
            m_phi_eb[amrlev].define(amrex::convert(m_grids[amrlev][0],IntVect(1)),
                                    m_dmap[amrlev][0], 1, 1);
            m_phi_eb[amrlev].setVal(0.0);
            auto const& flags = factory->getMultiEBCellFlagFab();
            auto const& levset = factory->getLevelSet();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(m_phi_eb[amrlev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& ndbx = mfi.growntilebox();
                const auto& flag = flags[mfi];
                if (flag.getType() != FabType::regular) {
                    Array4<Real const> const lstarr = levset.const_array(mfi);
                    Array4<Real> const& phi = m_phi_eb[amrlev].array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(ndbx, i, j, k,
                    {
                        if (lstarr(i,j,k) >= Real(0.0)) {
                            phi(i,j,k) = f(AMREX_D_DECL(problo[0]+i*cellsize[0],
                                                        problo[1]+j*cellsize[1],
                                                        problo[2]+k*cellsize[2]));
                        }
                    });
                }
            }
        }
    }
}

#endif

}

#endif
