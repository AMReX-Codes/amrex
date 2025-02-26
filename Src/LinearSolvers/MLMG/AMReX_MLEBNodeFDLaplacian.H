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
// del dot (sigma grad phi) - alpha/r^2 phi = rhs, for RZ where alpha is a
// scalar constant that is zero by default.
//
// New feature: for non-RZ, sigma can also be a single-component
// cell-centered multifab.

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

    ~MLEBNodeFDLaplacian () override = default;

    MLEBNodeFDLaplacian (const MLEBNodeFDLaplacian&) = delete;
    MLEBNodeFDLaplacian (MLEBNodeFDLaplacian&&) = delete;
    MLEBNodeFDLaplacian& operator= (const MLEBNodeFDLaplacian&) = delete;
    MLEBNodeFDLaplacian& operator= (MLEBNodeFDLaplacian&&) = delete;

    void setSigma (Array<Real,AMREX_SPACEDIM> const& a_sigma) noexcept;

    void setSigma (int amrlev, MultiFab const& a_sigma);

    void setRZ (bool flag);

    void setAlpha (Real a_alpha);

#ifdef AMREX_USE_EB

    // Phi on EB
    void setEBDirichlet (Real a_phi_eb);
    //
    template <typename F>
    std::enable_if_t<IsCallableR<Real,F,AMREX_D_DECL(Real,Real,Real)>::value>
    setEBDirichlet (F const& f);

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

    [[nodiscard]] std::unique_ptr<FabFactory<FArrayBox> > makeFactory (int amrlev, int mglev) const final;

    [[nodiscard]] bool scaleRHS (int amrlev, MultiFab* rhs) const final;

#endif

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info);

    [[nodiscard]] std::string name () const override { return std::string("MLEBNodeFDLaplacian"); }

    void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;

    void prepareForSolve () final;
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final;
    void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

    [[nodiscard]] bool isSingular (int) const final { return false; }
    [[nodiscard]] bool isBottomSingular () const final { return false; }

    void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                   MultiFab& sol, Location /*loc*/) const override;

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    void fillIJMatrix (MFIter const& mfi,
                       Array4<HypreNodeLap::AtomicInt const> const& gid,
                       Array4<int const> const& lid,
                       HypreNodeLap::Int* ncols,
                       HypreNodeLap::Int* cols,
                       Real* mat) const override;

    void fillRHS (MFIter const& mfi,
                  Array4<int const> const& lid,
                  Real* rhs,
                  Array4<Real const> const& bfab) const override;
#endif

    void postSolve (Vector<MultiFab*> const& sol) const override;

private:
    GpuArray<Real,AMREX_SPACEDIM> m_sigma{{AMREX_D_DECL(1_rt,1_rt,1_rt)}};
    Vector<Vector<std::unique_ptr<MultiFab>>> m_sigma_mf;
    bool m_has_sigma_mf = false;
    Real m_s_phi_eb = std::numeric_limits<Real>::lowest();
    Vector<MultiFab> m_phi_eb;
    int m_rz = false;
    Real m_rz_alpha = 0._rt;
};

#ifdef AMREX_USE_EB

template <typename F>
std::enable_if_t<IsCallableR<Real,F,AMREX_D_DECL(Real,Real,Real)>::value>
MLEBNodeFDLaplacian::setEBDirichlet (F const& f)
{
    m_phi_eb.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        auto const* factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
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
