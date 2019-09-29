#ifndef AMREX_MAC_PROJECTOR_H_
#define AMREX_MAC_PROJECTOR_H_

#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

#ifdef AMREX_USE_EB
#include <AMReX_MLEBABecLap.H>
#endif

namespace amrex {

class MacProjector
{
public:

    MacProjector (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  const Vector<Array<MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  const Vector<Geometry>& a_geom,
                  const LPInfo& a_lpinfo,  
                  const Vector<MultiFab const*>& a_divu = {});

    MacProjector (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  const Vector<Array<MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  const Vector<Geometry>& a_geom,
                  const Vector<MultiFab const*>& a_divu = {})
     : MacProjector(a_umac, a_beta, a_geom, LPInfo(), a_divu) {};

    void setDomainBC (const Array<LinOpBCType,AMREX_SPACEDIM>& lobc,
                      const Array<LinOpBCType,AMREX_SPACEDIM>& hibc);

#ifndef AMREX_USE_EB
    void project (const Vector<MultiFab*>& phi_in, Real reltol, Real atol)
     { project(phi_in, reltol, atol, MLMG::Location::FaceCenter); }
    void project (Real reltol, Real atol)
     { project(reltol, atol, MLMG::Location::FaceCenter); }
#endif

    void project (const Vector<MultiFab*>& phi_in, Real reltol, Real atol, MLMG::Location);
    void project (Real reltol, Real atol, MLMG::Location);

    void setVerbose            (int  v) noexcept { m_verbose = v; }
    void setCGVerbose          (int  v) noexcept { m_cg_verbose = v; }
    void setBottomVerbose      (int  v) noexcept { m_cg_verbose = v; }
    void setBottomTolerance    (Real t) noexcept { m_bottom_reltol = t; }
    void setBottomToleranceAbs (Real t) noexcept { m_bottom_abstol = t; }
    void setMaxIter (int v) noexcept { m_maxiter = v; }
    void setCGMaxIter (int v) noexcept { m_cg_maxiter = v; }

    void setBottomSolver (MLMG::BottomSolver bottom_solver_type_in) noexcept
     { bottom_solver_type = bottom_solver_type_in;}

    MLLinOp& getLinOp () noexcept { return *m_linop; }

private:

    std::unique_ptr<MLABecLaplacian> m_abeclap;
#ifdef AMREX_USE_EB
    std::unique_ptr<MLEBABecLap> m_eb_abeclap;
    Vector<EBFArrayBoxFactory const*> m_eb_factory;
#endif
    MLLinOp* m_linop = nullptr;

    std::unique_ptr<MLMG> m_mlmg;

    Vector<Array<MultiFab*,AMREX_SPACEDIM> > m_umac;
    Vector<MultiFab> m_rhs;
    Vector<MultiFab> m_phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > m_fluxes;

    Vector<Geometry> m_geom;

    int m_verbose = 0;
    int m_cg_verbose = 0;
    int m_maxiter = 200;
    int m_cg_maxiter = 200;
    int m_max_coarsening_level = 100;

    Real m_bottom_reltol = 1.e-4;
    Real m_bottom_abstol = -1.0;

    MLMG::BottomSolver bottom_solver_type = MLMG::BottomSolver::bicgstab;
};

}

#endif
