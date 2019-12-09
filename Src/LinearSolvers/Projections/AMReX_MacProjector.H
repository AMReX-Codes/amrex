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

    void setLevelBC  (int amrlev, const MultiFab* levelbcdata);

#ifndef AMREX_USE_EB
    void project (const Vector<MultiFab*>& phi_in, Real reltol, Real atol)
     { project(phi_in, reltol, atol, MLMG::Location::FaceCenter); }
    void project (Real reltol, Real atol)
     { project(reltol, atol, MLMG::Location::FaceCenter); }
#endif

    void project (const Vector<MultiFab*>& phi_in, Real reltol, Real atol, MLMG::Location);
    void project (Real reltol, Real atol, MLMG::Location);

    // Methods to set verbosity
    void setVerbose            (int  v) noexcept { m_verbose = v; }
    void setBottomVerbose      (int  v) noexcept { m_bottom_verbose = v; }

    // Methods to control MLMG behavior
    void setBottomSolver  (std::string b) noexcept { m_bottom_solver = b;}
    void setMaxIter               (int i) noexcept { m_maxiter = i; }
    void setBottomMaxIter         (int i) noexcept { m_bottom_maxiter = i; }
    void setBottomTolerance      (Real t) noexcept { m_bottom_rtol = t; }
    void setBottomToleranceAbs   (Real t) noexcept { m_bottom_atol = t; }


    MLLinOp& getLinOp () noexcept { return *m_linop; }

private:

    void readParameters ();
    void setup ();

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

    // Verbosity
    int m_verbose = 0;
    int m_bottom_verbose = 0;

    // Control MLMG behavior
    std::string m_bottom_solver = "bicg";
    int  m_maxiter        = 200;
    int  m_bottom_maxiter = 200;
    Real m_bottom_rtol    = 1.0e-4;
    Real m_bottom_atol    = -1.0;

    int m_max_coarsening_level = 100;
    bool m_needs_domain_bcs = true;

};

}

#endif
