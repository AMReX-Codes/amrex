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

    //
    // Constructors
    //
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

    //
    // Methods to set BCs and coarse/fine values
    //
    // These methods are wrappers of the linop methods of the same name
    // However, use of these is preferred to make sure operations
    // are performed in the correct order
    //
    void setDomainBC (const Array<LinOpBCType,AMREX_SPACEDIM>& lobc,
                      const Array<LinOpBCType,AMREX_SPACEDIM>& hibc);

    void setLevelBC  (int amrlev, const MultiFab* levelbcdata);

    void setCoarseFineBC (const MultiFab* crse, int crse_ratio)
        { m_linop->setCoarseFineBC(crse, crse_ratio);}

    //
    // Methods to perform projection
    //
#ifndef AMREX_USE_EB
    void project (const Vector<MultiFab*>& phi_in, Real reltol, Real atol)
     { project(phi_in, reltol, atol, MLMG::Location::FaceCenter); }
    void project (Real reltol, Real atol)
     { project(reltol, atol, MLMG::Location::FaceCenter); }
#endif

    void project (const Vector<MultiFab*>& phi_in, Real reltol, Real atol, MLMG::Location);
    void project (Real reltol, Real atol, MLMG::Location);

    //
    // Setters and getters
    //
    void setVerbose            (int  v) noexcept { m_verbose = v; }

    // Methods to get underlying objects
    // Use these to modify properties of MLMG and linear operator
    MLLinOp& getLinOp () noexcept { return *m_linop; }
    MLMG&    getMLMG  () noexcept { return *m_mlmg;  }

private:

    void setOptions ();

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

    bool m_needs_domain_bcs = true;

};

}

#endif
