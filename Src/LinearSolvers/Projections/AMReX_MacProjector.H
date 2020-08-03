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
                  MLMG::Location a_umac_loc,
                  const Vector<Array<MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  MLMG::Location a_beta_loc,
                  MLMG::Location  a_phi_loc,
                  const Vector<Geometry>& a_geom,
                  const LPInfo& a_lpinfo,
                  const Vector<MultiFab const*>& a_divu,
                  MLMG::Location a_divu_loc,
                  const Vector<iMultiFab const*>& a_overset_mask = {});

    MacProjector (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  MLMG::Location a_umac_loc,
                  const Vector<Array<MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  MLMG::Location a_beta_loc,
                  MLMG::Location  a_phi_loc,
                  const Vector<Geometry>& a_geom,
                  const LPInfo& a_lpinfo,
                  const Vector<MultiFab const*>& a_divu = {})
     : MacProjector(a_umac, a_umac_loc, a_beta, a_beta_loc, a_phi_loc,
                    a_geom, a_lpinfo, a_divu, MLMG::Location::CellCenter) {}

#ifndef AMREX_USE_EB
    MacProjector (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  const Vector<Array<MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  const Vector<Geometry>& a_geom,
                  const LPInfo& a_lpinfo,
                  const Vector<MultiFab const*>& a_divu = {}) 
     : MacProjector(a_umac, MLMG::Location::FaceCenter,
                    a_beta, MLMG::Location::FaceCenter, MLMG::Location::CellCenter,
                    a_geom, a_lpinfo, a_divu, MLMG::Location::CellCenter) {}

    MacProjector (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_umac,
                  const Vector<Array<MultiFab const*,AMREX_SPACEDIM> >& a_beta,
                  const Vector<Geometry>& a_geom,
                  const Vector<MultiFab const*>& a_divu = {})
     : MacProjector(a_umac, MLMG::Location::FaceCenter,
                    a_beta, MLMG::Location::FaceCenter, MLMG::Location::CellCenter,
                    a_geom, LPInfo(), a_divu, MLMG::Location::CellCenter) {}
#endif

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
    void project (const Vector<MultiFab*>& phi_in, Real reltol, Real atol);
    void project (Real reltol, Real atol);

    //
    // Setters and getters
    //
    void setVerbose            (int  v) noexcept 
       { m_verbose = v;
         m_mlmg->setVerbose(m_verbose); }

    // Methods to get underlying objects
    // Use these to modify properties of MLMG and linear operator
    MLLinOp& getLinOp () noexcept { return *m_linop; }
    MLMG&    getMLMG  () noexcept { return *m_mlmg;  }

private:

    void setOptions ();

    void averageDownVelocity ();

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
    Vector<MultiFab> m_divu;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > m_fluxes;

    Vector<Geometry> m_geom;

    int m_verbose = 0;

    bool m_needs_domain_bcs = true;

    // Location of umac -- face center vs face centroid
    MLMG::Location m_umac_loc;

    // Location of divu (RHS -- optional) -- cell center vs cell centroid
    MLMG::Location m_divu_loc;

};

}

#endif
