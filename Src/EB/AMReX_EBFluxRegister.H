#ifndef AMREX_EBFLUXREGISTER_H_
#define AMREX_EBFLUXREGISTER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>

#include <array>

namespace amrex {

class EBFluxRegister
{
public:
    
    EBFluxRegister () {;}

    EBFluxRegister (const BoxArray& fba, const BoxArray& cba,
                    const DistributionMapping& fdm, const DistributionMapping& cdm,
                    const Geometry& fgeom, const Geometry& cgeom,
                    const IntVect& ref_ratio, int fine_lev, int nvar);

    void define (const BoxArray& fba, const BoxArray& cba,
                 const DistributionMapping& fdm, const DistributionMapping& cdm,
                 const Geometry& fgeom, const Geometry& cgeom,
                 const IntVect& ref_ratio, int fine_lev, int nvar);

    void reset ();

    void CrseAdd (const MFIter& mfi, const std::array<FArrayBox,AMREX_SPACEDIM>& flux,
                  const Real* dx, Real dt);

    void FineAdd (const MFIter& mfi, const std::array<FArrayBox,AMREX_SPACEDIM>& flux,
                  const Real* dx, Real dt);

    void Reflux (MultiFab& state);

private:

    enum {
        crse_cell = 0, crse_fine_boundary_cell, fine_cell
    };

    MultiFab m_crse_data;
    iMultiFab m_crse_flag;
    Array<int> m_crse_fab_flag;
    
    MultiFab m_cfpatch;                   // This is built on crse/fine patches
    Array<Array<FArrayBox*> > m_cfp_fab;  // The size of this is (# of local fine grids (# of crse/fine patches for that grid))

    Geometry m_fine_geom;
    Geometry m_crse_geom;

    IntVect m_ratio;
    int m_fine_level;
    int m_ncomp;
};

}

#endif
