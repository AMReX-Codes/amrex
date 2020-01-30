#ifndef AMREX_FLASH_FLUX_REGISTER_H_
#define AMREX_FLASH_FLUX_REGISTER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>

#include <map>

namespace amrex {

class FlashFluxRegister
{
public:

    FlashFluxRegister () = default;

    FlashFluxRegister (const BoxArray& fba, const BoxArray& cba,
                       const DistributionMapping& fdm, const DistributionMapping& cdm,
                       const Geometry& fgeom, const Geometry& cgeom,
                       IntVect const& ref_ratio, int nvar);

    void define (const BoxArray& fba, const BoxArray& cba,
                 const DistributionMapping& fdm, const DistributionMapping& cdm,
                 const Geometry& fgeom, const Geometry& cgeom,
                 IntVect const& ref_ratio, int nvar);

    void store (int fine_global_index, int dir, FArrayBox const& fine_flux, Real scaling_factor);

    void communicate ();

    void load (int coarse_global_index, int dir, FArrayBox& crse_flux, Real scaling_factor) const;

protected:

    BoxArray m_fine_grids;
    BoxArray m_crse_grids;

    DistributionMapping m_fine_dmap;
    DistributionMapping m_crse_dmap;

    Geometry m_fine_geom;
    Geometry m_crse_geom;

    int m_ncomp;

    //
    Array<MultiFab,AMREX_SPACEDIM> m_fine_fluxes;
    Array<MultiFab,AMREX_SPACEDIM> m_crse_fluxes;
    std::map<int,Array<FArrayBox*,AMREX_SPACEDIM> > m_fine_map;
    std::map<int,Array<FArrayBox*,2*AMREX_SPACEDIM> > m_crse_map;
    
};

}

#endif
