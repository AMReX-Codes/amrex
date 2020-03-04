#ifndef AMREX_FLASH_FLUX_REGISTER_H_
#define AMREX_FLASH_FLUX_REGISTER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>

#include <map>

//
// Currently ref_ratio must be 2.  If needed, we can support arbitrary factors.
//

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

    // flux_in_register = scaling_factor * \sum{fine_flux} / (ref_ratio[0]*ref_ratio[1]*ref_ratio[2])
    void store (int fine_global_index, int dir, FArrayBox const& fine_flux, Real scaling_factor);

    // flux_in_register = scaling_factor * \sum{fine_flux * area}
    void store (int fine_global_index, int dir, FArrayBox const& fine_flux, FArrayBox const& area,
                Real scaling_factor);

    // flux_in_register = scaling_factor * \sum{fine_flux * area}, if the component is flux density
    //                    scaling_factor * \sum{fine_flux}       , otherwise
    void store (int fine_global_index, int dir, FArrayBox const& fine_flux, FArrayBox const& area,
                const int* isFluxDensity, Real scaling_factor);

    void communicate ();

    // crse_flux = flux_in_register * scaling_factor
    void load (int coarse_global_index, int dir, FArrayBox& crse_flux, Real scaling_factor) const;

    // crse_flux = flux_in_register * sf_f + cflux * sf_c
    void load (int coarse_global_index, int dir, FArrayBox& crse_flux, FArrayBox const& cflux,
               Real sf_f, Real sf_c) const;

    // crse_flux = flux_in_register / area
    void load (int coarse_global_index, int dir, FArrayBox& crse_flux, FArrayBox const& area) const;

    // crse_flux = flux_in_register/area * sf_f + cflux * sf_c
    void load (int coarse_global_index, int dir, FArrayBox& crse_flux, FArrayBox const& cflux,
               FArrayBox const& area, Real sf_f, Real sf_c) const;

    // crse_flux = flux_in_register/area * sf_f + cflux * sf_c, if the component is flux density
    //             flux_in_register      * sf_f + cflux * sf_c, otherwise
    void load (int coarse_global_index, int dir, FArrayBox& crse_flux, FArrayBox const& cflux,
               FArrayBox const& area, const int* isFluxDensity, Real sf_f, Real sf_c) const;

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

    mutable Vector<Gpu::PinnedVector<int> > m_h_ifd;
    mutable Vector<Gpu::DeviceVector<int> > m_d_ifd;
};

}

#endif
