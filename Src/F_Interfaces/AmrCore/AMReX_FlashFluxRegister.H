#ifndef AMREX_FLASH_FLUX_REGISTER_H_
#define AMREX_FLASH_FLUX_REGISTER_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_OpenMP.H>

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
    void store (int fine_global_index, int dir, FArrayBox const& fine_flux, Real sf);

    // flux_in_register = scaling_factor * \sum{fine_flux * area}
    void store (int fine_global_index, int dir, FArrayBox const& fine_flux, FArrayBox const& area,
                Real sf);

    // flux_in_register = scaling_factor * \sum{fine_flux * area}, if the component is flux density
    //                    scaling_factor * \sum{fine_flux}       , otherwise
    void store (int fine_global_index, int dir, FArrayBox const& fine_flux, FArrayBox const& area,
                const int* isFluxDensity, Real sf);

    // flux_in_register += scaling_factor * \sum{fine_flux} / (ref_ratio[0]*ref_ratio[1]*ref_ratio[2])
    void add (int fine_global_index, int dir, FArrayBox const& fine_flux, Real sf);

    // flux_in_register += scaling_factor * \sum{fine_flux * area}
    void add (int fine_global_index, int dir, FArrayBox const& fine_flux, FArrayBox const& area,
                Real sf);

    // flux_in_register += scaling_factor * \sum{fine_flux * area}, if the component is flux density
    //                     scaling_factor * \sum{fine_flux}       , otherwise
    void add (int fine_global_index, int dir, FArrayBox const& fine_flux, FArrayBox const& area,
                const int* isFluxDensity, Real sf);

    void communicate ();

    // crse_flux = flux_in_register * scaling_factor
    void load (int crse_global_index, int dir, FArrayBox& crse_flux, Real sf) const;

    // crse_flux = flux_in_register * sf_f + cflux * sf_c
    void load (int crse_global_index, int dir, FArrayBox& crse_flux, FArrayBox const& cflux,
               Real sf_f, Real sf_c) const;

    // crse_flux = flux_in_register / area
    void load (int crse_global_index, int dir, FArrayBox& crse_flux, FArrayBox const& area) const;

    // crse_flux = flux_in_register/area * sf_f + cflux * sf_c
    void load (int crse_global_index, int dir, FArrayBox& crse_flux, FArrayBox const& cflux,
               FArrayBox const& area, Real sf_f, Real sf_c) const;

    // crse_flux = flux_in_register/area * sf_f + cflux * sf_c, if the component is flux density
    //             flux_in_register      * sf_f + cflux * sf_c, otherwise
    void load (int crse_global_index, int dir, FArrayBox& crse_flux, FArrayBox const& cflux,
               FArrayBox const& area, const int* isFluxDensity, Real sf_f, Real sf_c) const;

    enum struct OpType { Store, Add };

    template <OpType op>
    void store_or_add (int fine_global_index, int dir, FArrayBox const& fine_flux, Real sf);

    template <OpType op>
    void store_or_add (int fine_global_index, int dir, FArrayBox const& fine_flux, FArrayBox const& area,
                       Real sf);

    template <OpType op>
    void store_or_add (int fine_global_index, int dir, FArrayBox const& fine_flux, FArrayBox const& area,
                       const int* isFluxDensity, Real sf);


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

template <FlashFluxRegister::OpType op>
void FlashFluxRegister::store_or_add (int fine_global_index, int dir,
                                      FArrayBox const& fine_flux, Real sf)
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_fine_map.find(fine_global_index);
    if (found != m_fine_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,AMREX_SPACEDIM> const& fab_a = found->second;
        if (fab_a[dir]) {
            Box const& b = fab_a[dir]->box();
            Array4<Real> const& dest = fab_a[dir]->array();
            Array4<Real const> const& src = fine_flux.const_array();
            if (dir == 0) {
#if (AMREX_SPACEDIM == 1)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    amrex::ignore_unused(j,k,dest);
                    auto rhs = src(2*i,0,0,n)*sf;
                    if constexpr (op == OpType::Store) {
                        dest(i,0,0,n) = rhs;
                    } else {
                        dest(i,0,0,n) += rhs;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    amrex::ignore_unused(k,dest);
                    auto rhs = (src(2*i,2*j  ,0,n) +
                                src(2*i,2*j+1,0,n)) * (Real(0.5)*sf);
                    if constexpr (op == OpType::Store) {
                        dest(i,j,0,n) = rhs;
                    } else {
                        dest(i,j,0,n) += rhs;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    auto rhs = (src(2*i,2*j  ,2*k  ,n) +
                                src(2*i,2*j+1,2*k  ,n) +
                                src(2*i,2*j  ,2*k+1,n) +
                                src(2*i,2*j+1,2*k+1,n)) * (Real(0.25)*sf);
                    amrex::ignore_unused(dest); // for cuda
                    if constexpr (op == OpType::Store) {
                        dest(i,j,k,n) = rhs;
                    } else {
                        dest(i,j,k,n) += rhs;
                    }
                });
#endif
            }
#if (AMREX_SPACEDIM >= 2)
            else if (dir == 1) {
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    amrex::ignore_unused(k,dest);
                    auto rhs = (src(2*i  ,2*j,0,n) +
                                src(2*i+1,2*j,0,n)) * (Real(0.5)*sf);
                    if constexpr (op == OpType::Store) {
                        dest(i,j,0,n) = rhs;
                    } else {
                        dest(i,j,0,n) += rhs;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    auto rhs = (src(2*i  ,2*j,2*k  ,n) +
                                src(2*i+1,2*j,2*k  ,n) +
                                src(2*i  ,2*j,2*k+1,n) +
                                src(2*i+1,2*j,2*k+1,n)) * (Real(0.25)*sf);
                    amrex::ignore_unused(dest); // for cuda
                    if constexpr (op == OpType::Store) {
                        dest(i,j,k,n) = rhs;
                    } else {
                        dest(i,j,k,n) += rhs;
                    }
                });
#endif
            }
#if (AMREX_SPACEDIM == 3)
            else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    auto rhs = (src(2*i  ,2*j  ,2*k,n) +
                                src(2*i+1,2*j  ,2*k,n) +
                                src(2*i  ,2*j+1,2*k,n) +
                                src(2*i+1,2*j+1,2*k,n)) * (Real(0.25)*sf);
                    amrex::ignore_unused(dest); // for cuda
                    if constexpr (op == OpType::Store) {
                        dest(i,j,k,n) = rhs;
                    } else {
                        dest(i,j,k,n) += rhs;
                    }
                });
            }
#endif
#endif
        }
    }
}

template <FlashFluxRegister::OpType op>
void FlashFluxRegister::store_or_add (int fine_global_index, int dir,
                                      FArrayBox const& fine_flux,
                                      FArrayBox const& fine_area, Real sf)
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_fine_map.find(fine_global_index);
    if (found != m_fine_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,AMREX_SPACEDIM> const& fab_a = found->second;
        if (fab_a[dir]) {
            Box const& b = fab_a[dir]->box();
            Array4<Real> const& dest = fab_a[dir]->array();
            Array4<Real const> const& src = fine_flux.const_array();
            Array4<Real const> const& area = fine_area.const_array();
            if (dir == 0) {
#if (AMREX_SPACEDIM == 1)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    amrex::ignore_unused(j,k,dest);
                    auto rhs = src(2*i,0,0,n)*area(2*i,0,0)*sf;
                    if constexpr (op == OpType::Store) {
                        dest(i,0,0,n) = rhs;
                    } else {
                        dest(i,0,0,n) += rhs;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    amrex::ignore_unused(k,dest);
                    auto rhs = (src(2*i,2*j  ,0,n)*area(2*i,2*j  ,0) +
                                src(2*i,2*j+1,0,n)*area(2*i,2*j+1,0)) * sf;
                    if constexpr (op == OpType::Store) {
                        dest(i,j,0,n) = rhs;
                    } else {
                        dest(i,j,0,n) += rhs;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    auto rhs = (src(2*i,2*j  ,2*k  ,n)*area(2*i,2*j  ,2*k  ) +
                                src(2*i,2*j+1,2*k  ,n)*area(2*i,2*j+1,2*k  ) +
                                src(2*i,2*j  ,2*k+1,n)*area(2*i,2*j  ,2*k+1) +
                                src(2*i,2*j+1,2*k+1,n)*area(2*i,2*j+1,2*k+1)) * sf;
                    amrex::ignore_unused(dest); // for cuda
                    if constexpr (op == OpType::Store) {
                        dest(i,j,k,n) = rhs;
                    } else {
                        dest(i,j,k,n) += rhs;
                    }
                });
#endif
            }
#if (AMREX_SPACEDIM >= 2)
            else if (dir == 1) {
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    amrex::ignore_unused(k,dest);
                    auto rhs = (src(2*i  ,2*j,0,n)*area(2*i  ,2*j,0) +
                                src(2*i+1,2*j,0,n)*area(2*i+1,2*j,0)) * sf;
                    if constexpr (op == OpType::Store) {
                        dest(i,j,0,n) = rhs;
                    } else {
                        dest(i,j,0,n) += rhs;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    auto rhs = (src(2*i  ,2*j,2*k  ,n)*area(2*i  ,2*j,2*k  ) +
                                src(2*i+1,2*j,2*k  ,n)*area(2*i+1,2*j,2*k  ) +
                                src(2*i  ,2*j,2*k+1,n)*area(2*i  ,2*j,2*k+1) +
                                src(2*i+1,2*j,2*k+1,n)*area(2*i+1,2*j,2*k+1)) * sf;
                    amrex::ignore_unused(dest); // for cuda
                    if constexpr (op == OpType::Store) {
                        dest(i,j,k,n) = rhs;
                    } else {
                        dest(i,j,k,n) += rhs;
                    }
                });
#endif
            }
#if (AMREX_SPACEDIM == 3)
            else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    auto rhs = (src(2*i  ,2*j  ,2*k,n)*area(2*i  ,2*j  ,2*k) +
                                src(2*i+1,2*j  ,2*k,n)*area(2*i+1,2*j  ,2*k) +
                                src(2*i  ,2*j+1,2*k,n)*area(2*i  ,2*j+1,2*k) +
                                src(2*i+1,2*j+1,2*k,n)*area(2*i+1,2*j+1,2*k)) * sf;
                    amrex::ignore_unused(dest); // for cuda
                    if constexpr (op == OpType::Store) {
                        dest(i,j,k,n) = rhs;
                    } else {
                        dest(i,j,k,n) += rhs;
                    }
                });
            }
#endif
#endif
        }
    }
}

template <FlashFluxRegister::OpType op>
void FlashFluxRegister::store_or_add (int fine_global_index, int dir,
                                      FArrayBox const& fine_flux,
                                      FArrayBox const& fine_area,
                                      const int* isFluxDensity, Real sf)
{
    auto& h_ifd = m_h_ifd[OpenMP::get_thread_num()];
    auto& d_ifd = m_d_ifd[OpenMP::get_thread_num()];

    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_fine_map.find(fine_global_index);
    if (found != m_fine_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,AMREX_SPACEDIM> const& fab_a = found->second;
        if (fab_a[dir]) {
            bool allsame = true;
            for (int n = 0; n < m_ncomp; ++n) {
                if (h_ifd[n] != isFluxDensity[n]) {
                    allsame = false;
                    h_ifd[n] = isFluxDensity[n];
                }
            }
            if (d_ifd.empty()) {
                allsame = false;
                d_ifd.resize(m_ncomp);
            }
            if (! allsame) {
                Gpu::copyAsync(Gpu::HostToDevice(), h_ifd.begin(), h_ifd.end(), d_ifd.begin());
            }

            Box const& b = fab_a[dir]->box();
            Array4<Real> const& dest = fab_a[dir]->array();
            Array4<Real const> const& src = fine_flux.const_array();
            Array4<Real const> const& area = fine_area.const_array();
            const int* ifd = d_ifd.data();
            if (dir == 0) {
#if (AMREX_SPACEDIM == 1)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    amrex::ignore_unused(j,k,dest);
                    Real rhs;
                    if (ifd[n]) {
                        rhs = src(2*i,0,0,n)*area(2*i,0,0)*sf;
                    } else {
                        rhs = src(2*i,0,0,n)*sf;
                    }
                    if constexpr (op == OpType::Store) {
                        dest(i,0,0,n) = rhs;
                    } else {
                        dest(i,0,0,n) += rhs;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    amrex::ignore_unused(k,dest);
                    Real rhs;
                    if (ifd[n]) {
                        rhs = (src(2*i,2*j  ,0,n)*area(2*i,2*j  ,0) +
                               src(2*i,2*j+1,0,n)*area(2*i,2*j+1,0)) * sf;
                    } else {
                        rhs = (src(2*i,2*j  ,0,n) +
                               src(2*i,2*j+1,0,n)) * sf;
                    }
                    if constexpr (op == OpType::Store) {
                        dest(i,j,0,n) = rhs;
                    } else {
                        dest(i,j,0,n) += rhs;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    Real rhs;
                    if (ifd[n]) {
                        rhs = (src(2*i,2*j  ,2*k  ,n)*area(2*i,2*j  ,2*k  ) +
                               src(2*i,2*j+1,2*k  ,n)*area(2*i,2*j+1,2*k  ) +
                               src(2*i,2*j  ,2*k+1,n)*area(2*i,2*j  ,2*k+1) +
                               src(2*i,2*j+1,2*k+1,n)*area(2*i,2*j+1,2*k+1)) * sf;
                    } else {
                        rhs = (src(2*i,2*j  ,2*k  ,n) +
                               src(2*i,2*j+1,2*k  ,n) +
                               src(2*i,2*j  ,2*k+1,n) +
                               src(2*i,2*j+1,2*k+1,n)) * sf;
                    }
                    amrex::ignore_unused(dest); // for cuda
                    if constexpr (op == OpType::Store) {
                        dest(i,j,k,n) = rhs;
                    } else {
                        dest(i,j,k,n) += rhs;
                    }
                });
#endif
            }
#if (AMREX_SPACEDIM >= 2)
            else if (dir == 1) {
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    amrex::ignore_unused(k,dest);
                    Real rhs;
                    if (ifd[n]) {
                        rhs = (src(2*i  ,2*j,0,n)*area(2*i  ,2*j,0) +
                               src(2*i+1,2*j,0,n)*area(2*i+1,2*j,0)) * sf;
                    } else {
                        rhs = (src(2*i  ,2*j,0,n) +
                               src(2*i+1,2*j,0,n)) * sf;
                    }
                    if constexpr (op == OpType::Store) {
                        dest(i,j,0,n) = rhs;
                    } else {
                        dest(i,j,0,n) += rhs;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    Real rhs;
                    if (ifd[n]) {
                        rhs = (src(2*i  ,2*j,2*k  ,n)*area(2*i  ,2*j,2*k  ) +
                               src(2*i+1,2*j,2*k  ,n)*area(2*i+1,2*j,2*k  ) +
                               src(2*i  ,2*j,2*k+1,n)*area(2*i  ,2*j,2*k+1) +
                               src(2*i+1,2*j,2*k+1,n)*area(2*i+1,2*j,2*k+1)) * sf;
                    } else {
                        rhs = (src(2*i  ,2*j,2*k  ,n) +
                               src(2*i+1,2*j,2*k  ,n) +
                               src(2*i  ,2*j,2*k+1,n) +
                               src(2*i+1,2*j,2*k+1,n)) * sf;
                    }
                    amrex::ignore_unused(dest); // for cuda
                    if constexpr (op == OpType::Store) {
                        dest(i,j,k,n) = rhs;
                    } else {
                        dest(i,j,k,n) += rhs;
                    }
                });
#endif
            }
#if (AMREX_SPACEDIM == 3)
            else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    Real rhs;
                    if (ifd[n]) {
                        rhs = (src(2*i  ,2*j  ,2*k,n)*area(2*i  ,2*j  ,2*k) +
                               src(2*i+1,2*j  ,2*k,n)*area(2*i+1,2*j  ,2*k) +
                               src(2*i  ,2*j+1,2*k,n)*area(2*i  ,2*j+1,2*k) +
                               src(2*i+1,2*j+1,2*k,n)*area(2*i+1,2*j+1,2*k)) * sf;
                    } else {
                        rhs = (src(2*i  ,2*j  ,2*k,n) +
                               src(2*i+1,2*j  ,2*k,n) +
                               src(2*i  ,2*j+1,2*k,n) +
                               src(2*i+1,2*j+1,2*k,n)) * sf;
                    }
                    amrex::ignore_unused(dest); // for cuda
                    if constexpr (op == OpType::Store) {
                        dest(i,j,k,n) = rhs;
                    } else {
                        dest(i,j,k,n) += rhs;
                    }
                });
            }
#endif
#endif
        }
    }
}

}

#endif
