#ifndef AMREX_EB_DATA_H_
#define AMREX_EB_DATA_H_
#include <AMReX_Config.H>

#include <AMReX_EBCellFlag.H>
#include <AMReX_Random.H>

namespace amrex
{

enum struct EBData_t : int
{
    levelset,     // level set
    volfrac,      // volume fraction
    centroid,     // volume centroid
    bndrycent,    // boundary centroid
    bndrynorm,    // boundary normal
    bndryarea,    // boundary area
    AMREX_D_DECL(apx, apy, apz), // area fraction
    AMREX_D_DECL(fcx, fcy, fcz), // face centroid
    AMREX_D_DECL(ecx, ecy, ecz), // edge centroid
    cellflag      // EBCellFlag
};

struct EBData
{
    template <EBData_t T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto get (int i, int j, int k) const noexcept
    {
        if constexpr (T == EBData_t::cellflag) {
            return (*m_cell_flag)(i,j,k);
        } else {
            return m_real_data[static_cast<int>(T)](i,j,k);
        }
    }

    template <EBData_t T, std::enable_if_t<   T == EBData_t::centroid
                                           || T == EBData_t::bndrycent
                                           || T == EBData_t::bndrynorm
                                           AMREX_D_TERM(|| T==EBData_t::fcx,
                                                        || T==EBData_t::fcy,
                                                        || T==EBData_t::fcz), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto get (int i, int j, int k, int n) const noexcept
    {
        return m_real_data[static_cast<int>(T)](i,j,k,n);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuArray<Real,AMREX_SPACEDIM>
    randomPointOnEB (int i, int j, int k, RandomEngine const& engine) const
    {
        Real nx = this->get<EBData_t::bndrynorm>(i,j,k,0);
        Real ny = this->get<EBData_t::bndrynorm>(i,j,k,1);
        Real bcx = this->get<EBData_t::bndrycent>(i,j,k,0);
        Real bcy = this->get<EBData_t::bndrycent>(i,j,k,1);
        int dir = (std::abs(nx) >= std::abs(ny)) ? 0 : 1;
#if (AMREX_SPACEDIM == 2)
        auto f = [&] (Real n0, Real n1, Real bc0, Real bc1)
        {
            Real rn = amrex::Random(engine);
            if (n1 == 0) {
                return amrex::makeTuple(bc0, rn-Real(0.5));
            } else {
                Real nn = n0/n1; // Note that we have n0 >= n1. So nn != 0.
                Real ym = bc1+nn*(bc0-Real(-0.5)); // where x=-0.5 and EB intersects
                Real yp = bc1+nn*(bc0-Real( 0.5)); // where x= 0.5 and EB intersects
                Real ymin = std::min(ym,yp);
                Real ymax = std::max(ym,yp);
                ymin = std::max(ymin, Real(-0.5));
                ymax = std::min(ymax, Real( 0.5));
                Real y = rn/(ymax-ymin) + ymin;
                Real x = bc0 - (y-bc1)*n1/n0;
                return amrex::makeTuple(x,y);
            }
        };

        if (dir == 0) {
            auto [x,y] = f( nx, ny,
                           bcx,bcy);
            return GpuArray<Real,2>{x,y};
        } else {
            auto [y,x] = f( ny, nx,
                           bcy,bcx);
            return GpuArray<Real,2>{x,y};
        }
#else
        Real nz = this->get<EBData_t::bndrynorm>(i,j,k,2);
        Real bcz = this->get<EBData_t::bndrycent>(i,j,k,2);
        if (std::abs(nz) > std::abs(nx) && std::abs(nz) > std::abs(ny)) {
            dir = 2;
        }
        auto f = [&] (Real n0, Real n1, Real n2, Real bc0, Real bc1, Real bc2)
        { // Note that n0 >= n1 >= n2;
            if (n1 == 0 && n2 == 0) {
                return amrex::makeTuple(bc0,
                                        amrex::Random(engine)-Real(0.5),
                                        amrex::Random(engine)-Real(0.5));
            } else if (n2 == 0) {
                Real nn = n0/n1;
                Real ym = bc1+nn*(bc0-Real(-0.5));
                Real yp = bc1+nn*(bc0-Real( 0.5));
                Real ymin = std::min(ym,yp);
                Real ymax = std::max(ym,yp);
                ymin = std::max(ymin, Real(-0.5));
                ymax = std::min(ymax, Real( 0.5));
                Real y = amrex::Random(engine)/(ymax-ymin) + ymin;
                Real z = amrex::Random(engine) - Real(0.5);
                Real x = bc0 - ((y-bc1)*n1+(z-bc2)*n2)/n0;
                return amrex::makeTuple(x,y,z);
            } else {
                Real y0 = bc1 - ((Real(-0.5)-bc0)*n0+(Real(-0.5)-bc2)*n2)/n1;
                Real y1 = bc1 - ((Real( 0.5)-bc0)*n0+(Real(-0.5)-bc2)*n2)/n1;
                Real y2 = bc1 - ((Real(-0.5)-bc0)*n0+(Real( 0.5)-bc2)*n2)/n1;
                Real y3 = bc1 - ((Real( 0.5)-bc0)*n0+(Real( 0.5)-bc2)*n2)/n1;
                Real ymin = std::min({y0,y1,y2,y3});
                Real ymax = std::max({y0,y1,y2,y3});
                ymin = std::max(ymin, Real(-0.5));
                ymax = std::min(ymax, Real( 0.5));
                Real z0 = bc2 - ((Real(-0.5)-bc0)*n0+(Real(-0.5)-bc1)*n1)/n2;
                Real z1 = bc2 - ((Real( 0.5)-bc0)*n0+(Real(-0.5)-bc1)*n1)/n2;
                Real z2 = bc2 - ((Real(-0.5)-bc0)*n0+(Real( 0.5)-bc1)*n1)/n2;
                Real z3 = bc2 - ((Real( 0.5)-bc0)*n0+(Real( 0.5)-bc1)*n1)/n2;
                Real zmin = std::min({z0,z1,z2,z3});
                Real zmax = std::max({z0,z1,z2,z3});
                zmin = std::max(zmin, Real(-0.5));
                zmax = std::min(zmax, Real( 0.5));
                Real x, y, z;
                do {
                    y = amrex::Random(engine)/(ymax-ymin) + ymin;
                    z = amrex::Random(engine)/(zmax-zmin) + zmin;
                    x = bc0 - ((y-bc1)*n1+(z-bc2)*n2)/n0;
                } while (x > Real(0.5) || x < Real(-0.5));
                return amrex::makeTuple(x,y,z);
            }
        };
        if (dir == 0) {
            if (std::abs(ny) >= std::abs(nz)) {
                auto [x,y,z] = f( nx, ny, nz,
                                 bcx,bcy,bcz);
                return GpuArray<Real,3>{x, y, z};
            } else {
                auto [x,z,y] = f( nx, nz, ny,
                                 bcx,bcz,bcy);
                return GpuArray<Real,3>{x, y, z};
            }
        } else if (dir == 1) {
            if (std::abs(nx) >= std::abs(nz)) {
                auto [y,x,z] = f( ny, nx, nz,
                                 bcy,bcx,bcz);
                return GpuArray<Real,3>{x, y, z};
            } else {
                auto [y,z,x] = f( ny, nz, nx,
                                 bcy,bcz,bcx);
                return GpuArray<Real,3>{x, y, z};
            }
        } else {
            if (std::abs(nx) >= std::abs(ny)) {
                auto [z,x,y] = f( nz, nx, ny,
                                 bcz,bcx,bcy);
                return GpuArray<Real,3>{x, y, z};
            } else {
                auto [z,y,x] = f( nz, ny, nx,
                                 bcz,bcy,bcx);
                return GpuArray<Real,3>{x, y, z};
            }
        }
#endif
    }

    static constexpr int real_data_size = static_cast<int>(EBData_t::cellflag);

    Array4<EBCellFlag const> const* m_cell_flag = nullptr;
    Array4<Real const> const* m_real_data = nullptr;
};

}
#endif
