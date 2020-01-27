#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_CYLINDRICAL_YEE_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_CYLINDRICAL_YEE_H_

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>

struct CylindricalYeeAlgorithm {

    static void InitializeStencilCoefficients(
        std::array<amrex::Real,3>& cell_size,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_r,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_z ) {

        // Store the inverse cell size along each direction in the coefficients
        stencil_coefs_r.resize(1);
        stencil_coefs_r[0] = 1./cell_size[0];  // 1./dr
        stencil_coefs_z.resize(1);
        stencil_coefs_z[0] = 1./cell_size[2];  // 1./dz
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDr(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        amrex::Real inv_dr = coefs_r[0];
        return inv_dr*( F(i+1,j,k,comp) - F(i,j,k,comp) );
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDr(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        amrex::Real inv_dr = coefs_r[0];
        return inv_dr*( F(i,j,k,comp) - F(i-1,j,k,comp) );
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDz(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_z, int const n_coefs_z,
        int const i, int const j, int const k, int const comp ) {

        amrex::Real inv_dz = coefs_z[0];
        return inv_dz*( F(i,j+1,k,comp) - F(i,j,k,comp) );
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDz(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_z, int const n_coefs_z,
        int const i, int const j, int const k, int const comp ) {

        amrex::Real inv_dz = coefs_z[0];
        return inv_dz*( F(i,j,k,comp) - F(i,j-1,k,comp) );
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DivideByR(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const r, int const dr, int const m,
        int const i, int const j, int const k, int const comp) {

        if (r != 0) {
            return F(i,j,k,comp)/r;
        } else { // r==0 ; singularity when dividing by r
            if (m==1) {
                // For m==1, F is linear in r, for small r
                // Therefore, the formula below regularizes the singularity
                return F(i+1,j,k,comp)/dr;
            } else {
                return 0;
            }
        }
    };

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CYLINDRICAL_YEE_H_
