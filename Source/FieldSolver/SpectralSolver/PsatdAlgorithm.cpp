#include <PsatdAlgorithm.H>
#include <WarpXConst.H>
#include <cmath>

using namespace amrex;

PsatdAlgorithm::PsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                         const DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const Real dt)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate the 1D vectors
    modified_kx_vec = SpectralKVector( ba, dm );
    modified_ky_vec = SpectralKVector( ba, dm );
    modified_kz_vec = SpectralKVector( ba, dm );
    // Allocate and fill them by computing the modified vector
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){
        Box bx = ba[mfi];
        ComputeModifiedKVector(
            modified_kx_vec[mfi], spectral_kspace.kx_vec[mfi],
            bx, spectral_kspace.dx[0], norder_x );
        ComputeModifiedKVector(
            modified_ky_vec[mfi], spectral_kspace.ky_vec[mfi],
            bx, spectral_kspace.dx[1], norder_y );
        ComputeModifiedKVector(
            modified_kz_vec[mfi], spectral_kspace.kz_vec[mfi],
            bx, spectral_kspace.dx[2], norder_z );
    }

    // Allocate the arrays of coefficients
    C_coef = SpectralCoefficients( ba, dm, 1, 0 );
    S_ck_coef = SpectralCoefficients( ba, dm, 1, 0 );
    X1_coef = SpectralCoefficients( ba, dm, 1, 0 );
    X2_coef = SpectralCoefficients( ba, dm, 1, 0 );
    X3_coef = SpectralCoefficients( ba, dm, 1, 0 );

    // Fill them with the right values:
    // Loop over boxes
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){

        const Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const Real* modified_kx = modified_kx_vec[mfi].dataPtr();
        const Real* modified_ky = modified_ky_vec[mfi].dataPtr();
        const Real* modified_kz = modified_kz_vec[mfi].dataPtr();
        // Extract arrays for the coefficients
        Array4<Real> C = C_coef[mfi].array();
        Array4<Real> S_ck = S_ck_coef[mfi].array();
        Array4<Real> X1 = X1_coef[mfi].array();
        Array4<Real> X2 = X2_coef[mfi].array();
        Array4<Real> X3 = X3_coef[mfi].array();

        // Loop over indices within one box
        ParallelFor( bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of vector
            const Real k_norm = std::sqrt(
                std::pow( modified_kx[i], 2 ) +
                std::pow( modified_ky[j], 2 ) +
                std::pow( modified_kz[k], 2 ) );

            // Calculate coefficients
            constexpr Real c = PhysConst::c;
            constexpr Real ep0 = PhysConst::ep0;
            if ( k_norm != 0 ){
                C(i,j,k) = std::cos( c*k_norm*dt );
                S_ck(i,j,k) = std::sin( c*k_norm*dt )/( c*k_norm );
                X1(i,j,k) = (1. - C(i,j,k))/(ep0 * c*c * k_norm*k_norm);
                X2(i,j,k) = (1. - S_ck(i,j,k)/dt )/(ep0 * k_norm*k_norm);
                X3(i,j,k) = (C(i,j,k) - S_ck(i,j,k)/dt )/(ep0 * k_norm*k_norm);
            } else { // Handle k_norm = 0, by using the analytical limit
                C(i,j,k) = 1.;
                S_ck(i,j,k) = dt;
                X1(i,j,k) = 0.5 * dt*dt / ep0;
                X2(i,j,k) = c*c * dt*dt / (6.*ep0);
                X3(i,j,k) = - c*c * dt*dt / (3.*ep0);
            }
        });
    }
};

void
PsatdAlgorithm::pushSpectralFields( SpectralData& f ) const{

    // Loop over boxes
    for ( MFIter mfi(f.Ex); mfi.isValid(); ++mfi ){

        const Box& bx = f.Ex[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> Ex_arr = f.Ex[mfi].array();
        Array4<Complex> Ey_arr = f.Ey[mfi].array();
        Array4<Complex> Ez_arr = f.Ez[mfi].array();
        Array4<Complex> Bx_arr = f.Bx[mfi].array();
        Array4<Complex> By_arr = f.By[mfi].array();
        Array4<Complex> Bz_arr = f.Bz[mfi].array();
        // Extract arrays for J and rho
        Array4<const Complex> Jx_arr = f.Jx[mfi].array();
        Array4<const Complex> Jy_arr = f.Jy[mfi].array();
        Array4<const Complex> Jz_arr = f.Jz[mfi].array();
        Array4<const Complex> rho_old_arr = f.rho_old[mfi].array();
        Array4<const Complex> rho_new_arr = f.rho_new[mfi].array();
        // Extract arrays for the coefficients
        Array4<const Real> C_arr = C_coef[mfi].array();
        Array4<const Real> S_ck_arr = S_ck_coef[mfi].array();
        Array4<const Real> X1_arr = X1_coef[mfi].array();
        Array4<const Real> X2_arr = X2_coef[mfi].array();
        Array4<const Real> X3_arr = X3_coef[mfi].array();
        // Extract pointers for the k vectors
        const Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
        const Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
        const Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor( bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Record old values of the fields to be updated
            const Complex Ex_old = Ex_arr(i,j,k);
            const Complex Ey_old = Ey_arr(i,j,k);
            const Complex Ez_old = Ez_arr(i,j,k);
            const Complex Bx_old = Bx_arr(i,j,k);
            const Complex By_old = By_arr(i,j,k);
            const Complex Bz_old = Bz_arr(i,j,k);
            // Shortcut for the values of J and rho
            const Complex Jx = Jx_arr(i,j,k);
            const Complex Jy = Jy_arr(i,j,k);
            const Complex Jz = Jz_arr(i,j,k);
            const Complex rho_old = rho_old_arr(i,j,k);
            const Complex rho_new = rho_new_arr(i,j,k);
            // k vector values, and coefficients
            const Real kx = modified_kx_arr[i];
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
            constexpr Real c2 = PhysConst::c*PhysConst::c;
            constexpr Real inv_ep0 = 1./PhysConst::ep0;
            constexpr Complex I = Complex{0,1};
            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);
            const Real X1 = X1_arr(i,j,k);
            const Real X2 = X2_arr(i,j,k);
            const Real X3 = X3_arr(i,j,k);

            // Update E (see WarpX online documentation: theory section)
            Ex_arr(i,j,k) = C*Ex_old
                        + S_ck*( c2*I*(ky*Bz_old - kz*By_old) - inv_ep0*Jx )
                        - I*( X2*rho_new - X3*rho_old )*kx;
            Ey_arr(i,j,k) = C*Ey_old
                        + S_ck*( c2*I*(kz*Bx_old - kx*Bz_old) - inv_ep0*Jy )
                        - I*( X2*rho_new - X3*rho_old )*ky;
            Ez_arr(i,j,k) = C*Ez_old
                        + S_ck*( c2*I*(kx*By_old - ky*Bx_old) - inv_ep0*Jz )
                        - I*( X2*rho_new - X3*rho_old )*kz;
            // Update B (see WarpX online documentation: theory section)
            Bx_arr(i,j,k) = C*Bx_old
                        - S_ck*I*(ky*Ez_old - kz*Ey_old)
                        +   X1*I*(ky*Jz     - kz*Jy );
            By_arr(i,j,k) = C*By_old
                        - S_ck*I*(kz*Ex_old - kx*Ez_old)
                        +   X1*I*(kz*Jx     - kx*Jz );
            Bz_arr(i,j,k) = C*Bz_old
                        - S_ck*I*(kx*Ey_old - ky*Ex_old)
                        +   X1*I*(kx*Jy     - ky*Jx );
        });
    }
};
