#include <PsatdAlgorithm.H>
#include <WarpXConst.H>
#include <cmath>

using namespace amrex;

/* \brief Initialize coefficients for the update equation */
PsatdAlgorithm::PsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                         const DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal, const Real dt)
     // Initialize members of base class
     : SpectralBaseAlgorithm( spectral_kspace, dm,
                              norder_x, norder_y, norder_z, nodal )
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate the arrays of coefficients
    C_coef = SpectralCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralCoefficients(ba, dm, 1, 0);
    X1_coef = SpectralCoefficients(ba, dm, 1, 0);
    X2_coef = SpectralCoefficients(ba, dm, 1, 0);
    X3_coef = SpectralCoefficients(ba, dm, 1, 0);

    // Fill them with the right values:
    // Loop over boxes and allocate the corresponding coefficients
    // for each box owned by the local MPI proc
    for (MFIter mfi(ba, dm); mfi.isValid(); ++mfi){

        const Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const Real* modified_kx = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz = modified_kz_vec[mfi].dataPtr();
        // Extract arrays for the coefficients
        Array4<Real> C = C_coef[mfi].array();
        Array4<Real> S_ck = S_ck_coef[mfi].array();
        Array4<Real> X1 = X1_coef[mfi].array();
        Array4<Real> X2 = X2_coef[mfi].array();
        Array4<Real> X3 = X3_coef[mfi].array();

        // Loop over indices within one box
        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of vector
            const Real k_norm = std::sqrt(
                std::pow(modified_kx[i], 2) +
#if (AMREX_SPACEDIM==3)
                std::pow(modified_ky[j], 2) +
                std::pow(modified_kz[k], 2));
#else
                std::pow(modified_kz[j], 2));
#endif

            // Calculate coefficients
            constexpr Real c = PhysConst::c;
            constexpr Real ep0 = PhysConst::ep0;
            if (k_norm != 0){
                C(i,j,k) = std::cos(c*k_norm*dt);
                S_ck(i,j,k) = std::sin(c*k_norm*dt)/(c*k_norm);
                X1(i,j,k) = (1. - C(i,j,k))/(ep0 * c*c * k_norm*k_norm);
                X2(i,j,k) = (1. - S_ck(i,j,k)/dt)/(ep0 * k_norm*k_norm);
                X3(i,j,k) = (C(i,j,k) - S_ck(i,j,k)/dt)/(ep0 * k_norm*k_norm);
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

/* Advance the E and B field in spectral space (stored in `f`)
 * over one time step */
void
PsatdAlgorithm::pushSpectralFields(SpectralFieldData& f) const{

    // Loop over boxes
    for (MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        const Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = f.fields[mfi].array();
        // Extract arrays for the coefficients
        Array4<const Real> C_arr = C_coef[mfi].array();
        Array4<const Real> S_ck_arr = S_ck_coef[mfi].array();
        Array4<const Real> X1_arr = X1_coef[mfi].array();
        Array4<const Real> X2_arr = X2_coef[mfi].array();
        Array4<const Real> X3_arr = X3_coef[mfi].array();
        // Extract pointers for the k vectors
        const Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Record old values of the fields to be updated
            using Idx = SpectralFieldIndex;
            const Complex Ex_old = fields(i,j,k,Idx::Ex);
            const Complex Ey_old = fields(i,j,k,Idx::Ey);
            const Complex Ez_old = fields(i,j,k,Idx::Ez);
            const Complex Bx_old = fields(i,j,k,Idx::Bx);
            const Complex By_old = fields(i,j,k,Idx::By);
            const Complex Bz_old = fields(i,j,k,Idx::Bz);
            // Shortcut for the values of J and rho
            const Complex Jx = fields(i,j,k,Idx::Jx);
            const Complex Jy = fields(i,j,k,Idx::Jy);
            const Complex Jz = fields(i,j,k,Idx::Jz);
            const Complex rho_old = fields(i,j,k,Idx::rho_old);
            const Complex rho_new = fields(i,j,k,Idx::rho_new);
            // k vector values, and coefficients
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0;
            const Real kz = modified_kz_arr[j];
#endif
            constexpr Real c2 = PhysConst::c*PhysConst::c;
            constexpr Real inv_ep0 = 1./PhysConst::ep0;
            constexpr Complex I = Complex{0,1};
            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);
            const Real X1 = X1_arr(i,j,k);
            const Real X2 = X2_arr(i,j,k);
            const Real X3 = X3_arr(i,j,k);

            // Update E (see WarpX online documentation: theory section)
            fields(i,j,k,Idx::Ex) = C*Ex_old
                        + S_ck*(c2*I*(ky*Bz_old - kz*By_old) - inv_ep0*Jx)
                        - I*(X2*rho_new - X3*rho_old)*kx;
            fields(i,j,k,Idx::Ey) = C*Ey_old
                        + S_ck*(c2*I*(kz*Bx_old - kx*Bz_old) - inv_ep0*Jy)
                        - I*(X2*rho_new - X3*rho_old)*ky;
            fields(i,j,k,Idx::Ez) = C*Ez_old
                        + S_ck*(c2*I*(kx*By_old - ky*Bx_old) - inv_ep0*Jz)
                        - I*(X2*rho_new - X3*rho_old)*kz;
            // Update B (see WarpX online documentation: theory section)
            fields(i,j,k,Idx::Bx) = C*Bx_old
                        - S_ck*I*(ky*Ez_old - kz*Ey_old)
                        +   X1*I*(ky*Jz     - kz*Jy);
            fields(i,j,k,Idx::By) = C*By_old
                        - S_ck*I*(kz*Ex_old - kx*Ez_old)
                        +   X1*I*(kz*Jx     - kx*Jz);
            fields(i,j,k,Idx::Bz) = C*Bz_old
                        - S_ck*I*(kx*Ey_old - ky*Ex_old)
                        +   X1*I*(kx*Jy     - ky*Jx);
        });
    }
};
