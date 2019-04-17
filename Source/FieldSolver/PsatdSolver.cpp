#include <WarpXConst.H>
#include <cmath>

using namespace amrex;
using namespace Gpu;


class PsatdSolver
{
    using SpectralCoefficients = FabArray<BaseFab<Real>>
    using SpectralVector = LayoutData<ManagedVector<Real>>

    public:
        PsatdSolver( const BoxArray& ba, const DistributionMapping& dm, const Real* dx );
        void pushSpectralFields( SpectralData& f ) const;

    private:
        SpectralVector kx, ky, kz;
        SpectralCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
};

/*
 * ba: BoxArray for spectral space
 * dm: DistributionMapping for spectral space
 */
PsatdSolver::PsatdSolver( const BoxArray& ba, const DistributionMapping& dm,
                          const Real* dx, const Real dt )
{
    // Allocate the 1D vectors
    kx = SpectralVector( ba, dm );
    ky = SpectralVector( ba, dm );
    kz = SpectralVector( ba, dm );
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){
        Box bx = ba[mfi];
        AllocateAndFillKvector( kx[mfi], bx, dx, 0 )
        AllocateAndFillKvector( ky[mfi], bx, dx, 1 )
        AllocateAndFillKvector( kz[mfi], bx, dx, 2 )
    }

    // Allocate the arrays of coefficients
    C_coef = SpectralMatrix( ba, dm, 1, 0 );
    S_ck_coef = SpectralMatrix( ba, dm, 1, 0 );
    X1_coef = SpectralMatrix( ba, dm, 1, 0 );
    X2_coef = SpectralMatrix( ba, dm, 1, 0 );
    X3_coef = SpectralMatrix( ba, dm, 1, 0 );

    // Fill them with the right values:
    // Loop over boxes
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){

        const Box& bx = mfi.box();

        // Extract pointers for the k vectors
        const Real* kx = kx[mfi].dataPtr();
        const Real* ky = ky[mfi].dataPtr();
        const Real* kz = kz[mfi].dataPtr();
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
            const Real k_norm = std::sqrt( kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k] );

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
}

void
PsatdSolver::pushSpectralFields( SpectralFields& f ) const{

    // Loop over boxes
    for ( MFIter mfi(f.Ex); mfi.isValid(); ++mfi ){

        const Box& bx = mfi.box();

        // Extract arrays for the fields to be updated
        Array4<Complex> Ex_arr = f.Ex[mfi].array();
        Array4<Complex> Ey_arr = f.Ey[mfi].array();
        Array4<Complex> Ez_arr = f.Ez[mfi].array();
        Array4<Complex> Bx_arr = f.Bx[mfi].array();
        Array4<Complex> By_arr = f.By[mfi].array();
        Array4<Complex> Bz_arr = f.Bz[mfi].array();
        // Extract arrays for J
        const Array4<Complex> Jx_arr = f.Jx[mfi].array();
        const Array4<Complex> Jy_arr = f.Jy[mfi].array();
        const Array4<Complex> Jz_arr = f.Jz[mfi].array();
        const Array4<Complex> rho_old_arr = f.rho_old[mfi].array();
        const Array4<Complex> rho_new_arr = f.rho_new[mfi].array();
        // Extract arrays for the coefficients
        const Array4<Real> C_arr = C_coef[mfi].array();
        const Array4<Real> S_ck_arr = S_ck_coef[mfi].array();
        const Array4<Real> inv_k2_arr =
        // Extract pointers for the k vectors
        const Real* kx_arr = kx[mfi].dataPtr();
        const Real* ky_arr = ky[mfi].dataPtr();
        const Real* kz_arr = kz[mfi].dataPtr();

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
            // k vector values, and coefficients
            const Real kx = kx_arr[i];
            const Real ky = ky_arr[j];
            const Real kz = kz_arr[k];
            constexpr Real c2 = PhysConst::c*PhysConst::c;
            constexpr Real inv_ep0 = 1./PhysConst::ep0;
            constexpr Complex I = Complex{0,1};
            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);
            const Real X1 = X1_arr(i,j,k);
            const Real X2 = X2_arr(i,j,k);
            const Real X3 = X3_arr(i,j,k);
            // Short cut for the values of J and rho
            const Complex Jx = Jx_arr(i,j,k);
            const Complex Jy = Jy_arr(i,j,k);
            const Complex Jz = Jz_arr(i,j,k);

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
                        +   X1*I*(ky*Jz_old - kz*Jy_old);
            By_arr(i,j,k) = C*By_old
                        - S_ck*I*(kz*Ex_old - kx*Ez_old)
                        +   X1*I*(kz*Jx_old - kx*Jz_old);
            Bz_arr(i,j,k) = C*Bz_old
                        - S_ck*I*(kx*Ey_old - ky*Ex_old)
                        +   X1*I*(kx*Jy_old - ky*Jx_old);
        });
    }
}

AllocateAndFillKvector( ManagedVector<Real>& k, const Box& bx, const Real* dx, const int i_dim )
{
    // Alllocate k to the right size
    int N = bx.length( i_dim );
    k.resize( N );

    // Fill the k vector
    const Real PI = std::atan(1.0)*4;
    const Real dk = 2*PI/(N*dx[i_dim]);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( bx.smallEnd(i_dim) == 0,
        "Expected box to start at 0, in spectral space.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( bx.bigEnd(i_dim) == N-1,
        "Expected different box end index in spectral space.");
    // Fill positive values of k (FFT conventions: first half is positive)
    for (int i=0; i<(N+1)/2; i++ ){
        k[i] = i*dk;
    }
    // Fill negative values of k (FFT conventions: second half is negative)
    for (int i=(N+1)/2, i<N; i++){
        k[i] = (N-i)*dk;
    }
    // TODO: This should be quite different for the hybrid spectral code:
    // In that case we should take into consideration the actual indices of the box
    // and distinguish the size of the local box and that of the global FFT
    // TODO: For real-to-complex,

}
