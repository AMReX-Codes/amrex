
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
        SpectralCoefficients C_coef, S_coef;
        Real dt;
};

/*
 * ba: BoxArray for spectral space
 * dm: DistributionMapping for spectral space
 */
PsatdSolver::PsatdSolver( const BoxArray& ba, const DistributionMapping& dm,
                          const Real* dx, const Real dt ) : dt(dt)
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
    S_coef = SpectralMatrix( ba, dm, 1, 0 );
    // Fill them with the right values
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){
        FillCoefficients( C_coef[mfi], S_coef[mfi], kx[mfi], ky[mfi], kz[mfi] );
    }
}

void
PsatdSolver::pushSpectralFields( SpectralFields& f ) const{


    // Loop over boxes
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){

        const Box& bx = mfi.box();

        // Extract arrays for the fields to be updated
        Array4<Complex> Ex = f.Ex[mfi].array();
        Array4<Complex> Ey = f.Ey[mfi].array();
        Array4<Complex> Ez = f.Ez[mfi].array();
        Array4<Complex> Bx = f.Bx[mfi].array();
        Array4<Complex> By = f.By[mfi].array();
        Array4<Complex> Bz = f.Bz[mfi].array();
        // Extract arrays for J
        const Array4<Complex> Jx_arr = f.Jx[mfi].array();
        const Array4<Complex> Jy_arr = f.Jy[mfi].array();
        const Array4<Complex> Jz_arr = f.Jz[mfi].array();
        const Array4<Complex> rho_old = f.rho_old[mfi].array();
        const Array4<Complex> rho_new = f.rho_new[mfi].array();
        // Extract arrays for the coefficients
        const Array4<Real> C_arr = C_coef[mfi].array();
        const Array4<Real> S_over_k_arr = S_over_k_coef[mfi].array();
        const Array4<Real> inv_k2 =
        // Extract pointers for the k vectors
        const Real* kx_arr = kx[mfi].dataPtr();
        const Real* ky_arr = ky[mfi].dataPtr();
        const Real* kz_arr = kz[mfi].dataPtr();

        ParallelFor( bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // Record old values of the fields to be updated
            const Complex Ex_old = Ex(i,j,k);
            const Complex Ey_old = Ey(i,j,k);
            const Complex Ez_old = Ez(i,j,k);
            const Complex Bx_old = Bx(i,j,k);
            const Complex By_old = By(i,j,k);
            const Complex Bz_old = Bz(i,j,k);
            // k vector values, and coefficients
            const Real kx = kx_arr[i];
            const Real ky = ky_arr[j];
            const Real kz = kz_arr[k];
            const Real C = C_arr(i,j,k);
            const Real S_over_k = S_over_k_arr(i,j,k);
            const Real inv_k2 = inv_k2_arr(i,j,k);
            const I = Complex{0,1};
            // Calculate divE and divJ
            const k_dot_E = (kx*Ex_old + ky*Ey_old + kz*Ez_old);
            const k_dot_J = (kx*Jx + ky*Jy + kz*Jz);

            // Update E
            Ex(i,j,k) = C*Ex_old + I*S_over_k*( ky*Bz_old - kz*By_old )
                        - S_over_k*Jx
                        + (1.-C)*inv_k2*k_dot_E*kx
                        + (S_over_k - dt)*inv_k2*k_dot_J*kx;



            // Update B


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
