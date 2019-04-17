
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
        SpectralCoefficients C, S;
};

/*
 * ba: BoxArray for spectral space
 * dm: DistributionMapping for spectral space
 */
PsatdSolver::PsatdSolver( const BoxArray& ba, const DistributionMapping& dm, const Real* dx )
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
    C = SpectralMatrix( ba, dm, 1, 0 );
    S = SpectralMatrix( ba, dm, 1, 0 );
    // Fill them with the right values
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){
        FillCoefficients( C[mfi], S[mfi], kx[mfi], ky[mfi], kz[mfi] );
    }
}

void
PsatdSolver::pushSpectralFields( SpectralFields& f ) const{

    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){


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
