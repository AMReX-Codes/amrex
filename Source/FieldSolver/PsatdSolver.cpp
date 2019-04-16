
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
 * ba: BoxArray for spectral space ( index i corresponds to k=2*pi*i/(n*dx) )
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
    int n = bx.length( i_dim );
    k.resize( n );

    // Fill with the right values
    for (int i=bx.smallEnd(i_dim); i<=bx.bigEnd(i_dim); i++ ){
        k[i] = ...
    }
}
