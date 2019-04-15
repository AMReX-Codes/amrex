
using namespace amrex;
using namespace Gpu;


class PsatdSolver
{
    using SpectralMatrix = FabArray<BaseFab<Real>>
    using SpectralVector = LayoutData<ManagedVector<Real>>

    public:
        PsatdSolver( BoxArray ba, DistributionMapping dm, Real* dx );
        pushSpectralFields( SpectralData data );

    private:
        SpectralVector kx, ky, kz;
        SpectralMatrix C, S;
}

/*
 * ba: BoxArray for spectral space ( index i corresponds to k=pi*i/dx )
 * dm: DistributionMapping for spectral space
 */
PsatdSolver::PsatSolver( BoxArray ba, DistributionMapping dm, Real* dx )
{

    // Allocate the 1D vectors
    kx = SpectralVector( ba, dm );
    ky = SpectralVector( ba, dm );
    kz = SpectralVector( ba, dm );
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){
        kx[mfi].resize( ba[mfi].length(0) );
        ky[mfi].resize( ba[mfi].length(1) );
        kz[mfi].resize( ba[mfi].length(2) );
    }

    // Allocate the matrix of coefficients
    C = SpectralMatrix( ba, dm, 1, 0 );
    S = SpectralMatrix( ba, dm, 1, 0 );

    // Fill vector and matrices with the right values
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){
        kx[mfi] = ...
        ky[mfi] = ...
        kz[mfi] = ...
        C[mfi] = ...
        S[mfi] = ...
    }
}
