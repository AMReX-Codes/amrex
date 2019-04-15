
using namespace amrex;
using namespace Gpu;


class PsatdSolver
{
    using SpectralMatrix = FabArray<BaseFab<Real>>
    using SpectralVector = LayoutData<ManagedVector>

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
    kx = SpectralVector( ba, dm );
    ky = SpectralVector( ba, dm );
    kz = SpectralVector( ba, dm );

    C = SpectralMatrix( ba, dm, 1, 0 );
    S = SpectralMatrix( ba, dm, 1, 0 );

    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){
        kx[mfi] = ...
        ky[mfi] = ...
        kz[mfi] = ...
        C[mfi] = ...
        S[mfi] = ...
    }
}
