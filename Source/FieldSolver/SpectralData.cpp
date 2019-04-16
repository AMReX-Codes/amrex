#include<WarpXComplex.H>
#include<AMReX_FabArray.H>

using namespace amrex;
using Complex = std::complex<Real>;

/* \brief Class that stores the fields in spectral space
 *  and performs the spectral transforms to/from real space
 */
class SpectralData
{
    using SpectralField = FabArray<BaseFab<Complex>>;
#ifdef AMREX_USE_GPU
    // Add cuFFT-specific code
#else
    using FFTplans = LayoutData<fftw_plan>;
#endif

    public:
        SpectralData( const BoxArray& ba, const DistributionMapping& dm );
        ~SpectralData();
        void ForwardTransform( const MultiFab& mf, const std::string& field_name );
        void InverseTransform( MultiFab& mf, const std::string& field_name );

    private:
        SpectralField Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new;
        SpectralField tmpRealField, tmpSpectralField; // Store fields before/after transform
        FFTplans forward_plan, inverse_plan;
};


SpectralData::SpectralData( const BoxArray& ba, const DistributionMapping& dm )
{
    // Allocate the arrays that contain the fields in spectral space
    Ex = SpectralField(ba, dm, 1, 0);
    Ey = SpectralField(ba, dm, 1, 0);
    Ez = SpectralField(ba, dm, 1, 0);
    Bx = SpectralField(ba, dm, 1, 0);
    By = SpectralField(ba, dm, 1, 0);
    Bz = SpectralField(ba, dm, 1, 0);
    Jx = SpectralField(ba, dm, 1, 0);
    Jy = SpectralField(ba, dm, 1, 0);
    Jz = SpectralField(ba, dm, 1, 0);
    rho_old = SpectralField(ba, dm, 1, 0);
    rho_new = SpectralField(ba, dm, 1, 0);

    // Allocate temporary arrays
    tmpRealField = SpectralField(ba, dm, 1, 0);
    tmpSpectralField = SpectralField(ba, dm, 1, 0);

    // Allocate and initialize the FFT plans
    forward_plan = FFTplans(ba, dm);
    inverse_plan = FFTplans(ba, dm);
    for ( MFIter mfi(ba, dm); mfi.isValid(); ++mfi ){
        Box bx = ba[mfi];
#ifdef AMREX_USE_GPU
        // Add cuFFT-specific code
#else
        // Create FFTW plans
        forward_plan[mfi] = fftw_plan_dft_3d(
            // Swap dimensions: AMReX data is Fortran-order, but FFTW is C-order
            bx.length(2), bx.length(1), bx.length(0),
            reinterpret_cast<fftw_complex*>( tmpRealField[mfi].dataPtr() ),
            reinterpret_cast<fftw_complex*>( tmpSpectralField[mfi].dataPtr() ),
            FFTW_FORWARD, FFTW_ESTIMATE );
        inverse_plan[mfi] = fftw_plan_dft_3d(
            // Swap dimensions: AMReX data is Fortran-order, but FFTW is C-order
            bx.length(2), bx.length(1), bx.length(0),
            reinterpret_cast<fftw_complex*>( tmpSpectralField[mfi].dataPtr() ),
            reinterpret_cast<fftw_complex*>( tmpRealField[mfi].dataPtr() ),
            FFTW_BACKWARD, FFTW_ESTIMATE );
        // TODO: Add 2D code
        // TODO: Do real-to-complex transform instead of complex-to-complex
        // TODO: Let the user decide whether to use FFTW_ESTIMATE or FFTW_MEASURE
#endif
    }
}

SpectralData::~SpectralData()
{
    for ( MFIter mfi(tmpRealField); mfi.isValid(); ++mfi ){
#ifdef AMREX_USE_GPU
        // Add cuFFT-specific code
#else
        // Destroy FFTW plans
        fftw_destroy_plan( forward_plan[mfi] );
        fftw_destroy_plan( inverse_plan[mfi] );
#endif
    }
}
