#include<WarpXComplex.H>
#include<AMReX_FabArray.H>

using namespace amrex;

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
        SpectralData( const BoxArray& realspace_ba,
                      const BoxArray& spectralspace_ba,
                      const DistributionMapping& dm );
        ~SpectralData();
        void ForwardTransform( const MultiFab& mf, const std::string& field_name );
        void InverseTransform( MultiFab& mf, const std::string& field_name );

    private:
        SpectralField Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new;
        SpectralField tmpRealField, tmpSpectralField; // Store fields before/after transform
        FFTplans forward_plan, inverse_plan;
};


SpectralData::SpectralData( const BoxArray& realspace_ba,
                            const BoxArray& spectralspace_ba,
                            const DistributionMapping& dm )
{
    // Allocate the arrays that contain the fields in spectral space
    Ex = SpectralField(spectralspace_ba, dm, 1, 0);
    Ey = SpectralField(spectralspace_ba, dm, 1, 0);
    Ez = SpectralField(spectralspace_ba, dm, 1, 0);
    Bx = SpectralField(spectralspace_ba, dm, 1, 0);
    By = SpectralField(spectralspace_ba, dm, 1, 0);
    Bz = SpectralField(spectralspace_ba, dm, 1, 0);
    Jx = SpectralField(spectralspace_ba, dm, 1, 0);
    Jy = SpectralField(spectralspace_ba, dm, 1, 0);
    Jz = SpectralField(spectralspace_ba, dm, 1, 0);
    rho_old = SpectralField(spectralspace_ba, dm, 1, 0);
    rho_new = SpectralField(spectralspace_ba, dm, 1, 0);

    // Allocate temporary arrays - over different boxarrays
    tmpRealField = SpectralField(realspace_ba, dm, 1, 0);
    tmpSpectralField = SpectralField(spectralspace_ba, dm, 1, 0);

    // Allocate and initialize the FFT plans
    forward_plan = FFTplans(spectralspace_ba, dm);
    inverse_plan = FFTplans(spectralspace_ba, dm);
    for ( MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi ){
        Box bx = spectralspace_ba[mfi];
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


void
ForwardTransform( const MultiFab& mf, const std::string& field_name )
{
    // Loop over boxes
    for ( MFIter mfi(mf); mfi.isValid(); ++mfi ){

        // Copy the real-space field `mf` to the temporary field `tmpRealField`
        // This ensures that all fields have the same number of points
        // before the Fourier transform.
        // As a consequence, the copy discards the *last* point of `mf`
        // in any direction that has *nodal* index type.
        const Box realspace_bx = mf[mfi].box().enclosedCells(); // discards last point in each nodal direction
        AMREX_ALWAYS_ASSERT( realspace_bx == tmpRealField[mfi].box() );
        const Array4<Real> mf_arr = mf[mfi].array();
        Array4<Complex> tmp_arr = tmpRealField[mfi].array();
        ParallelFor( realspace_bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tmp_arr(i,j,k) = mf_arr(i,j,k);
        });

        // Perform Fourier transform from `tmpRealField` to `tmpSpectralField`
#ifdef AMREX_USE_GPU
        // Add cuFFT-specific code ; make sure that this is done on the same
        // GPU stream as the above copy
#else
        fftw_execute( forward_plan[mfi] );
#endif

        // Copy the spectral-space field `tmpSpectralField` to the appropriate field
        const Box spectralspace_bx = tmpSpectralField[mfi].box();
        const Array4<Complex> tmp_arr = tmpSpectralField[mfi].array();
        Array4<Complex> field_arr = Ex[mfi].array(); // TODO: Pick the right field
        ParallelFor( spectralspace_bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            field_arr(i,j,k) = tmp_arr(i,j,k);
        });
    }
}
