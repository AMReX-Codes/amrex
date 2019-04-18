#include <WarpXConst.H>
#include <SpectralKSpace.H>

SpectralKSpace::SpectralKSpace( const BoxArray& realspace_ba,
                                const DistributionMapping& dm,
                                const Real* realspace_dx )
{
    // Create the box array that corresponds to spectral space
    BoxList spectral_bl; // Create empty box list
    // Loop over boxes and fill the box list
    for (int i=0; i < realspace_ba.size(); i++ ) {
        // For local FFTs, each box in spectral space starts at 0 in each direction
        // and has the same number of points as the real space box (including guard cells)
        Box realspace_bx = realspace_ba[i];
        Box bx = Box( IntVect::TheZeroVector(), realspace_bx.bigEnd() - realspace_bx.smallEnd() );
        spectral_bl.push_back( bx );
    }
    spectralspace_ba.define( spectral_bl );

    // Allocate the 1D vectors
    kx_vec = SpectralKVector( spectralspace_ba, dm );
    ky_vec = SpectralKVector( spectralspace_ba, dm );
    kz_vec = SpectralKVector( spectralspace_ba, dm );
    // Initialize the values on each box
    for ( MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi ){
        Box bx = spectralspace_ba[mfi];
        AllocateAndFillKvector( kx_vec[mfi], bx, dx, 0 );
        AllocateAndFillKvector( ky_vec[mfi], bx, dx, 1 );
        AllocateAndFillKvector( kz_vec[mfi], bx, dx, 2 );
    }

    // Store the cell size
    dx = realspace_dx;
}

void
AllocateAndFillKvector( ManagedVector<Real>& k, const Box& bx, const Real* dx, const int i_dim )
{
    // Alllocate k to the right size
    int N = bx.length( i_dim );
    k.resize( N );

    // Fill the k vector
    const Real dk = 2*MathConst::pi/(N*dx[i_dim]);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( bx.smallEnd(i_dim) == 0,
        "Expected box to start at 0, in spectral space.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( bx.bigEnd(i_dim) == N-1,
        "Expected different box end index in spectral space.");
    // Fill positive values of k (FFT conventions: first half is positive)
    for (int i=0; i<(N+1)/2; i++ ){
        k[i] = i*dk;
    }
    // Fill negative values of k (FFT conventions: second half is negative)
    for (int i=(N+1)/2; i<N; i++){
        k[i] = (N-i)*dk;
    }
    // TODO: This should be quite different for the hybrid spectral code:
    // In that case we should take into consideration the actual indices of the box
    // and distinguish the size of the local box and that of the global FFT
    // TODO: For real-to-complex,

}

void
ComputeModifiedKVector( ManagedVector<Real>& modified_k,
                        const ManagedVector<Real>& k,
                        const Box& bx, const Real dx, const int norder )
{
    // Allocate modified_k to the right size
    int N = k.size();
    modified_k.resize( N );

    // For now, this simply copies the infinite order k
    for (int i=0; i<N; i++ ){
        modified_k[i] = k[i];
    }

}
