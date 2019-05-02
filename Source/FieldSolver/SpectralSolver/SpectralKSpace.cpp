#include <WarpXConst.H>
#include <SpectralKSpace.H>
#include <cmath>

using namespace amrex;
using namespace Gpu;

/* \brief Initialize k space object.
 *
 * \param realspace_ba Box array that corresponds to the decomposition
 * of the fields in real space (cell-centered ; includes guard cells)
 * \param dm Indicates which MPI proc owns which box, in realspace_ba.
 * \param realspace_dx Cell size of the grid in real space
 */
SpectralKSpace::SpectralKSpace( const BoxArray& realspace_ba,
                                const DistributionMapping& dm,
                                const RealVect realspace_dx )
    : dx(realspace_dx)  // Store the cell size as member `dx`
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        realspace_ba.ixType()==IndexType::TheCellType(),
        "SpectralKSpace expects a cell-centered box.");

    // Create the box array that corresponds to spectral space
    BoxList spectral_bl; // Create empty box list
    // Loop over boxes and fill the box list
    for (int i=0; i < realspace_ba.size(); i++ ) {
        // For local FFTs, boxes in spectral space start at 0 in
        // each direction and have the same number of points as the
        // (cell-centered) real space box
        // TODO: this will be different for the hybrid FFT scheme
        Box realspace_bx = realspace_ba[i];
        IntVect fft_size = realspace_bx.length();
        // Because the spectral solver uses real-to-complex FFTs, we only
        // need the positive k values along the fastest axis
        // (first axis for AMReX Fortran-order arrays) in spectral space.
        // This effectively reduces the size of the spectral space by half
        // see e.g. the FFTW documentation for real-to-complex FFTs
        IntVect spectral_bx_size = fft_size;
        spectral_bx_size[0] = fft_size[0]/2 + 1;
        // Define the corresponding box
        Box spectral_bx = Box( IntVect::TheZeroVector(),
                               spectral_bx_size - IntVect::TheUnitVector() );
        spectral_bl.push_back( spectral_bx );
    }
    spectralspace_ba.define( spectral_bl );

    // Allocate the components of the k vector: kx, ky (only in 3D), kz
    bool only_positive_k;
    for (int i_dim=0; i_dim<AMREX_SPACEDIM; i_dim++) {
        if (i_dim==0) {
            // Real-to-complex FFTs: first axis contains only the positive k
            only_positive_k = true;
        } else {
            only_positive_k = false;
        }
        k_vec[i_dim] = getKComponent(dm, realspace_ba, i_dim, only_positive_k);
    }
}

/* For each box, in `spectralspace_ba`, which is owned by the local MPI rank
 * (as indicated by the argument `dm`), compute the values of the
 * corresponding k coordinate along the dimension specified by `i_dim`
 */
KVectorComponent
SpectralKSpace::getKComponent( const DistributionMapping& dm,
                               const BoxArray& realspace_ba,
                               const int i_dim,
                               const bool only_positive_k ) const
{
    // Initialize an empty ManagedVector in each box
    KVectorComponent k_comp(spectralspace_ba, dm);
    // Loop over boxes and allocate the corresponding ManagedVector
    // for each box owned by the local MPI proc
    for ( MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi ){
        Box bx = spectralspace_ba[mfi];
        ManagedVector<Real>& k = k_comp[mfi];

        // Allocate k to the right size
        int N = bx.length( i_dim );
        k.resize( N );

        // Fill the k vector
        IntVect fft_size = realspace_ba[mfi].length();
        const Real dk = 2*MathConst::pi/(fft_size[i_dim]*dx[i_dim]);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( bx.smallEnd(i_dim) == 0,
            "Expected box to start at 0, in spectral space.");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( bx.bigEnd(i_dim) == N-1,
            "Expected different box end index in spectral space.");
        if (only_positive_k){
            // Fill the full axis with positive k values
            // (typically: first axis, in a real-to-complex FFT)
            for (int i=0; i<N; i++ ){
                k[i] = i*dk;
            }
        } else {
            const int mid_point = (N+1)/2;
            // Fill positive values of k
            // (FFT conventions: first half is positive)
            for (int i=0; i<mid_point; i++ ){
                k[i] = i*dk;
            }
            // Fill negative values of k
            // (FFT conventions: second half is negative)
            for (int i=mid_point; i<N; i++){
                k[i] = (i-N)*dk;
            }
        }
        // TODO: this will be different for the hybrid FFT scheme
    }
    return k_comp;
}

/* For each box, in `spectralspace_ba`, which is owned by the local MPI rank
 * (as indicated by the argument `dm`), compute the values of the
 * corresponding correcting "shift" factor, along the dimension
 * specified by `i_dim`.
 *
 * (By default, we assume the FFT is done from/to a nodal grid in real space
 * It the FFT is performed from/to a cell-centered grid in real space,
 * a correcting "shift" factor must be applied in spectral space.)
 */
SpectralShiftFactor
SpectralKSpace::getSpectralShiftFactor( const DistributionMapping& dm,
                                        const int i_dim,
                                        const int shift_type ) const
{
    // Initialize an empty ManagedVector in each box
    SpectralShiftFactor shift_factor( spectralspace_ba, dm );
    // Loop over boxes and allocate the corresponding ManagedVector
    // for each box owned by the local MPI proc
    for ( MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi ){
        const ManagedVector<Real>& k = k_vec[i_dim][mfi];
        ManagedVector<Complex>& shift = shift_factor[mfi];

        // Allocate shift coefficients
        shift.resize( k.size() );

        // Fill the shift coefficients
        Real sign = 0;
        switch (shift_type){
            case ShiftType::TransformFromCellCentered: sign = -1.; break;
            case ShiftType::TransformToCellCentered: sign = 1.;
        }
        constexpr Complex I{0,1};
        for (int i=0; i<k.size(); i++ ){
            shift[i] = std::exp( I*sign*k[i]*0.5*dx[i_dim] );
        }
    }
    return shift_factor;
}

/* \brief For each box, in `spectralspace_ba`, which is owned by the local MPI
 * rank (as indicated by the argument `dm`), compute the values of the
 * corresponding finite-order modified k vector, along the
 * dimension specified by `i_dim`
 *
 * The finite-order modified k vector is the spectral-space representation
 * of a finite-order stencil in real space.
 *
 * \param n_order Order of accuracy of the stencil, in discretizing
 *                a spatial derivative
 * \param nodal Whether the stencil is to be applied to a nodal or
                staggered set of fields
 */
KVectorComponent
SpectralKSpace::getModifiedKComponent( const DistributionMapping& dm,
                                       const int i_dim,
                                       const int n_order,
                                       const bool nodal ) const
{
    // Initialize an empty ManagedVector in each box
    KVectorComponent modified_k_comp(spectralspace_ba, dm);

    // Compute real-space stencil coefficients
    Vector<Real> stencil_coef = getFonbergStencilCoefficients(n_order, nodal);

    // Loop over boxes and allocate the corresponding ManagedVector
    // for each box owned by the local MPI proc
    for ( MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi ){
        Real delta_x = dx[i_dim];
        const ManagedVector<Real>& k = k_vec[i_dim][mfi];
        ManagedVector<Real>& modified_k = modified_k_comp[mfi];

        // Allocate modified_k to the same size as k
        modified_k.resize( k.size() );

        // Fill the modified k vector
        for (int i=0; i<k.size(); i++ ){
            modified_k[i] = 0;
            for (int n=1; n<stencil_coef.size(); n++){
                if (nodal){
                    modified_k[i] += stencil_coef[n]* \
                        std::sin( k[i]*n*delta_x )/( n*delta_x );
                } else {
                    modified_k[i] += stencil_coef[n]* \
                        std::sin( k[i]*(n-0.5)*delta_x )/( (n-0.5)*delta_x );
                }
            }
        }
    }
    return modified_k_comp;
}

/* Returns an array of coefficients, corresponding to the weight
 * of each point in a finite-difference approximation (to order `n_order`)
 * of a derivative.
 *
 * `nodal` indicates whether this finite-difference approximation is
 * taken on a nodal grid or a staggered grid.
 */
Vector<Real>
getFonbergStencilCoefficients( const int n_order, const bool nodal )
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( n_order%2 == 0,
                                      "n_order should be even.");
    const int m = n_order/2;
    Vector<Real> coefs;
    coefs.resize( m+1 );

    // Note: there are closed-form formula for these coefficients,
    // but they result in an overflow when evaluated numerically.
    // One way to avoid the overflow is to calculate the coefficients
    // by recurrence.

    // Coefficients for nodal (a.k.a. centered) finite-difference
    if (nodal == true) {
        coefs[0] = -2.; // First coefficient
        for (int n=1; n<m+1; n++){ // Get the other coefficients by recurrence
            coefs[n] = - (m+1-n)*1./(m+n)*coefs[n-1];
        }
    }
    // Coefficients for staggered finite-difference
    else {
        Real prod = 1.;
        for (int k=1; k<m+1; k++){
            prod *= (m+k)*1./(4*k);
        }
        coefs[0] = 4*m*prod*prod; // First coefficient
        for (int n=1; n<m+1; n++){ // Get the other coefficients by recurrence
            coefs[n] = - ((2*n-3)*(m+1-n))*1./((2*n-1)*(m-1+n))*coefs[n-1];
        }
    }
    return coefs;
}
