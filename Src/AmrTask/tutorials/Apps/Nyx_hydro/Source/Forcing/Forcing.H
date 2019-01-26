#ifndef _Forcing_H_
#define _Forcing_H_

#include <fstream>

#include <AMReX_AmrLevel.H>

#define MAX_DIMENSION 3

/* Type of spectral profile */
typedef int spect_profile_type;
const spect_profile_type
  None       = 0,
  Plane      = 1,
  Band       = 2,
  Parabolic  = 3;

/**************************************************************************
 *
 *  STOCHASTIC FORCING CLASS
 *
 *  written by: Wolfram Schmidt
 *  date:       May, 2005
 *  modified1:  May, 2017: ported to Nyx
 *
 *  PURPOSE: composes and evolves a stochastic driving force in Fourier space
 *
 **************************************************************************/

class StochasticForcing {

  private:

    static int verbose;
    static int SpectralRank;            // rank of the spectrum
    int NumModes;                       // number of Fourier modes in the spectral domain
    int NumNonZeroModes;                // number of non-zero Fourier modes
    int i1, i2, j1, j2, k1, k2;         // index boundaries specifying the range of wave numbers
    int decay;                          // if set to non-zero value, the force field decays
    int seed;                           // random seed
    amrex::Real AmpltThresh;            // threshold value for the normalised amplitude of modes counted as non-zero
//
//  Integral scale parameters (these are generally NOT the expectation values but characteristic quantities!)
//
    int alpha[MAX_DIMENSION];                 // ratio of domain size to integral scale
    amrex::Real BandWidth[MAX_DIMENSION];     // bandwidth of the forcing spectrum in units of 2*alpha
    amrex::Real IntgrVelocity[MAX_DIMENSION]; // integral velocity V in the statistically stationary regime
    amrex::Real IntgrTime[MAX_DIMENSION];     // integral (large-eddy turn-over) time T=L/V of the flow evolution
    amrex::Real IntgrLength[MAX_DIMENSION];   // integral length L of the force field
    amrex::Real WaveNumber[MAX_DIMENSION];    // characteristic wave number = 2 pi/L
//  
//  Stochastic process parameters (in Fourier space)
// 
    spect_profile_type SpectProfile;          // profile of forcing spectrum
    amrex::Real AutoCorrlTime[MAX_DIMENSION]; // autocorrelation time of the stochastic process
    amrex::Real SolenoidalWeight;             // determines weight of solenoidal relative to dilatational components
    amrex::Real DecayInitTime;                // time at which the decay of the force field is initiated
//
//  Spectral data
//
    int *mask;                                 // flags identifying non-zero modes
    amrex::Real *Amplitude[MAX_DIMENSION];     // amplitudes of forcing modes
    amrex::Real *InjectionEven[MAX_DIMENSION]; // random increments (cos modes)
    amrex::Real *InjectionOdd[MAX_DIMENSION];  // random increments (sin modes)
    amrex::Real *SpectrumEven[MAX_DIMENSION];  // forcing cos modes
    amrex::Real *SpectrumOdd[MAX_DIMENSION];   // forcing sin modes

  public:

//
// Constructor (nullifies data)
//
    StochasticForcing();
//
// Destructor
//
    ~StochasticForcing();
//
// Set forcing amplitudes and initialise spectrum
//
    void init(int rank, const amrex::Real* prob_lo, const amrex::Real* prob_hi);
//
// Evolve the force field over a time step which is small compared to the autocorrelation time
//
    void evolve(amrex::Real dt);
//
// Distribute spectrum
//
    void distribute(void);
//
// Calculate the instantaneous RMS force
//
    amrex::Real rms(void);
//
// Write/read forcing spectrum to/from output file
//
    int ReadSpectrum(char *fname);
    int WriteSpectrum(char *fname);
// 
// Initiate the decay of the force field
//
    void set_decay(void);
//
// Set the weighting parameter
//
    void set_SolenoidalWeight(int my_soln_weight);
//
// Get boundary indices and range of spectral domain
//
    int get_SpectralRank(void);
    int get_LeftBoundary(int dim);
    int get_RightBoundary(int dim);
    int get_Range(int dim);
//
// Get number of forcing modes
//
    int get_NumModes(void);
    int get_NumNonZeroModes(void);
//
// Get spectral profile
//
    spect_profile_type get_SpectProfile(void);
//
// Get the characteristic wave number
//
    amrex::Real get_WaveNumber(int dim);
//
// Get the integral length scale
//
    amrex::Real get_IntgrLength(int dim);
//
// Get the integral time scale
//
    amrex::Real get_IntgrTime(int dim);
//
// Copy flags identifying non-zero modes
//
    void copy_mask(int* target);
//
// Copy non-zero modes of the forcing spectrum
//
    void copy_SpectrumEven(int dim, amrex::Real* target);
    void copy_SpectrumOdd(int dim, amrex::Real* target);
//
// Copy complete forcing spectrum into a single array for subsequent processing
//
    void copy_ExpandedSpectrum(int dim, amrex::Real* target);
//
// Read spectrum from file handle
//
    void read_Spectrum(std::ifstream& input);
//
// Write spectrum to file handle
//
    void write_Spectrum(std::ofstream& output);

  private:

//
// Compute random injection
//
    void inject(void);
//
// Random number generator
//
    void gauss_deviate(amrex::Real amplt, amrex::Real *x, amrex::Real *y);
//
// Read forcing parameters
//
    void read_params(void);
};


inline void StochasticForcing::set_decay(void) 
{
    decay = 1;
}

inline void StochasticForcing::set_SolenoidalWeight(int my_soln_weight)
{
    if (my_soln_weight >= 1.0) SolenoidalWeight = 1.0;
    if (my_soln_weight  < 1.0) SolenoidalWeight = my_soln_weight;
    if (my_soln_weight <= 0.0) SolenoidalWeight = 0.0;
}

inline int StochasticForcing::get_SpectralRank(void)
{
    return SpectralRank;
}

inline int StochasticForcing::get_LeftBoundary(int dim)
{
    switch (dim) {
	case 1:
	    return i1;
	case 2:
	    return j1;
	case 3:
	    return k1;
	default:
	    return 0;
    }
}

inline int StochasticForcing::get_RightBoundary(int dim)
{
    switch (dim) {
	case 1:
	    return i2;
	case 2:
	    return j2;
	case 3:
	    return k2;
	default:
	    return 0;
    }
}

inline int StochasticForcing::get_NumModes(void)
{
    return NumModes;
}

inline int StochasticForcing::get_NumNonZeroModes(void)
{
    return NumNonZeroModes;
}

inline spect_profile_type StochasticForcing::get_SpectProfile(void)
{
    return SpectProfile;
}

inline amrex::Real StochasticForcing::get_WaveNumber(int dim)
{
    return WaveNumber[dim];
}

inline amrex::Real StochasticForcing::get_IntgrLength(int dim)
{
    return IntgrLength[dim];
}

inline amrex::Real StochasticForcing::get_IntgrTime(int dim)
{
    return IntgrTime[dim];
}

inline void StochasticForcing::copy_mask(int* target)
{
    for (int n = 0; n < NumModes; n++) 
	target[n] = mask[n];
}

inline void StochasticForcing::copy_SpectrumOdd(int dim, amrex::Real* target)
{
    for (int m = 0; m < NumNonZeroModes; m++) 
	target[m] = SpectrumOdd[dim][m];
}

inline void StochasticForcing::copy_SpectrumEven(int dim, amrex::Real* target)
{
    for (int m = 0; m < NumNonZeroModes; m++) 
	target[m] = SpectrumEven[dim][m];
}

inline void StochasticForcing::copy_ExpandedSpectrum(int dim, amrex::Real* target)
{
    for (int n = 0, m = NumNonZeroModes-1; n < NumModes; n++) 
	target[n] = mask[NumModes-1-n] ? SpectrumOdd[dim][m--] : 0.0;

    target[NumModes] = 0.0;

    for (int n = NumModes+1, m = 0; n <= 2*NumModes; n++) 
	target[n] = mask[n-NumModes-1] ? SpectrumEven[dim][m++] : 0.0;
}

inline void StochasticForcing::read_Spectrum(std::ifstream& input)
{
    for (int dim = 0; dim < SpectralRank; dim++) {
	for (int m = 0; m < NumNonZeroModes; m++) 
            input >> SpectrumEven[dim][m];
	for (int m = 0; m < NumNonZeroModes; m++) 
            input >> SpectrumOdd[dim][m];
    }
}

inline void StochasticForcing::write_Spectrum(std::ofstream& output)
{
    for (int dim = 0; dim < SpectralRank; dim++) {
	for (int m = 0; m < NumNonZeroModes; m++) 
            output << SpectrumEven[dim][m] << '\n';
	for (int m = 0; m < NumNonZeroModes; m++) 
            output << SpectrumOdd[dim][m] << '\n';
    }
}
#endif

