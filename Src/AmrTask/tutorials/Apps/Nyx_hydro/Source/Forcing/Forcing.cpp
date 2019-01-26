#include <Nyx.H>
#include <Nyx_F.H>

#include "Forcing.H"

using namespace amrex;

unsigned long int mt_random();

int StochasticForcing::verbose      = 0;
int StochasticForcing::SpectralRank = 3;

extern "C"
{void fort_get_grav_const(Real* Gconst);}

//
//  Default constructor
//
StochasticForcing::StochasticForcing() 
{
    i1 = i2 = 0;
    j1 = j2 = 0;
    k1 = k2 = 0;
    NumModes = 0;
    NumNonZeroModes = 0;
    decay = 0; 
    seed = 27011974;

    SpectProfile = Parabolic;

    AmpltThresh = 1.0e-4;
    SolenoidalWeight = 1.0;
    DecayInitTime = 0.0;

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
	alpha[dim]         = 2;         
	BandWidth[dim]     = 1.0;
	IntgrVelocity[dim] = 0.0; 
	IntgrLength[dim]   = 0.0;
	WaveNumber[dim]    = 0.0;
	IntgrTime[dim]     = 0.0;
	AutoCorrlTime[dim] = 1.0;

	Amplitude[dim]     = NULL;
	InjectionEven[dim] = NULL;
	InjectionOdd[dim]  = NULL;
	SpectrumEven[dim]  = NULL;
	SpectrumOdd[dim]   = NULL;
    }
    mask = NULL;
}

//
//  Default destructor
//
StochasticForcing::~StochasticForcing() 
{
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
	delete [] Amplitude[dim];
	delete [] InjectionEven[dim];
	delete [] InjectionOdd[dim];
	delete [] SpectrumEven[dim];
	delete [] SpectrumOdd[dim];
    }
    delete [] mask;
}

/***********************************************************************
 *
 *  STOCHASTIC FORCING CLASS: evolve
 *
 *  written by: Wolfram Schmidt
 *  date:       May, 2005
 *  modified1:  Oct, 2014: updated to support Enzo 2.4 // P. Grete
 *  modified2:  May, 2017: ported to Nyx
 *
 *  PURPOSE: evolves the random forcing spectrum in the fashion of
 *           a multi-dimensional Ornstein-Uhlenbeck process
 *           
 *           Parameters:
 *           dt -- time step (small compared to AutoCorrlTime)
 *
 *  AUXILIARIES: inject, gauss_deviate, distribute, rms
 *
 ***********************************************************************/

void StochasticForcing::evolve(Real dt)
{
    if (ParallelDescriptor::IOProcessor()) {

	Real DriftCoeff[MAX_DIMENSION], DiffCoeff[MAX_DIMENSION];

	if (decay == 0) {

	    inject();
		    
	    /* Increment forcing spectrum (drift and random diffusion) 
	     * For general properties of Ornstein-Uhlenbeck process, see e.g.
	     * Turbulent Flows by Pope (2000) Appendix J with 
	     * drift and diffusion coefficients given eq (J.41)
	     */

	    for (int dim = 0; dim < SpectralRank; dim++) {
		DriftCoeff[dim] = exp(-dt/AutoCorrlTime[dim]);
		DiffCoeff[dim]  = sqrt(1 - DriftCoeff[dim]*DriftCoeff[dim]);
		for (int n = 0, m = 0; n < NumModes; n++)
		    if (mask[n]) {
			SpectrumEven[dim][m] = DriftCoeff[dim] * SpectrumEven[dim][m] +
			                       DiffCoeff [dim] * InjectionEven[dim][n];
			SpectrumOdd [dim][m] = DriftCoeff[dim] * SpectrumOdd [dim][m] + 
			                       DiffCoeff [dim] * InjectionOdd [dim][n];
			++m;
		    }
	    }

	} else {

	    /* increment forcing spectrum (drift only) */

	    for (int dim = 0; dim < SpectralRank; dim++) {
		DriftCoeff[dim] = exp(-dt/AutoCorrlTime[dim]);
		for (int m = 0; m < NumNonZeroModes; m++) {
		    SpectrumEven[dim][m] = DriftCoeff[dim] * SpectrumEven[dim][m];
		    SpectrumOdd [dim][m] = DriftCoeff[dim] * SpectrumOdd [dim][m];
		}
	    }
	}
    }

    /* communicate spectrum among processors */

    distribute();
}

//
// Compute new random injection
//
void StochasticForcing::inject(void)
{
    if (ParallelDescriptor::IOProcessor()) {

	int i, j, k, n, dim;
	Real a, b, contr, div;

	/* compute Gaussian deviates */

	for (dim = 0; dim < SpectralRank; dim++)
	    for (n = 0; n < NumModes; n++) {
		if (mask[n]) {
		    gauss_deviate(Amplitude[dim][n], &a, &b);
		} else {
		    a = 0.0; b = 0.0;
		}
		InjectionEven[dim][n] = a;
		InjectionOdd[dim][n]  = b;
	    }

	/* project modes 
	 * see eq (8) in Schmidt et al., A&A (2009)
	 * http://dx.doi.org/10.1051/0004-6361:200809967 */

	for (i = 0; i < i2; i++) { // wave vectors in positive x-direction
	    InjectionEven[0][i] = (1.0 - SolenoidalWeight) * InjectionEven[0][i];
	    InjectionOdd[0][i]  = (1.0 - SolenoidalWeight) * InjectionOdd[0][i];
	}

	if (SpectralRank > 1) {
	    
	    for (n = 0; n < i2; n++) { // wave vectors in positive x-direction
		InjectionEven[1][n] = SolenoidalWeight * InjectionEven[1][n];
		InjectionOdd [1][n] = SolenoidalWeight * InjectionOdd [1][n];
	    }
	    
	    n = i2;
	    for (j = 1; j <= j2; j++) { // wave vectors in xy-plane
		for (i = i1; i <= i2; i++) {
		    contr = (1.0 - 2.0 * SolenoidalWeight) * 
			(i*InjectionEven[0][n] + 
			 j*InjectionEven[1][n]) / Real(i*i + j*j);
		    InjectionEven[0][n] = SolenoidalWeight * InjectionEven[0][n] + i*contr;
		    InjectionEven[1][n] = SolenoidalWeight * InjectionEven[1][n] + j*contr;
		    contr = (1.0 - 2.0 * SolenoidalWeight) * 
			(i*InjectionOdd[0][n] + 
			 j*InjectionOdd[1][n]) / Real(i*i + j*j);
		    InjectionOdd[0][n] = SolenoidalWeight * InjectionOdd[0][n] + i*contr;
		    InjectionOdd[1][n] = SolenoidalWeight * InjectionOdd[1][n] + j*contr;
		    ++n;
		}
	    }
	    
	    if (SpectralRank > 2) {
		
		for (n = 0; n < i2 + j2*(i2-i1+1); n++) { // wave vectors in xy-plane
		    InjectionEven[2][n] = SolenoidalWeight * InjectionEven[2][n];
		    InjectionOdd[2][n]  = SolenoidalWeight * InjectionOdd [2][n];
		}
		
		for (k = 1; k <= k2; k++) { // wave vectors not aligned to xy-plane
		    for (j = j1; j <= j2; j++) {
			for (i = i1; i <= i2; i++) {
			    contr = (1.0 - 2.0 * SolenoidalWeight) * 
				(i*InjectionEven[0][n] + 
				 j*InjectionEven[1][n] + 
				 k*InjectionEven[2][n] ) / Real(i*i + j*j + k*k);
			    InjectionEven[0][n] = SolenoidalWeight * InjectionEven[0][n] + i*contr;
			    InjectionEven[1][n] = SolenoidalWeight * InjectionEven[1][n] + j*contr;
			    InjectionEven[2][n] = SolenoidalWeight * InjectionEven[2][n] + k*contr;
			    contr = (1.0 - 2.0 * SolenoidalWeight) * 
				(i*InjectionOdd[0][n] + 
				 j*InjectionOdd[1][n] +
				 k*InjectionOdd[2][n]) / Real(i*i + j*j + k*k);
			    InjectionOdd[0][n] = SolenoidalWeight * InjectionOdd[0][n] + i*contr;
			    InjectionOdd[1][n] = SolenoidalWeight * InjectionOdd[1][n] + j*contr;
			    InjectionOdd[2][n] = SolenoidalWeight * InjectionOdd[2][n] + k*contr;
			    ++n;

			}
		    }
		}
	    }
	}

    }
}

//
// Generate couple of normally distributed random deviates (Box-Muller-Algorithm)
//
void StochasticForcing::gauss_deviate(Real amplt, Real *x, Real *y)
{
	Real v_sqr, v1, v2;
	Real norm;

	do {
	    v1 = 2.0* (Real)(mt_random()%2147483563)/(2147483563.0) - 1.0;
	    v2 = 2.0* (Real)(mt_random()%2147483563)/(2147483563.0) - 1.0;
	    v_sqr = v1*v1+v2*v2;
	} while (v_sqr >= 1.0 || v_sqr == 0.0);
	
	norm = amplt * sqrt(-2.0*log(v_sqr)/v_sqr);

	*x = norm * v1; *y = norm * v2;
}

//
// Distribute the spectrum
//
void StochasticForcing::distribute(void)
{
    /* communicate spectrum among processors */

    for (int dim = 0; dim < SpectralRank; dim++) {
	ParallelDescriptor::Bcast(SpectrumEven[dim], NumNonZeroModes, ParallelDescriptor::IOProcessorNumber());
	ParallelDescriptor::Bcast(SpectrumOdd[dim],  NumNonZeroModes, ParallelDescriptor::IOProcessorNumber());
    }

    /* copy sepctrum to forcing_spect_module */

    for (int dim = 0; dim < SpectralRank; dim++)
        fort_set_modes(SpectrumEven[dim], SpectrumOdd[dim], &NumNonZeroModes, &dim);
}

//
// Compute RMS magnitude
// 
Real StochasticForcing::rms(void)
{
    int m;
    Real sum_even = 0.0, sum_odd = 0.0;

    for (int dim = 0; dim < SpectralRank; dim++) {
        for (m = 0; m < NumNonZeroModes; m++)
            sum_even += SpectrumEven[dim][m] * SpectrumEven[dim][m];
        for (m = 0; m < NumNonZeroModes; m++)
            sum_odd  += SpectrumOdd[dim][m]  * SpectrumOdd[dim][m];
    }

    return sqrt(sum_even + sum_odd);
}

