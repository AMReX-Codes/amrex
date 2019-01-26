#include <iostream>
#include <iomanip>
#include <float.h>

#include <AMReX_ParmParse.H>
#include <Nyx.H>
#include <Nyx_F.H>

#include "Forcing.H"

using namespace amrex;

void mt_init(unsigned int seed);

void mt_read(std::ifstream& input);

/***********************************************************************
/
/  STOCHASTIC FORCING CLASS METHOD: init
/
/  written by: Wolfram Schmidt
/  date:       May, 2005
/  modified1:  Oct, 2014: updated to support Enzo 2.4 // P. Grete
/  modified2:  May, 2017: ported to Nyx
/
/  PURPOSE: initializes StochasticForcing object with given parameters; 
/           de facto, this is a parametrized constructor;
/           since the parameters are not known prior to the declaration
/           of the object in Nyx, however, it is defined as a
/           regular method
/
/           Parameters:
/           rank -- dimension of Fouries space
/           prob_lo -- lower domain boundaries
/           prob_hi -- higher domain boundaries
/
/  AUXILIARIES: read_param()
/           reads from inputs:
/           forcing.seed -- random seed
/           forcing.profile -- shape of forcing power spectrum
/                              (1: plane waves, 2: band, 3: parabolic) 
/           forcing.alpha -- ratio of domain length to integral length
/                            for each dimension (L = X/alpha)
/           forcing.band_width -- determines band width of the forcing 
/                                 spectrum relative to alpha (maximal 
/                                 value = 1)
/           forcing.intgr_vel -- characteristic velocity scale for each 
/                                dimension (charcteristic force per unit 
/                                mass F = V*V/L)
/           forcing.auto_corrl -- determines autocorrelation time of the
/                                 stochastic force in units of the integral
/                                 time scale T = L/V
/           forcing.soln_weight -- determines weight of solenoidal relative
/                                  to dilatational modes (1 = purely 
/                                  solenoidal, 0 = purely dilatational)
/           read_Spectrum()
/
************************************************************************/

void StochasticForcing::init(int rank, const Real* prob_lo, const Real* prob_hi)
{
    BL_PROFILE("StochasticForcing::init()");
    int i, j, k, m, n;

    /* set parameters */

    Real domain_length[SpectralRank];
    for (int i = 0; i < SpectralRank; i++) {
	domain_length[i] = prob_hi[i] - prob_lo[i];
	IntgrLength[i] = domain_length[i]; // store intermediate result
    }

    SpectralRank = rank;

    read_params();

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nStochastic forcing field:"
                  << "\n   spectral profile = " << SpectProfile
                  << "\n   weight of solenoidal component = " << SolenoidalWeight  
                  << "\n   alpha = ("
                  << alpha[0] << " " << alpha[1] << " " << alpha[2] << ")"
                  << "\n   band width = ("
                  << BandWidth[0] << " " << BandWidth[1] << " " << BandWidth[2]<< ")"
                  << "\n   integral velocity = ("
                  << IntgrVelocity[0] << " " << IntgrVelocity[1] << " " << IntgrVelocity[2] << ")"
                  << "\n   integral length = ("
                  << IntgrLength[0] << " " << IntgrLength[1] << " " << IntgrLength[2] << ")"
                  << "\n   integral time = ("
                  << IntgrTime[0] << " " << IntgrTime[1] << " " << IntgrTime[2] << ")"
                  << "\n   autocorrelation time = ("
                  << AutoCorrlTime[0] << " " << AutoCorrlTime[1] << " " << AutoCorrlTime[2] << ")" << std::endl;
    }

    /* determine boundary indices of spectral domain */

    if (SpectralRank > 0) {
	i1 = -2*alpha[0]+1; if (i1 > 0) i1 = 0;
	i2 =  2*alpha[0]-1; if (i2 < 0) i2 = 0;
    }
    if (SpectralRank > 1) {
	j1 = -2*alpha[1]+1; if (j1 > 0) j1 = 0;
	j2 =  2*alpha[1]-1; if (j2 < 0) j2 = 0;
    }
    if (SpectralRank > 2) {
	k1 = -2*alpha[2]+1; if (k1 > 0) k1 = 0;
	k2 =  2*alpha[2]-1; if (k2 < 0) k2 = 0;
    }

    /* determine number of linearly independent modes */
    
    NumModes = i2 + j2*(i2-i1+1) + k2*(j2-j1+1)*(i2-i1+1);
    
    mask = new int[NumModes];

    if (ParallelDescriptor::IOProcessor()) {
	if (verbose) 
           std::cout << "Total number of stochastic forcing modes = " << NumModes << std::endl;

	/* determine forcing field amplitudes */

	for (int dim = 0; dim < SpectralRank; dim++) {
	    Amplitude[dim] = new Real[NumModes];
	}

	if (SpectProfile == Plane) {

            if (verbose > 1)
	       std::cout << "\ni1 = " << i1 << ", i2 = " << i2 
	                 << "\nj1 = " << j1 << ", j2 = " << j2 
	                 << "\nk1 = " << k1 << ", k2 = " << k2 << std::endl;

	    for (n = 0; n < NumModes; n++)
		Amplitude[0][n] = 0.0;

	    for (i = 1; i <= i2; i++) {
            	if (i == alpha[0])
	           Amplitude[0][i-1] = 1.0;
		if (verbose > 1) 
 	           std::cout <<   "i = " << i 
	                     << ", n = " << i-1 << ", Amplitude = " << Amplitude[0][i-1] << std::endl;
	    }

	    if (SpectralRank > 1) {

		n = i2;
		for (j = 1; j <= j2; j++)
		    for (i = i1; i <= i2; i++) {
	                if (i == 0 && j == alpha[1])
                           Amplitude[0][n] = 1.0; 
	                if (verbose > 1)
	                   std::cout <<   "i = " << i <<  ", j = " << j 
	                             << ", n = " << n << ", Amplitude = " << Amplitude[0][n] << std::endl;
			++n;
	            }

		if (SpectralRank > 2) {

	            for (k = 1; k <= k2; k++)
	                for (j = j1; j <= j2; j++)
	                    for (i = i1; i <= i2; i++) {
	                        if (i == 0 && j == 0 && k == alpha[2])
                                   Amplitude[0][n] = 1.0;
	                        if (verbose > 1)
	                           std::cout <<   "i = " << i <<  ", j = " << j << ", k = " << k 
	                                     << ", n = " << n << ", Amplitude = " << Amplitude[0][n] << std::endl;	                       
	                        ++n;
	                    }
	        }
	    }

	} else {
	    
	    Real x1, x2;
	    Real a, a1 = 1.0, a2 = 1.0;
 
	    a1 = 1.0-BandWidth[0];
	    a2 = 1.0+BandWidth[0];
	    if (verbose > 1) 
	        std::cout << "i1 = " << i1 << ", i2 = " << i2 
	                  << "\na1 = " << a1 << ", a2 = " << a2 << std::endl;

	    /* compute amplitude factors for wave numbers within the interval [a1, a2] */

	    for (i = 1; i <= i2; i++) {
		a = 1.0; a = Real(i) / Real(alpha[0]);
		x1 = (a - a1); 
		x2 = (a2 - a);
		if (SpectProfile ==  Parabolic) {
		    Amplitude[0][i-1] = x1*x2;
		    if (Amplitude[0][i-1] <= 0.0) Amplitude[0][i-1] = 0.0;
		    Amplitude[0][i-1] *= Amplitude[0][i-1];
		} else if (SpectProfile == Band) {
		    if ((x1 >= 0.0) && (x2 >= 0.0)) Amplitude[0][i-1] = 1.0;
		}

 		if (verbose > 1) 
		    std::cout << "i = " << i << ", a = " << a 
		              << ", x1 = " << x1 << ", x2 = " << x2
		              << ", n = " << i-1 << ", Amplitude = " << Amplitude[0][i-1] << std::endl;
	    }

	    if (SpectralRank > 1) {
		
		Real b, b1 = 1.0, b2 = 1.0;
		Real f1 = 0.0, f2 = 0.0;
		Real g1 = 1.0, g2 = 1.0;

		b1 = (a1 > 0.0) ? (1.0-BandWidth[1]) : 0.0; 
		b2 = (a1 > 0.0) ? (1.0+BandWidth[1]) : 2.0;
		if (b1 > 0.0) {
		    f1 = a1;
		    g1 = a1/b1; g1 *= g1;
		}
		if (b2 > 0.0) {
		    f2 = a2;
		    g2 = a2/b2; g2 *= g2;
		}
        
		if (verbose > 1) 
		   std::cout <<   "j1 = " << j1 << ", j2 = " << j2
		             << "\nb1 = " << b1 << ", b2 = " << b2
		             << "\nf1 = " << f1 << ", f2 = " << f2
		             << "\ng1 = " << g1 << ", g2 = " << g2 << std::endl;

		/* compute amplitude factors for wave numbers bounded by the
		   ellipses with semi axes a1, b1 and a2, b2, respectively */
		
		n = i2;
		for (j = 1; j <= j2; j++) {
		    b = 0.0; if (alpha[1] > 0.0) b = Real(j) / Real(alpha[1]); b *= b;
		    for (i = i1; i <= i2; i++) {
			a = 0.0; if (alpha[0] > 0.0) a = Real(i) / Real(alpha[0]); a *= a;
			x1 = sqrt(a + g1*b) - f1;
			x2 = f2 - sqrt(a + g2*b);
			if (SpectProfile ==  Parabolic) {
			    Amplitude[0][n] = x1*x2;
			    if (Amplitude[0][n] <= 0.0) Amplitude[0][n] = 0.0;
			    Amplitude[0][n] *= Amplitude[0][n];
			} else if (SpectProfile == Band) {
			    if ((x1 >= 0.0) && (x2 >= 0.0)) Amplitude[0][n] = 1.0;
			}
            
			if (verbose > 1) 
			    std::cout <<   "i = " << i << ", a = " << a 
			              << ", j = " << j << ", b = " << b
			              << ", x1 = " << x1 << ", x2 = " << x2
			              << ", n = " << n << ", Amplitude = " << Amplitude[0][n] << std::endl;
           
			++n;
		    }
		}

		if (SpectralRank > 2) {

		    Real c, c1 = 1.0, c2 = 1.0;
		    Real h1 = 1.0, h2 = 1.0;
		    
		    c1 = (b1 > 0.0) ? (1.0-BandWidth[2]) : 0.0;
		    c2 = (b1 > 0.0) ? (1.0+BandWidth[2]) : 0.0;
		    if (c1 > 0.0) {
			h1 = a1/c1; h1 *= h1;
		    }
		    if (c2 > 0.0) {
			h2 = a2/c2; h2 *= h2;
		    }
            
		    if (verbose > 1)
		        std::cout <<   "c1 = " << c1 << ", c2 = " << c2
		                  << "\nh1 = " << h1 << ", h2 = " << h2 << std::endl;
		    
		    /* compute amplitude factors for wave numbers bounded by the
		       ellipsoids with semi axes a1, b1, c1 and a2, b2, c2, respectively */

		    for (k = 1; k <= k2; k++) {
			c = 0.0; if (alpha[2] > 0.0) c = Real(k) / Real(alpha[2]); c *= c;
			for (j = j1; j <= j2; j++) {
			    b = 0.0; if (alpha[1] > 0.0) b = Real(j) / Real(alpha[1]); b *= b;
			    for (i = i1; i <= i2; i++) {
				a = 0.0; if (alpha[0] > 0.0) a = Real(i) / Real(alpha[0]); a*= a;
				x1 = sqrt(a + g1*b + h1*c) - f1;
				x2 = f2 - sqrt(a + g2*b + h2*c);
				if (SpectProfile ==  Parabolic) {
				    Amplitude[0][n] = x1*x2;
				    if (Amplitude[0][n] < 0.0) Amplitude[0][n] = 0.0;
				    Amplitude[0][n] *= Amplitude[0][n];
				} else if (SpectProfile == Band) {
				    if ((x1 >= 0.0) && (x2 >= 0.0)) Amplitude[0][n] = 1.0;
				}
                
				if (verbose > 1)
				    std::cout <<   "i = " << i << ", a = " << a 
				              << ", j = " << j << ", b = " << b
				              << ", k = " << k << ", c = " << c
				              << ", x1 = " << x1 << ", x2 = " << x2
				              << ", n = " << n << ", Amplitude = " << Amplitude[0][n] << std::endl;
                
				++n;
			    }
			}
		    }
		}
	    }
	}

	/* normalise amplitude factors and 
           set flags for modes with amplitude larger than the threshold */

	Real norm = 0.0;
	for (n = 0; n < NumModes; n++) norm += 2*Amplitude[0][n]*Amplitude[0][n];
	norm = 1/sqrt(norm);

	NumNonZeroModes = 0;
	for (n = 0; n < NumModes; n++) {
	    Amplitude[0][n] *= norm;
	    if (Amplitude[0][n] > AmpltThresh) {
		mask[n] = 1;
		++NumNonZeroModes;
	    } else {
		mask[n] = 0;
	    }
	    if (SpectralRank > 1) Amplitude[1][n] = Amplitude[0][n];
	    if (SpectralRank > 2) Amplitude[2][n] = Amplitude[1][n];
	}
	if (verbose)	    	    
	    std::cout << "Number of non-zero stochastic forcing modes = " << NumNonZeroModes << std::endl;
	
	if (NumNonZeroModes == 0)        
            amrex::Error("Nyx::vanishing number of forcing modes");
	
	for (int dim = 0; dim < SpectralRank; dim++) {

	    InjectionEven[dim] = new Real[NumModes];
	    InjectionOdd[dim]  = new Real[NumModes];

            // normalization such that RMS acceleration is V/T for isotropic characteristic scales 
	    norm = (IntgrVelocity[dim]/IntgrTime[dim]) /
	           sqrt(1.0 - 2.0*SolenoidalWeight + 3.0*SolenoidalWeight*SolenoidalWeight);

	    if (verbose) 
	        std::cout << "Normalization factor[" << dim << "] = " << norm << std::endl;

	    for (n = 0; n < NumModes; n++) {
		Amplitude[dim][n] *= norm;
		InjectionEven[dim][n] = 0.0;
		InjectionOdd [dim][n] = 0.0;
	    }
        }

	/* initialise new sequence of random numbers */

	mt_init((unsigned int) seed);

	/* compute initial set of random deviates */

	inject();
    }

    /* communicate mask among processors */

    ParallelDescriptor::Bcast(&NumNonZeroModes, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(mask, NumModes, ParallelDescriptor::IOProcessorNumber());

    /* allocate memory for the forcing spectrum and set inital random deviates */

    for (int dim = 0; dim < SpectralRank; dim++) {

	SpectrumEven[dim] = new Real[NumNonZeroModes];
	SpectrumOdd[dim]  = new Real[NumNonZeroModes];
	
	if (ParallelDescriptor::IOProcessor()) {  
	    if (verbose > 1)
	       std::cout << "Master #" << ParallelDescriptor::MyProc() << " spectrum:\n";
	    for (n = 0, m = 0; n < NumModes; n++) {
		if (mask[n]) { // set only non-zero modes
		    SpectrumEven[dim][m] = InjectionEven[dim][n];
		    SpectrumOdd [dim][m] = InjectionOdd [dim][n];
	            if (verbose > 1)
		       std::cout << "   " << dim << ", mode " << n << " = (" 
	                         << SpectrumEven[dim][m] << "," << SpectrumOdd [dim][m] << ")\n";
		    ++m;
		}
	    }
            std::cout << std::endl;
	}
    }

    /* communicate spectrum among processors */

    for (int dim = 0; dim < SpectralRank; dim++) {
	ParallelDescriptor::Bcast(SpectrumEven[dim], NumNonZeroModes, ParallelDescriptor::IOProcessorNumber());
	ParallelDescriptor::Bcast(SpectrumOdd[dim],  NumNonZeroModes, ParallelDescriptor::IOProcessorNumber());
    }

    /* set wave vectors and copy sepctrum to forcing_spect_module */

    fort_alloc_spect(&NumNonZeroModes);

    int kvect[SpectralRank];
 
    m = 0;
    for (i = 1; i <= i2; i++)
        if (mask[i-1]) {
	   kvect[0] = i; kvect[1] = 0; kvect[2] = 0;
           fort_set_wavevector(kvect, &m);
           ++m;
    }
    if (SpectralRank > 1) {
        
       n = i2;
       for (j = 1; j <= j2; j++) 
           for (i = i1; i <= i2; i++)
               if (mask[n++]) {
                  kvect[0] = i; kvect[1] = j; kvect[2] = 0;
                  fort_set_wavevector(kvect, &m);
                  ++m;
               }

       if (SpectralRank > 2) {
       
          for (k = 1; k <= k2; k++)
              for (j = j1; j <= j2; j++)
                  for (i = i1; i <= i2; i++)
                      if (mask[n++]) {
                         kvect[0] = i; kvect[1] = j; kvect[2] = k;
                         fort_set_wavevector(kvect, &m);
                         ++m;
                      }
       }
    }

    for (int dim = 0; dim < SpectralRank; dim++)
        fort_set_modes(SpectrumEven[dim], SpectrumOdd[dim], &NumNonZeroModes, &dim);

    return;
}

//
// Read forcing parameters from inputs
//
void StochasticForcing::read_params()
{
    BL_PROFILE("StochasticForcing::read_params()");
    static bool done = false;

    if (!done)
    {
        const Real strt = ParallelDescriptor::second();

        ParmParse pp("forcing");
 
        pp.query("v", verbose);

        pp.query("seed", seed);

	int profile = 3;
	pp.query("profile", profile);
        switch(profile) {
           case 0: SpectProfile = None; break;
           case 1: SpectProfile = Plane; break;
           case 2: SpectProfile = Band; break;
           case 3: SpectProfile = Parabolic; break;
           default: if (ParallelDescriptor::IOProcessorNumber())
                       std::cerr << "Spectral profile must be 0, 1, 2, or 3." << std::endl;
                    amrex::Abort("StochasticForcing::invalid parameter");
        }

        pp.query("soln_weight", SolenoidalWeight);
        if (SolenoidalWeight > 1.0 || SolenoidalWeight < 0.0) {
           if (ParallelDescriptor::IOProcessorNumber())
              std::cerr << "Solenoidal weight must be between 0 and 1 " << std::endl;
           amrex::Abort("StochasticForcing::invalid parameter");
        }
        if (SpectralRank == 1) SolenoidalWeight = 0.0;

        Vector<int> input_int(SpectralRank);

        pp.getarr("alpha", input_int, 0, SpectralRank);
        for (int i = 0; i < SpectralRank; i++) {
            alpha[i] = input_int[i];
            if (alpha[i] > 0) 
               IntgrLength[i] /= alpha[i];
            else {
               if (ParallelDescriptor::IOProcessorNumber())
                  std::cerr << "alpha must be > 0." << std::endl;
               amrex::Abort("StochasticForcing::invalid parameter");
            }
        }

        Vector<Real> input_real(SpectralRank);

        pp.getarr("band_width", input_real, 0, SpectralRank);
        for (int i = 0; i < SpectralRank; i++)
            BandWidth[i] = (input_real[i] > 0.0) ? input_real[i] : 0.0;

        pp.getarr("intgr_vel", input_real, 0, SpectralRank);
        for (int i = 0; i < SpectralRank; i++) {
            IntgrVelocity[i] = (input_real[i] > 0.0) ? input_real[i] : 0.0;
            IntgrTime[i] = (IntgrVelocity[i] > 0.0) ? (IntgrLength[i] / IntgrVelocity[i]) : FLT_MAX;
        }

        pp.getarr("auto_corrl", input_real, 0, SpectralRank);
        for (int i = 0; i < SpectralRank; i++)
            AutoCorrlTime[i] = input_real[i]*IntgrTime[i];

        done = true;
    }
}

//
// Read the spectrum from checkpoint file
//
void
Nyx::forcing_post_restart (const std::string& restart_file)
{
    BL_PROFILE("StochasticForcing::forcing_post_restart()");
    if (level > 0)
        return;

    if (ParallelDescriptor::IOProcessor())
    {
        std::string FileName = restart_file + "/forcing";
        std::ifstream File;
        File.open(FileName.c_str(), std::ios::in);
        if (!File.good())
            amrex::FileOpenFailed(FileName);
        forcing->read_Spectrum(File);
        File.close();

        FileName = restart_file + "/mt";
        File.open(FileName.c_str(), std::ios::in);
        if (!File.good())
            amrex::FileOpenFailed(FileName);
        mt_read(File);
    }

    forcing->distribute();
}

