#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <new>
#include <unistd.h>

#include <AMReX_Box.H>
#include <AMReX_DataServices.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include "rfftw_mpi.h" // Machine specific option (i.e. Jagwar)
//#include "rfftw_mpi.h"

#include "AmrDeriveSpectrum.H"

//
// Usage note:
//
// FFTW uses a slab decomposition and decides how to arrange the data
// Suppose there are N cells in the k direction and P processors
// The domain is decomposed into P slabs with N/P cells
// It is therefore sensible to choose P such that P divides N
// P cannot exceed N, and my experience suggests P=N/4 is a good choice
//
// Takes an input file to define the following variables:
//
// verbose = 0 or 1
// infile = [list of input plotfiles to use]
// vars = [list of variables to load]
// div_free = 0 or 1
// density_weighting = 0 or 1
// use_cutoff_density = 0 or 1
// cutoff_density = density below which to zero velocities.
// transpose_dp = 0 or 1
// density
// do_filter = 0 or 1
// filterWN = list of filter wavenumbers
// finestLevel = finest level to use
//

int main (int argc, char* argv[])
{
    //
    // Initialize
    //
  
    amrex::Initialize(argc,argv);
  
    nProcs = ParallelDescriptor::NProcs();
    myProc = ParallelDescriptor::MyProc();
    IOProc = ParallelDescriptor::IOProcessorNumber();
  
    //
    // Read input data
    //
    ParmParse pp;

    if (ParallelDescriptor::IOProcessor()) 
      std::cout << "getting started" << std::endl;

    verbose=0;
    pp.query("verbose",verbose);
    if (ParallelDescriptor::IOProcessor()) 
      std::cout << "setting verbose = " << verbose << std::endl;

    int nPlotFiles(pp.countval("infile"));
    if(nPlotFiles <= 0) {
      std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
      std::cerr << "Exiting." << std::endl;
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    // Make an array of srings containing paths of input plot files
    Vector<std::string> plotFileNames(nPlotFiles);
    for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
      pp.get("infile", plotFileNames[iPlot], iPlot);
    }
    if (ParallelDescriptor::IOProcessor()){ 
      std::cout << "number of plotfiles  = " << nPlotFiles << std::endl;
      std::cout << "first infile = " << plotFileNames[0] << std::endl;
    }

    nVars=pp.countval("vars");
    if (nVars==0)
        amrex::Abort("Must specify vars to load");
    if (ParallelDescriptor::IOProcessor())
        std::cout << "nVars = " << nVars << std::endl;

    div_free=0;
    pp.query("div_free",div_free);
    if (div_free && nVars!=3) 
      amrex::Abort("Must specify three vars (assumed to be x y z velocities) if using div_free");

    transpose_dp=1;
    pp.query("transpose_dp",transpose_dp);

	use_cutoff_density=0;
	pp.query("use_cutoff_density",use_cutoff_density);

	cutoff_density=0.0;
	pp.query("cutoff_density",cutoff_density);

    density_weighting=0;
    pp.query("density_weighting",density_weighting);
    if (div_free && density_weighting) 
      amrex::Abort("Density weighting and div free at the same time doesn't make sense (yet)");
    // Has to be exactly one to get nVars right
    if (density_weighting) density_weighting=1; else density_weighting=0;

    whichVar.resize(nVars+density_weighting);
    if (ParallelDescriptor::IOProcessor())
	std::cout << "vars = ";

    for (int i=0; i<nVars; i++) {
	  pp.get("vars",whichVar[i],i);
	  if (ParallelDescriptor::IOProcessor())
	    std::cout << " " << whichVar[i];
    }

    if (density_weighting) {
	  pp.get("density",whichVar[nVars]);
	  if (ParallelDescriptor::IOProcessor())
	    std::cout << " " << whichVar[nVars];
    }

    if (ParallelDescriptor::IOProcessor())
	  std::cout << std::endl;

    Vector<int> destFills(nVars+density_weighting);
    for (int c=0; c<nVars+density_weighting; c++ ) destFills[c] = c;

    //size arrays for holding data
    sum  = (Real*) malloc(sizeof(Real)*nVars);
    sum2 = (Real*) malloc(sizeof(Real)*nVars);

    spectrum.resize(nVars+1);
    if (div_free) {
      spectrumS.resize(nVars);
      spectrumC.resize(nVars);
    }

    Qx.resize(nVars);
    Qy.resize(nVars);
    Qz.resize(nVars);

    local_data.resize(nVars);
    local_data_c.resize(nVars);

    if (ParallelDescriptor::IOProcessor())
	  std::cout << std::endl;

    pp.query("do_filter",do_filter);
    if (do_filter) {
      nFilters=pp.countval("filterWN");
      filterWN.resize(nFilters);
      for (int i=0; i<nFilters; i++)
		pp.get("filterWN",filterWN[i],i);
      if (ParallelDescriptor::IOProcessor())
		std::cout << "Filtering on wavenumber(s)..." << std::endl;
    }
    
    //
    // Read plot file info
    //
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {

      // initialize sum to zero
      for (int iVar=0; iVar<nVars; iVar++)
		sum[iVar]=sum2[iVar]=0.;

      infile=plotFileNames[iPlot];
      if (ParallelDescriptor::IOProcessor())
		std::cout << "working on " << plotFileNames[iPlot] << std::endl;

      /* from Spherical polars
	  DataServices dataServices(infile, fileType);

	  AmrData& amrData = dataServices.AmrDataRef();
	  */

      DataServices *dataServices = new DataServices(infile, fileType);
      if( ! dataServices->AmrDataOk())
      	DataServices::Dispatch(DataServices::ExitRequest, NULL);

      AmrData amrData(dataServices->AmrDataRef());

      Time      = amrData.Time();
      timeSteps = amrData.LevelSteps()[0];

      int finestLevel = amrData.FinestLevel();
      int finestLevelIn(-1);
      pp.query("finestLevel",finestLevelIn);
      if (finestLevelIn>=0 && finestLevelIn<finestLevel) {
		finestLevel=finestLevelIn;
      }

	  if (ParallelDescriptor::IOProcessor())
		std::cout << "Using finestLevel = " << finestLevel << std::endl;

	  Box probDomain(amrData.ProbDomain()[finestLevel]);

	  // Set AMReX and FourierTranform array sizes
	  // Note this defaults to a transposition
	  BLix = FTkx = probDomain.length(XDIR);
	  BLjx = FTjx = probDomain.length(YDIR);
	  BLkx = FTix = probDomain.length(ZDIR);

	  // Figure out the maximum length and wavenumber scaling factors for non-cubic domains
	  FTmx = FTix;
	  if (FTjx>FTmx) FTmx=FTjx;
	  if (FTkx>FTmx) FTmx=FTkx;

	  FTis = FTmx/FTix;
	  FTjs = FTmx/FTjx;
	  FTks = FTmx/FTkx;

	  // Half kx+1 - accounts for the fftw padding
	  FThkxpo = FTkx/2+1;

	  // Number of wavenumbers in the spectra
	  wavenumbers = FTmx/2;

	  // Size of correlation functions
	  Qix = BLix/2;
	  Qjx = BLjx/2;
	  Qkx = BLkx/2;

	  // Declare memory for spectra (plus one for counting hits)
	  for (int iVar=0; iVar<=nVars; iVar++) {
		spectrum[iVar]=(Real*)malloc(sizeof(Real)*wavenumbers);
		for (int wn=0; wn<wavenumbers; wn++)
		  spectrum[iVar][wn] = 0.0;
	  }

	  if (div_free) {
		for (int iVar=0; iVar<nVars; iVar++) {
		  spectrumS[iVar]=(Real*)malloc(sizeof(Real)*wavenumbers);
		  spectrumC[iVar]=(Real*)malloc(sizeof(Real)*wavenumbers);
		  for (int wn=0; wn<wavenumbers; wn++) {
			spectrumS[iVar][wn] = 0.0;
			spectrumC[iVar][wn] = 0.0;
		  }
		}
	  }

	  // Declare memory for correlation functions
	  // Qx, Qy and Qz are the correlations in the three directions
	  for (int iVar=0; iVar<nVars; iVar++) {
		Qx[iVar]=(Real*)malloc(sizeof(Real)*Qix);    for (int i=0; i<Qix; i++)    Qx[iVar][i] = 0.0;
		Qy[iVar]=(Real*)malloc(sizeof(Real)*Qjx);    for (int j=0; j<Qjx; j++)    Qy[iVar][j] = 0.0;
		Qz[iVar]=(Real*)malloc(sizeof(Real)*Qkx);    for (int k=0; k<Qkx; k++)    Qz[iVar][k] = 0.0;
	  }

	  // Other AMReX stuff
	  probLo=amrData.ProbLo();
	  probHi=amrData.ProbHi();

	  Lx = probHi[0]-probLo[0];
	  Ly = probHi[0]-probLo[0];
	  Lz = probHi[0]-probLo[0];

	  dx = Lx/(Real)BLix;
	  dy = Ly/(Real)BLjx;
	  dz = Lz/(Real)BLkx;

	  Real dxyz=dx*dy*dz;

	  //
	  // Plan ffts and make boxes, distribution etc.
	  //
	  BoxArray   domainBoxArray(nProcs);
	  Vector<int> pmap(nProcs+1);

	  plan_ffts(probDomain,domainBoxArray,pmap);

	  DistributionMapping domainDistMap(pmap);

	  //
	  // Load plot file into prescribed data structure
	  // 
	  int ngrow(0);
	  MultiFab mf;
	  mf.define(domainBoxArray, domainDistMap, nVars+density_weighting, ngrow, MFInfo().SetAlloc(true));

	  Real timer_start = ParallelDescriptor::second();

	  amrData.FillVar(mf, finestLevel, whichVar, destFills);

	  //
	  // Zero velocity components if density < density_lower_bound
	  //
	  if (use_cutoff_density) {

		for (MFIter ntmfi(mf); ntmfi.isValid(); ++ntmfi) {
		  Real denTemp;
		  FArrayBox &myFab = mf[ntmfi];
		  int XVEL(0), YVEL(1), ZVEL(2), DEN(3);

		  for (int ni(XVEL); ni <= ZVEL; ++ni) {

			for (int n(0); n < myFab.box().numPts(); ++n) {
			  denTemp = myFab.dataPtr(DEN)[n];
			  if (denTemp < cutoff_density)
				myFab.dataPtr(ni)[n] = 0.0;
			}
		  }
		}
	  }

	  for (int n=0; n<nVars+density_weighting; n++)
		amrData.FlushGrids(amrData.StateNumber(whichVar[n]));

	  Real timer_stop = ParallelDescriptor::second();

	  if (ParallelDescriptor::IOProcessor())
		std::cout << "   ...done (" << timer_stop-timer_start << "s)." << std::endl;

	  //
	  // Allocate memory for fftw data
	  //
      for (int iVar=0; iVar<nVars; iVar++) {
          local_data[iVar]   = (fftw_real*)    malloc(sizeof(fftw_real) * total_local_size);
          if (local_data[iVar] == NULL) amrex::Abort("Malloc fail (local_data)");
          local_data_c[iVar] = (fftw_complex*) local_data[iVar];
          for (int i=0; i<total_local_size; i++)
              local_data[iVar][i] = 0.;
      }
    
	  //
	  // Evaluate fft
	  //
	  Spectra(mf, probDomain);

	  if (ParallelDescriptor::IOProcessor()) {
		std::string suffix;
		suffix = "";
		if (div_free) {
		  if (transpose_dp)
			suffix += "_df_tdp";
		  else
			suffix += "_df";
		}
		if (density_weighting)
		  suffix += "_dw";
		suffix += ".dat";

		std::cout << "Outputting to file..." << std::endl;
		for (int iVar=0; iVar<=nVars; iVar++) {
		  std::string outfile;
		  if (iVar==nVars)
			outfile = infile + "/spectrum_count" + suffix;
		  else
			outfile = infile + "/" + whichVar[iVar] + "_spectrum" + suffix;
		  FILE* file=fopen(outfile.c_str(),"w");
		  for (int wn=0; wn<wavenumbers; wn++)
			if (div_free)
			  fprintf(file,"%i %e %e %e\n",wn,spectrum[iVar][wn],spectrumS[iVar][wn],spectrumC[iVar][wn]);
			else
			  fprintf(file,"%i %e\n",wn,spectrum[iVar][wn]);
		  fclose(file);
		}
		for (int iVar=0; iVar<nVars; iVar++) {
		  std::string outfile;
		  FILE* file;

		  // Integrals
		  outfile = infile + "/" + whichVar[iVar] + "_Int" + suffix;
		  file=fopen(outfile.c_str(),"w");
		  fprintf(file,"%e %e %e\n",Time,sum[iVar]*dxyz,sum2[iVar]*dxyz);
		  fclose(file);

		  // Qx
		  outfile = infile + "/" + whichVar[iVar] + "_Qx" + suffix;
		  file=fopen(outfile.c_str(),"w");
		  for (int i=0; i<Qix; i++)
			fprintf(file,"%e %e\n",dx*(0.5+(Real)i),Qx[iVar][i]);
		  fclose(file);

		  // Qy
		  outfile = infile + "/" + whichVar[iVar] + "_Qy" + suffix;
		  file=fopen(outfile.c_str(),"w");
		  for (int i=0; i<Qjx; i++)
			fprintf(file,"%e %e\n",dy*(0.5+(Real)i),Qy[iVar][i]);
		  fclose(file);

		  // Qz
		  outfile = infile + "/" + whichVar[iVar] + "_Qz" + suffix;
		  file=fopen(outfile.c_str(),"w");
		  for (int i=0; i<Qkx; i++)
			fprintf(file,"%e %e\n",dz*(0.5+(Real)i),Qz[iVar][i]);
		  fclose(file);
		}
		std::cout << "   ...done." << std::endl;
	  }

	  rfftwnd_mpi_destroy_plan(plan_real2cplx);
	  rfftwnd_mpi_destroy_plan(plan_cplx2real);

	  for (int iVar=0; iVar<nVars; iVar++) {
		free(local_data[iVar]);
		free(Qx[iVar]);
		free(Qy[iVar]);
		free(Qz[iVar]);
		free(spectrum[iVar]);
	  }
	  if (div_free) {
		for (int iVar=0; iVar<nVars; iVar++) {
		  free(spectrumS[iVar]);
		  free(spectrumC[iVar]);
		}
	  }

    } // iPlot

    amrex::Finalize();
}





void Spectra(MultiFab &mf, Box &probDomain)
{  
  //
  // Populate fft data
  //
  if (ParallelDescriptor::IOProcessor())
	std::cout << "Populating fft data..." << std::endl;

  Real timer_start = ParallelDescriptor::second();

  for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	const FArrayBox &myFab = mf[mfi];
	
	const int  *dlo = myFab.loVect();
	const int  *dhi = myFab.hiVect();

	const Box&  vBox = mfi.validbox();
	const int  *lo   = vBox.loVect();
	const int  *hi   = vBox.hiVect();
	const int   mfix = hi[0] - lo[0] + 1;
	const int   mfjx = hi[1] - lo[1] + 1;
	const int   mfkx = hi[2] - lo[2] + 1;

	if (verbose>1) {
	  for (int iProc=0; iProc<nProcs; iProc++) {
		if (iProc==myProc) {
		  std::cout << "--" << '\n'
					<< "Proc " << iProc << '\n'
					<< "BLix BLjx BLkx " << BLix << " " << BLjx << " " << BLkx << " " << '\n'
					<< "FTix FTjx FTkx " << FTix << " " << FTjx << " " << FTkx << " " << '\n'
					<< "mfix mfjx mfkx " << mfix << " " << mfjx << " " << mfkx << " " << '\n'
					<< "FTmx " << FTmx << '\n'
					<< "FTis FTjs FTks " << FTis << " " << FTjs << " " << FTks << " " << '\n'
					<< "FThkxpo " << FThkxpo << '\n'
					<< "local_ix local_i_start " << local_ix << " " << local_i_start << '\n'
					<< "local_jx_after_transpose local_j_start_after_transpose " << local_jx_after_transpose << " " << local_j_start_after_transpose << '\n'
					<< std::endl;
		  std::cout.flush();
		}
		ParallelDescriptor::Barrier();
	  }
	}

	const Real* density_data;
	if (density_weighting)
	  density_data = myFab.dataPtr(nVars);

	for (int iVar=0; iVar<nVars; iVar++) {

	  const Real* mf_data = myFab.dataPtr(iVar);

	  for (int i=0; i<mfix; i++) {
		int FTk=i;
		for (int j=0; j<mfjx; j++) {
		  int FTj=j;
		  for (int k=0; k<mfkx; k++) {
			int FTi=k;
			int dat_cell=(k*mfjx+j)*mfix+i;
			int fft_cell=(FTi*FTjx+FTj)*(2*FThkxpo)+FTk;
			Real val = mf_data[dat_cell];
			if (density_weighting)
			  val *= pow(density_data[dat_cell],(1.0/3.0));
			local_data[iVar][fft_cell] = val;
			sum[iVar]  += val;
			sum2[iVar] += val*val;
		  }
		}
	  }

	} // iVar

  } // mfi

  ParallelDescriptor::ReduceRealSum(sum,nVars,IOProc);
  ParallelDescriptor::ReduceRealSum(sum2,nVars,IOProc);

  Real timer_stop = ParallelDescriptor::second();

  if (ParallelDescriptor::IOProcessor())
	std::cout << "   ...done (" << timer_stop-timer_start << "s)." << std::endl;

  //
  // Perform transforms
  //
  if (ParallelDescriptor::IOProcessor())
	std::cout << "Performing real to complex transform..." << std::endl;

  timer_start = ParallelDescriptor::second();

  for (int iVar=0; iVar<nVars; iVar++)
	rfftwnd_mpi(plan_real2cplx, 1, local_data[iVar], NULL, FFTW_TRANSPOSED_ORDER);

  timer_stop = ParallelDescriptor::second();

  if (ParallelDescriptor::IOProcessor())
	std::cout << "   ...done (" << timer_stop-timer_start << "s)." << std::endl;

  if (do_filter) {
	//
	// Filter data
	//
	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Filtering..." << std::endl;
      
	timer_start = ParallelDescriptor::second();
      
	filter(mf, probDomain);
      
	timer_stop = ParallelDescriptor::second();
      
	if (ParallelDescriptor::IOProcessor())
	  std::cout << "   ...done (" << timer_stop-timer_start << "s)." << std::endl;
  }

  //
  // Integrate spectra
  //
  if (ParallelDescriptor::IOProcessor())
	std::cout << "Evaluating energy spectrum..." << std::endl;
    
  timer_start = ParallelDescriptor::second();

  // Divisor to normalise transform
  Real div = ((Real)FTix)*((Real)FTjx)*((Real)FTkx);
    
  for (int j=0; j<local_jx_after_transpose; j++) {
	int jp = j+local_j_start_after_transpose;
	int jj = FTjx-jp; if (jp<jj) jj = jp; else jj = -jj;
	jj *= FTjs; // Scale if non-cubic
	for (int i=0; i<FTix; i++) {
	  int ii = FTix-i; if (i<ii) ii = i; else ii = -ii;
	  ii *= FTis; // Scale if non-cubic
	  for (int k=0; k<FThkxpo; k++) {
		int kk = FTkx-k; if (k<kk) kk = k; else kk = -kk;
		kk *= FTks; // Scale if non-cubic

        // for binning
		int wn = (int) (0.5+sqrt((Real)(ii*ii+jj*jj+kk*kk)));

		int ccell = (j*FTix+i)*FThkxpo+k;

		if (div_free && wn<wavenumbers) {
		  // Extract components of velocity
		  Real uir = local_data_c[0][ccell].re/div;
		  Real uii = local_data_c[0][ccell].im/div;
		  Real ujr = local_data_c[1][ccell].re/div;
		  Real uji = local_data_c[1][ccell].im/div;
		  Real ukr = local_data_c[2][ccell].re/div;
		  Real uki = local_data_c[2][ccell].im/div;
		  // Construct dot product / k2
		  Real dpii, dpjj, dpkk;
		  if (transpose_dp) {
		    dpii = (Real) kk;
		    dpjj = (Real) jj;
		    dpkk = (Real) ii;
		  } else {
		    dpii = (Real) ii;
		    dpjj = (Real) jj;
		    dpkk = (Real) kk;
		  }
		  Real dpr(0.), dpi(0.);
		  Real mks = dpii*dpii + dpjj*dpjj + dpkk*dpkk;
		  if (mks>0.) {
		    dpr = ( (uir*dpii) + (ujr*dpjj) + (ukr*dpkk) ) / mks;
		    dpi = ( (uii*dpii) + (uji*dpjj) + (uki*dpkk) ) / mks;
		  }
		  // Make v the compressible part
		  Real vir = dpr * dpii;
		  Real vii = dpi * dpii;
		  Real vjr = dpr * dpjj;
		  Real vji = dpi * dpjj;
		  Real vkr = dpr * dpkk;
		  Real vki = dpi * dpkk;
		  // Subtract this from u to make u the div free part (total will get done below)
		  uir -= vir;
		  uii -= vii;
		  ujr -= vjr;
		  uji -= vji;
		  ukr -= vkr;
		  uki -= vki;
		  // Make the spectra
		  spectrumS[0][wn] += 0.5*(uir*uir + uii*uii);
		  spectrumS[1][wn] += 0.5*(ujr*ujr + uji*uji);
		  spectrumS[2][wn] += 0.5*(ukr*ukr + uki*uki);
		  spectrumC[0][wn] += 0.5*(vir*vir + vii*vii);
		  spectrumC[1][wn] += 0.5*(vjr*vjr + vji*vji);
		  spectrumC[2][wn] += 0.5*(vkr*vkr + vki*vki);
		}

		for (int iVar=0; iVar<nVars; iVar++) {
		  Real re = local_data_c[iVar][ccell].re/div;
		  Real im = local_data_c[iVar][ccell].im/div;
		  Real sq = re*re + im*im;
		  if (wn<wavenumbers)
			spectrum[iVar][wn] += 0.5 * sq;

		  re = local_data_c[iVar][ccell].re = sq;
		  im = local_data_c[iVar][ccell].im = 0.;
		}
		// Let's count the number of hits
		if (wn<wavenumbers)
		  spectrum[nVars][wn] += 1;
	  }
	}
  }

  for (int iVar=0; iVar<=nVars; iVar++) {
	ParallelDescriptor::ReduceRealSum(spectrum[iVar],wavenumbers,IOProc);
	if (div_free) {
	  ParallelDescriptor::ReduceRealSum(spectrumS[iVar],wavenumbers,IOProc);
	  ParallelDescriptor::ReduceRealSum(spectrumC[iVar],wavenumbers,IOProc);
	}
  }

  timer_stop = ParallelDescriptor::second();

  if (ParallelDescriptor::IOProcessor())
	std::cout << "   ...done (" << timer_stop-timer_start << "s)." << std::endl;

  //
  // Invert and calculate correlation functions for integral length scale
  //

  if (ParallelDescriptor::IOProcessor())
	std::cout << "Inverting for correlation tensor..." << std::endl;
    
  timer_start = ParallelDescriptor::second();

  for (int iVar=0; iVar<nVars; iVar++)
	rfftwnd_mpi(plan_cplx2real, 1, local_data[iVar], NULL, FFTW_TRANSPOSED_ORDER);

  timer_stop = ParallelDescriptor::second();

  if (ParallelDescriptor::IOProcessor())
	std::cout << "   ...done (" << timer_stop-timer_start << "s)." << std::endl;

  //
  // Evaluate x and y correlations on bottom plane, and z correlation on the part of domain we have
  //
  // Access fft data through:   int fft_cell=(FTi*FTjx+FTj)*(2*FThkxpo)+FTk;
  //
  for (int iVar=0; iVar<nVars; iVar++) {
	if (local_i_start==0) {
	  // This is the process with the bottom plane (probably proc 0)
	  // Correlation function in x direction (remember it's the transpose)
	  for (int i=0; i<Qix; i++)
		Qx[iVar][i] = local_data[iVar][i];
	  // Correlation function in y direction
	  for (int j=0; j<Qjx; j++)
		Qy[iVar][j] = local_data[iVar][j*(2*FThkxpo)];
	}
	// Correlation function in z direction (remember it's the transpose)
	// Only do the bit on this processor
	for (int k=0; ( (k<local_ix) && ((k+local_i_start)<Qkx) ); k++)
	  Qz[iVar][k+local_i_start] = local_data[iVar][k*FTjx*(2*FThkxpo)];
	
	// And reduce
	ParallelDescriptor::ReduceRealSum(Qx[iVar],Qix,IOProc);
	ParallelDescriptor::ReduceRealSum(Qy[iVar],Qjx,IOProc);
	ParallelDescriptor::ReduceRealSum(Qz[iVar],Qkx,IOProc);
  }
    
}





void plan_ffts(Box &probDomain, BoxArray &domainBoxArray, Vector<int> &pmap)
{
  plan_real2cplx = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
										   FTix, FTjx, FTkx,
										   FFTW_REAL_TO_COMPLEX,
										   FFTW_ESTIMATE);
  
  plan_cplx2real = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
										   FTix, FTjx, FTkx,
										   FFTW_COMPLEX_TO_REAL,
										   FFTW_ESTIMATE);
  //
  // FFTW prescribes the data structure
  //
  rfftwnd_mpi_local_sizes(plan_real2cplx,
						  &local_ix, &local_i_start,
						  &local_jx_after_transpose,
						  &local_j_start_after_transpose,
						  &total_local_size);
  
  if (ParallelDescriptor::IOProcessor())
	std::cout << "Total_local_size = " << total_local_size << std::endl;

  //
  // Collect the gridding
  // 
  int local_ix_array[nProcs];
  int local_i_start_array[nProcs];

  for (int iProc=0; iProc<nProcs; iProc++) {
	if (iProc==myProc) {
	  local_ix_array[iProc]      = local_ix;
	  local_i_start_array[iProc] = local_i_start;
	} else {
	  local_ix_array[iProc]      = 0;
	  local_i_start_array[iProc] = 0;
	}
  }
  
  ParallelDescriptor::ReduceIntSum(local_ix_array,nProcs,IOProc);
  ParallelDescriptor::ReduceIntSum(local_i_start_array,nProcs,IOProc); 
  
  if (ParallelDescriptor::IOProcessor()) {
	local_xlo.resize(nProcs);
	local_xhi.resize(nProcs);
	for (int iProc=0; iProc<nProcs; iProc++) {
	  local_xlo[iProc] = dx*(Real)local_i_start_array[iProc];
	  local_xhi[iProc] = local_xlo[iProc]+dx*(Real)local_ix_array[iProc];
	}

	if (verbose>1) {
	  for (int iProc=0; iProc<nProcs; iProc++) {
		std::cout << "Proc " << iProc
				  << " : local_ix = " << local_ix_array[iProc]
				  << " : local_i_start = " << local_i_start_array[iProc]
				  << " : xlo = " << local_xlo[iProc]
				  << " : xhi = " << local_xhi[iProc]
				  << std::endl;
	  }
	}
  }

  // Need to make sure box is ok, can have issues
  // e.g. for 1024^3 on 48 procs, slabs of 21 are insufficient
  //      but slabs of 22 fit on 47 procs, so the last is empty
  if (local_ix==0)
	amrex::Abort("Number of processors doesn't divide domain size (see usage note)");
 
  // Now do the slabs
  Box        tempBox(probDomain);    
  Vector<int> tempBoxSmall(nProcs,0);
  Vector<int> tempBoxBig(nProcs,0);
      
  // When using "AmrDeriveSubdomains" probLo is no longer 0, so let's add it in...
  tempBoxSmall[myProc] = probDomain.smallEnd(2) + local_i_start;
  tempBoxBig[myProc]   = probDomain.smallEnd(2) + local_i_start + local_ix - 1;
  

  for (int iProc=0; iProc<nProcs; iProc++) {
	ParallelDescriptor::ReduceIntSum(tempBoxSmall[iProc]);
	ParallelDescriptor::ReduceIntSum(tempBoxBig[iProc]);
  }

  for (int iProc=0; iProc<nProcs; iProc++) {
	tempBox.setSmall(ZDIR, tempBoxSmall[iProc]);
	tempBox.setBig(ZDIR, tempBoxBig[iProc]); 
	domainBoxArray.set(iProc, tempBox);
  }
  
  // And now the distibution mapping
  for (int iProc=0; iProc<nProcs; iProc++)
	pmap[iProc] = iProc;

  pmap[nProcs] = myProc;

}





void filter(MultiFab &mf, Box &probDomain)
{
  //
  // Allocate memory for a copy of fftw data
  // Let's do 2 copies for gt and lt threshold wavenumber
  //
  int nFiltVars(2*nVars);
  Vector<fftw_real*>    filtered_data(nFiltVars);
  Vector<fftw_complex*> filtered_data_c(nFiltVars);

  for (int iVar=0; iVar<nFiltVars; iVar++) {
    filtered_data[iVar]   = (fftw_real*)    malloc(sizeof(fftw_real) * total_local_size);   if (filtered_data[iVar] == NULL) amrex::Abort("Malloc fail (filtered_data)");
    filtered_data_c[iVar] = (fftw_complex*) filtered_data[iVar];
  }

  //
  // Allocate memory for output data
  //
  int nOutVars(nFiltVars);
  int ngrow(0);
  MultiFab mfOut;
  mfOut.define(mf.boxArray(), mf.DistributionMap(), nOutVars, ngrow, MFInfo().SetAlloc(true));
  
  //
  // Loop over filter wavenumbers
  //
  for (int iFilt=0; iFilt<nFilters; iFilt++) {
    
    if (ParallelDescriptor::IOProcessor())
      std::cout << "   Copying for filter width " << filterWN[iFilt] << "..." << std::endl;
    
    Real timer_start = ParallelDescriptor::second();
    
    // Reset data
    for (int iVar=0; iVar<nFiltVars; iVar++) {
      for (int i=0; i<total_local_size; i++)
		filtered_data[iVar][i] = 0.;
    }
    mf.setVal(0.);
    //
    // Filter
    //
    for (int iVar=0; iVar<nVars; iVar++) {

      for (int j=0; j<local_jx_after_transpose; j++) {
		int jp = j+local_j_start_after_transpose;
		int jj = FTjx-jp; if (jp<jj) jj = jp; else jj = -jj;
		jj *= FTjs; // Scale if non-cubic

		for (int i=0; i<FTix; i++) {
		  int ii = FTix-i; if (i<ii) ii = i; else ii = -ii;
		  ii *= FTis; // Scale if non-cubic

		  for (int k=0; k<FThkxpo; k++) {
			int kk = FTkx-k; if (k<kk) kk = k; else kk = -kk;
			kk *= FTks; // Scale if non-cubic
	    
			int wn = (int) (0.5+sqrt((Real)(ii*ii+jj*jj+kk*kk)));
	    
			int ccell = (j*FTix+i)*FThkxpo+k;
	    
			// This filters *spherically* based on wavenumber filterWN
			if (wn < filterWN[iFilt]) {
			  filtered_data_c[iVar][ccell].re = local_data_c[iVar][ccell].re;
			  filtered_data_c[iVar][ccell].im = local_data_c[iVar][ccell].im;
			  filtered_data_c[iVar+nVars][ccell].re = 0.;
			  filtered_data_c[iVar+nVars][ccell].im = 0.;
			} else {
			  filtered_data_c[iVar][ccell].re = 0.;
			  filtered_data_c[iVar][ccell].im = 0.;
			  filtered_data_c[iVar+nVars][ccell].re = local_data_c[iVar][ccell].re;
			  filtered_data_c[iVar+nVars][ccell].im = local_data_c[iVar][ccell].im;
			}
	    
		  }
		}
      }
      
    } // iVar
    
    Real timer_stop = ParallelDescriptor::second();
    if (ParallelDescriptor::IOProcessor())
      std::cout << "      ...done (" << timer_stop-timer_start << "s)." << std::endl;
    
    //
    // Invert
    //
    if (ParallelDescriptor::IOProcessor())
      std::cout << "   Inverting..." << std::endl;
    
    timer_start = ParallelDescriptor::second();
    
    for (int iVar=0; iVar<nFiltVars; iVar++)
      rfftwnd_mpi(plan_cplx2real, 1, filtered_data[iVar], NULL, FFTW_TRANSPOSED_ORDER);
    
    timer_stop = ParallelDescriptor::second();
    
    if (ParallelDescriptor::IOProcessor())
      std::cout << "      ...done (" << timer_stop-timer_start << "s)." << std::endl;
    
    //
    // And unload the data into the multifab
    //
    if (ParallelDescriptor::IOProcessor())
      std::cout << "   Unloading..." << std::endl;  
    timer_start = ParallelDescriptor::second();
    
    Real div = ((Real)FTix)*((Real)FTjx)*((Real)FTkx);
    
    for(MFIter mfi(mfOut); mfi.isValid(); ++mfi) {
      FArrayBox &myFab = mfOut[mfi];
      
      const int  *dlo = myFab.loVect();
      const int  *dhi = myFab.hiVect();

	  const Box&  vBox = mfi.validbox();	  
      const int  *lo   = vBox.loVect();
      const int  *hi   = vBox.hiVect();
      const int   mfix = hi[0] - lo[0] + 1;
      const int   mfjx = hi[1] - lo[1] + 1;
      const int   mfkx = hi[2] - lo[2] + 1;
      
      for (int iVar=0; iVar<nOutVars; iVar++) {
		Real* mf_data = myFab.dataPtr(iVar);
	
		for (int i=0; i<mfix; i++) {
		  int FTk=i;

		  for (int j=0; j<mfjx; j++) {
			int FTj=j;

			for (int k=0; k<mfkx; k++) {
			  int FTi=k;
			  int dat_cell=(k*mfjx+j)*mfix+i;
			  int fft_cell=(FTi*FTjx+FTj)*(2*FThkxpo)+FTk;
			  // Load into multi fab for output to plot file
			  mf_data[dat_cell] = filtered_data[iVar][fft_cell] / div;
			}
		  }
		}
	
      } // iVar
      
    } // mfi
    
    timer_stop = ParallelDescriptor::second();
    if (ParallelDescriptor::IOProcessor())
      std::cout << "      ...done (" << timer_stop-timer_start << "s)." << std::endl;
    
    //
    // And write out the plot files
    //
    if (ParallelDescriptor::IOProcessor())
      std::cout << "   Outputting data..." << std::endl;
    
    string pltfile = infile;
    for (int iVar=0; iVar<nVars; iVar++) {
      pltfile += "_" + whichVar[iVar];
    }
    if (density_weighting)
      pltfile += "_dw";
    char suffix[64];
    sprintf(suffix,"_filtered_%i",filterWN[iFilt]);
    pltfile += suffix;
    
    if (ParallelDescriptor::IOProcessor())
      if (!amrex::UtilCreateDirectory(pltfile, 0755))
		amrex::CreateDirectoryFailed(pltfile);
    ParallelDescriptor::Barrier();
    
    std::string HeaderFileName = pltfile + "/Header";
    
    static const std::string the_plot_file_type("NavierStokes-V1.1");
    
    std::ofstream os;
    
    if (ParallelDescriptor::IOProcessor()) {
      os.open(HeaderFileName.c_str());
      
      int old_prec = os.precision(15);
      
      // The plot file type
      os << the_plot_file_type << '\n';

      // The number of variables
      os << nOutVars << '\n';

      // The variable names
      for (int iVar=0; iVar<nVars; iVar++)
		os << whichVar[iVar] << "gt" << filterWN[iFilt] << '\n';

      for (int iVar=0; iVar<nVars; iVar++)
		os << whichVar[iVar] << "lt" << filterWN[iFilt] << '\n';

      // The number of space dimensions
      os << BL_SPACEDIM << '\n'; 

      // Time
      os << Time << '\n'; 

      // Finest level
      os << "0" << '\n'; 

      // Domain
      for (int i=0; i<BL_SPACEDIM; i++)
		os << probLo[i] << ' ';

      os << '\n';

      for (int i=0; i<BL_SPACEDIM; i++)
		os << probHi[i] << ' ';

      os << '\n';

      // Refinement ratios
      os << '\n';

      // Cell sizes
      os << probDomain << ' ';
      os << '\n';

      // Time steps
      os << timeSteps << '\n';

      // dx
      os << dx << ' ' << dy << ' ' << dz << '\n';

      // CoordSys & bndry
      os << "0\n0\n";
    }

    //
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    std::string Level = "Level_0";
    std::string FullPath = pltfile;

    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
      FullPath += '/';

    FullPath += Level;

    if (ParallelDescriptor::IOProcessor())
      if (!amrex::UtilCreateDirectory(FullPath, 0755))
		amrex::CreateDirectoryFailed(FullPath);

    ParallelDescriptor::Barrier();
    
    if (ParallelDescriptor::IOProcessor())
      {
		std::cout << "      Grids = " << mf.boxArray().size() << std::endl;
	
		os << 0 << ' ' << mf.boxArray().size() << ' ' << Time << '\n';
		os << timeSteps << '\n';
	
		for (int i = 0; i < mf.boxArray().size(); ++i)
		  {
			os << probLo[0]    << ' ' << probHi[0] << ' '
			   << probLo[1]    << ' ' << probHi[1] << ' '
			   << local_xlo[i] << ' ' << local_xhi[i] << '\n';
		  }
		std::string PathNameInHeader = Level;
		PathNameInHeader += BaseName;
		os << PathNameInHeader << '\n';
      }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(mfOut,TheFullPath);
    
    timer_stop = ParallelDescriptor::second();
    if (ParallelDescriptor::IOProcessor())
      std::cout << "      ...done (" << timer_stop-timer_start << "s)." << std::endl;

  } // iFilt

}
