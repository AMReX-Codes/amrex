#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BLFort.H>
#include <AMReX_Print.H>

// declare a fortran subroutine
extern "C"
{
  void work(const int* lo, const int* hi, BL_FORT_FAB_ARG(dfab));
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::Print() << "AMReX version " << amrex::Version() << "\n";

    {
	auto tstart = amrex::ParallelDescriptor::second();

	// define the lower and upper corner of a 3D domain
	amrex::IntVect domain_lo(0 , 0, 0); 
	amrex::IntVect domain_hi(63,63,63); 
	
	// build a box for the domain
	amrex::Box domain(domain_lo, domain_hi);
	
	// build a box array from the 64^3 domain box
	amrex::BoxArray ba(domain);
	// break the box array into 32^3 boxes
	ba.maxSize(32);

	// build a distribution map for this box array
	amrex::DistributionMapping dm(ba);
		
	// build a multifab on the box array with 1 component, 0 ghost cells
	amrex::MultiFab data(ba, dm, 1, 0);  

	// loop over boxes and do some work
	for (amrex::MFIter mfi(data); mfi.isValid(); ++mfi)
	{
	    const amrex::Box& bx = mfi.validbox();  // box for this fab
	    
	    // call a fortran subroutine
	    work(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(data[mfi]));
	}
	
	// print some information for checking
	/* I got
	   min      = 7.945268926e-11
	   max      = 1
	   max norm = 1
	   L1  norm = 8680.319857
	   L2  norm = 56.24354515
	*/
	amrex::AllPrint()
	    << " Proc. " << amrex::ParallelDescriptor::MyProc() << "\n"
	    << "    min      = " << data.min(0)  << "\n"
	    << "    max      = " << data.max(0)  << "\n"
	    << "    max norm = " << data.norm0() << "\n"
	    << "    L1  norm = " << data.norm1() << "\n"
	    << "    L2  norm = " << data.norm2() << "\n";

	auto tend = amrex::ParallelDescriptor::second();
	amrex::Print() << "Run time: " << tend-tstart << "\n";
    }

    amrex::Finalize();
}

