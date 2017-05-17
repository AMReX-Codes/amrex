#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BLFort.H>

using namespace amrex;

// declare a fortran subroutine
extern "C"
{
  void work(const int& flag, const int* lo, const int* hi, BL_FORT_FAB_ARG(dfab));
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    // define the lower and upper corner of a 3D domain
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(63,63,63); 
 
    // build a box for the domain
    Box domain(domain_lo, domain_hi);

    // build a box array from the 64^3 domain box
    BoxArray ba(domain);
    // break the box array into 32^3 boxes
    ba.maxSize(32);

    DistributionMapping dm{ba};

    // build a multifab on the box array with 1 component, 0 ghost cells
    MultiFab data(ba, dm, 1, 0);  

    // loop over boxes and initialize the data
	
    bool tiling = true;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(data,tiling); mfi.isValid(); ++mfi)
    { // Tiling with default tile size in FabArray.cpp

	const Box& bx = mfi.tilebox();  // box for this tile

	// call a fortran subroutine
	work(0, bx.loVect(), bx.hiVect(),
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
    std::cout << std::endl;
    std::cout << "min      = " << data.min(0)  << std::endl;
    std::cout << "max      = " << data.max(0)  << std::endl;
    std::cout << "max norm = " << data.norm0() << std::endl;
    std::cout << "L1  norm = " << data.norm1() << std::endl;
    std::cout << "L2  norm = " << data.norm2() << std::endl;

    MultiFab data2(ba, dm, 1, 0);
    MultiFab::Copy(data2, data, 0, 0, 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(data2); mfi.isValid(); ++mfi) 
    { // Without explicitly turn on tiling, this becomes OMP over fabs.

	const Box& bx = mfi.validbox();  // valid box for this fab

	// call a fortran subroutine
	work(1, bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(data2[mfi]));
    }


    // This time we will explicitly control tile size.
    IntVect tilesize(AMREX_D_DECL(10240,8,32));

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(data2,tilesize); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	// call a fortran subroutine
	work(2, bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(data2[mfi]));
    }

    std::cout << "\nThe second MultiFab should be the same up to roundoff errors." << std::endl;
    std::cout << "min      = " << data2.min(0)  << std::endl;
    std::cout << "max      = " << data2.max(0)  << std::endl;
    std::cout << "max norm = " << data2.norm0() << std::endl;
    std::cout << "L1  norm = " << data2.norm1() << std::endl;
    std::cout << "L2  norm = " << data2.norm2() << std::endl;
    
    MultiFab::Subtract(data2, data, 0, 0, 1, 0);
    std::cout << "\nThe difference should be almost zero." << std::endl;
    std::cout << "min      = " << data2.min(0)  << std::endl;
    std::cout << "max      = " << data2.max(0)  << std::endl;
    std::cout << "max norm = " << data2.norm0() << std::endl;
    std::cout << "L1  norm = " << data2.norm1() << std::endl;
    std::cout << "L2  norm = " << data2.norm2() << std::endl;    

    amrex::Finalize();
}

