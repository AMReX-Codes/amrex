#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <BLFort.H>

// declare a fortran subroutine
BL_FORT_PROC_DECL(WORK, work) (const int* lo, const int* hi, BL_FORT_FAB_ARG(dfab));

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    // define the lower and upper corner of a 3D domain
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(63,63,63); 
 
    // build a box for the domain
    Box domain(domain_lo, domain_hi);

    // build a box array from the 64^3 domain box
    BoxArray ba(domain);
    // break the box array into 32^3 boxes
    ba.maxSize(32);

    // build a multifab on the box array with 1 component, 0 ghost cells
    MultiFab data(ba, 1, 0);  

    // loop over boxes and do some work
    for (MFIter mfi(data); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.validbox();  // box for this fab

	// call a fortran subroutine
	BL_FORT_PROC_CALL(WORK, work)
	    (bx.loVect(), bx.hiVect(),
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
    std::cout << "min      = " << data.min(0)  << std::endl;
    std::cout << "max      = " << data.max(0)  << std::endl;
    std::cout << "max norm = " << data.norm0() << std::endl;
    std::cout << "L1  norm = " << data.norm1() << std::endl;
    std::cout << "L2  norm = " << data.norm2() << std::endl;

    BoxLib::Finalize();
}

