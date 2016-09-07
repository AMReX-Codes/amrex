#include <iostream>

#include <BoxLib.H>
#include <MultiFabUtil.H>
#include <MultiFab_C_F.H>

extern "C" {
    void ff ();
}

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    // define the lower and upper corner of a 3D domain
    IntVect domain_lo(0,0,0); 
    IntVect domain_hi(7,7,7); 
 
    // build a box for the domain
    Box domain(domain_lo, domain_hi);

    // build a box array from the 64^3 domain box
    BoxArray ba(domain);
    // break the box array into 32^3 boxes
    ba.maxSize(4);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Number of boxes: " << ba.size() << std::endl;
    }    

    // This defines the physical size of the box.  Right now the box is [-1,1] in each direction.
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
	real_box.setLo(n,-1.0);
	real_box.setHi(n, 1.0);
    }

    // This says we are using Cartesian coordinates
    int coord = 0;
    
    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1; 
    
    // This defines a Geometry object which is useful for writing the plotfiles  
    Geometry geom(domain, &real_box, coord, is_per);

    // build a multifab on the box array with 2 component, 1 ghost cells
    int nc = 2, ng=1;
    MultiFab mf(ba, nc, ng);

    const DistributionMapping& dmap = mf.DistributionMap();

    MultiFab_C_to_F mfctof(geom, dmap, ba);

    mfctof.share(mf, "data");

    ff(); // a Fortran function

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "\n" << " In C++\n";
	std::cout << "     " << "ncomp  = " << mf.nComp() << "\n";
	std::cout << "     " << "nghost = " << mf.nGrow() << "\n";
    }

    mf.setVal(1.0);
    Real norm = mf.norm1(nc-1,ng);
    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "     " << "1-norm   " << norm << "\n";
    }

    mf.setVal(3.0, 0); // not setting ghost cells
    norm = mf.norm1(nc-1,ng);
    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "     " << "1-norm   " << norm << "\n";
    }

    mf.FillBoundary(geom.periodicity());
    norm = mf.norm1(nc-1,ng);
    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "     " << "1-norm   " << norm << "\n";
	std::cout << std::endl;
    }

    BoxLib::Finalize();
}

