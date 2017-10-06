#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

//***********************************************
// This can be moved into a C++ header file and #included here instead

#include <AMReX_BLFort.H>

using namespace amrex;

extern "C"
{
  void work_on_data(Real* data, const int* ng, const int* nc, 
		    const int* lo, const int* hi);
}

//***********************************************

void main_main ()
{

  // define a box with bounds (0,0) to (15,15) in 2D
  // or a box with bounds (0,0,0) to (15,15,15) in 3D
#if (BL_SPACEDIM == 2)
  IntVect lo(0,0), hi(15,15);
#elif (BL_SPACEDIM == 3)
  IntVect lo(0,0,0), hi(15,15,15);
#endif
  Box bx(lo,hi);

  BoxArray ba(bx);  // ba has one 16^2 box
  ba.maxSize(8);    // ba now has four 8^2 boxes

  // This stores the physical coordinates of the problem domain
  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++)
  {
    real_box.setLo(n,-1.0);
    real_box.setHi(n, 1.0);
  }

  // This says we are using Cartesian coordinates
  int coord = 0;

  // This says the domain is periodic in each direction
  int is_per[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    is_per[i] = 1;
  }

  // Define a Geometry object to store problem domain physical coordinates,
  // incides, coordinate system, and periodicity
  Geometry geom;
  geom.define(bx, &real_box, coord, is_per);

  // define a multifab with 2 components and 6 ghost cells
  int Ncomp = 2 , Nghost = 6;

  // the "1" means only 1 MultiFab; you can declare more if needed
  Vector <std::unique_ptr<MultiFab> > data(1);

  DistributionMapping dm(ba);

  // build the "0th" MultiFab
  data[0].reset(new MultiFab(ba, dm, Ncomp, Nghost));

  // MFIter is a "MultiFab Iterator" that loops over 
  // the grids in the MultiFab
  for ( MFIter mfi(*data[0]); mfi.isValid(); ++mfi)
  {
    // get the box associated with the current grid
    const Box& bx = mfi.validbox();
    work_on_data((*data[0])[mfi].dataPtr(), &Nghost, &Ncomp, 
		 bx.loVect(), bx.hiVect());
  }

  // fill periodic domain boundary and neighboring grid ghost cells
  data[0]->FillBoundary(geom.periodicity());
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}
