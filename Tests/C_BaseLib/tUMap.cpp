

#include <fstream>

#include <AMReX_BLFort.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BaseUMAP.H>

//BL_FORT_PROC_DECL(FILLFAB,fillfab)(Real* d, const int* nx, const int* ny);

using namespace amrex;

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    //
    // This in only for 2D.
    //

    const int NY = 4;
    const int NX = 3;
    const int NZ = 2;

    Box bx(IntVect::TheZeroVector(),IntVect(D_DECL(NX-1,NY-1,NZ-1)));

    FArrayBox fab(bx,2);
    fab.setVal(1.e200);

    BaseUmap<Real> umap(bx, 2);

    Real v;
    v = 2.0;
    umap.setVal( v, IntVect(2,1,0), 0, 1);

    for (IntVect iv=bx.smallEnd(); iv<=bx.bigEnd(); bx.next(iv)) {
        v = umap.getVal(iv, 0, 1);
        std::cout << iv << " has value " << v <<std::endl;
//
    }

    std::cout << "hi";

    amrex::Finalize();

    return 0;
}
