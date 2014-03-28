
#include <winstd.H>

#include <fstream>

#include <BLFort.H>
#include <Utility.H>
#include <MultiFab.H>

BL_FORT_PROC_DECL(FILLFAB,fillfab)(Real* d, const int* nx, const int* ny);

int
main (int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);
    //
    // This in only for 2D.
    //
#if BL_SPACEDIM==2

    const int NY = 1720;
    const int NX = 2000;

    Box bx(IntVect::TheZeroVector(),IntVect(NX-1,NY-1));

    FArrayBox fab(bx,2);

    fab.setVal(1.e200);

    std::ofstream ofs;

    ofs.open("out.fab", std::ios::out|std::ios::trunc);

    if (!ofs.good())
        BoxLib::FileOpenFailed("out.fab");

    BL_FORT_PROC_CALL(FILLFAB,fillfab)(fab.dataPtr(), &NX, &NY);

    fab.writeOn(ofs);

    ofs.close();

    if (!ofs.good())
        BoxLib::Error("Write failed");

#endif

    BoxLib::Finalize();

    return 0;
}
