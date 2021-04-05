#include "MyTest.H"

void MyTest::initializeLinearData(int ilev) {
#if (AMREX_SPACEDIM == 2)
    initializeLinearDataFor2D(ilev);
#else
    initializeLinearDataFor3D(ilev);
#endif
}
