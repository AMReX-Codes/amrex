#include "MyTest.H"

void MyTest::initializePoiseuilleData(int ilev) {
#if (AMREX_SPACEDIM == 2)
    initializePoiseuilleDataFor2D(ilev);
#else
    initializePoiseuilleDataFor3D(ilev);
#endif
}
