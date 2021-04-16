#include "MyTest.H"

void MyTest::initializeTriangleWaveData(int ilev) {
#if (AMREX_SPACEDIM == 2)
    initializeTriangleWaveDataFor2D(ilev);
#endif
}
