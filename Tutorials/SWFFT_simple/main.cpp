#include <AMReX.H>
#include <SWFFT_Test.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    {
        BL_PROFILE("main()");
        SWFFT_Test sw_test;
        sw_test.computeFFT();
    }
    
    amrex::Finalize();
    return 0;
}
