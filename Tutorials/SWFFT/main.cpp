#include <AMReX.H>
#include <SWFFT_Solver.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    {
        BL_PROFILE("main()");
        SWFFT_Solver sw_solver;
        sw_solver.solve();
    }
    
    amrex::Finalize();
    return 0;
}
