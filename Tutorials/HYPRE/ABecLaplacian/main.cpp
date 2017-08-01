#include <AMReX.H>
#include <ABL.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    {
        BL_PROFILE("main()");
        ABL abl;
        abl.solve();
    }
    
    amrex::Finalize();
    return 0;
}

