#include <AMReX.H>
#include <AMReX_Print.H>

void init();
void work();
void init2();
void work2();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        init();
        work();

        init2();
        work2();
    }
    amrex::Finalize();
}
