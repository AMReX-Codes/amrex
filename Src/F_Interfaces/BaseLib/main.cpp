
#include <BoxLib.H>

extern "C" { void fmain(); }

int main (int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    fmain();

    BoxLib::Finalize();
}
