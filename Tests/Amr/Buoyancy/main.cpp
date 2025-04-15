#include "AMReX.H"

int MyMain();

int main(int argc, char** argv) {
    amrex::Initialize(argc, argv);
    const int ret = MyMain();
    amrex::Finalize();
    return ret;
}

int MyMain()
{
    return 0;
}
