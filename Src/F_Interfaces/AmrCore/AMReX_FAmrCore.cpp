
#include <AMReX_FAmrCore.H>


amrex::FAmrCore::FAmrCore ()
{
    std::cout << "constructor of FAmrCore\n";
}

amrex::FAmrCore::~FAmrCore ()
{
    std::cout << "destructor of FAmrCore\n";
}

void
amrex::FAmrCore::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    std::cout << "in ErrorEst\n";
}
