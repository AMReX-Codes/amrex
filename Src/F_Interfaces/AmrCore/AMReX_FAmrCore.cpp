
#include <AMReX_FAmrCore.H>


amrex::FAmrCore::FAmrCore ()
    : amrex::AmrCore()
{
}

amrex::FAmrCore::~FAmrCore ()
{
}

void
amrex::FAmrCore::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    std::cout << "in ErrorEst\n";
}
