
#include <CNS.H>

using namespace amrex;

void
CNS::variableSetUp ()
{
    read_params();

    // xxxxx set up StateData

    ErrorSetUp();
}

void
CNS::variableCleanUp ()
{
    desc_lst.clear();
}
