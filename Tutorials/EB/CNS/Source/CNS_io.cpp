
#include <CNS.H>

using namespace amrex;

void
CNS::restart (Amr& papa, std::istream& is, bool bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    // xxxxx need to do flux register stuff

    buildMetrics();

    //  xxxxx initialize EB
}

void 
CNS::checkPoint (const std::string& dir, std::ostream& os, VisMF::How how, bool dump_old) 
{
    AmrLevel::checkPoint(dir, os, how, dump_old);
}

void
CNS::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
    AmrLevel::writePlotFile(dir, os, how);
}

