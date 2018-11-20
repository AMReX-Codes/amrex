
#include <CNS.H>

using namespace amrex;

void
CNS::restart (Amr& papa, std::istream& is, bool bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    if (do_reflux && level > 0) {
//        flux_reg.define(grids, papa.boxArray(level-1),
//                        dmap, papa.DistributionMap(level-1),
//                        geom, papa.Geom(level-1),
//                        papa.refRatio(level-1), level, NUM_STATE);
    }

    buildMetrics();
}

void 
CNS::checkPoint (const std::string& dir, std::ostream& os, VisMF::How how, bool dump_old) 
{
    AmrLevel::checkPoint(dir, os, how, dump_old);
}

void
CNS::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
    BL_PROFILE("CNS::writePlotFile()");
    AmrLevel::writePlotFile(dir, os, how);
}
