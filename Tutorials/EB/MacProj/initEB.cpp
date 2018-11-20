#include <AMReX_ParmParse.H>
#include <AMReX_EB2_IF.H>
#include "initEB.H"

using namespace amrex;

void initEB (const Geometry& geom)
{
    EB2::Build(geom, 0, 30);
}


