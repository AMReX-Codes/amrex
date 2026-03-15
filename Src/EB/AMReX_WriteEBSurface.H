#ifndef AMREX_WRITE_EBSURFACE_H
#define AMREX_WRITE_EBSURFACE_H
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H> // IWYU pragma: export
#include <AMReX_Geometry.H> // IWYU pragma: export

namespace amrex {

class EBFArrayBoxFactory;

void WriteEBSurface (const amrex::BoxArray & ba, const amrex::DistributionMapping & dmap, const amrex::Geometry & geom,
                     const amrex::EBFArrayBoxFactory * ebf);

}

#endif
