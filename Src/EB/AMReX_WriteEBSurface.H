#ifndef AMREX_WRITE_EBSURFACE_H
#define AMREX_WRITE_EBSURFACE_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

namespace amrex {

class EBFArrayBoxFactory;

void WriteEBSurface (const amrex::BoxArray & ba, const amrex::DistributionMapping & dmap, const amrex::Geometry & geom,
                     const amrex::EBFArrayBoxFactory * ebf); 

}

#endif

