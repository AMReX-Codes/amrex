#ifndef AMREX_WRITE_EBSURFACE_H
#define AMREX_WRITE_EBSURFACE_H
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

#include <string>

/**
 * \file AMReX_WriteEBSurface.H
 *
 * Writes EB surface geometry for diagnostic visualization.
 */

namespace amrex {

class EBFArrayBoxFactory;

/**
 * \brief Dump EB surfaces described by \p ebf to pvd format file.
 */
void WriteEBSurface (const amrex::BoxArray & ba, const amrex::DistributionMapping & dmap, const amrex::Geometry & geom,
                     const amrex::EBFArrayBoxFactory * ebf);

/**
 * \brief Triangulate and write the planar EB reconstruction to an ASCII STL file.
 */
void WriteEBSurfaceSTL (const amrex::BoxArray& ba,
                        const amrex::DistributionMapping& dmap,
                        const amrex::Geometry& geom,
                        const amrex::EBFArrayBoxFactory* ebf,
                        std::string const& filename);

}

#endif
