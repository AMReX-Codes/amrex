#ifndef AMREX_CONDUIT_BLUEPRINT_H_
#define AMREX_CONDUIT_BLUEPRINT_H_

///////////////////////////////////////////////////////////////////////////////
//
// Conduit Mesh Blueprint Support for AMReX
//
// This file is only included when USE_CONDUIT = TRUE
//
///////////////////////////////////////////////////////////////////////////////

#include <string>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#endif

///////////////////////////////////////////////////////////////////////////////
// These functions generate a Conduit Blueprint Mesh representation of
// AMReX AMR meshes and particle containers. 
// 
// Data from MultiFABs and ParticleContainers are zero copied and annotated
// in Conduit Node Trees. This is the format used to pass meshes to ALPINE
// Ascent for in situ rendering and analysis.
//
// For more details, please see:
//
// Conduit:
//  Repo: https://github.com/llnl/conduit
//  Docs http://llnl-conduit.readthedocs.io/en/latest/
//  Blueprint Docs: http://llnl-conduit.readthedocs.io/en/latest/blueprint.html
//
// Ascent:
//  Ascent Repo: http://github.com/alpine-dav/ascent
//  Ascent Docs: http://ascent.readthedocs.io/en/latest/
//
//
// Or ping Cyrus Harrison <cyrush@llnl.gov> or Matt Larsen <<larsen30@llnl.gov>
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///  Current Support:
///  * 2D + 3D
///  * single + multi-level (w/o nesting)
///  * ghosts (indicator field created using `grow`)
///  * particles
///
/// TODO:
///  * AMR nesting
///
///////////////////////////////////////////////////////////////////////////////

#include <conduit/conduit.hpp>

namespace amrex
{

    // Creates an in-memory Mesh Blueprint representation of a single level
    // AMReX mesh.
    //
    // This method's signature mirrors WriteSingleLevelPlotfile(), sans the
    // file name parameters.
    void SingleLevelToBlueprint (const MultiFab &mf,
                                 const Vector<std::string> &varnames,
                                 const Geometry &geom,
                                 Real time_value,
                                 int level_step,
                                 conduit::Node &bp_mesh);

    // Creates an in-memory Mesh Blueprint representation of a multi-level
    // AMReX mesh.
    //
    // This method's signature mirrors WriteMultiLevelPlotfile(), sans the
    // file name parameters.
    void MultiLevelToBlueprint (int n_levels,
                                const Vector<const MultiFab*> &mfs,
                                const Vector<std::string> &varnames,
                                const Vector<Geometry> &geoms,
                                Real time_value,
                                const Vector<int> &level_steps,
                                const Vector<IntVect> &ref_ratio,
                                conduit::Node &bp_mesh);

#ifdef AMREX_PARTICLES
    // Creates an in-memory Mesh Blueprint representation of an AMReX
    // particle container.
    //
    //  This method's signature is similar to ParticleContainer::Checkpoint()
    //   real_comp_names are first mapped in order to the aos then soa reals
    //   int_comp_names are first mapped in order to the aos then soa ints
    //
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    void ParticleContainerToBlueprint (const ParticleContainer<NStructReal,
                                                               NStructInt,
                                                               NArrayReal,
                                                               NArrayInt> &pc,
                                       const Vector<std::string> &real_comp_names,
                                       const Vector<std::string> &int_comp_names,
                                       conduit::Node &bp_mesh);
#endif

    // Writes a Mesh Blueprint representation to a set of files that
    // can be visualized in VisIt using the Blueprint plugin.
    //
    // Note: The default protocol "json", is ASCII-based, which is great
    //       for debugging, but can be quite slow for large data sets.
    //       For binary output, you can use:
    //         "conduit_bin"
    //         "hdf5" (if conduit was built with hdf5 support)
    //
    void WriteBlueprintFiles (const conduit::Node &bp_mesh,
                              const std::string &fname_base = "blueprint_out_",
                              int   step  = 0,
                              const std::string &protocol = "json");

}

#ifdef AMREX_PARTICLES
    #include "AMReX_Conduit_Blueprint_ParticlesI.H"
#endif


#endif
