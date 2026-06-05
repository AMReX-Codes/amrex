/* Copyright 2026 The AMReX Community
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef AMREX_PARTICLEHEADER_H_
#define AMREX_PARTICLEHEADER_H_
#include <AMReX_Config.H>

#include <AMReX_INT.H>
#include <AMReX_Vector.H>

#include <iosfwd>
#include <string>

namespace amrex {

/**
 * \brief Metadata prefix of an AMReX particle "Header" file.
 *
 * This is the single source of truth for parsing the leading metadata of a
 * particle plotfile/checkpoint "Header" file, as written by
 * WriteBinaryParticleData. It is used both by ParticleContainer::Restart (to
 * read a file into a container) and standalone (e.g. from language bindings) to
 * discover the on-disk layout - number and names of the real/int components,
 * the precision, and whether the file is a checkpoint - before constructing or
 * configuring a container to read it.
 *
 * The component counts follow the same convention as the on-disk format: for
 * pure Struct-of-Arrays particles the AMREX_SPACEDIM position components are
 * implicit and \p num_real counts only the remaining (runtime) real components.
 */
struct ParticleHeader
{
    //! raw version string, e.g. "Version_Two_Dot_One_double"
    std::string version;
    //! precision the particle data was written in: "single" or "double"
    std::string how;
    //! whether particle ids need to be converted (Version_Two_Dot_One and later)
    bool convert_ids = false;
    //! AMREX_SPACEDIM the file was written with
    int dim = 0;
    //! number of real components (pure SoA: excludes the AMREX_SPACEDIM positions)
    int num_real = 0;
    //! names of the real components, size() == num_real
    Vector<std::string> real_comp_names;
    //! number of integer components
    int num_int = 0;
    //! names of the integer components, size() == num_int
    Vector<std::string> int_comp_names;
    //! true if the file is a checkpoint (full restart), false for a plotfile
    bool is_checkpoint = false;
    //! total number of particles in the file
    Long num_particles = 0;
    //! the next particle id to hand out (maxnextid)
    Long next_id = 0;
    //! finest level present in the file
    int finest_level = 0;

    /**
     * \brief Parse the metadata prefix from an input stream.
     *
     * Consumes exactly the leading metadata block (version through finest
     * level), leaving the stream positioned right after the finest-level entry
     * and before the per-level grid table. This lets ParticleContainer::Restart
     * continue reading the grid offsets from the same stream.
     *
     * \param is input stream positioned at the start of a particle Header
     */
    void parse (std::istream& is);

    /**
     * \brief Read and parse the "Header" file of a particle plotfile/checkpoint.
     *
     * Collective: the file is read on the I/O process and broadcast.
     *
     * \param dir   plotfile/checkpoint directory
     * \param file  particle sub-directory name (e.g. "particle0")
     * \return the parsed header metadata
     */
    static ParticleHeader read (const std::string& dir, const std::string& file);
};

}

#endif
