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
 *
 * After the metadata prefix, the "Header" file contains a per-level grid
 * table that locates each grid's binary particle data on disk. It can be
 * parsed with parse_grid_table(); read() returns a fully parsed header
 * including the table.
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
     * \brief One entry of the per-level grid table.
     *
     * Locates the binary particle data of a single grid: the particles of a
     * grid are stored contiguously in the data file "DATA_XXXXX" (with XXXXX
     * = \p which), starting \p where bytes into that file.
     */
    struct GridEntry
    {
        //! index of the binary data file (DATA_XXXXX) holding this grid's particles
        int which = 0;
        //! number of particles stored for this grid
        int count = 0;
        //! byte offset of this grid's particle data within the data file
        Long where = 0;
    };

    /**
     * \brief Per level and grid: where each grid's particle data is stored.
     *
     * Filled by parse_grid_table(); grids.size() == finest_level + 1 and
     * grids[lev].size() is the number of grids of level \p lev. Grids without
     * particles have a zero \p count.
     */
    Vector<Vector<GridEntry>> grids;

    /**
     * \brief Parse the metadata prefix from an input stream.
     *
     * Consumes exactly the leading metadata block (version through finest
     * level), leaving the stream positioned right after the finest-level entry
     * and before the per-level grid table. The table can then be parsed from
     * the same stream with parse_grid_table().
     *
     * \param is input stream positioned at the start of a particle Header
     */
    void parse (std::istream& is);

    /**
     * \brief Parse the per-level grid table into \p grids.
     *
     * Must be called on a stream positioned right after the metadata prefix,
     * i.e. directly following a call to parse().
     *
     * \param is input stream positioned at the start of the grid table
     */
    void parse_grid_table (std::istream& is);

    /**
     * \brief Read and parse the "Header" file of a particle plotfile/checkpoint.
     *
     * Parses the complete header: the metadata prefix and the per-level grid
     * table. Collective: the file is read on the I/O process and broadcast.
     *
     * \param dir   plotfile/checkpoint directory
     * \param file  particle sub-directory name (e.g. "particle0")
     * \return the parsed header metadata
     */
    static ParticleHeader read (const std::string& dir, const std::string& file);
};

}

#endif
