/* Copyright 2026 The AMReX Community
 *
 * License: BSD-3-Clause-LBNL
 */
#include <AMReX_ParticleHeader.H>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <istream>
#include <sstream>

namespace amrex {

void
ParticleHeader::parse (std::istream& is)
{
    is >> version;
    AMREX_ASSERT(!version.empty());

    // What do our version strings mean?
    // "Version_One_Dot_Zero" -- hard-wired to write out in double precision.
    // "Version_One_Dot_One" -- can write out either as either single or double precision.
    // Appended to the latter version string are either "_single" or "_double" to
    // indicate how the particles were written.
    // "Version_Two_Dot_Zero" -- this is the AMReX particle file format
    // "Version_Two_Dot_One" -- expanded particle ids to allow for 2**39-1 per proc
    convert_ids = false;
    if (version.find("Version_Two_Dot_One") != std::string::npos) {
        convert_ids = true;
    }
    if (version.find("Version_One_Dot_Zero") != std::string::npos) {
        how = "double";
    }
    else if (version.find("Version_One_Dot_One")  != std::string::npos ||
             version.find("Version_Two_Dot_Zero") != std::string::npos ||
             version.find("Version_Two_Dot_One") != std::string::npos) {
        if (version.find("_single") != std::string::npos) {
            how = "single";
        }
        else if (version.find("_double") != std::string::npos) {
            how = "double";
        }
        else {
            std::string msg("ParticleHeader::parse(): bad version string: ");
            msg += version;
            amrex::Error(msg.c_str());
        }
    }
    else {
        std::string msg("ParticleHeader::parse(): unknown version string: ");
        msg += version;
        amrex::Abort(msg.c_str());
    }

    is >> dim;

    is >> num_real;
    real_comp_names.resize(num_real);
    for (int i = 0; i < num_real; ++i) {
        is >> real_comp_names[i];
    }

    is >> num_int;
    int_comp_names.resize(num_int);
    for (int i = 0; i < num_int; ++i) {
        is >> int_comp_names[i];
    }

    is >> is_checkpoint;

    is >> num_particles;

    is >> next_id;

    is >> finest_level;
}

ParticleHeader
ParticleHeader::read (const std::string& dir, const std::string& file)
{
    AMREX_ASSERT(!dir.empty());
    AMREX_ASSERT(!file.empty());

    std::string fullname = dir;
    if (!fullname.empty() && fullname.back() != '/') {
        fullname += '/';
    }
    fullname += file;

    std::string HdrFileName = fullname;
    if (!HdrFileName.empty() && HdrFileName.back() != '/') {
        HdrFileName += '/';
    }
    HdrFileName += "Header";

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(HdrFileName, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream HdrFile(fileCharPtrString, std::istringstream::in);

    ParticleHeader header;
    header.parse(HdrFile);
    return header;
}

}
