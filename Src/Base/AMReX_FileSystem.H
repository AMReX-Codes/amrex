#ifndef AMREX_FILE_SYSTEM_H_
#define AMREX_FILE_SYSTEM_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_FileSystem.H
 * \brief File-system utility functions (directory creation, path query, removal).
 */

#include <string>

#ifdef _WIN32
using mode_t = unsigned short;
#else
#include <sys/types.h>  // for mode_t
#endif

namespace amrex::FileSystem {

/// \ingroup amrex_utilities
//! Create all directories along \p path with permissions \p mode; returns true on success. When \p verbose is true diagnostic messages are printed.
bool
CreateDirectories (std::string const& path, mode_t mode, bool verbose = false);

/// \ingroup amrex_utilities
//! Return the current working directory as an absolute path string.
std::string
CurrentPath ();

/// \ingroup amrex_utilities
//! Return true if \p filename exists in the file system.
bool
Exists (std::string const& filename);

/// \ingroup amrex_utilities
//! Remove the file or empty directory at \p filename; returns true on success.
bool
Remove (std::string const& filename);

/// \ingroup amrex_utilities
//! Recursively remove \p p and all its contents; returns true on success.
bool
RemoveAll (std::string const& p); // recursive remove

}

#endif
