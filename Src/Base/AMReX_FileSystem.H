#ifndef AMREX_FILE_SYSTEM_H_
#define AMREX_FILE_SYSTEM_H_
#include <AMReX_Config.H>

#include <string>

#ifdef _WIN32
typedef unsigned short mode_t;
#else
#include <sys/types.h>  // for mode_t
#endif

namespace amrex::FileSystem {

bool
CreateDirectories (std::string const& path, mode_t mode, bool verbose = false);

std::string
CurrentPath ();

bool
Exists (std::string const& filename);

bool
Remove (std::string const& filename);

bool
RemoveAll (std::string const& p); // recursive remove

}

#endif
