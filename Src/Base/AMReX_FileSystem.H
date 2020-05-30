#ifndef AMREX_FILE_SYSTEM_H_
#define AMREX_FILE_SYSTEM_H_

#include <string>

namespace amrex {
namespace FileSystem {

std::string
CurrentPath ();

bool
Exists (std::string const& filename);

bool
Remove (std::string const& filename);

bool
RemoveAll (std::string const& p); // recursive remove

}}

#endif
