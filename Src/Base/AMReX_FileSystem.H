#ifndef AMREX_FILE_SYSTEM_H_
#define AMREX_FILE_SYSTEM_H_
#include <AMReX_Config.H>

#include <string>

#ifdef _WIN32
using mode_t = unsigned short;
#else
#include <sys/types.h>  // for mode_t
#endif

namespace amrex::FileSystem {

/**
 * \brief Create a directory and any missing parent directories.
 *
 * @param path directory path to create
 * @param mode requested permissions for newly-created directories on POSIX
 *        platforms; ignored on Windows
 * @param verbose print a diagnostic message when directory creation fails
 * @return true if the directory already existed or was created without error,
 *         false if directory creation failed. This differs from the return
 *         value of `std::filesystem::create_directories`, which is true only
 *         when a directory was actually created and false both when the
 *         directory already existed and when creation failed with an error.
 */
bool
CreateDirectories (std::string const& path, mode_t mode, bool verbose = false);

/**
 * \brief Return the current working directory.
 *
 * @return current working directory path, or an empty string if it cannot be
 *         determined
 */
std::string
CurrentPath ();

/**
 * \brief Check if the given path exists.
 * Symbolic links are checked without following the link, so a dangling
 * symbolic link is considered to exist.
 *
 * @param filename path to check
 * @return true if the path exists, false otherwise
 */
bool
Exists (std::string const& filename);

/**
 * \brief Remove a file, symbolic link, or empty directory.
 *
 * @param filename path to remove
 * @return true if a path was removed, false otherwise
 */
bool
Remove (std::string const& filename);

/**
 * \brief Recursively remove a path.
 *
 * @param p path to remove
 * @return true if removal completed without error, including when the path did
 *         not exist; false otherwise
 */
bool
RemoveAll (std::string const& p);

}

#endif
