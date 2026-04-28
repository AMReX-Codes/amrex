#include <AMReX_FileSystem.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX.H>

#include <filesystem>
#include <system_error>

#if !defined(_WIN32)
#include <cstdio>
#include <cstddef>
#include <cstring>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif

namespace amrex::FileSystem {

bool
Exists (std::string const& filename)
{
    std::error_code ec;
    auto const status = std::filesystem::symlink_status(std::filesystem::path{filename}, ec);
    bool const r = std::filesystem::exists(status);
    if (ec && (status.type() != std::filesystem::file_type::not_found) && amrex::Verbose() > 0) {
        amrex::AllPrint() << "amrex::FileSystem::Exists failed. " << ec.message() << '\n';
    }
    return r;
}

std::string
CurrentPath ()
{
    std::error_code ec;
    auto path = std::filesystem::current_path(ec);
    if (ec && amrex::Verbose() > 0) {
        amrex::AllPrint() << "amrex::FileSystem::CurrentPath failed. " << ec.message() << '\n';
    }
    return path.string();
}

bool
Remove (std::string const& filename)
{
    std::error_code ec;
    bool r = std::filesystem::remove(std::filesystem::path{filename},ec);
    if (ec && amrex::Verbose() > 0) {
        amrex::AllPrint() << "amrex::FileSystem::Remove failed. " << ec.message() << '\n';
    }
    return r;
}

bool
RemoveAll (std::string const& p)
{
    std::error_code ec;
    std::filesystem::remove_all(std::filesystem::path{p},ec);
    if (ec) {
        amrex::Error("amrex::FileSystem::RemoveAll failed to remove " + p
                     + ": " + ec.message());
        return false;
    } else {
        return true;
    }
}

#if defined(_WIN32) // || __cplusplus >= 201703L

bool
CreateDirectories (std::string const& p, mode_t /*mode*/, bool verbose)
{
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path{p}, ec);
    if (ec && verbose) {
        amrex::AllPrint() << "amrex::UtilCreateDirectory failed to create "
                          << p << ": " << ec.message() << '\n';
    }
    return !ec;
}

#else

bool
CreateDirectories (std::string const& path, mode_t mode, bool verbose)
{
    bool retVal = false;
    Vector<std::pair<std::string, int> > pathError;

    const char* path_sep_str = "/";

    if (path.empty() || path == path_sep_str) {
        return true;
    }

    errno = 0;

    if(std::strchr(path.c_str(), *path_sep_str) == nullptr) {
        //
        // No slashes in the path.
        //
        errno = 0;
        if(mkdir(path.c_str(), mode) < 0 && errno != EEXIST) {
            retVal = false;
        } else {
            retVal = true;
        }
        pathError.push_back(std::make_pair(path, errno));
    } else {
        //
        // Make copy of the directory pathname so we can write to it.
        //
        char *dir = new char[path.length() + 1];
        (void) std::strncpy(dir, path.c_str(), path.length()+1);

        char *slash = std::strchr(dir, *path_sep_str);

        if(dir[0] == *path_sep_str) {  // full pathname.
            do {
                if(*(slash+1) == 0) {
                    break;
                }
                if((slash = std::strchr(slash+1, *path_sep_str)) != nullptr) { // NOLINT(bugprone-assignment-in-if-condition)
                    *slash = 0;
                }
                errno = 0;
                if(mkdir(dir, mode) < 0 && errno != EEXIST) {
                    retVal = false;
                } else {
                    retVal = true;
                }
                pathError.push_back(std::make_pair(dir, errno));
                if(slash) {
                    *slash = *path_sep_str;
                }
            } while(slash);

        } else {  // relative pathname.

            do {
                *slash = 0;
                errno = 0;
                if(mkdir(dir, mode) < 0 && errno != EEXIST) {
                    retVal = false;
                } else {
                    retVal = true;
                }
                pathError.push_back(std::make_pair(dir, errno));
                *slash = *path_sep_str;
            } while((slash = std::strchr(slash+1, *path_sep_str)) != nullptr); // NOLINT(bugprone-assignment-in-if-condition)

            errno = 0;
            if(mkdir(dir, mode) < 0 && errno != EEXIST) {
                retVal = false;
            } // No `else` because the retVal has been set in the do-while loop above.
            pathError.push_back(std::make_pair(dir, errno));
        }

        delete [] dir;
    }

    if(retVal == false  || verbose == true) {
      for(auto & i : pathError) {
          amrex::AllPrint()<< "amrex::UtilCreateDirectory:: path errno:  "
                           << i.first << " :: "
                           << strerror(i.second)
                           << '\n';
      }
    }

    return retVal;
}

#endif

}
