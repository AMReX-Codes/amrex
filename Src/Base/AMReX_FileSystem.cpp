#include <AMReX_FileSystem.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX.H>

#if defined(_WIN32) // || __cplusplus >= 201703L

#include <filesystem>
#include <system_error>

namespace amrex {
namespace FileSystem {

bool
CreateDirectories (std::string const& p, mode_t /*mode*/, bool verbose)
{
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path{p}, ec);
    if (ec and verbose) {
        amrex::AllPrint() << "amrex::UtilCreateDirectory failed to create "
                          << p << ": " << ec.message() << std::endl;
    }
    return !ec;    
}

bool
Exists (std::string const& filename)
{
    std::error_code ec;
    bool r = std::filesystem::exists(std::filesystem::path{filename}, ec);
    if (ec and amrex::Verbose() > 0) {
        amrex::AllPrint() << "amrex::FileSystem::Exists failed. " << ec.message() << std::endl;
    }
    return r;
}

std::string
CurrentPath ()
{
    std::error_code ec;
    auto path = std::filesystem::current_path(ec);
    if (ec and amrex::Verbose() > 0) {
        amrex::AllPrint() << "amrex::FileSystem::CurrentPath failed. " << ec.message() << std::endl;
    }
    return path.string();
}

bool
Remove (std::string const& filename)
{
    std::error_code ec;
    bool r = std::filesystem::remove(std::filesystem::path{filename},ec);
    return !ec;
}

bool
RemoveAll (std::string const& p)
{
    std::error_code ec;
    std::filesystem::remove_all(std::filesystem::path{p},ec);
    return !ec;
}

}}

#else

#include <cstdio>
#include <cstddef>
#include <cstring>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>

namespace amrex {
namespace FileSystem {

bool
CreateDirectories (std::string const& path, mode_t mode, bool verbose)
{
    bool retVal(false);
    Vector<std::pair<std::string, int> > pathError;

    const char* path_sep_str = "/";

    if (path.length() == 0 || path == path_sep_str) {
        return true;
    }

    errno = 0;

    if(std::strchr(path.c_str(), *path_sep_str) == 0) {
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
        (void) strcpy(dir, path.c_str());

        char *slash = std::strchr(dir, *path_sep_str);

        if(dir[0] == *path_sep_str) {  // full pathname.
            do {
                if(*(slash+1) == 0) {
                    break;
                }
                if((slash = std::strchr(slash+1, *path_sep_str)) != 0) {
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
            } while((slash = std::strchr(slash+1, *path_sep_str)) != 0);

            errno = 0;
            if(mkdir(dir, mode) < 0 && errno != EEXIST) {
                retVal = false;
            } else {
                retVal = true;
            }
            pathError.push_back(std::make_pair(dir, errno));
        }

        delete [] dir;
    }

    if(retVal == false  || verbose == true) {
      for(int i(0); i < pathError.size(); ++i) {
          amrex::AllPrint()<< "amrex::UtilCreateDirectory:: path errno:  "
                           << pathError[i].first << " :: "
                           << strerror(pathError[i].second)
                           << std::endl;
      }
    }

    return retVal;
}

bool
Exists (std::string const& filename)
{
    struct stat statbuff;
    return (lstat(filename.c_str(), &statbuff) != -1);
}

std::string
CurrentPath ()
{
    constexpr int bufSize = 1024;
    char temp[bufSize];
    char *rCheck = getcwd(temp, bufSize);
    if(rCheck == 0) {
        amrex::Abort("**** Error:  getcwd buffer too small.");
    }
    return std::string(rCheck);
}

bool
Remove (std::string const& filename)
{
    return unlink(filename.c_str());
}

bool
RemoveAll (std::string const& p)
{
    if (p.size() >= 1990) {
        amrex::Error("FileSystem::RemoveAll: Path name too long");
        return false;
    }
    char command[2000];
    std::snprintf(command, 2000, "\\rm -rf %s", p.c_str());;
    int retVal = std::system(command);
    if (retVal == -1 || WEXITSTATUS(retVal) != 0) {
        amrex::Error("Removing old directory failed.");
        return false;
    }
    return true;
}

}}

#endif
