#include <AMReX_FileSystem.H>
#include <AMReX_Print.H>
#include <AMReX.H>

#if defined(_WIN32) || __cplusplus >= 201703L

#include <filesystem>
#include <system_error>

namespace amrex {
namespace FileSystem {

bool
Exists (std::string const& filename)
{
    std::error_code ec;
    return std::filesystem::exists(std::filesystem::path{filename});
}

std::string
CurrentPath ()
{
    std::error_code ec;
    auto path = std::filesystem::current_path(ec);
    if (ec) {
        amrex::AllPrint() << "amrex::FileSystem::CurrentPath failed. " << ec.message() << std::endl;
    }
    return path.string();
}

bool
Remove (std::string const& filename)
{
    std::error_code ec;
    bool r = std::filesystem::remove(std::filesystem::path(filename),ec);
    return (!ec and r);
}

}}

#else

#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>

namespace amrex {
namespace FileSystem {

bool
Exists (std::string const& filename)
{
    struct stat statbuff;
    return(lstat(filename.c_str(), &statbuff) != -1);
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

}}

#endif
