#include <AMReX_Machine.H>
#include <AMReX.H>
#include <AMReX_String.H>

#include <cstdlib>

namespace amrex::Machine
{

namespace {
    std::string s_name;
}

void Initialize ()
{
    // Known machines:
    //   nersc.perlmutter: NERSC_HOST=perlmutter
    //                     LMOD_SITE_NAME=perlmutter
    //   olcf.frontier   : LMOD_SITE_NAME=OLCF
    //                     LMOD_SYSTEM_NAME=frontier

    auto const* env_nersc_host = std::getenv("NERSC_HOST");
    auto const* env_lmod_site_name = std::getenv("LMOD_SITE_NAME");
    auto const* env_lmod_system_name = std::getenv("LMOD_SYSTEM_NAME");
    auto const* env_slurm_cluster_name = std::getenv("SLURM_CLUSTER_NAME");

    if (env_nersc_host && env_lmod_system_name) {
        s_name = std::string("nersc.");
        s_name.append(env_lmod_system_name);
    } else if (env_lmod_site_name && env_lmod_system_name) {
        s_name = std::string(env_lmod_site_name);
        s_name.append(".").append(env_lmod_system_name);
    } else if (env_slurm_cluster_name) {
        s_name = std::string(env_slurm_cluster_name);
    }

    if ( ! s_name.empty()) {
        s_name = amrex::toLower(std::move(s_name));
    }

    amrex::ExecOnFinalize(Machine::Finalize);
}

void Finalize () {}

std::string const& name ()
{
    return s_name;
}

}
