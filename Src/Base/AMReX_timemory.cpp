
#include "AMReX_timemory.H"
#include <timemory/timemory.hpp>
#include <timemory/utility/argparse.hpp>

#include <unordered_map>
#include <memory>
#include <string>
#include <set>

#if defined(AMREX_USE_HDF5) && defined(TIMEMORY_USE_GOTCHA)
#define AMREX_HDF5_WRAPPERS
#include "hdf5.h"
#endif

namespace amrex
{
//
#if defined(AMREX_HDF5_WRAPPERS)
struct hdf5_wrappers
{};
using hdf5_bundle_t = tim::component_tuple<AMREX_TIMEMORY_COMPONENTS>;
//
// this component will wrap all hdf5 functions with 'hdf5_bundle_t'
using hdf5_gotcha_t = tim::component::gotcha<20, hdf5_bundle_t, hdf5_wrappers>;
//
using hdf5_handle = tim::lightweight_tuple<hdf5_gotcha_t>;
//
static void
configure_hdf5_wrappers ()
{
    hdf5_gotcha_t::get_initializer () = [] () {
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 0, H5Aclose);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 1, H5Awrite);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 2, H5Dclose);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 3, H5Fclose);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 4, H5Gclose);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 5, H5Pclose);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 6, H5Pset_alignment);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 7, H5Pset_all_coll_metadata_ops);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 8, H5Pset_coll_metadata_write);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 9, H5Pset_dxpl_async);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 10, H5Pset_dxpl_mpio);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 11, H5Pset_fapl_mpio);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 12, H5Pset_vol_async);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 13, H5Sclose);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 14, H5Sselect_hyperslab);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 15, H5Tclose);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 16, H5Tinsert);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 17, H5Tset_size);
        TIMEMORY_C_GOTCHA (hdf5_gotcha_t, 18, H5Tset_strpad);
    };
}
//
#else
using hdf5_handle = tim::lightweight_tuple<>;
#endif
//
namespace
{
static hdf5_handle&
get_hdf5_wrapper ()
{
    static hdf5_handle _instance{};
    return _instance;
}
//
static bool&
get_use_hdf5_wrappers ()
{
    static auto _instance = tim::get_env<bool> ("BL_PROFILE_HDF5", false);
    return _instance;
}
//
}  // namespace
//
static auto
configure_hdf5_wrappers ()
{
    fprintf (stderr, "Warning! Profiling HDF5 wrappers is not supported! "
                     "Requires GOTCHA support in timemory.\n");
}
//
//--------------------------------------------------------------------------------------//
//
BLProfiler_timemory_region_map_t&
BL_timemory_get_regions ()
{
    static thread_local BLProfiler_timemory_region_map_t _instance{};
    return _instance;
}
//
//--------------------------------------------------------------------------------------//
//
void
BL_timemory_configure (int argc, char** argv)
{
    using parser_t = tim::argparse::argument_parser;

    if (argc < 1 || argv == nullptr)
    {
        char* _argv = new char[1];
        _argv[0]    = '\0';
        tim::timemory_init (1, &_argv);
        return;
    }

    tim::timemory_init (argc, argv);
    auto parser = parser_t (argv[0]);

    parser.enable_help ();
    parser.add_argument ({ "--blt-verbose" }, "Set the verbosity of profiler")
        .max_count (1)
        .action ([] (parser_t& p) {
            if (p.get<std::string> ("blt-verbose").empty ())
                tim::settings::verbose () = 1;
            else
                tim::settings::verbose () = p.get<int> ("blt-verbose");
        });
    parser
        .add_argument ({ "--blt-debug" },
                       "Enable debug output for profiler (also affected by verbosity")
        .max_count (1)
        .action ([] (parser_t& p) {
            tim::settings::debug () = true;
            // note: this intentionally sets verbose
            if (p.get<std::string> ("blt-debug").empty ())
                tim::settings::verbose () = 0;
            else
                tim::settings::verbose () = p.get<int> ("blt-debug");
        });
    parser.add_argument ({ "--blt-output" }, "Output directory for BL_timemory")
        .max_count (1)
        .action ([] (parser_t& p) {
            tim::settings::output_path () = p.get<std::string> ("blt-output");
        });
    auto fmt =
        parser
            .add_argument (
                { "--blt-format" },
                "Output formats for BL_timemory (see `timemory-avail -S | grep OUTPUT`)")
            .choices ({
                "text",
                "cout",
                "plot",
                "json",
                "diff",
                "time",
                "dart",
                "flamegraph",
            })
            .action ([] (parser_t& p) {
                auto choices = p.get<std::set<std::string>> ("blt-format");
                auto has     = [&] (const char* c) {
                    return choices.find (c) != choices.end ();
                };
                tim::settings::text_output ()       = has ("text");
                tim::settings::cout_output ()       = has ("cout");
                tim::settings::plot_output ()       = has ("plot");
                tim::settings::json_output ()       = has ("json");
                tim::settings::diff_output ()       = has ("diff");
                tim::settings::diff_output ()       = has ("time");
                tim::settings::dart_output ()       = has ("dart");
                tim::settings::flamegraph_output () = has ("flamegraph");
            });
    parser
        .add_argument ({ "--blt-input" },
                       "Input directory for BL_timemory (enable run comparison mode)")
        .max_count (1)
        .action ([] (parser_t& p) {
            tim::settings::input_path ()  = p.get<std::string> ("blt-input");
            tim::settings::time_output () = true;  // time-stamped folders for this run
            tim::settings::diff_output () = true;  // difference between input folder
                                                   // and time-stamp folder for run
        });
    parser
        .add_argument (
            { "--blt-per-process" },
            "Output individual results for each process (i.e. rank) instead of "
            "reporting the aggregation")
        .count (0)
        .action ([] (parser_t&) { tim::settings::collapse_processes () = false; });
    parser
        .add_argument (
            { "--blt-per-process" },
            "Output individual results for each process (i.e. rank) instead of "
            "reporting the aggregation")
        .count (0)
        .action ([] (parser_t&) { tim::settings::collapse_processes () = false; });
    parser
        .add_argument ({ "--blt-per-thread" },
                       "Output individual results for each thread instead of "
                       "reporting the aggregation")
        .count (0)
        .action ([] (parser_t&) { tim::settings::collapse_threads () = false; });
    parser
        .add_argument ({ "--blt-flat-profile" },
                       "All call-stack entries will have a depth of zero")
        .count (0)
        .action ([] (parser_t&) { tim::settings::flat_profile () = true; });
    parser
        .add_argument ({ "--blt-timeline-profile" },
                       "All call-stack entries will be unique")
        .count (0)
        .action ([] (parser_t&) { tim::settings::flat_profile () = true; });
    parser
        .add_argument ({ "--blt-components" },
                       "Components to profile (see `timemory-avail -Cs` for valid string "
                       "identifier) ")
        .action ([] (parser_t& p) {
            tim::settings::flat_profile () = true;
            auto vec = p.get<std::vector<std::string>> ("blt-components");
            if (!vec.empty ())
                tim::configure<BL_timemory_bundle> (vec);
        });
    parser
        .add_argument ({ "--blt-papi" },
                       "Enable collecting CPU HW counters via PAPI (see "
                       "`timemory-avail -H` and `papi_native_avail`)")
        .action ([] (parser_t& p) {
            auto vec = p.get<std::vector<std::string>> ("blt-papi");
            if (!vec.empty ())
            {
                auto& _evts = tim::settings::papi_events ();
                for (const auto& itr : vec)
                    _evts = TIMEMORY_JOIN (",", _evts, itr);
                BL_timemory_bundle::configure<tim::component::papi_vector> ();
            }
        });
    parser
        .add_argument ({ "--blt-cupti" },
                       "Enable collecting NVIDIA GPU HW counters via CUPTI (see "
                       "`timemory-avail -H`)")
        .action ([] (parser_t& p) {
            auto vec = p.get<std::vector<std::string>> ("blt-cupti");
            if (!vec.empty ())
            {
                auto& _evts = tim::settings::cupti_events ();
                for (const auto& itr : vec)
                    _evts = TIMEMORY_JOIN (",", _evts, itr);
                BL_timemory_bundle::configure<tim::component::cupti_counters> ();
            }
        });
    parser.add_argument ({ "--blt-hdf5" }, "Profile HDF5 functions")
        .count (0)
        .action ([] (parser_t&) { get_use_hdf5_wrappers () = true; });

    // if CI=true (i.e. continuous integration) in set env
    // (Travis sets this, for example) enable echoing the
    // performance metrics for dashboard
    if (tim::get_env<bool> ("CI", false))
        fmt.set_default (std::string ("cout, dart"));

    auto err = parser.parse (argc, argv);
    if (err)
    {
        std::cerr << err << std::endl;
        parser.print_help ();
        exit (EXIT_FAILURE);
    }
}
//
//--------------------------------------------------------------------------------------//
//
void
BL_timemory_initialize (int argc, char** argv, bool parse_args)
{
    std::string default_components = "wall_clock, ";
#if defined(_OPENMP)
    default_components += "ompt_handle, ";
#endif
    auto _env = tim::get_env ("BL_PROFILE_COMPONENTS", default_components) + ", " +
                tim::settings::global_components ();
    tim::configure<BL_timemory_bundle> (tim::delimit (_env, ", \t:;"));

    if (argc > 0 && argv && parse_args)
    {
        BL_timemory_configure (argc, argv);
    }
    else if (argc < 1 || argv == nullptr)
    {
        char* _argv = new char[1];
        _argv[0]    = '\0';
        tim::timemory_init (1, &_argv);
    }
    else
    {
        tim::timemory_init (argc, argv);
    }

    if (get_use_hdf5_wrappers ())
    {
        // generate wrappers
        configure_hdf5_wrappers ();
        // start the wrappers
        get_hdf5_wrapper ().start ();
    }
}
//
void
BL_timemory_finalize ()
{
    if (get_use_hdf5_wrappers ())
    {
        // stop the wrappers
        get_hdf5_wrapper ().stop ();
    }
    BL_timemory_get_regions ().clear ();
    tim::timemory_finalize ();
}
//
TimemoryProfiler::TimemoryProfiler (std::string&& funcname, bool _start)
    : m_handle (std::forward<std::string> (funcname))
{
    if (_start)
        m_handle.start ();
}
//
TimemoryProfiler::TimemoryProfiler (const char* funcname, bool _start)
    : m_handle (funcname)
{
    if (_start)
        m_handle.start ();
}
//
TimemoryProfiler::TimemoryProfiler (size_t _hash, bool _start)
    : m_handle (_hash)
{
    if (_start)
        m_handle.start ();
}
//
TimemoryProfiler::TimemoryProfiler (const tim::source_location::captured& _loc,
                                    bool                                  _start)
    : m_handle (_loc)
{
    if (_start)
        m_handle.start ();
}
//
TimemoryProfiler::~TimemoryProfiler () { m_handle.stop (); }
//
void
TimemoryProfiler::start ()
{
    m_handle.start ();
}
//
void
TimemoryProfiler::stop ()
{
    m_handle.stop ();
}
//
TimemoryTplProfiler::value_type*
TimemoryTplProfiler::get ()
{
    static value_type _instance{};
    return &_instance;
}
//
void
TimemoryTplProfiler::start ()
{
    get ()->start ();
}
//
void
TimemoryTplProfiler::stop ()
{
    get ()->stop ();
}
//
}  // namespace amrex

//--------------------------------------------------------------------------------------//
//
//  Everything after this is "unused"
//
//--------------------------------------------------------------------------------------//

template class tim::component_tuple<AMREX_TIMEMORY_COMPONENTS>;
template class tim::auto_tuple<AMREX_TIMEMORY_COMPONENTS>;

namespace
{
bool
bl_timemory_load()
{
    std::string default_components = "wall_clock, ";
#if defined(_OPENMP)
    default_components += "ompt_bundle, ";
#endif

    auto functor = [=]() {
        return tim::get_env("BL_PROFILE_COMPONENTS", default_components) + ", " +
               tim::settings::global_components();
    };

    constexpr auto Idx = amrex::BL_timemory_bundle_idx;
    tim::env::get_user_bundle_variables ()[Idx].push_back (functor);

    // default settings
    tim::settings::mpi_init()       = false;
    tim::settings::mpi_finalize()   = false;
    tim::settings::upcxx_init()     = false;
    tim::settings::upcxx_finalize() = false;
    tim::settings::cout_output ()   = false;

    // if CI=true (i.e. continuous integration) is set in env
    // enable cout the performance metrics for logs
    // and dart for dashboard
    if(tim::get_env<bool>("CI", false))
    {
        tim::settings::cout_output() = true;
        tim::settings::dart_output() = true;
    }

    // allow env overrides
    tim::settings::parse();

    // the following only currently works on linux
    // via the /proc/<PID>/command_line procfs file
    tim::config::read_command_line (amrex::BL_timemory_initialize);

    return true;
}

// appends the configurations when the library is loaded
auto bl_timemory_is_loaded = bl_timemory_load();
}  // namespace
