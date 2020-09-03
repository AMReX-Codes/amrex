
#include "AMReX_timemory.H"
#include "AMReX_ParmParse.H"

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

extern "C"
{
#if defined(TIMEMORY_USE_MPIP_LIBRARY)
    void timemory_register_mpip ();
    void timemory_deregister_mpip ();
#endif
//
#if defined(TIMEMORY_USE_OMPT_LIBRARY)
    void timemory_register_ompt ();
    void timemory_deregister_ompt ();
#endif
}

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
static bool&
get_use_mpip_wrappers ()
{
#if defined(TIMEMORY_USE_MPIP_LIBRARY)
    static auto _instance = tim::get_env<bool> ("BL_PROFILE_MPI", true);
    return _instance;
#else
    static auto _instance = false;
    return _instance;
#endif
}
//
static bool&
get_use_ompt_wrappers ()
{
#if defined(TIMEMORY_USE_OMPT_LIBRARY)
    static auto _instance = tim::get_env<bool> ("BL_PROFILE_OMP", true);
    return _instance;
#else
    static auto _instance = false;
    return _instance;
#endif
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
namespace
{
template <size_t Idx>
void
insert_avail_impl (std::vector<std::string>& _vec)
{
    using type           = typename tim::component::enumerator<Idx>::type;
    constexpr bool value = tim::component::enumerator<Idx>::value;
    if (value)
        _vec.push_back (tim::component::properties<type>::id ());
}
//
template <size_t... Idx>
auto print_avail_impl (std::index_sequence<Idx...>)
{
    std::vector<std::string> _vec;
    TIMEMORY_FOLD_EXPRESSION (insert_avail_impl<Idx> (_vec));
    return _vec;
}
}  // namespace
//
//--------------------------------------------------------------------------------------//
//
void
BL_timemory_configure (int argc, char** argv)
{
    if (argc < 1 || argv == nullptr)
    {
        char* _argv = new char[1];
        _argv[0]    = '\0';
        tim::timemory_init (1, &_argv);
        return;
    }

    tim::timemory_init (argc, argv);
}
//
//--------------------------------------------------------------------------------------//
//
void
BL_timemory_initialize (int argc, char** argv, bool parse_args)
{
    std::string default_components = "wall_clock";
    auto        _env = tim::get_env ("BL_PROFILE_COMPONENTS", default_components) + ", " +
                tim::settings::global_components ();

    std::vector<std::string> components   = tim::delimit (_env, ", \t:;");
    std::vector<std::string> papi_events  = {};
    std::vector<std::string> cupti_events = {};
    bool                     print_avail  = false;
    bool                     per_thread   = !tim::settings::collapse_threads ();
    bool                     per_process  = !tim::settings::collapse_processes ();

    ParmParse pp ("timemory");
    pp.query ("enabled", tim::settings::enabled ());
    // diagnostic
    pp.query ("verbose", tim::settings::verbose ());
    pp.query ("debug", tim::settings::debug ());
    pp.query ("print_available", print_avail);
    // misc
    pp.query ("output_path", tim::settings::output_path ());
    pp.query ("input_path", tim::settings::input_path ());
    // output formats
    pp.query ("text_output", tim::settings::text_output ());
    pp.query ("cout_output", tim::settings::cout_output ());
    pp.query ("plot_output", tim::settings::plot_output ());
    pp.query ("json_output", tim::settings::json_output ());
    pp.query ("diff_output", tim::settings::diff_output ());
    pp.query ("time_output", tim::settings::time_output ());
    pp.query ("dart_output", tim::settings::dart_output ());
    pp.query ("flamegraph_output", tim::settings::flamegraph_output ());
    // collection schemes
    pp.query ("flat_profile", tim::settings::flat_profile ());
    pp.query ("timeline_profile", tim::settings::timeline_profile ());
    pp.query ("per_process", per_process);
    pp.query ("per_thread", per_thread);
    // component collection configuration
    pp.queryarr ("components", components);
    pp.queryarr ("cpu_hw_counters", papi_events);
    pp.queryarr ("gpu_hw_counters", cupti_events);
    // output formatting
    pp.query ("precision", tim::settings::precision ());
    pp.query ("scientific", tim::settings::scientific ());
    pp.query ("timing_units", tim::settings::timing_units ());
    pp.query ("timing_precision", tim::settings::timing_precision ());
    pp.query ("timing_scientific", tim::settings::timing_scientific ());
    pp.query ("memory_units", tim::settings::memory_units ());
    pp.query ("memory_precision", tim::settings::memory_precision ());
    pp.query ("memory_scientific", tim::settings::memory_scientific ());
    // empirical roofline toolkit (ERT) settings
    pp.query ("ert_max_data_size", tim::settings::ert_max_data_size ());
    pp.query ("ert_min_working_size", tim::settings::ert_min_working_size ());
    pp.query ("ert_num_threads", tim::settings::ert_num_threads ());

#if defined(AMREX_HDF5_WRAPPERS)
    pp.query ("profile_hdf5", get_use_hdf5_wrappers ());
#endif

#if defined(TIMEMORY_USE_MPIP_LIBRARY)
    pp.query ("profile_mpi", get_use_mpip_wrappers ());
#endif

#if defined(TIMEMORY_USE_OMPT_LIBRARY)
    pp.query ("profile_omp", get_use_ompt_wrappers ());
#endif

    // only print on rank 0 for MPI and/or UPC++
    if (print_avail && tim::dmp::rank () == 0)
    {
        auto available =
            print_avail_impl (std::make_index_sequence<TIMEMORY_COMPONENTS_END>{});

        std::stringstream sout;
        sout << "Available timemory components:\n";
        for (const auto& itr : available)
            sout << "    - " << itr << '\n';
        std::cout << sout.str () << std::flush;
    }

    tim::settings::collapse_processes () = !per_process;
    tim::settings::collapse_threads ()   = !per_thread;

    if (!components.empty ())
    {
        tim::configure<BL_timemory_bundle> (components);

#if defined(TIMEMORY_USE_MPIP_LIBRARY)
        tim::configure<tim::component::user_mpip_bundle> (components);
#endif
#if defined(TIMEMORY_USE_OMPT_LIBRARY)
        tim::configure<tim::component::user_ompt_bundle> (components);
#endif
    }
    else
    {
        // disable the component
        tim::trait::runtime_enabled<BL_timemory_bundle>::set (false);
    }

    if (!papi_events.empty ())
    {
        auto& _evts = tim::settings::papi_events ();
        for (const auto& itr : papi_events)
            _evts = TIMEMORY_JOIN (",", _evts, itr);
        BL_timemory_bundle::configure<tim::component::papi_vector> ();
    }

    if (!cupti_events.empty ())
    {
        auto& _evts = tim::settings::cupti_events ();
        for (const auto& itr : cupti_events)
            _evts = TIMEMORY_JOIN (",", _evts, itr);
        BL_timemory_bundle::configure<tim::component::cupti_counters> ();
    }

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

#if defined(TIMEMORY_USE_MPIP_LIBRARY)
    if (get_use_mpip_wrappers ())
        timemory_register_mpip ();
#endif

#if defined(TIMEMORY_USE_OMPT_LIBRARY)
    if (get_use_ompt_wrappers ())
        timemory_register_ompt ();
#endif
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

#if defined(TIMEMORY_USE_MPIP_LIBRARY)
    if (get_use_mpip_wrappers ())
        timemory_deregister_mpip ();
#endif

#if defined(TIMEMORY_USE_OMPT_LIBRARY)
    if (get_use_ompt_wrappers ())
        timemory_deregister_ompt ();
#endif

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
