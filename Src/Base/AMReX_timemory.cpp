
#if defined(AMREX_USE_GPU) && defined(AMREX_USE_OMP) && defined(TIMEMORY_USE_OMPT_LIBRARY)
#   undef TIMEMORY_USE_OMPT_LIBRARY
#endif

#include "AMReX_timemory.H"
#include "AMReX_ParmParse.H"

#include <timemory/timemory.hpp>

using api_t = amrex::BL_timemory_tag;

#if defined(TIMEMORY_USE_MPI) && defined(TIMEMORY_USE_GOTCHA) && !defined(TIMEMORY_MPI_GOTCHA)
#    define TIMEMORY_MPI_GOTCHA
#endif

#if defined(TIMEMORY_MPI_GOTCHA)
#include <timemory/components/gotcha/mpip.hpp>

// this component will track communication sizes
// via the mpi_data_tracker_t component below
TIMEMORY_DECLARE_COMPONENT(mpi_comm_data)

using mpi_data_tracker_t = tim::component::data_tracker<float, api_t>;

// TIMEMORY_STATISTICS_TYPE(mpi_data_tracker_t, float)
TIMEMORY_DEFINE_CONCRETE_TRAIT(uses_memory_units, mpi_data_tracker_t, true_type)
TIMEMORY_DEFINE_CONCRETE_TRAIT(is_memory_category, mpi_data_tracker_t, true_type)
TIMEMORY_DEFINE_CONCRETE_TRAIT(base_has_last, mpi_data_tracker_t, true_type)
TIMEMORY_DEFINE_CONCRETE_TRAIT(base_has_accum, mpi_data_tracker_t, true_type)

using mpi_comm_data_t    = tim::component::mpi_comm_data;
using mpip_bundle_t      = tim::component_tuple<mpi_comm_data_t,
						AMREX_TIMEMORY_COMPONENTS>;

#endif

#if defined(TIMEMORY_USE_OMPT_LIBRARY)
#include <timemory/components/ompt.hpp>

extern "C" void timemory_register_ompt();
extern "C" void timemory_deregister_ompt();
#endif

#include <unordered_map>
#include <memory>
#include <string>
#include <set>

namespace amrex
{
//
namespace
{
//
#if defined(TIMEMORY_MPI_GOTCHA)
static bool&
get_use_mpip_wrappers ()
{
    static auto _instance = tim::get_env<bool> ("BL_PROFILE_MPI", true);
    return _instance;
}
#endif
//
#if defined(TIMEMORY_USE_OMPT_LIBRARY)
static bool&
get_use_ompt_wrappers ()
{
    static auto _instance = tim::get_env<bool> ("BL_PROFILE_OMP", true);
    return _instance;
}
#endif
//
}  // namespace
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
BL_timemory_initialize (int argc, char** argv)
{
    static bool _once = false;
    if (_once)
        return;
    _once = true;

    static auto _manager  = tim::manager::instance();
    static auto _settings = tim::settings::instance();

    // default settings
    tim::settings::mpi_init()       = false;
    tim::settings::mpi_finalize()   = false;
    tim::settings::upcxx_init()     = false;
    tim::settings::upcxx_finalize() = false;

    // if CI=true (i.e. continuous integration) is set in env
    // enable cout the performance metrics for logs
    // and dart for dashboard
    if(tim::get_env<bool>("CI", false))
    {
        tim::settings::cout_output() = true;
        tim::settings::dart_output() = true;
    }

    if (argc < 1 || argv == nullptr)
    {
        char* _argv = new char[1];
        _argv[0]    = '\0';
        tim::timemory_init (1, &_argv);
        return;
    }

    tim::timemory_init (argc, argv);
    
    tim::settings::cout_output () = false;
    tim::settings::plot_output () = false;

    tim::consume_parameters(_manager, _settings);
}
//
//--------------------------------------------------------------------------------------//
//
void
BL_timemory_configure ()
{
    std::string default_components = "wall_clock";
    auto        _env = tim::get_env ("BL_PROFILE_COMPONENTS", default_components) + ", " +
                tim::settings::global_components ();

    std::vector<std::string> components   = tim::delimit (_env, ", \t:;");
    std::vector<std::string> papi_events  = {};
    std::vector<std::string> cupti_events = {};
    std::vector<std::string> mpip_permit  = {};
    std::vector<std::string> mpip_reject  = {};
    bool                     print_avail  = false;
    bool                     mpip_comm    = true;
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

#if defined(TIMEMORY_MPI_GOTCHA)
    pp.query ("profile_mpi", get_use_mpip_wrappers ());
    pp.query ("profile_mpi_comm_data", mpip_comm);
    pp.queryarr ("profile_mpi_permit", mpip_permit);
    pp.queryarr ("profile_mpi_reject", mpip_reject);
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

    auto vec_to_set = [](const auto& _vec)
    {
	using type       = std::decay_t<decltype(_vec)>;
	using value_type = typename type::value_type;
	std::set<value_type> _ret{};
	for(const auto& itr : _vec)
	    _ret.insert(itr);
	return _ret;
    };
    
#if defined(TIMEMORY_MPI_GOTCHA)
    if (get_use_mpip_wrappers ())
    {
	auto _permit = vec_to_set(mpip_permit);
	auto _reject = vec_to_set(mpip_reject);
	tim::component::configure_mpip<mpip_bundle_t, api_t>(_permit, _reject);
	tim::component::activate_mpip<mpip_bundle_t, api_t>();
	tim::trait::runtime_enabled<mpi_comm_data_t>::set(mpip_comm);
    }
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
#if defined(TIMEMORY_MPI_GOTCHA)
    if (get_use_mpip_wrappers ())
    {
	tim::component::deactivate_mpip<mpip_bundle_t, api_t>(0);
        tim::trait::runtime_enabled<mpi_comm_data_t>::set(false);
    }
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

template class tim::component_tuple<AMREX_TIMEMORY_COMPONENTS>;
template class tim::auto_tuple<AMREX_TIMEMORY_COMPONENTS>;

#if defined(TIMEMORY_MPI_GOTCHA)
//
namespace tim
{
namespace component
{
//
struct mpi_comm_data : base<mpi_comm_data, void>
{
    using value_type = void;
    using this_type  = mpi_comm_data;
    using base_type  = base<this_type, value_type>;
    using tracker_t  = tim::auto_tuple<mpi_data_tracker_t>;
    using data_type  = float;

    TIMEMORY_DEFAULT_OBJECT(mpi_comm_data)

    static void preinit()
    {
        mpi_data_tracker_t::label()       = "mpi_comm_data";
        mpi_data_tracker_t::description() = "Tracks MPI communication data";
    }

    // MPI_Send
    void audit(const std::string& _name, const void*, int count, MPI_Datatype datatype,
               int dst, int, MPI_Comm)
    {
        int size = 0;
        MPI_Type_size(datatype, &size);
	if(count * size == 0)
	    return;
	PRINT_HERE("Logging %s", _name.c_str());
        tracker_t _t(_name);
        add(_t, count * size);
        add_secondary(_t, TIMEMORY_JOIN("_", _name, "dst", dst), count * size);
    }

    // MPI_Recv
    void audit(const std::string& _name, void*, int count, MPI_Datatype datatype, int dst,
               int, MPI_Comm, MPI_Status*)
    {
        int size = 0;
        MPI_Type_size(datatype, &size);
	if(count * size == 0)
	    return;
	PRINT_HERE("Logging %s", _name.c_str());
        tracker_t _t(_name);
        add(_t, count * size);
        add_secondary(_t, TIMEMORY_JOIN("_", _name, "dst", dst), count * size);
    }

    // MPI_Bcast
    void audit(const std::string& _name, void*, int count, MPI_Datatype datatype,
               int root, MPI_Comm)
    {
        int size = 0;
        MPI_Type_size(datatype, &size);
	if(count * size == 0)
	    return;
	PRINT_HERE("Logging %s", _name.c_str());
        add(_name, count * size, TIMEMORY_JOIN("_", _name, "root", root));
    }

    // MPI_Allreduce
    void audit(const std::string& _name, const void*, void*, int count,
               MPI_Datatype datatype, MPI_Op, MPI_Comm)
    {
        int size = 0;
        MPI_Type_size(datatype, &size);
	if(count * size == 0)
	    return;
	PRINT_HERE("Logging %s", _name.c_str());
        add(_name, count * size);
    }

    // MPI_Sendrecv
    void audit(const std::string& _name, const void*, int sendcount,
               MPI_Datatype sendtype, int, int sendtag, void*, int recvcount,
               MPI_Datatype recvtype, int, int recvtag, MPI_Comm, MPI_Status*)
    {
        int send_size = 0;
        int recv_size = 0;
        MPI_Type_size(sendtype, &send_size);
        MPI_Type_size(recvtype, &recv_size);
	if(sendcount * send_size == 0 || recvcount * recv_size == 0)
	    return;
	PRINT_HERE("Logging %s", _name.c_str());
        tracker_t _t(_name);
        add(_t, sendcount * send_size + recvcount * recv_size);
        add_secondary(_t, TIMEMORY_JOIN("_", _name, "send"), sendcount * send_size,
                      TIMEMORY_JOIN("_", _name, "send", "tag", sendtag));
        add_secondary(_t, TIMEMORY_JOIN("_", _name, "recv"), recvcount * recv_size,
                      TIMEMORY_JOIN("_", _name, "recv", "tag", recvtag));
    }

    // MPI_Gather
    void audit(const std::string& _name, const void*, int sendcount,
               MPI_Datatype sendtype, void*, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm)
    {
        int send_size = 0;
        int recv_size = 0;
        MPI_Type_size(sendtype, &send_size);
        MPI_Type_size(recvtype, &recv_size);
	if(sendcount * send_size == 0 || recvcount * recv_size == 0)
	    return;
	PRINT_HERE("Logging %s", _name.c_str());
        tracker_t _t(_name);
        add(_t, sendcount * send_size + recvcount * recv_size);
        tracker_t _r(TIMEMORY_JOIN("_", _name, "root", root));
        add(_r, sendcount * send_size + recvcount * recv_size);
        add_secondary(_r, TIMEMORY_JOIN("_", _name, "root", root, "send"),
                      sendcount * send_size);
        add_secondary(_r, TIMEMORY_JOIN("_", _name, "root", root, "recv"),
                      recvcount * recv_size);
    }

    // MPI_Scatter
    void audit(const std::string& _name, void*, int sendcount, MPI_Datatype sendtype,
               void*, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm)
    {
        int send_size = 0;
        int recv_size = 0;
        MPI_Type_size(sendtype, &send_size);
        MPI_Type_size(recvtype, &recv_size);
	if(sendcount * send_size == 0 || recvcount * recv_size == 0)
	    return;
	PRINT_HERE("Logging %s", _name.c_str());
        tracker_t _t(_name);
        add(_t, sendcount * send_size + recvcount * recv_size);
        tracker_t _r(TIMEMORY_JOIN("_", _name, "root", root));
        add(_r, sendcount * send_size + recvcount * recv_size);
        add_secondary(_r, TIMEMORY_JOIN("_", _name, "root", root, "send"),
                      sendcount * send_size);
        add_secondary(_r, TIMEMORY_JOIN("_", _name, "root", root, "recv"),
                      recvcount * recv_size);
    }

    // MPI_Alltoall
    void audit(const std::string& _name, void*, int sendcount, MPI_Datatype sendtype,
               void*, int recvcount, MPI_Datatype recvtype, MPI_Comm)
    {
        int send_size = 0;
        int recv_size = 0;
        MPI_Type_size(sendtype, &send_size);
        MPI_Type_size(recvtype, &recv_size);
	if(sendcount * send_size == 0 || recvcount * recv_size == 0)
	    return;
	PRINT_HERE("Logging %s", _name.c_str());
        tracker_t _t(_name);
        add(_t, sendcount * send_size + recvcount * recv_size);
        add_secondary(_t, TIMEMORY_JOIN("_", _name, "send"), sendcount * send_size);
        add_secondary(_t, TIMEMORY_JOIN("_", _name, "recv"), recvcount * recv_size);
    }

private:
    template <typename... Args>
    void add(tracker_t& _t, data_type value, Args&&... args)
    {
	if(value <= 0.0)
	    return;
        _t.store(std::plus<data_type>{}, value);
        TIMEMORY_FOLD_EXPRESSION(add_secondary(_t, std::forward<Args>(args), value));
    }

    template <typename... Args>
    void add(const std::string& _name, data_type value, Args&&... args)
    {
	if(value <= 0.0)
	    return;
        tracker_t _t(_name);
        add(_t, value, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void add_secondary(tracker_t&, const std::string& _name, data_type value,
                       Args&&... args)
    {
        if(tim::settings::add_secondary() && value > 0.0)
        {
            tracker_t _s(_name);
            add(_s, value, std::forward<Args>(args)...);
        }
    }
};
}  // namespace component
}  // namespace tim
//
TIMEMORY_STORAGE_INITIALIZER(mpi_comm_data, mpi_comm_data)
TIMEMORY_STORAGE_INITIALIZER(mpi_data_tracker_t, mpi_data_tracker_t)
//
#endif
