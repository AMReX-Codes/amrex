
#ifndef BL_TIMEMORY_H_
#define BL_TIMEMORY_H_

#include <timemory/hash.hpp>
#include <timemory/components/timing/wall_clock.hpp>
#include <timemory/components/user_bundle.hpp>
#include <timemory/general/source_location.hpp>
#include <timemory/variadic/auto_tuple.hpp>
#include <timemory/variadic/component_tuple.hpp>

//
//  These declarations are used by multiple profilers
//
namespace amrex
{
// this is used to make the user_bundle specific to AMReX
struct BL_timemory_tag
{};
}  // namespace amrex

TIMEMORY_BUNDLE_INDEX(BL_bundle_idx, 0)
TIMEMORY_COMPONENT_ALIAS(BL_timemory_bundle,
                         user_bundle<BL_bundle_idx, amrex::BL_timemory_tag>)

// This just makes the extern template declarations easier
#if !defined(AMREX_TIMEMORY_COMPONENTS)
#define AMREX_TIMEMORY_COMPONENTS tim::component::BL_timemory_bundle
#endif

extern template class tim::component_tuple<AMREX_TIMEMORY_COMPONENTS>;
extern template class tim::auto_tuple<AMREX_TIMEMORY_COMPONENTS>;

namespace amrex
{
static constexpr auto BL_timemory_bundle_idx = tim::component::BL_bundle_idx;

/// this is a generic wrapper for components which can be selected at runtime.
/// the wrapper does not support passing in arguments to the components.
/// if you want a component which can receive input from the application
/// and store/track something, list it explicitly in \ref AMREX_TIMEMORY_COMPONENTS
/// and create an overload in TimemoryProfiler to accept the arguments
using BL_timemory_bundle = tim::component::BL_timemory_bundle;

/// this handle requires start and stop to be explicitly called
using BL_timemory_handle = tim::component_tuple<AMREX_TIMEMORY_COMPONENTS>;

/// this handle scopes the start and stop via RAII
using BL_timemory_auto_handle = tim::auto_tuple<AMREX_TIMEMORY_COMPONENTS>;

/// \fn void BL_timemory_configure(int argc, char** argv)
/// \param argc number of arguments
/// \param argv arguments
///
/// \brief This function support command line arguments processing
void
BL_timemory_configure (int argc, char** argv);

void
BL_timemory_initialize (int argc = 0, char** argv = nullptr, bool parse_args = true);

void
BL_timemory_finalize ();

//! A simple profiler that returns basic performance information (e.g. min, max, and
//! average running time)
class TimemoryProfiler
{
public:
    using value_type   = BL_timemory_handle;
    using region_map_t = std::unordered_map<std::string, std::unique_ptr<value_type>>;

public:
    explicit TimemoryProfiler (std::string&& funcname, bool _start = true);
    explicit TimemoryProfiler (const char* funcname, bool _start = true);
    explicit TimemoryProfiler (size_t _hash, bool _start = true);
    explicit TimemoryProfiler (const tim::source_location::captured& _hash,
                               bool                                  _start = true);
    ~TimemoryProfiler ();

    TimemoryProfiler (const TimemoryProfiler&)     = delete;
    TimemoryProfiler (TimemoryProfiler&&) noexcept = default;
    TimemoryProfiler& operator= (const TimemoryProfiler&) = delete;
    TimemoryProfiler& operator= (TimemoryProfiler&&) noexcept = default;

    void start ();
    void stop ();

private:
    BL_timemory_handle m_handle;
};

class TimemoryTplProfiler
{
public:
    using value_type = tim::lightweight_tuple<
        tim::component::allinea_map, tim::component::craypat_record,
        tim::component::cuda_profiler, tim::component::vtune_profiler,
        tim::component::gperftools_cpu_profiler,
        tim::component::gperftools_heap_profiler>;

    static void        start ();
    static void        stop ();
    static value_type* get ();
};

using BLProfiler_timemory_region_map_t =
    std::unordered_map<std::string, std::unique_ptr<TimemoryProfiler>>;

BLProfiler_timemory_region_map_t&
BL_timemory_get_regions ();

}  // namespace amrex

#endif
