// We only support BL_PROFILE, BL_PROFILE_VAR, BL_PROFILE_VAR_STOP, BL_PROFILE_VAR_START,
// BL_PROFILE_VAR_NS, and BL_PROFILE_REGION.

#include <AMReX_TinyProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#ifdef AMREX_USE_GPU
#include <AMReX_GpuDevice.H>
#endif
#include <AMReX_Print.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#ifdef AMREX_USE_CUDA
#if __has_include(<nvtx3/nvToolsExt.h>)
#  include <nvtx3/nvToolsExt.h>
#else
#  include <nvToolsExt.h>
#endif
#endif

#if defined(AMREX_USE_HIP) && defined(AMREX_USE_ROCTX)
#include <roctracer/roctx.h>
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <set>

namespace amrex {

std::deque<const TinyProfiler*> TinyProfiler::mem_stack;
#ifdef AMREX_USE_OMP
std::vector<TinyProfiler::aligned_deque> TinyProfiler::mem_stack_thread_private;
#endif
std::vector<std::map<std::string, MemStat>*> TinyProfiler::all_memstats;
std::vector<std::string> TinyProfiler::all_memnames;

std::vector<std::string>          TinyProfiler::regionstack;
std::deque<std::tuple<double,double,std::string*> > TinyProfiler::ttstack;
std::map<std::string,std::map<std::string, TinyProfiler::Stats> > TinyProfiler::statsmap;
double TinyProfiler::t_init = std::numeric_limits<double>::max();
bool TinyProfiler::device_synchronize_around_region = false;
int TinyProfiler::n_print_tabs = 0;
int TinyProfiler::verbose = 0;
double TinyProfiler::print_threshold = 1.;
bool TinyProfiler::enabled = true;
bool TinyProfiler::memprof_enabled = true;
std::string TinyProfiler::output_file;

namespace {
    constexpr char mainregion[] = "main";
}

TinyProfiler::TinyProfiler (std::string funcname) noexcept
    : fname(std::move(funcname))
{
    start();
}

TinyProfiler::TinyProfiler (std::string funcname, bool start_) noexcept
    : fname(std::move(funcname))
{
    if (start_) { start(); }
}

TinyProfiler::TinyProfiler (const char* funcname) noexcept
    : fname(funcname)
{
    start();
}

TinyProfiler::TinyProfiler (const char* funcname, bool start_) noexcept
    : fname(funcname)
{
    if (start_) { start(); }
}

TinyProfiler::~TinyProfiler ()
{
    stop();
}

void
TinyProfiler::start () noexcept
{
    if (!enabled) { return; }

    memory_start();

#ifdef AMREX_USE_OMP
#pragma omp master
#endif
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(stats.empty(), "TinyProfiler cannot be started twice");
    }

#ifdef AMREX_USE_OMP
#pragma omp master
#endif
    if (!regionstack.empty()) {

#ifdef AMREX_USE_GPU
        if (device_synchronize_around_region) {
            amrex::Gpu::streamSynchronize();
        }
#endif

        const double t = amrex::second();

        ttstack.emplace_back(t, 0.0, &fname);
        global_depth = static_cast<int>(ttstack.size());
#ifdef AMREX_USE_OMP
        in_parallel_region = omp_in_parallel();
#else
        in_parallel_region = false;
#endif

#ifdef AMREX_USE_CUDA
        nvtxRangePush(fname.c_str());
#elif defined(AMREX_USE_HIP) && defined(AMREX_USE_ROCTX)
        roctxRangePush(fname.c_str());
#endif

        for (auto const& region : regionstack)
        {
            Stats& st = statsmap[region][fname];
            ++st.depth;
            stats.push_back(&st);
        }

        if (verbose) {
            ++n_print_tabs;
            std::string whitespace;
            for (int itab = 0; itab < n_print_tabs; ++itab) {
                whitespace += "  ";
            }
            // If we try to print to output_file here, it may not be thread
            // safe. Also note that this is controlled by verbose already.
            amrex::Print() << whitespace << "TP: Entering " << fname << '\n';
        }
    }
}

void
TinyProfiler::stop () noexcept
{
    if (!enabled) { return; }

    memory_stop();

#ifdef AMREX_USE_OMP
#pragma omp master
#endif
    if (!stats.empty()) {

#ifdef AMREX_USE_GPU
        if (device_synchronize_around_region) {
            amrex::Gpu::streamSynchronize();
        }
#endif

        const double t = amrex::second();

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<int>(ttstack.size()) == global_depth,
            "TinyProfiler sections must be nested with respect to each other");
#ifdef AMREX_USE_OMP
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(in_parallel_region == omp_in_parallel(),
            "TinyProfiler sections must be nested with respect to parallel regions");
#endif

        {
            const std::tuple<double,double,std::string*>& tt = ttstack.back();

            // first: wall time when the pair is pushed into the stack
            // second: accumulated dt of children
            double dtin = t - std::get<0>(tt); // elapsed time since start() is called.
            double dtex = dtin - std::get<1>(tt);

            for (Stats* st : stats)
            {
                --(st->depth);
                ++(st->n);
                if (st->depth == 0) {
                    st->dtin += dtin;
                }
                st->dtex += dtex;
            }

            ttstack.pop_back();
            if (!ttstack.empty()) {
                std::tuple<double,double,std::string*>& parent = ttstack.back();
                std::get<1>(parent) += dtin;
            }

#ifdef AMREX_USE_CUDA
            nvtxRangePop();
#elif defined(AMREX_USE_HIP) && defined(AMREX_USE_ROCTX)
            roctxRangePop();
#endif
        }

        stats.clear();

        if (verbose) {
            std::string whitespace;
            for (int itab = 0; itab < n_print_tabs; ++itab) {
                whitespace += "  ";
            }
            --n_print_tabs;
            // If we try to print to output_file here, it may not be thread
            // safe. Also note that this is controlled by verbose already.
            amrex::Print() << whitespace << "TP: Leaving  " << fname << '\n';
        }
    }
}

void
TinyProfiler::memory_start () const noexcept
{
    if (!memprof_enabled) { return; }

    // multiple omp threads may share the same TinyProfiler object so this function must be const
    // it is NOT allowed to double start a section
#ifdef AMREX_USE_OMP
    if (omp_in_parallel()) {
        mem_stack_thread_private[omp_get_thread_num()].deque.push_back(this);
    } else
#endif
    {
        mem_stack.push_back(this);
    }
}

void
TinyProfiler::memory_stop () const noexcept
{
    if (!memprof_enabled) { return; }

    // multiple omp threads may share the same TinyProfiler object so this function must be const
    // it IS allowed to double stop a section
#ifdef AMREX_USE_OMP
    if (omp_in_parallel()) {
        if (!mem_stack_thread_private[omp_get_thread_num()].deque.empty() &&
            mem_stack_thread_private[omp_get_thread_num()].deque.back() == this) {
            mem_stack_thread_private[omp_get_thread_num()].deque.pop_back();
        }
    } else
#endif
    {
        if (!mem_stack.empty() && mem_stack.back() == this) {
            mem_stack.pop_back();
        }
    }
}

MemStat*
TinyProfiler::memory_alloc (std::size_t nbytes, std::map<std::string, MemStat>& memstats) noexcept
{
    if (!memprof_enabled) { return nullptr; }

    // this function is not thread safe for the same memstats
    // the caller of this function (CArena::alloc) has a mutex
    MemStat* stat = nullptr;
#ifdef AMREX_USE_OMP
    if (omp_in_parallel() && !mem_stack_thread_private[omp_get_thread_num()].deque.empty()) {
        stat = &memstats[
            mem_stack_thread_private[omp_get_thread_num()].deque.back()->fname
        ];
    } else
#endif
    if (!mem_stack.empty()) {
        stat = &memstats[mem_stack.back()->fname];
    } else {
        stat = &memstats["Unprofiled"];
    }

    ++stat->nalloc;
    stat->currentmem += static_cast<Long>(nbytes);
    stat->avgmem -= static_cast<double>(nbytes) * amrex::second();
    stat->maxmem = std::max(stat->maxmem, stat->currentmem);

    return stat;
}

void
TinyProfiler::memory_free (std::size_t nbytes, MemStat* stat) noexcept
{
    if (!memprof_enabled) { return; }

    // this function is not thread safe for the same stat
    // the caller of this function (CArena::free) has a mutex
    if (stat) {
        ++stat->nfree;
        stat->avgmem += static_cast<double>(nbytes) * amrex::second();
        stat->currentmem -= static_cast<Long>(nbytes);
    }
}


void
TinyProfiler::Initialize () noexcept
{
    {
        amrex::ParmParse pp("tiny_profiler");
        pp.queryAdd("device_synchronize_around_region", device_synchronize_around_region);
        if (! pp.query("verbose", "v", verbose)) {
            pp.add("verbose", verbose);
        }
        // Specify the maximum percentage of inclusive time
        // that the "Other" section in the output can have (default 1%)
        pp.queryAdd("print_threshold", print_threshold);

        pp.queryAdd("enabled", enabled);
        pp.queryAdd("output_file", output_file);
    }

    if (!enabled) { return; }

    if (ParallelDescriptor::IOProcessor()) {
        static bool first = true;
        if (first && !output_file.empty() && output_file != "/dev/null") {
            if (FileSystem::Exists(output_file)) {
                FileSystem::Remove(output_file);
            }
            first = false;
        }
    }

    regionstack.emplace_back(mainregion);
    t_init = amrex::second();
}

void
TinyProfiler::MemoryInitialize () noexcept
{
    {
        amrex::ParmParse pp("tiny_profiler");
        pp.queryAdd("enabled", enabled);
        pp.queryAdd("memprof_enabled", memprof_enabled);
        memprof_enabled = memprof_enabled && enabled;
    }

    if (!memprof_enabled) { return; }

#ifdef AMREX_USE_OMP
    mem_stack_thread_private.resize(omp_get_max_threads());
#endif
}

void
TinyProfiler::Finalize (bool bFlushing) noexcept
{
    if (!enabled) { return; }

    static bool finalized = false;
    if (!bFlushing) {                // If flushing, don't make this the last time!
        if (finalized) {
            return;
        } else {
            finalized = true;
        }
    }

    double t_final = amrex::second();

    // make a local copy so that any functions call after this will not be recorded in the local copy.
    auto lstatsmap = statsmap;

    int nprocs = ParallelDescriptor::NProcs();
    int ioproc = ParallelDescriptor::IOProcessorNumber();

    double dt_max = t_final - t_init;
    ParallelReduce::Max(dt_max, ioproc, ParallelDescriptor::Communicator());
    double dt_min = t_final - t_init;
    ParallelReduce::Min(dt_min, ioproc, ParallelDescriptor::Communicator());
    double dt_avg = t_final - t_init;
    ParallelReduce::Sum(dt_avg, ioproc, ParallelDescriptor::Communicator());
    dt_avg /= double(nprocs);

    std::ofstream ofs;
    std::ostream* os = nullptr;
    std::streamsize oldprec = 0;
    if (ParallelDescriptor::IOProcessor()) {
        if (output_file.empty()) {
            os = &(amrex::OutStream());
        } else if (output_file != "/dev/null") {
            ofs.open(output_file, std::ios_base::app);
            if (!ofs.is_open()) {
                amrex::Error("TinyProfiler failed to open "+output_file);
            }
            os = static_cast<std::ostream*>(&ofs);
        }
    }

    if (os)
    {
        os->precision(4);
        *os << "\n\nTinyProfiler total time across processes [min...avg...max]: "
            << dt_min << " ... " << dt_avg << " ... " << dt_max << "\n";
    }

    // make sure the set of regions is the same on all processes.
    {
        Vector<std::string> localRegions, syncedRegions;
        bool alreadySynced;

        for (auto const& kv : lstatsmap) {
            localRegions.push_back(kv.first);
        }

        amrex::SyncStrings(localRegions, syncedRegions, alreadySynced);

        if (!alreadySynced) {
            for (auto const& s : syncedRegions) {
                if (lstatsmap.find(s) == lstatsmap.end()) {
                    lstatsmap.insert(std::make_pair(s,std::map<std::string,Stats>()));
                }
            }
        }
    }

    PrintStats(lstatsmap[mainregion], dt_max, os);
    for (auto& kv : lstatsmap) {
        if (kv.first != mainregion) {
            if (os) {
                *os << "\n\nBEGIN REGION " << kv.first << "\n";
            }
            PrintStats(kv.second, dt_max, os);
            if (os) {
                *os << "END REGION " << kv.first << "\n";
            }
        }
    }

    if(os) { os->precision(oldprec); }
}

void
TinyProfiler::MemoryFinalize (bool bFlushing) noexcept
{
    if (!memprof_enabled) { return; }

    // This function must be called BEFORE the profiled arenas are deleted

    static bool finalized = false;
    if (!bFlushing) {                // If flushing, don't make this the last time!
        if (finalized) {
            return;
        } else {
            finalized = true;
        }
    }

    double t_final = amrex::second();
    double dt_max = t_final - t_init;
    int ioproc = ParallelDescriptor::IOProcessorNumber();
    ParallelReduce::Max(dt_max, ioproc, ParallelDescriptor::Communicator());

    std::ofstream ofs;
    std::ostream* os = nullptr;
    std::streamsize oldprec = 0;
    if (ParallelDescriptor::IOProcessor()) {
        if (output_file.empty()) {
            os = &(amrex::OutStream());
        } else if (output_file != "/dev/null") {
            ofs.open(output_file, std::ios_base::app);
            if (!ofs.is_open()) {
                amrex::Error("TinyProfiler failed to open "+output_file);
            }
            os = static_cast<std::ostream*>(&ofs);
        }
    }

    for (std::size_t i = 0; i < all_memstats.size(); ++i) {
        PrintMemStats(*(all_memstats[i]), all_memnames[i], dt_max, t_final, os);
    }

    if (!bFlushing) {
        all_memstats.clear();
        all_memnames.clear();
    }

    if(os) { os->precision(oldprec); }
}

void
TinyProfiler::RegisterArena (const std::string& memory_name,
                             std::map<std::string, MemStat>& memstats) noexcept
{
    if (!memprof_enabled) { return; }

    all_memstats.push_back(&memstats);
    all_memnames.push_back(memory_name);
}

void
TinyProfiler::DeregisterArena (std::map<std::string, MemStat>& memstats) noexcept
{
    if (!memprof_enabled) { return; }

    for (std::size_t i = 0; i < all_memstats.size();) {
        if (all_memstats[i] == &memstats) {
            all_memstats.erase(all_memstats.begin() + i); // NOLINT
            all_memnames.erase(all_memnames.begin() + i); // NOLINT
        } else {
            ++i;
        }
    }
}

void
TinyProfiler::PrintStats (std::map<std::string,Stats>& regstats, double dt_max,
                          std::ostream* os)
{
    // make sure the set of profiled functions is the same on all processes
    {
        Vector<std::string> localStrings, syncedStrings;
        bool alreadySynced;

        for(auto const& kv : regstats) {
            localStrings.push_back(kv.first);
        }

        amrex::SyncStrings(localStrings, syncedStrings, alreadySynced);

        if (! alreadySynced) {  // add the new name
            for (auto const& s : syncedStrings) {
                if (regstats.find(s) == regstats.end()) {
                    regstats.insert(std::make_pair(s, Stats()));
                }
            }
        }
    }

    if (regstats.empty()) { return; }

    int nprocs = ParallelDescriptor::NProcs();
    int ioproc = ParallelDescriptor::IOProcessorNumber();

    std::vector<ProcStats> allprocstats;
    int maxfnamelen = 0;
    Long maxncalls = 0;

    // now collect global data onto the ioproc
    for (const auto & regstat : regstats)
    {
        Long n = regstat.second.n;
        double dts[2] = {regstat.second.dtin, regstat.second.dtex};

        std::vector<Long> ncalls(nprocs);
        std::vector<double> dtdt(2*nprocs);

        if (ParallelDescriptor::NProcs() == 1)
        {
            ncalls[0] = n;
            dtdt[0] = dts[0];
            dtdt[1] = dts[1];
        } else
        {
            ParallelDescriptor::Gather(&n, 1, ncalls.data(), 1, ioproc);
            ParallelDescriptor::Gather(dts, 2, dtdt.data(), 2, ioproc);
        }

        if (ParallelDescriptor::IOProcessor()) {
            ProcStats pst;
            for (int i = 0; i < nprocs; ++i) {
                pst.nmin  = std::min(pst.nmin, ncalls[i]);
                pst.navg +=                    ncalls[i];
                pst.nmax  = std::max(pst.nmax, ncalls[i]);
                pst.dtinmin  = std::min(pst.dtinmin, dtdt[2*i]);
                pst.dtinavg +=                       dtdt[2*i];
                pst.dtinmax  = std::max(pst.dtinmax, dtdt[2*i]);
                pst.dtexmin  = std::min(pst.dtexmin, dtdt[2*i+1]);
                pst.dtexavg +=                       dtdt[2*i+1];
                pst.dtexmax  = std::max(pst.dtexmax, dtdt[2*i+1]);
            }
            pst.navg /= nprocs;
            pst.dtinavg /= nprocs;
            pst.dtexavg /= nprocs;
            pst.fname = regstat.first;
            allprocstats.push_back(pst);
            maxfnamelen = std::max(maxfnamelen, int(pst.fname.size()));
            maxncalls = std::max(maxncalls, pst.nmax);
        }
    }

    if (ParallelDescriptor::IOProcessor() && os)
    {
        *os << std::setfill(' ') << std::setprecision(4);
        int wt = 9;

        int wnc = (int) std::log10 ((double) maxncalls) + 1;
        wnc = std::max(wnc, int(std::string("NCalls").size()));
        wt  = std::max(wt,  int(std::string("Excl. Min").size()));
        int wp = 6;
        wp  = std::max(wp,  int(std::string("Max %").size()));

        const std::string hline(maxfnamelen+wnc+2+(wt+2)*3+wp+2,'-');

        ProcStats other_procstat;
        bool print_other_procstat = false;

        // try to combine low-performance impact functions into "Other" to clean up the output
        if (print_threshold > 0.) {
            // initialize other_procstat to zero
            other_procstat.nmin = 0;
            other_procstat.dtinmin = 0.;
            other_procstat.dtexmin = 0.;
            other_procstat.fname = "Other";
            int num_procstats_in_other = 0;

            // sort by exclusive time and iterate backwards over the profiled functions
            std::sort(allprocstats.begin(), allprocstats.end(), ProcStats::compin);
            for (Long i = static_cast<Long>(allprocstats.size())-1; i >= 0; --i) {
                // include function in "Other" if together they are below the threshold
                if ((other_procstat.dtinmax + allprocstats[i].dtinmax)*(100.0/dt_max)
                        < print_threshold) {
                    allprocstats[i].do_print = false;
                    ++num_procstats_in_other;

                    // add time for function to "Other"
                    // for min and max this is not exact but produces an upper limit
                    other_procstat.nmin += allprocstats[i].nmin;
                    other_procstat.navg += allprocstats[i].navg;
                    other_procstat.nmax += allprocstats[i].nmax;

                    other_procstat.dtinmin += allprocstats[i].dtinmin;
                    other_procstat.dtinavg += allprocstats[i].dtinavg;
                    other_procstat.dtinmax += allprocstats[i].dtinmax;

                    other_procstat.dtexmin += allprocstats[i].dtexmin;
                    other_procstat.dtexavg += allprocstats[i].dtexavg;
                    other_procstat.dtexmax += allprocstats[i].dtexmax;
                } else {
                    break;
                }
            }

            if (num_procstats_in_other == 1) {
                // if only one function would be included in "Other"
                // the output would not get shorter
                allprocstats.back().do_print = true;
            } else if (num_procstats_in_other >= 2) {
                print_other_procstat = true;
            }
        }

        // Exclusive time
        std::sort(allprocstats.begin(), allprocstats.end(), ProcStats::compex);
        if (print_other_procstat) {
            // make sure "Other" is printed at the end of the list
            allprocstats.push_back(other_procstat);
        }
        *os << "\n" << hline << "\n";
        *os << std::left
            << std::setw(maxfnamelen) << "Name"
            << std::right
            << std::setw(wnc+2) << "NCalls"
            << std::setw(wt+2) << "Excl. Min"
            << std::setw(wt+2) << "Excl. Avg"
            << std::setw(wt+2) << "Excl. Max"
            << std::setw(wp+2)  << "Max %"
            << "\n" << hline << "\n";
        for (const auto & allprocstat : allprocstats)
        {
            if (!allprocstat.do_print) {
                continue;
            }
            *os << std::setprecision(4) << std::left
                << std::setw(maxfnamelen) << allprocstat.fname
                << std::right
                << std::setw(wnc+2) << allprocstat.navg
                << std::setw(wt+2) << allprocstat.dtexmin
                << std::setw(wt+2) << allprocstat.dtexavg
                << std::setw(wt+2) << allprocstat.dtexmax
                << std::setprecision(2) << std::setw(wp+1) << std::fixed
                << allprocstat.dtexmax*(100.0/dt_max) << "%";
            os->unsetf(std::ios_base::fixed);
            *os << "\n";
        }
        *os << hline << "\n";
        if (print_other_procstat) {
            allprocstats.pop_back();
        }

        // Inclusive time
        std::sort(allprocstats.begin(), allprocstats.end(), ProcStats::compin);
        if (print_other_procstat) {
            // make sure "Other" is printed at the end of the list
            allprocstats.push_back(other_procstat);
        }
        *os << "\n" << hline << "\n";
        *os << std::left
            << std::setw(maxfnamelen) << "Name"
            << std::right
            << std::setw(wnc+2) << "NCalls"
            << std::setw(wt+2) << "Incl. Min"
            << std::setw(wt+2) << "Incl. Avg"
            << std::setw(wt+2) << "Incl. Max"
            << std::setw(wp+2)  << "Max %"
            << "\n" << hline << "\n";
        for (const auto & allprocstat : allprocstats)
        {
            if (!allprocstat.do_print) {
                continue;
            }
            *os << std::setprecision(4) << std::left
                << std::setw(maxfnamelen) << allprocstat.fname
                << std::right
                << std::setw(wnc+2) << allprocstat.navg
                << std::setw(wt+2) << allprocstat.dtinmin
                << std::setw(wt+2) << allprocstat.dtinavg
                << std::setw(wt+2) << allprocstat.dtinmax
                << std::setprecision(2) << std::setw(wp+1) << std::fixed
                << allprocstat.dtinmax*(100.0/dt_max) << "%";
            os->unsetf(std::ios_base::fixed);
            *os << "\n";
        }
        *os << hline << "\n\n";
    }
}

void
TinyProfiler::PrintMemStats (std::map<std::string, MemStat>& memstats,
                             std::string const& memname, double dt_max,
                             double t_final, std::ostream* os)
{
    // make sure the set of profiled functions is the same on all processes
    {
        Vector<std::string> localStrings, syncedStrings;
        bool alreadySynced;

        for(auto const& kv : memstats) {
            localStrings.push_back(kv.first);
        }

        amrex::SyncStrings(localStrings, syncedStrings, alreadySynced);

        if (! alreadySynced) {  // add the new name
            for (auto const& s : syncedStrings) {
                if (memstats.find(s) == memstats.end()) {
                    memstats[s]; // insert
                }
            }
        }
    }

    if (memstats.empty()) { return; }

    const int nprocs = ParallelDescriptor::NProcs();
    const int ioproc = ParallelDescriptor::IOProcessorNumber();

    std::vector<MemProcStats> allprocstats;

    // now collect global data onto the ioproc
    for (const auto & it : memstats)
    {
        Long nalloc = it.second.nalloc;
        Long nfree = it.second.nfree;
        // simulate the freeing of remaining memory currentmem for the avgmem metric
        Long avgmem = static_cast<Long>(
            (it.second.avgmem + static_cast<double>(it.second.currentmem) * t_final) / dt_max);
        Long maxmem = it.second.maxmem;

        std::vector<Long> nalloc_vec(nprocs);
        std::vector<Long> nfree_vec(nprocs);
        std::vector<Long> avgmem_vec(nprocs);
        std::vector<Long> maxmem_vec(nprocs);

        if (nprocs == 1)
        {
            nalloc_vec[0] = nalloc;
            nfree_vec[0] = nfree;
            avgmem_vec[0] = avgmem;
            maxmem_vec[0] = maxmem;
        } else
        {
            ParallelDescriptor::Gather(&nalloc, 1, nalloc_vec.data(), 1, ioproc);
            ParallelDescriptor::Gather(&nfree , 1,  nfree_vec.data(), 1, ioproc);
            ParallelDescriptor::Gather(&maxmem, 1, maxmem_vec.data(), 1, ioproc);
            ParallelDescriptor::Gather(&avgmem, 1, avgmem_vec.data(), 1, ioproc);
        }

        if (ParallelDescriptor::IOProcessor()) {
            MemProcStats pst;
            for (int i = 0; i < nprocs; ++i) {

                pst.nalloc += nalloc_vec[i];
                pst.nfree += nfree_vec[i];
                pst.avgmem_min = std::min(pst.avgmem_min, avgmem_vec[i]);
                pst.avgmem_avg += avgmem_vec[i];
                pst.avgmem_max = std::max(pst.avgmem_max, avgmem_vec[i]);
                pst.maxmem_min = std::min(pst.maxmem_min, maxmem_vec[i]);
                pst.maxmem_avg += maxmem_vec[i];
                pst.maxmem_max = std::max(pst.maxmem_max, maxmem_vec[i]);
            }
            pst.avgmem_avg /= nprocs;
            pst.maxmem_avg /= nprocs;
            pst.fname = it.first;
            allprocstats.push_back(pst);
        }
    }

    std::sort(allprocstats.begin(), allprocstats.end(), MemProcStats::compmem);

    std::vector<std::vector<std::string>> allstatsstr;

    if (nprocs == 1) {
        allstatsstr.push_back({"Name", "Nalloc", "Nfree", "AvgMem", "MaxMem"});
    } else {
        allstatsstr.push_back({"Name", "Nalloc", "Nfree",
                               "AvgMem min", "AvgMem avg", "AvgMem max",
                               "MaxMem min", "MaxMem avg", "MaxMem max"});
    }

    auto mem_to_string = [] (Long nbytes) {
        std::string unit = "   B";
        if (nbytes >= 10000) {
            nbytes /= 1024;
            unit = " KiB";
        }
        if (nbytes >= 10000) {
            nbytes /= 1024;
            unit = " MiB";
        }
        if (nbytes >= 10000) {
            nbytes /= 1024;
            unit = " GiB";
        }
        if (nbytes >= 10000) {
            nbytes /= 1024;
            unit = " TiB";
        }
        return std::to_string(nbytes) + unit;
    };

    for (auto& stat : allprocstats) {
        if (stat.nalloc != 0 || stat.nfree != 0 || stat.maxmem_max != 0) {
            if (nprocs == 1) {
                allstatsstr.push_back({stat.fname,
                                    std::to_string(stat.nalloc),
                                    std::to_string(stat.nfree),
                                    mem_to_string(stat.avgmem_max),
                                    mem_to_string(stat.maxmem_max)});
            } else {
                allstatsstr.push_back({stat.fname,
                                    std::to_string(stat.nalloc),
                                    std::to_string(stat.nfree),
                                    mem_to_string(stat.avgmem_min),
                                    mem_to_string(stat.avgmem_avg),
                                    mem_to_string(stat.avgmem_max),
                                    mem_to_string(stat.maxmem_min),
                                    mem_to_string(stat.maxmem_avg),
                                    mem_to_string(stat.maxmem_max)});
            }
        }
    }

    std::vector<int> maxlen(allstatsstr[0].size(), 0);
    for (auto& strvec : allstatsstr) {
        for (std::size_t i=0; i<maxlen.size(); ++i) {
            maxlen[i] = std::max(maxlen[i], static_cast<int>(strvec[i].size()));
        }
    }

    for (std::size_t i=1; i<maxlen.size(); ++i) {
        maxlen[i] += 2;
    }

    if (allstatsstr.size() == 1 || !os) { return; }

    int lenhline = 0;
    for (auto i : maxlen) {
        lenhline += i;
    }
    const std::string hline(lenhline, '-');

    *os << memname << " Usage:\n";
    *os << hline << "\n";
    for (std::size_t i=0; i<allstatsstr.size(); ++i) {
        *os << std::left << std::setw(maxlen[0]) << allstatsstr[i][0];
        for (std::size_t j=1; j<maxlen.size(); ++j) {
            *os << std::right << std::setw(maxlen[j]) << allstatsstr[i][j];
        }
        *os << '\n';
        if (i==0) {
            *os << hline << "\n";
        }
    }
    *os << hline << "\n\n";
}

void
TinyProfiler::StartRegion (std::string regname) noexcept
{
    if (!enabled) { return; }

    if (std::find(regionstack.begin(), regionstack.end(), regname) == regionstack.end()) {
        regionstack.emplace_back(std::move(regname));
    }
}

void
TinyProfiler::StopRegion (const std::string& regname) noexcept
{
    if (!enabled) { return; }

    if (regname == regionstack.back()) {
        regionstack.pop_back();
    }
}

TinyProfileRegion::TinyProfileRegion (std::string a_regname) noexcept
    : regname(std::move(a_regname)),
      tprof(std::string("REG::")+regname, false)
{
    TinyProfiler::StartRegion(regname);
    tprof.start();
}

TinyProfileRegion::TinyProfileRegion (const char* a_regname) noexcept
    : regname(a_regname),
      tprof(std::string("REG::")+std::string(a_regname), false)
{
    TinyProfiler::StartRegion(a_regname);
    tprof.start();
}

TinyProfileRegion::~TinyProfileRegion ()
{
    tprof.stop();
    TinyProfiler::StopRegion(regname);
}

void
TinyProfiler::PrintCallStack (std::ostream& os)
{
    if (!enabled) { return; }

    os << "===== TinyProfilers ======\n";
    for (auto const& x : ttstack) {
        os << *(std::get<2>(x)) << "\n";
    }
}

}
