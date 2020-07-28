// We only support BL_PROFILE, BL_PROFILE_VAR, BL_PROFILE_VAR_STOP, BL_PROFILE_VAR_START,
// BL_PROFILE_VAR_NS, and BL_PROFILE_REGION.

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>

#include <AMReX_TinyProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>

#ifdef AMREX_USE_CUPTI
#include <AMReX_CuptiTrace.H>
#include <cupti.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


namespace amrex {

std::vector<std::string>          TinyProfiler::regionstack;
std::deque<std::tuple<double,double,std::string*> > TinyProfiler::ttstack;
std::map<std::string,std::map<std::string, TinyProfiler::Stats> > TinyProfiler::statsmap;
double TinyProfiler::t_init = std::numeric_limits<double>::max();

namespace {
    std::set<std::string> improperly_nested_timers;
    static constexpr char mainregion[] = "main";
}

TinyProfiler::TinyProfiler (std::string funcname) noexcept
    : fname(std::move(funcname)), uCUPTI(false)
{
    start();
}

TinyProfiler::TinyProfiler (std::string funcname, bool start_, bool useCUPTI) noexcept
    : fname(std::move(funcname)), uCUPTI(useCUPTI)
{
    if (start_) start();
}

TinyProfiler::TinyProfiler (const char* funcname) noexcept
    : fname(funcname), uCUPTI(false)
{
    start();
}

TinyProfiler::TinyProfiler (const char* funcname, bool start_, bool useCUPTI) noexcept
    : fname(funcname), uCUPTI(useCUPTI)
{
    if (start_) start();
}

TinyProfiler::~TinyProfiler ()
{
    stop();
}

void
TinyProfiler::start () noexcept
{
#ifdef _OPENMP
#pragma omp master
#endif
    if (stats.empty() && !regionstack.empty())
    {
        double t;
	if (!uCUPTI) {
	    t = amrex::second();
	} else {
#ifdef AMREX_USE_CUPTI
	    cudaDeviceSynchronize();
	    cuptiActivityFlushAll(0);
	    activityRecordUserdata.clear();
	    t = amrex::second();
#endif
	}

	ttstack.emplace_back(std::make_tuple(t, 0.0, &fname));
	global_depth = ttstack.size();

#ifdef AMREX_USE_CUDA
	nvtxRangePush(fname.c_str());
#endif

        for (auto const& region : regionstack)
        {
            Stats& st = statsmap[region][fname];
            ++st.depth;
            stats.push_back(&st);
        }
    }
}

void
TinyProfiler::stop () noexcept
{
#ifdef _OPENMP
#pragma omp master
#endif
    if (!stats.empty()) 
    {
        double t;
	int nKernelCalls = 0;
#ifdef AMREX_USE_CUPTI
        if (uCUPTI) {
            cudaDeviceSynchronize();
            cuptiActivityFlushAll(0);
            t = computeElapsedTimeUserdata(activityRecordUserdata);
            nKernelCalls = activityRecordUserdata.size();
        } else
#endif
        {
	    t = amrex::second();
        }

	while (static_cast<int>(ttstack.size()) > global_depth) {
	    ttstack.pop_back();
	};

	if (static_cast<int>(ttstack.size()) == global_depth)
	{
	    const std::tuple<double,double,std::string*>& tt = ttstack.back();

	    // first: wall time when the pair is pushed into the stack
	    // second: accumulated dt of children
	    double dtin;
	    double dtex;
	    if (!uCUPTI) {
	        dtin = t - std::get<0>(tt); // elapsed time since start() is called.
	        dtex = dtin - std::get<1>(tt);
	    } else {
	        dtin = t;
	        dtex = dtin - std::get<1>(tt); 
	    }

        for (Stats* st : stats)
        {
            --(st->depth);
            ++(st->n);
            if (st->depth == 0) {
                st->dtin += dtin;
            }
            st->dtex += dtex;
            st->usesCUPTI = uCUPTI;
            if (uCUPTI) {
                st->nk += nKernelCalls;
            }
        }

        ttstack.pop_back();
        if (!ttstack.empty()) {
            std::tuple<double,double,std::string*>& parent = ttstack.back();
            std::get<1>(parent) += dtin;
        }

#ifdef AMREX_USE_CUDA
        nvtxRangePop();
#endif
	} else {
	    improperly_nested_timers.insert(fname);
	}

    stats.clear();
    }
}

#ifdef AMREX_USE_CUPTI
void
TinyProfiler::stop (unsigned boxUintID) noexcept
{
#ifdef _OPENMP
#pragma omp master
#endif
    if (!stats.empty()) 
    {
        double t;
        cudaDeviceSynchronize();
        cuptiActivityFlushAll(0);
        t = computeElapsedTimeUserdata(activityRecordUserdata);
        int nKernelCalls = activityRecordUserdata.size();

        for (auto& record : activityRecordUserdata) 
        {
            record->setUintID(boxUintID);
        }

        while (static_cast<int>(ttstack.size()) > global_depth) 
        {
            ttstack.pop_back();
        };

        if (static_cast<int>(ttstack.size()) == global_depth)
        {
            const std::tuple<double,double,std::string*>& tt = ttstack.back();

            // first: wall time when the pair is pushed into the stack
            // second: accumulated dt of children
            double dtin;
            double dtex;

            dtin = t;
            dtex = dtin - std::get<1>(tt); 

            for (Stats* st : stats)
            {
                --(st->depth);
                ++(st->n);
                if (st->depth == 0) 
                {
                    st->dtin += dtin;
                }
                st->dtex += dtex;
                st->usesCUPTI = uCUPTI;
                st->nk += nKernelCalls;		
            }

            ttstack.pop_back();
            if (!ttstack.empty()) 
            {
                std::tuple<double,double,std::string*>& parent = ttstack.back();
                std::get<1>(parent) += dtin;
            }

#ifdef AMREX_USE_CUDA
            nvtxRangePop();
#endif
        } else 
        {
            improperly_nested_timers.insert(fname);
        }

        stats.clear();
    }
}
#endif

void
TinyProfiler::Initialize () noexcept
{
    regionstack.push_back(mainregion);
    t_init = amrex::second();
}

void
TinyProfiler::Finalize (bool bFlushing) noexcept
{
    static bool finalized = false;
    if (!bFlushing) {		// If flushing, don't make this the last time!
        if (finalized) {
            return;
        } else {
            finalized = true;
        }
    }

    double t_final = amrex::second();

    // make a local copy so that any functions call after this will not be recorded in the local copy.
    auto lstatsmap = statsmap;

    bool properly_nested = improperly_nested_timers.size() == 0;
    ParallelDescriptor::ReduceBoolAnd(properly_nested);
    if (!properly_nested) {
        Vector<std::string> local_imp, sync_imp;
        bool synced;
        for (std::set<std::string>::const_iterator it = improperly_nested_timers.begin();
             it != improperly_nested_timers.end(); ++it)
        {
            local_imp.push_back(*it);
        }

        amrex::SyncStrings(local_imp, sync_imp, synced);

        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "\nWARNING: TinyProfilers not properly nested!!!\n";
            for (int i = 0; i < sync_imp.size(); ++i) {
                amrex::Print() << "     " << sync_imp[i] << "\n";
            }
            amrex::Print() << std::endl;
        }
    }

    int nprocs = ParallelDescriptor::NProcs();
    int ioproc = ParallelDescriptor::IOProcessorNumber();

    double dt_max = t_final - t_init;
    ParallelReduce::Max(dt_max, ioproc, ParallelDescriptor::Communicator());
    double dt_min = t_final - t_init;
    ParallelReduce::Min(dt_min, ioproc, ParallelDescriptor::Communicator());
    double dt_avg = t_final - t_init;
    ParallelReduce::Sum(dt_avg, ioproc, ParallelDescriptor::Communicator());
    dt_avg /= double(nprocs);

    if  (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "\n\n";
        amrex::Print().SetPrecision(4)
            <<"TinyProfiler total time across processes [min...avg...max]: " 
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

    PrintStats(lstatsmap[mainregion], dt_max);
    for (auto& kv : lstatsmap) {
        if (kv.first != mainregion) {
            amrex::Print() << "\n\nBEGIN REGION " << kv.first << "\n";
            PrintStats(kv.second, dt_max);
            amrex::Print() << "END REGION " << kv.first << "\n";
        }
    }
}

void
TinyProfiler::PrintStats (std::map<std::string,Stats>& regstats, double dt_max)
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

    if (regstats.empty()) return;

    int nprocs = ParallelDescriptor::NProcs();
    int ioproc = ParallelDescriptor::IOProcessorNumber();

    std::vector<ProcStats> allprocstats;
    int maxfnamelen = 0;
    Long maxncalls = 0;

    // now collect global data onto the ioproc
    for (auto it = regstats.cbegin(); it != regstats.cend(); ++it)
    {
        Long n = it->second.n;
        double dts[2] = {it->second.dtin, it->second.dtex};

        std::vector<Long> ncalls(nprocs);
        std::vector<double> dtdt(2*nprocs);

        if (ParallelDescriptor::NProcs() == 1)
        {
            ncalls[0] = n;
            dtdt[0] = dts[0];
            dtdt[1] = dts[1];
        } else
        {
            ParallelDescriptor::Gather(&n, 1, &ncalls[0], 1, ioproc);
            ParallelDescriptor::Gather(dts, 2, &dtdt[0], 2, ioproc);
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
            pst.fname = it->first;
#ifdef AMREX_USE_CUPTI
            pst.usesCUPTI = it->second.usesCUPTI;
#endif
            allprocstats.push_back(pst);
            maxfnamelen = std::max(maxfnamelen, int(pst.fname.size()));
            maxncalls = std::max(maxncalls, pst.nmax);
        }
    }

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::OutStream() << std::setfill(' ') << std::setprecision(4);
        int wt = 9;

        int wnc = (int) std::log10 ((double) maxncalls) + 1;
        wnc = std::max(wnc, int(std::string("NCalls").size()));
        wt  = std::max(wt,  int(std::string("Excl. Min").size()));
        int wp = 6;
        wp  = std::max(wp,  int(std::string("Max %").size()));

        const std::string hline(maxfnamelen+wnc+2+(wt+2)*3+wp+2,'-');
#ifdef AMREX_USE_CUPTI
        const std::string hlinehlf((maxfnamelen+wnc+2+(wt+2)*3+wp+2)/2-12,'-');
#endif
        // Exclusive time
        std::sort(allprocstats.begin(), allprocstats.end(), ProcStats::compex);
        amrex::OutStream() << "\n" << hline << "\n";
        amrex::OutStream() << std::left
                           << std::setw(maxfnamelen) << "Name"
                           << std::right
                           << std::setw(wnc+2) << "NCalls"
                           << std::setw(wt+2) << "Excl. Min"
                           << std::setw(wt+2) << "Excl. Avg"
                           << std::setw(wt+2) << "Excl. Max"
                           << std::setw(wp+2)  << "Max %"
                           << "\n" << hline << "\n";
        for (auto it = allprocstats.cbegin(); it != allprocstats.cend(); ++it)
        {
#ifdef AMREX_USE_CUPTI
            if (it->usesCUPTI)
            {
                amrex::OutStream() << hlinehlf << "START CUPTI Trace Stats-"
                                   << hlinehlf << "\n";
            }
#endif
            amrex::OutStream() << std::setprecision(4) << std::left
                               << std::setw(maxfnamelen) << it->fname
                               << std::right
                               << std::setw(wnc+2) << it->navg
                               << std::setw(wt+2) << it->dtexmin
                               << std::setw(wt+2) << it->dtexavg
                               << std::setw(wt+2) << it->dtexmax
                               << std::setprecision(2) << std::setw(wp+1) << std::fixed 
                               << it->dtexmax*(100.0/dt_max) << "%";
            amrex::OutStream().unsetf(std::ios_base::fixed);
            amrex::OutStream() << "\n";	    
#ifdef AMREX_USE_CUPTI	    
            if (it->usesCUPTI)
            {
                amrex::OutStream() << std::setprecision(4) << std::left
                                   << std::setw(maxfnamelen) // it->fname
                                   << std::right
                                   << std::setw(wnc+2) // it->navg
                                   << std::setw(wt+2) // it->dtexmin
                                   << std::setw(wt+2) // it->dtexavg
                                   << std::setw(wt+2) // it->dtexmax
                                   << std::setprecision(2) << std::setw(wp+1) << std::fixed; // it->dtexmax*(100.0/dt_max)
                amrex::OutStream().unsetf(std::ios_base::fixed);
                amrex::OutStream();
                amrex::OutStream() << hlinehlf << "--END CUPTI Trace Stats-" << hlinehlf << "\n";
            }
#endif
        }
        amrex::OutStream() << hline << "\n";

        // Inclusive time
        std::sort(allprocstats.begin(), allprocstats.end(), ProcStats::compin);
        amrex::OutStream() << "\n" << hline << "\n";
        amrex::OutStream() << std::left
                           << std::setw(maxfnamelen) << "Name"
                           << std::right
                           << std::setw(wnc+2) << "NCalls"
                           << std::setw(wt+2) << "Incl. Min"
                           << std::setw(wt+2) << "Incl. Avg"
                           << std::setw(wt+2) << "Incl. Max"
                           << std::setw(wp+2)  << "Max %"
                           << "\n" << hline << "\n";
        for (auto it = allprocstats.cbegin(); it != allprocstats.cend(); ++it)
        {
#ifdef AMREX_USE_CUPTI
            if (it->usesCUPTI) 
            {
                amrex::OutStream() << hlinehlf << "START CUPTI Trace Stats-" << hlinehlf << "\n";
            }
#endif
            amrex::OutStream() << std::setprecision(4) << std::left
                               << std::setw(maxfnamelen) << it->fname
                               << std::right
                               << std::setw(wnc+2) << it->navg
                               << std::setw(wt+2) << it->dtinmin
                               << std::setw(wt+2) << it->dtinavg
                               << std::setw(wt+2) << it->dtinmax
                               << std::setprecision(2) << std::setw(wp+1) << std::fixed 
                               << it->dtinmax*(100.0/dt_max) << "%";
            amrex::OutStream().unsetf(std::ios_base::fixed);
            amrex::OutStream() << "\n";
#ifdef AMREX_USE_CUPTI
            if (it->usesCUPTI) 
            {
                amrex::OutStream() << std::setprecision(4) << std::left
                                   << std::setw(maxfnamelen) // it->fname
                                   << std::right
                                   << std::setw(wnc+2) // it->navg
                                   << std::setw(wt+2) // it->dtexmin
                                   << std::setw(wt+2) // it->dtexavg
                                   << std::setw(wt+2) // it->dtexmax
                                   << std::setprecision(2) << std::setw(wp+1) << std::fixed; // it->dtexmax*(100.0/dt_max)
                amrex::OutStream().unsetf(std::ios_base::fixed);
                amrex::OutStream();
                amrex::OutStream() << hlinehlf << "--END CUPTI Trace Stats-" << hlinehlf << "\n";
            }
#endif
        }
        amrex::OutStream() << hline << "\n";
        amrex::OutStream() << std::endl;
    }
}

void
TinyProfiler::StartRegion (std::string regname) noexcept
{
    if (std::find(regionstack.begin(), regionstack.end(), regname) == regionstack.end()) {
        regionstack.emplace_back(std::move(regname));
    }
}

void
TinyProfiler::StopRegion (const std::string& regname) noexcept
{
    if (regname == regionstack.back()) {
        regionstack.pop_back();
    }
}

TinyProfileRegion::TinyProfileRegion (std::string a_regname) noexcept
    : regname(std::move(a_regname)),
      tprof(std::string("REG::")+regname, false, false)
{
    TinyProfiler::StartRegion(regname);
    tprof.start();
}

TinyProfileRegion::TinyProfileRegion (const char* a_regname) noexcept
    : regname(a_regname),
      tprof(std::string("REG::")+std::string(a_regname), false, false)
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
    os << "===== TinyProfilers ======\n";
    for (auto const& x : ttstack) {
        os << *(std::get<2>(x)) << "\n";
    }
}

}
