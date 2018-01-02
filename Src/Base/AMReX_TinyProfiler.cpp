// We only support BL_PROFILE, BL_PROFILE_VAR, BL_PROFILE_VAR_STOP, BL_PROFILE_VAR_START,
// BL_PROFILE_VAR_NS, and BL_PROFILE_REGION.

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>

#include <AMReX_TinyProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

std::vector<std::string>          TinyProfiler::regionstack;
std::stack<std::pair<Real,Real> > TinyProfiler::ttstack;
std::map<std::string,std::map<std::string, TinyProfiler::Stats> > TinyProfiler::statsmap;
Real TinyProfiler::t_init = std::numeric_limits<Real>::max();

namespace {
    std::set<std::string> improperly_nested_timers;
    static constexpr char mainregion[] = "main";
}

TinyProfiler::TinyProfiler (std::string funcname)
    : fname(std::move(funcname))
{
    start();
}

TinyProfiler::TinyProfiler (std::string funcname, bool start_)
    : fname(std::move(funcname))
{
    if (start_) start();
}

TinyProfiler::TinyProfiler (const char* funcname)
    : fname(funcname)
{
    start();
}

TinyProfiler::TinyProfiler (const char* funcname, bool start_)
    : fname(funcname)
{
    if (start_) start();
}

TinyProfiler::~TinyProfiler ()
{
    stop();
}

void
TinyProfiler::start ()
{
#ifdef _OPENMP
#pragma omp master
#endif
    if (stats.empty())
    {
	Real t = amrex::second();

	ttstack.push(std::make_pair(t, 0.0));
	global_depth = ttstack.size();

        for (auto const& region : regionstack)
        {
            Stats& st = statsmap[region][fname];
            ++st.depth;
            stats.push_back(&st);
        }

#ifdef AMREX_USE_CUDA
	nvtx_id = nvtxRangeStartA(fname.c_str());
#endif
    }
}

void
TinyProfiler::stop ()
{
#ifdef _OPENMP
#pragma omp master
#endif
    if (!stats.empty()) 
    {
	Real t = amrex::second();

	while (static_cast<int>(ttstack.size()) > global_depth) {
	    ttstack.pop();
	};

	if (static_cast<int>(ttstack.size()) == global_depth)
	{
	    const std::pair<Real,Real>& tt = ttstack.top();
	    
	    // first: wall time when the pair is pushed into the stack
	    // second: accumulated dt of children
	    
	    Real dtin = t - tt.first; // elapsed time since start() is called.
	    Real dtex = dtin - tt.second;

            for (Stats* st : stats)
            {
                --(st->depth);
                ++(st->n);
                if (st->depth == 0) {
                    st->dtin += dtin;
                }
                st->dtex += dtex;
            }
                
	    ttstack.pop();
	    if (!ttstack.empty()) {
		std::pair<Real,Real>& parent = ttstack.top();
		parent.second += dtin;
	    }

#ifdef AMREX_USE_CUDA
	    nvtxRangeEnd(nvtx_id);
#endif
	} else {
	    improperly_nested_timers.insert(fname);
	} 

        stats.clear();
    }
}

void
TinyProfiler::Initialize ()
{
    regionstack.push_back(mainregion);
    t_init = amrex::second();
}

void
TinyProfiler::Finalize ()
{
    static bool finalized = false;
    if (finalized) {
	return;
    } else {
	finalized = true;
    }

    Real t_final = amrex::second();

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
	    std::cout << "\nWARNING: TinyProfilers not properly nested!!!\n";
	    for (int i = 0; i < sync_imp.size(); ++i) {
		std::cout << "     " << sync_imp[i] << "\n";
	    }
	    std::cout << std::endl;
	}
    }

    int nprocs = ParallelDescriptor::NProcs();
    int ioproc = ParallelDescriptor::IOProcessorNumber();

    Real dt_max = t_final - t_init;
    ParallelDescriptor::ReduceRealMax(dt_max, ioproc);
    Real dt_min = t_final - t_init;
    ParallelDescriptor::ReduceRealMin(dt_min, ioproc);
    Real dt_avg = t_final - t_init;
    ParallelDescriptor::ReduceRealSum(dt_avg, ioproc);
    dt_avg /= Real(nprocs);

    if  (ParallelDescriptor::IOProcessor())
    {
	std::cout << "\n\n";
	std::cout << std::setprecision(4) 
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
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "\n\nBEGIN REGION " << kv.first << "\n";
            }
            PrintStats(kv.second, dt_max);
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "END REGION " << kv.first << "\n";
            }
        }
    }
}

void
TinyProfiler::PrintStats (std::map<std::string,Stats>& regstats, Real dt_max)
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
    long maxncalls = 0;

    // now collect global data onto the ioproc
    for (auto it = regstats.cbegin(); it != regstats.cend(); ++it)
    {
	long n = it->second.n;
	Real dts[2] = {it->second.dtin, it->second.dtex};

	std::vector<long> ncalls(nprocs);
	std::vector<Real> dtdt(2*nprocs);

	if (ParallelDescriptor::NProcs() == 1) {
	    ncalls[0] = n;
	    dtdt[0] = dts[0];
	    dtdt[1] = dts[1];
	} else {
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
	    
	    allprocstats.push_back(pst);
	    maxfnamelen = std::max(maxfnamelen, int(pst.fname.size()));
	    maxncalls = std::max(maxncalls, pst.nmax);
	}
    }

    if (ParallelDescriptor::IOProcessor()) {

	std::cout << std::setfill(' ') << std::setprecision(4);
	int wt = 9;

	int wnc = (int) std::log10 ((double) maxncalls) + 1;
	wnc = std::max(wnc, int(std::string("NCalls").size()));
	wt  = std::max(wt,  int(std::string("Excl. Min").size()));
	int wp = 6;
	wp  = std::max(wp,  int(std::string("Max %").size()));

	const std::string hline(maxfnamelen+wnc+2+(wt+2)*3+wp+2,'-');

	// Exclusive time
	std::sort(allprocstats.begin(), allprocstats.end(), ProcStats::compex);
	std::cout << "\n" << hline << "\n";
	std::cout << std::left
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
	    std::cout << std::setprecision(4) << std::left
		      << std::setw(maxfnamelen) << it->fname
		      << std::right
		      << std::setw(wnc+2) << it->navg
		      << std::setw(wt+2) << it->dtexmin
		      << std::setw(wt+2) << it->dtexavg
		      << std::setw(wt+2) << it->dtexmax
		      << std::setprecision(2) << std::setw(wp+1) << std::fixed 
		      << it->dtexmax*(100.0/dt_max) << "%";
	    std::cout.unsetf(std::ios_base::fixed);
	    std::cout << "\n";
	}
	std::cout << hline << "\n";

	// Inclusive time
	std::sort(allprocstats.begin(), allprocstats.end(), ProcStats::compin);
	std::cout << "\n" << hline << "\n";
	std::cout << std::left
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
	    std::cout << std::setprecision(4) << std::left
		      << std::setw(maxfnamelen) << it->fname
		      << std::right
		      << std::setw(wnc+2) << it->navg
		      << std::setw(wt+2) << it->dtinmin
		      << std::setw(wt+2) << it->dtinavg
		      << std::setw(wt+2) << it->dtinmax
		      << std::setprecision(2) << std::setw(wp+1) << std::fixed 
		      << it->dtinmax*(100.0/dt_max) << "%";
	    std::cout.unsetf(std::ios_base::fixed);
	    std::cout << "\n";
	}
	std::cout << hline << "\n";

	std::cout << std::endl;
    }
}

void
TinyProfiler::StartRegion (std::string regname)
{
    if (std::find(regionstack.begin(), regionstack.end(), regname) == regionstack.end()) {
        regionstack.emplace_back(std::move(regname));
    }
}

void
TinyProfiler::StopRegion (const std::string& regname)
{
    if (regname == regionstack.back()) {
        regionstack.pop_back();
    }
}

TinyProfileRegion::TinyProfileRegion (std::string a_regname)
    : regname(std::move(a_regname)),
      tprof(regname, false)
{
    TinyProfiler::StartRegion(regname);
    tprof.start();
}

TinyProfileRegion::TinyProfileRegion (const char* a_regname)
    : regname(a_regname),
      tprof(a_regname, false)
{
    TinyProfiler::StartRegion(a_regname);
    tprof.start();
}

TinyProfileRegion::~TinyProfileRegion ()
{
    tprof.stop();
    TinyProfiler::StopRegion(regname);
}

}
