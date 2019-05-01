#include <algorithm>
#include <cstdio>
#include <list>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <AMReX_Geometry.H>
#include <AMReX_TagBox.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_CoordSys.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_Cluster.H>
#include <AMReX_LevelBld.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_PROB_AMR_F.H>
#include <AMReX_Amr.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabSet.H>
#include <AMReX_StateData.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#ifdef BL_USE_ARRAYVIEW
#include <DatasetClient.H>
#endif

#ifdef USE_PERILLA
#include <WorkerThread.H>
#include <Perilla.H>
#ifdef USE_PERILLA
    pthread_mutex_t teamFinishLock=PTHREAD_MUTEX_INITIALIZER;
#ifdef PERILLA_USE_UPCXX
extern struct rMsgMap_t{
    std::map< int, std::map< int,  std::list< Package* > > > map;
    volatile int size=0;
    pthread_mutex_t lock= PTHREAD_MUTEX_INITIALIZER;
}rMsgMap;
extern struct sMsgMap_t{
    std::map< int, std::map< int,  std::list< Package* > > > map;
    volatile int size=0;
    pthread_mutex_t lock= PTHREAD_MUTEX_INITIALIZER;
}sMsgMap;

#endif
#endif
#endif

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrInSituBridge.H>
#endif

namespace amrex {

//
// Static class members.  Set defaults in Initialize()!!!
//
std::list<std::string> Amr::state_plot_vars;
std::list<std::string> Amr::state_small_plot_vars;
std::list<std::string> Amr::derive_plot_vars;
std::list<std::string> Amr::derive_small_plot_vars;
bool                   Amr::first_plotfile;
bool                   Amr::first_smallplotfile;
Vector<BoxArray>       Amr::initial_ba;
Vector<BoxArray>       Amr::regrid_ba;
#ifdef BL_USE_SENSEI_INSITU
AmrInSituBridge*       Amr::insitu_bridge;
#endif

namespace
{
    const std::string CheckPointVersion("CheckPointVersion_1.0");

    bool initialized = false;
}

//Tan Nov 24, 2017 : I removed this anonymous namespace so I could access the inner variables from other source files 
//namespace   
//{
    //
    // These are all ParmParse'd in.  Set defaults in Initialize()!!!
    //
    int  plot_nfiles;
    int  mffile_nstreams;
    int  probinit_natonce;
    bool plot_files_output;
    int  checkpoint_nfiles;
    int  regrid_on_restart;
    int  use_efficient_regrid;
    int  plotfile_on_restart;
    int  insitu_on_restart;
    int  checkpoint_on_restart;
    bool checkpoint_files_output;
    int  compute_new_dt_on_regrid;
    bool precreateDirectories;
    bool prereadFAHeaders;
    VisMF::Header::Version plot_headerversion(VisMF::Header::Version_v1);
    VisMF::Header::Version checkpoint_headerversion(VisMF::Header::Version_v1);
//}



bool
Amr::UsingPrecreateDirectories () noexcept
{
    return precreateDirectories;
}

void
Amr::Initialize ()
{
    if (initialized) return;
    //
    // Set all defaults here!!!
    //
    Amr::first_plotfile      = true;
    Amr::first_smallplotfile = true;
    plot_nfiles              = 64;
    mffile_nstreams          = 1;
    probinit_natonce         = 512;
    plot_files_output        = true;
    checkpoint_nfiles        = 64;
    regrid_on_restart        = 0;
    use_efficient_regrid     = 0;
    plotfile_on_restart      = 0;
    insitu_on_restart        = 0;
    checkpoint_on_restart    = 0;
    checkpoint_files_output  = true;
    compute_new_dt_on_regrid = 0;
    precreateDirectories     = true;
    prereadFAHeaders         = true;
    plot_headerversion       = VisMF::Header::Version_v1;
    checkpoint_headerversion = VisMF::Header::Version_v1;
#ifdef BL_USE_SENSEI_INSITU
    insitu_bridge            = nullptr;
#endif
    amrex::ExecOnFinalize(Amr::Finalize);

    initialized = true;
}

void
Amr::Finalize ()
{
    Amr::state_plot_vars.clear();
    Amr::derive_plot_vars.clear();
    Amr::derive_small_plot_vars.clear();
    Amr::regrid_ba.clear();
    Amr::initial_ba.clear();
    Amr::finalizeInSitu();

    initialized = false;
}

bool Amr::Plot_Files_Output () { return plot_files_output; }

std::ostream&
Amr::DataLog (int i)
{
    return *datalog[i];
}

int
Amr::NumDataLogs () noexcept
{
    return datalog.size();
}

bool
Amr::RegridOnRestart () const noexcept
{
    return regrid_on_restart;
}

void
Amr::setDtMin (const Vector<Real>& dt_min_in) noexcept
{
    for (int i = 0; i <= finest_level; i++)
        dt_min[i] = dt_min_in[i];
}

Vector<std::unique_ptr<AmrLevel> >&
Amr::getAmrLevels () noexcept
{
    return amr_level;
}

long
Amr::cellCount (int lev) noexcept
{
    return amr_level[lev]->countCells();
}

int
Amr::numGrids (int lev) noexcept
{
    return amr_level[lev]->numGrids();
}

std::unique_ptr<MultiFab>
Amr::derive (const std::string& name,
             Real               time,
             int                lev,
             int                ngrow)
{
    return amr_level[lev]->derive(name,time,ngrow);
}

Amr::Amr ()
    :
    AmrCore()
{
    Initialize();
    InitAmr();
}

Amr::Amr (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord)
    :
    AmrCore(rb,max_level_in,n_cell_in,coord)
{
    Initialize();
    InitAmr();
}

void
Amr::InitAmr ()
{
    BL_PROFILE("Amr::InitAmr()");
    //
    // Determine physics class.
    //
    levelbld = getLevelBld();
    //
    // Global function that define state variables.
    //
    levelbld->variableSetUp();
    //
    // Set default values.
    //
    plot_int               = -1;
    small_plot_int         = -1;
    last_plotfile          = 0;
    last_smallplotfile     = -1;
    last_checkpoint        = 0;
    record_run_info        = false;
    record_grid_info       = false;
    file_name_digits       = 5;
    record_run_info_terse  = false;
    bUserStopRequest       = false;
    message_int            = 10;
#ifdef BL_USE_SENSEI_INSITU
    insitu_bridge          = nullptr;
#endif

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        isPeriodic[i] = false;

    ParmParse pp("amr");
    //
    // Check for command line flags.
    //
    pp.query("regrid_on_restart",regrid_on_restart);
    pp.query("use_efficient_regrid",use_efficient_regrid);
    pp.query("plotfile_on_restart",plotfile_on_restart);
    pp.query("insitu_on_restart",insitu_on_restart);
    pp.query("checkpoint_on_restart",checkpoint_on_restart);

    pp.query("compute_new_dt_on_regrid",compute_new_dt_on_regrid);

    pp.query("mffile_nstreams", mffile_nstreams);
    pp.query("probinit_natonce", probinit_natonce);

    probinit_natonce = std::max(1, std::min(ParallelDescriptor::NProcs(), probinit_natonce));

    pp.query("file_name_digits", file_name_digits);

    pp.query("initial_grid_file",initial_grids_file);
    pp.query("regrid_file"      , regrid_grids_file);

    pp.query("message_int", message_int);
    
    if (pp.contains("run_log"))
    {
        std::string log_file_name;
        pp.get("run_log",log_file_name);
        setRecordRunInfo(log_file_name);
    }
    if (pp.contains("run_log_terse"))
    {
        std::string log_file_name;
        pp.get("run_log_terse",log_file_name);
        setRecordRunInfoTerse(log_file_name);
    }
    if (pp.contains("grid_log"))
    {
        std::string grid_file_name;
        pp.get("grid_log",grid_file_name);
        setRecordGridInfo(grid_file_name);
    }

    if (pp.contains("data_log"))
    {
      int num_datalogs = pp.countval("data_log");
      datalog.resize(num_datalogs);
      datalogname.resize(num_datalogs);
      pp.queryarr("data_log",datalogname,0,num_datalogs);
      for (int i = 0; i < num_datalogs; i++) 
        setRecordDataInfo(i,datalogname[i]);
    }

    probin_file = "probin";  // Make "probin" the default

    if (pp.contains("probin_file"))
    {
        pp.get("probin_file",probin_file);
    }
    //
    // If set, then restart from checkpoint file.
    //
    pp.query("restart", restart_chkfile);
    //
    // If set, then restart from plotfile.
    //
    pp.query("restart_from_plotfile", restart_pltfile);

    int nlev     = max_level+1;
    dt_level.resize(nlev);
    level_steps.resize(nlev);
    level_count.resize(nlev);
    n_cycle.resize(nlev);
    dt_min.resize(nlev);
    amr_level.resize(nlev);
    //
    // Set bogus values.
    //
    for (int i = 0; i < nlev; i++)
    {
        dt_level[i]    = 1.e200; // Something nonzero so old & new will differ
        level_steps[i] = 0;
        level_count[i] = 0;
        n_cycle[i]     = 0;
        dt_min[i]      = 0.0;
    }

    // Make the default regrid_int = 1 for all levels.
    if (max_level > 0) 
    {
       regrid_int.resize(max_level);
       for (int i = 0; i < max_level; i++)
           regrid_int[i]  = 1;
    }
    
    //
    // Setup plot and checkpoint controls.
    //
    initPltAndChk();
    
    //
    // Setup insitu controls
    //
    initInSitu();

    //
    // Setup subcycling controls.
    //
    initSubcycle();

    //
    // Read in the regrid interval if max_level > 0.
    //
    if (max_level > 0) 
    {
       int numvals = pp.countval("regrid_int");
       if (numvals == 1)
       {
           //
           // Set all values to the single available value.
           //
           int the_regrid_int = 0;
           pp.query("regrid_int",the_regrid_int);
           for (int i = 0; i < max_level; i++)
           {
               regrid_int[i] = the_regrid_int;
           }
       }
       else if (numvals == 0)
       {
           if (verbose) {
               amrex::Print(amrex::ErrorStream()) << "Using default regrid_int = 1 at all levels!\n";
           }
       }
       else if (numvals < max_level)
       {
           amrex::Error("You did not specify enough values of regrid_int");
       }
       else 
       {
           //
           // Otherwise we expect a vector of max_level values
           //
           pp.queryarr("regrid_int",regrid_int,0,max_level);
       }
    }

    if (max_level > 0 && !initial_grids_file.empty())
    {
#define STRIP while( is.get() != '\n' ) {}
        std::ifstream is(initial_grids_file.c_str(),std::ios::in);

        if (!is.good())
            amrex::FileOpenFailed(initial_grids_file);

        int in_finest,ngrid;

        is >> in_finest;
        STRIP;
        initial_ba.resize(in_finest);

        use_fixed_upto_level = in_finest;
        if (in_finest > max_level)
           amrex::Error("You have fewer levels in your inputs file then in your grids file!");

        for (int lev = 1; lev <= in_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (int i = 0; i < ngrid; i++)
            {
                Box bx;
                is >> bx;
                STRIP;
                bx.refine(ref_ratio[lev-1]);
                bl.push_back(bx);
            }
            initial_ba[lev-1].define(bl);
        }
        is.close();
        if (verbose > 0) {
            amrex::Print() << "Read initial_ba. Size is " << initial_ba.size() << "\n";
        }

#undef STRIP
    }

    if (max_level > 0 && !regrid_grids_file.empty())
    {
#define STRIP while( is.get() != '\n' ) {}
        std::ifstream is(regrid_grids_file.c_str(),std::ios::in);

        if (!is.good())
            amrex::FileOpenFailed(regrid_grids_file);

        int in_finest,ngrid;

        is >> in_finest;
        STRIP;
        regrid_ba.resize(in_finest);
        for (int lev = 1; lev <= in_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (int i = 0; i < ngrid; i++)
            {
                Box bx;
                is >> bx;
                STRIP;
                 bx.refine(ref_ratio[lev-1]);
                 for (int idim = 0 ; idim < AMREX_SPACEDIM; ++idim)
                 {
                     if (bx.length(idim) > max_grid_size[lev][idim])
                     {
                         std::ostringstream ss;
                         ss << "Grid " << bx << " too large" << '\n';
                         amrex::Error(ss.str());
                     }
                 }
                 bl.push_back(bx);
            }
            regrid_ba[lev-1].define(bl);
        }
        is.close();
#undef STRIP
    }

    loadbalance_with_workestimates = 0;
    pp.query("loadbalance_with_workestimates", loadbalance_with_workestimates);

    loadbalance_level0_int = 2;
    pp.query("loadbalance_level0_int", loadbalance_level0_int);

    loadbalance_max_fac = 1.5;
    pp.query("loadbalance_max_fac", loadbalance_max_fac);
}

int
Amr::initInSitu()
{
#if defined(BL_USE_SENSEI_INSITU)
    insitu_bridge = new AmrInSituBridge;
    if (insitu_bridge->initialize())
    {
        amrex::ErrorStream() << "Amr::initInSitu : Failed to initialize." << std::endl;
        amrex::Abort();
    }
#endif
    return 0;
}

int
Amr::updateInSitu()
{
#if defined(BL_USE_SENSEI_INSITU)
    if (insitu_bridge && insitu_bridge->update(this))
    {
        amrex::ErrorStream() << "Amr::updateInSitu : Failed to update." << std::endl;
        amrex::Abort();
    }
#endif
    return 0;
}

int
Amr::finalizeInSitu()
{
#if defined(BL_USE_SENSEI_INSITU)
    if (insitu_bridge)
    {
        if (insitu_bridge->finalize())
            amrex::ErrorStream() << "Amr::finalizeInSitu : Failed to finalize." << std::endl;

        delete insitu_bridge;
        insitu_bridge = nullptr;
    }
#endif
    return 0;
}

bool
Amr::isStatePlotVar (const std::string& name)
{
    for (std::list<std::string>::const_iterator li = state_plot_vars.begin(), End = state_plot_vars.end();
         li != End;
         ++li)
    {
        if (*li == name) {
            return true;
	}
    }
    return false;
}

bool
Amr::isStateSmallPlotVar (const std::string& name)
{
    for (std::list<std::string>::const_iterator li = state_small_plot_vars.begin(), End = state_small_plot_vars.end();
         li != End;
         ++li)
    {
        if (*li == name)
            return true;
    }
    return false;
}

void
Amr::fillStatePlotVarList ()
{
    state_plot_vars.clear();
    const DescriptorList &desc_lst = AmrLevel::get_desc_lst();
    for (int typ(0); typ < desc_lst.size(); ++typ) {
        for (int comp(0); comp < desc_lst[typ].nComp(); ++comp) {
            if (desc_lst[typ].getType() == IndexType::TheCellType()) {
                state_plot_vars.push_back(desc_lst[typ].name(comp));
	    }
	}
    }
}

void
Amr::clearStatePlotVarList ()
{
    state_plot_vars.clear();
}

void
Amr::fillStateSmallPlotVarList ()
{
    state_small_plot_vars.clear();
    const DescriptorList &desc_lst = AmrLevel::get_desc_lst();
    for (int typ(0); typ < desc_lst.size(); ++typ) {
        for (int comp(0); comp < desc_lst[typ].nComp(); ++comp) {
            if (desc_lst[typ].getType() == IndexType::TheCellType()) {
                state_small_plot_vars.push_back(desc_lst[typ].name(comp));
	    }
	}
    }
}

void
Amr::clearStateSmallPlotVarList ()
{
    state_small_plot_vars.clear();
}

void
Amr::addStatePlotVar (const std::string& name)
{
    if ( ! isStatePlotVar(name)) {
        state_plot_vars.push_back(name);
    }
}

void
Amr::addStateSmallPlotVar (const std::string& name)
{
    if (!isStateSmallPlotVar(name))
        state_small_plot_vars.push_back(name);
}

void
Amr::deleteStatePlotVar (const std::string& name)
{
    if (isStatePlotVar(name)) {
        state_plot_vars.remove(name);
    }
}

bool
Amr::isDerivePlotVar (const std::string& name) noexcept
{
    for (std::list<std::string>::const_iterator li = derive_plot_vars.begin(), End = derive_plot_vars.end();
         li != End;
         ++li)
    {
        if (*li == name) {
            return true;
	}
    }

    return false;
}

bool
Amr::isDeriveSmallPlotVar (const std::string& name) noexcept
{
    for (std::list<std::string>::const_iterator li = derive_small_plot_vars.begin(), End = derive_small_plot_vars.end();
         li != End;
         ++li)
    {
        if (*li == name) {
            return true;
	}
    }

    return false;
}

void 
Amr::fillDerivePlotVarList ()
{
    derive_plot_vars.clear();
    DeriveList& derive_lst = AmrLevel::get_derive_lst();
    std::list<DeriveRec>& dlist = derive_lst.dlist();
    for (std::list<DeriveRec>::const_iterator it = dlist.begin(), End = dlist.end();
         it != End;
         ++it)
    {
        if (it->deriveType() == IndexType::TheCellType())
        {
            derive_plot_vars.push_back(it->name());
        }
    }
}

void 
Amr::fillDeriveSmallPlotVarList ()
{
    derive_small_plot_vars.clear();
    DeriveList& derive_lst = AmrLevel::get_derive_lst();
    std::list<DeriveRec>& dlist = derive_lst.dlist();
    for (std::list<DeriveRec>::const_iterator it = dlist.begin(), End = dlist.end();
         it != End;
         ++it)
    {
        if (it->deriveType() == IndexType::TheCellType())
        {
            derive_small_plot_vars.push_back(it->name());
        }
    }
}

void
Amr::clearDerivePlotVarList ()
{
    derive_plot_vars.clear();
}

void
Amr::clearDeriveSmallPlotVarList ()
{
    derive_small_plot_vars.clear();
}

void
Amr::addDerivePlotVar (const std::string& name)
{
    if (!isDerivePlotVar(name))
        derive_plot_vars.push_back(name);
}

void
Amr::addDeriveSmallPlotVar (const std::string& name)
{
    if (!isDeriveSmallPlotVar(name))
        derive_small_plot_vars.push_back(name);
}

void
Amr::deleteDerivePlotVar (const std::string& name)
{
    if (isDerivePlotVar(name))
        derive_plot_vars.remove(name);
}

void
Amr::deleteDeriveSmallPlotVar (const std::string& name)
{
    if (isDeriveSmallPlotVar(name))
        derive_small_plot_vars.remove(name);
}

Amr::~Amr ()
{
    levelbld->variableCleanUp();

    Amr::Finalize();
}

void
Amr::setRecordGridInfo (const std::string& filename)
{
    record_grid_info = true;
    if (ParallelDescriptor::IOProcessor())
    {
        gridlog.open(filename.c_str(),std::ios::out|std::ios::app);
        if (!gridlog.good())
            amrex::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordGridInfo");
}

void
Amr::setRecordRunInfo (const std::string& filename)
{
    record_run_info = true;
    if (ParallelDescriptor::IOProcessor())
    {
        runlog.open(filename.c_str(),std::ios::out|std::ios::app);
        if (!runlog.good())
            amrex::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordRunInfo");
}

void
Amr::setRecordRunInfoTerse (const std::string& filename)
{
    record_run_info_terse = true;
    if (ParallelDescriptor::IOProcessor())
    {
        runlog_terse.open(filename.c_str(),std::ios::out|std::ios::app);
        if (!runlog_terse.good())
            amrex::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordRunInfoTerse");
}

void
Amr::setRecordDataInfo (int i, const std::string& filename)
{
    if (ParallelDescriptor::IOProcessor())
    {
        datalog[i].reset(new std::fstream);
        datalog[i]->open(filename.c_str(),std::ios::out|std::ios::app);
        if (!datalog[i]->good())
            amrex::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordDataInfo");
}

void
Amr::setDtLevel (const Vector<Real>& dt_lev) noexcept
{
    for (int i = 0; i <= finest_level; i++)
        dt_level[i] = dt_lev[i];
}

void
Amr::setDtLevel (Real dt, int lev) noexcept
{
    dt_level[lev] = dt;
}

void
Amr::setNCycle (const Vector<int>& ns) noexcept
{
    for (int i = 0; i <= finest_level; i++)
        n_cycle[i] = ns[i];
}

long
Amr::cellCount () noexcept
{
    long cnt = 0;
    for (int i = 0; i <= finest_level; i++)
        cnt += amr_level[i]->countCells();
    return cnt;
}

int
Amr::numGrids () noexcept
{
    int cnt = 0;
    for (int i = 0; i <= finest_level; i++)
        cnt += amr_level[i]->numGrids();
    return cnt;
}

int
Amr::okToContinue () noexcept
{
    int ok = true;
    for (int i = 0; ok && (i <= finest_level); i++)
        ok = ok && amr_level[i]->okToContinue();
    if(bUserStopRequest) {
      ok = false;
    }
    return ok;
}

void
Amr::writePlotFile ()
{
    if ( ! Plot_Files_Output()) {
      return;
    }

    BL_PROFILE_REGION_START("Amr::writePlotFile()");
    BL_PROFILE("Amr::writePlotFile()");

    VisMF::SetNOutFiles(plot_nfiles);
    VisMF::Header::Version currentVersion(VisMF::GetHeaderVersion());
    VisMF::SetHeaderVersion(plot_headerversion);

    if (first_plotfile) {
        first_plotfile = false;
        amr_level[0]->setPlotVariables();
    }

    // Don't continue if we have no variables to plot.

    if (statePlotVars().size() == 0) {
      return;
    }

    Real dPlotFileTime0 = amrex::second();

    const std::string& pltfile = amrex::Concatenate(plot_file_root,level_steps[0],file_name_digits);

    if (verbose > 0) {
	amrex::Print() << "PLOTFILE: file = " << pltfile << '\n';
    }

    if (record_run_info && ParallelDescriptor::IOProcessor()) {
        runlog << "PLOTFILE: file = " << pltfile << '\n';
    }

  amrex::StreamRetry sretry(pltfile, abort_on_stream_retry_failure,
                             stream_max_tries);

  const std::string pltfileTemp(pltfile + ".temp");

  while(sretry.TryFileOutput()) {
    //
    //  if either the pltfile or pltfileTemp exists, rename them
    //  to move them out of the way.  then create pltfile
    //  with the temporary name, then rename it back when
    //  it is finished writing.  then stream retry can rename
    //  it to a bad suffix if there were stream errors.
    //

    if (precreateDirectories) {    // ---- make all directories at once
      amrex::UtilRenameDirectoryToOld(pltfile, false);      // dont call barrier
      amrex::UtilCreateCleanDirectory(pltfileTemp, false);  // dont call barrier
      for(int i(0); i <= finest_level; ++i) {
	amr_level[i]->CreateLevelDirectory(pltfileTemp);
      }
      ParallelDescriptor::Barrier("Amr::writePlotFile:PCD");

    } else {
      amrex::UtilRenameDirectoryToOld(pltfile, false);     // dont call barrier
      amrex::UtilCreateCleanDirectory(pltfileTemp, true);  // call barrier
    }

    std::string HeaderFileName(pltfileTemp + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec(0);

    if (ParallelDescriptor::IOProcessor()) {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out | std::ios::trunc |
	                                        std::ios::binary);
        if ( ! HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
	}
        old_prec = HeaderFile.precision(15);
    }

    for (int k(0); k <= finest_level; ++k) {
        amr_level[k]->writePlotFilePre(pltfileTemp, HeaderFile);
    }

    for (int k(0); k <= finest_level; ++k) {
        amr_level[k]->writePlotFile(pltfileTemp, HeaderFile);
    }

    for (int k(0); k <= finest_level; ++k) {
        amr_level[k]->writePlotFilePost(pltfileTemp, HeaderFile);
    }

    if (ParallelDescriptor::IOProcessor()) {
        HeaderFile.precision(old_prec);
        if ( ! HeaderFile.good()) {
            amrex::Error("Amr::writePlotFile() failed");
	}
    }

    last_plotfile = level_steps[0];

    if (verbose > 0) {
        const int IOProc        = ParallelDescriptor::IOProcessorNumber();
        Real      dPlotFileTime = amrex::second() - dPlotFileTime0;

        ParallelDescriptor::ReduceRealMax(dPlotFileTime,IOProc);

	amrex::Print() << "Write plotfile time = " << dPlotFileTime << "  seconds" << "\n\n";
    }
    ParallelDescriptor::Barrier("Amr::writePlotFile::end");

    if(ParallelDescriptor::IOProcessor()) {
      std::rename(pltfileTemp.c_str(), pltfile.c_str());
    }
    ParallelDescriptor::Barrier("Renaming temporary plotfile.");
    //
    // the plotfile file now has the regular name
    //

  }  // end while

  VisMF::SetHeaderVersion(currentVersion);
  
  BL_PROFILE_REGION_STOP("Amr::writePlotFile()");
}

void
Amr::writeSmallPlotFile ()
{
    if ( ! Plot_Files_Output()) {
      return;
    }

    BL_PROFILE_REGION_START("Amr::writeSmallPlotFile()");
    BL_PROFILE("Amr::writeSmallPlotFile()");

    VisMF::SetNOutFiles(plot_nfiles);
    VisMF::Header::Version currentVersion(VisMF::GetHeaderVersion());
    VisMF::SetHeaderVersion(plot_headerversion);

    if (first_smallplotfile) {
        first_smallplotfile = false;
        amr_level[0]->setSmallPlotVariables();
    }

    // Don't continue if we have no variables to plot.
    
    if (stateSmallPlotVars().size() == 0) {
      return;
    }

    Real dPlotFileTime0 = amrex::second();

    const std::string& pltfile = amrex::Concatenate(small_plot_file_root,
                                                     level_steps[0],
                                                     file_name_digits);

    if (verbose > 0) {
	amrex::Print() << "SMALL PLOTFILE: file = " << pltfile << '\n';
    }

    if (record_run_info && ParallelDescriptor::IOProcessor()) {
        runlog << "SMALL PLOTFILE: file = " << pltfile << '\n';
    }

  amrex::StreamRetry sretry(pltfile, abort_on_stream_retry_failure,
                             stream_max_tries);

  const std::string pltfileTemp(pltfile + ".temp");

  while(sretry.TryFileOutput()) {
    //
    //  if either the pltfile or pltfileTemp exists, rename them
    //  to move them out of the way.  then create pltfile
    //  with the temporary name, then rename it back when
    //  it is finished writing.  then stream retry can rename
    //  it to a bad suffix if there were stream errors.
    //
    if (precreateDirectories) {    // ---- make all directories at once
      amrex::UtilRenameDirectoryToOld(pltfile, false);      // dont call barrier
      amrex::UtilCreateCleanDirectory(pltfileTemp, false);  // dont call barrier
      for(int i(0); i <= finest_level; ++i) {
        amr_level[i]->CreateLevelDirectory(pltfileTemp);
      }
      ParallelDescriptor::Barrier("Amr::precreate smallplotfile Directories");
    } else {
      amrex::UtilRenameDirectoryToOld(pltfile, false);     // dont call barrier
      amrex::UtilCreateCleanDirectory(pltfileTemp, true);  // call barrier
    }


    std::string HeaderFileName(pltfileTemp + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec(0);

    if (ParallelDescriptor::IOProcessor()) {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out | std::ios::trunc |
	                                        std::ios::binary);
        if ( ! HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
	}
        old_prec = HeaderFile.precision(15);
    }

    for (int k(0); k <= finest_level; ++k) {
        amr_level[k]->writeSmallPlotFile(pltfileTemp, HeaderFile);
    }

    if (ParallelDescriptor::IOProcessor()) {
        HeaderFile.precision(old_prec);
        if ( ! HeaderFile.good()) {
            amrex::Error("Amr::writeSmallPlotFile() failed");
	}
    }

    last_smallplotfile = level_steps[0];

    if (verbose > 0) {
        const int IOProc        = ParallelDescriptor::IOProcessorNumber();
        Real      dPlotFileTime = amrex::second() - dPlotFileTime0;

        ParallelDescriptor::ReduceRealMax(dPlotFileTime,IOProc);

	amrex::Print() << "Write small plotfile time = " << dPlotFileTime << "  seconds" << "\n\n";
    }
    ParallelDescriptor::Barrier("Amr::writeSmallPlotFile::end");

    if(ParallelDescriptor::IOProcessor()) {
      std::rename(pltfileTemp.c_str(), pltfile.c_str());
    }
    ParallelDescriptor::Barrier("Renaming temporary plotfile.");
    //
    // the plotfile file now has the regular name
    //

  }  // end while

  VisMF::SetHeaderVersion(currentVersion);
  
  BL_PROFILE_REGION_STOP("Amr::writeSmallPlotFile()");
}

void
Amr::checkInput ()
{
    if (max_level < 0)
        amrex::Error("checkInput: max_level not set");
    //
    // Check that blocking_factor is a power of 2.
    //
    for (int i = 0; i < max_level; i++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            int k = blocking_factor[i][idim];
            while ( k > 0 && (k%2 == 0) )
                k /= 2;
            if (k != 1)
                amrex::Error("Amr::checkInput: blocking_factor not power of 2");
        }
    }
    //
    // Check level dependent values.
    //
    for (int i = 0; i < max_level; i++)
    {
        if (MaxRefRatio(i) < 2 || MaxRefRatio(i) > 12)
            amrex::Error("Amr::checkInput: bad ref_ratios");
    }
    const Box& domain = Geom(0).Domain();
    if (!domain.ok())
        amrex::Error("level 0 domain bad or not set");
    //
    // Check that domain size is a multiple of blocking_factor[0].
    //
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        int len = domain.length(i);
        if (len%blocking_factor[0][i] != 0)
            amrex::Error("domain size not divisible by blocking_factor");
    }
    //
    // Check that max_grid_size is even.
    //
    for (int i = 0; i < max_level; i++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (max_grid_size[i][idim]%2 != 0) {
                amrex::Error("max_grid_size is not even");
            }
        }
    }

    //
    // Check that max_grid_size is a multiple of blocking_factor at every level.
    //
    for (int i = 0; i < max_level; i++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (max_grid_size[i][idim]%blocking_factor[i][idim] != 0) {
                amrex::Error("max_grid_size not divisible by blocking_factor");
            }
        }
    }

    if( ! Geometry::ProbDomain().ok()) {
        amrex::Error("Amr::checkInput: bad physical problem size");
    }

    if(verbose > 0) {
	amrex::Print() << "Successfully read inputs file ... " << '\n';
    }
}

void
Amr::init (Real strt_time,
           Real stop_time)
{
    BL_PROFILE_REGION_START("Amr::init()");
    BL_PROFILE("Amr::init()");
    if( ! restart_chkfile.empty() && restart_chkfile != "init")
    {
        restart(restart_chkfile);
    }
    else
    {
        initialInit(strt_time,stop_time);
        checkPoint();

        if(plot_int > 0 || plot_per > 0 || plot_log_per > 0) {
            writePlotFile();
        }

        if (small_plot_int > 0 || small_plot_per > 0 || small_plot_log_per > 0)
	        writeSmallPlotFile();

        updateInSitu();
    }

#ifdef BL_COMM_PROFILING
    Vector<Box> probDomain(maxLevel()+1);
    for(int i(0); i < probDomain.size(); ++i) {
	probDomain[i] = Geom(i).Domain();
    }
    BL_COMM_PROFILE_INITAMR(finest_level, max_level, ref_ratio, probDomain);
#endif
    BL_PROFILE_REGION_STOP("Amr::init()");
}

void
Amr::readProbinFile (int& a_init)
{
    BL_PROFILE("Amr::readProbinFile()");
    //
    // Populate integer array with name of probin file.
    //
    int probin_file_length = probin_file.length();

    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
        probin_file_name[i] = probin_file[i];

    if (verbose > 0)
	amrex::Print() << "Starting to call amrex_probinit ... \n";

    const int nAtOnce = probinit_natonce;
    const int MyProc  = ParallelDescriptor::MyProc();
    const int NProcs  = ParallelDescriptor::NProcs();
    const int NSets   = (NProcs + (nAtOnce - 1)) / nAtOnce;
    const int MySet   = MyProc/nAtOnce;

    Real piStart = 0, piEnd = 0, piStartAll = amrex::second();

    for (int iSet = 0; iSet < NSets; ++iSet)
    {
        if (MySet == iSet)
        {
            //
            // Call the pesky probin reader.
            //
            piStart = amrex::second();

#ifdef AMREX_DIMENSION_AGNOSTIC

            amrex_probinit(&a_init,
			   probin_file_name.dataPtr(),
			   &probin_file_length,
			   AMREX_ZFILL(Geometry::ProbLo()),
			   AMREX_ZFILL(Geometry::ProbHi()));

#else

            amrex_probinit(&a_init,
			   probin_file_name.dataPtr(),
			   &probin_file_length,
			   Geometry::ProbLo(),
			   Geometry::ProbHi());
#endif

            piEnd = amrex::second();
            const int iBuff     = 0;
            const int wakeUpPID = (MyProc + nAtOnce);
            const int tag       = (MyProc % nAtOnce);
            if (wakeUpPID < NProcs)
                ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
        if (MySet == (iSet + 1))
        {
            //
            // Next set waits.
            //
            int iBuff;
            int waitForPID = (MyProc - nAtOnce);
            int tag        = (MyProc % nAtOnce);
            ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
        }
    }

    if (verbose > 1)
    {
        const int IOProc     = ParallelDescriptor::IOProcessorNumber();
        Real      piTotal    = piEnd - piStart;
        Real      piTotalAll = amrex::second() - piStartAll;

        ParallelDescriptor::ReduceRealMax(piTotal,    IOProc);
        ParallelDescriptor::ReduceRealMax(piTotalAll, IOProc);

	amrex::Print() << "amrex_probinit max time   = " << piTotal    << '\n'
		       << "amrex_probinit total time = " << piTotalAll << '\n';
    }

    if (verbose > 0)
	amrex::Print() << "Successfully run amrex_probinit\n";
}

void
Amr::initialInit (Real              strt_time,
                  Real              stop_time,
                  const BoxArray*   lev0_grids,
                  const Vector<int>* pmap)
{
    BL_PROFILE("Amr::initialInit()");
    InitializeInit(strt_time, stop_time, lev0_grids, pmap);

    // This is a subtlety, but in the case where we are initializing the data
    //   from a plotfile, we want to use the time read in from the plotfile as 
    //   the start time instead of using "strt_time".
    // The Amr data "cumtime" has been set in InitializeInit; if we are restarting 
    //   from a plotfile, then cumtime must be re-defined in that initialization routine. 
    //   Thus here we pass "cumtime" rather than "strt_time" to FinalizeInit.
    FinalizeInit  (cumtime, stop_time);
}

void
Amr::InitializeInit(Real              strt_time,
                    Real              stop_time,
                    const BoxArray*   lev0_grids,
                    const Vector<int>* pmap)
{
    BL_PROFILE("Amr::InitializeInit()");
    BL_COMM_PROFILE_NAMETAG("Amr::InitializeInit TOP");
    if (check_input) checkInput();
    //
    // Generate internal values from user-supplied values.
    //
    finest_level = 0;
    //
    // Init problem dependent data.
    //
    int linit = true;

    if (!probin_file.empty()) {
        readProbinFile(linit);
    }

    cumtime = strt_time;
    //
    // Define base level grids.  Note that if we are restarting from a plotfile, this
    //    routine will call the level 0 AmrLevel initialization which will overwrite cumtime.
    //
    defBaseLevel(strt_time, lev0_grids, pmap);
}

void
Amr::FinalizeInit (Real              strt_time,
                   Real              stop_time)
{
    BL_PROFILE("Amr::FinalizeInit()");
    //
    // Compute dt and set time levels of all grid data.
    //
    amr_level[0]->computeInitialDt(finest_level,
				   sub_cycle,
				   n_cycle,
				   ref_ratio,
				   dt_level,
				   stop_time);
    //
    // The following was added for multifluid.
    //
    Real dt0   = dt_level[0];
    dt_min[0]  = dt_level[0];
    n_cycle[0] = 1;

    if (max_level > 0)
        bldFineLevels(strt_time);

    for (int lev = 1; lev <= finest_level; lev++)
    {
        dt0           /= n_cycle[lev];
        dt_level[lev]  = dt0;
        dt_min[lev]    = dt_level[lev];
    }

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev]->setTimeLevel(strt_time,dt_level[lev],dt_level[lev]);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev]->post_regrid(0,finest_level);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        level_steps[lev] = 0;
        level_count[lev] = 0;
    }

    //
    // Perform any special post_initialization operations.
    //
    for(int lev(0); lev <= finest_level; ++lev) {
      amr_level[lev]->post_init(stop_time);
    }

    if (ParallelDescriptor::IOProcessor())
    {
       if (verbose > 1)
       {
           amrex::Print() << "INITIAL GRIDS \n";
           printGridInfo(amrex::OutStream(),0,finest_level);
       }
       else if (verbose > 0)
       { 
           amrex::Print() << "INITIAL GRIDS \n";
           printGridSummary(amrex::OutStream(),0,finest_level);
       }
    }

    if (record_grid_info && ParallelDescriptor::IOProcessor())
    {
        gridlog << "INITIAL GRIDS \n";
        printGridInfo(gridlog,0,finest_level);
    }

    BL_COMM_PROFILE_NAMETAG("Amr::initialInit BOTTOM");
}

void
Amr::restart (const std::string& filename)
{
    BL_PROFILE_REGION_START("Amr::restart()");
    BL_PROFILE("Amr::restart()");

    which_level_being_advanced = -1;

    Real dRestartTime0 = amrex::second();

    VisMF::SetMFFileInStreams(mffile_nstreams);

    if (verbose > 0) {
	amrex::Print() << "restarting calculation from file: " << filename << "\n";
    }

    if (record_run_info && ParallelDescriptor::IOProcessor()) {
        runlog << "RESTART from file = " << filename << '\n';
    }
    //
    // Init problem dependent data.
    //
    int linit = false;

    readProbinFile(linit);
    //
    // Start calculation from given restart file.
    //
    if (record_run_info && ParallelDescriptor::IOProcessor()) {
        runlog << "RESTART from file = " << filename << '\n';
    }

    // ---- preread and broadcast all FabArray headers if this file exists
    std::map<std::string, Vector<char> > faHeaderMap;
    if(prereadFAHeaders) {
      // ---- broadcast the file with the names of the fabarray headers
      std::string faHeaderFilesName(filename + "/FabArrayHeaders.txt");
      Vector<char> faHeaderFileChars;
      bool bExitOnError(false);  // ---- dont exit if this file does not exist
      ParallelDescriptor::ReadAndBcastFile(faHeaderFilesName, faHeaderFileChars,
                                           bExitOnError);
      if(faHeaderFileChars.size() > 0) {  // ---- headers were read
        std::string faFileCharPtrString(faHeaderFileChars.dataPtr());
        std::istringstream fais(faFileCharPtrString, std::istringstream::in);
        while ( ! fais.eof()) {  // ---- read and broadcast each header
          std::string faHeaderName;
          fais >> faHeaderName;
          if( ! fais.eof()) {
            std::string faHeaderFullName(filename + '/' + faHeaderName + "_H");
            Vector<char> &tempCharArray = faHeaderMap[faHeaderFullName];
            ParallelDescriptor::ReadAndBcastFile(faHeaderFullName, tempCharArray);
	    if(verbose > 2) {
		amrex::Print() 
		    << ":::: faHeaderName faHeaderFullName tempCharArray.size() = " << faHeaderName
		    << "  " << faHeaderFullName << "  " << tempCharArray.size() << "\n";
	    }
          }
        }
        StateData::SetFAHeaderMapPtr(&faHeaderMap);
      }
    }

    //
    // Open the checkpoint header file for reading.
    //
    std::string File(filename + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);
    //
    // Read global data.
    //
    // Attempt to differentiate between old and new CheckPointFiles.
    //
    int         spdim;
    bool        new_checkpoint_format = false;
    std::string first_line;

    std::getline(is,first_line);

    if (first_line == CheckPointVersion)
    {
        new_checkpoint_format = true;
        is >> spdim;
    }
    else
    {
        spdim = atoi(first_line.c_str());
    }

    if (spdim != AMREX_SPACEDIM)
    {
        amrex::ErrorStream() << "Amr::restart(): bad spacedim = " << spdim << '\n';
        amrex::Abort();
    }

    is >> cumtime;
    int mx_lev;
    is >> mx_lev;
    is >> finest_level;

    Vector<Box> inputs_domain(max_level+1);
    for (int lev = 0; lev <= max_level; ++lev)
    {
	Box bx(Geom(lev).Domain().smallEnd(),Geom(lev).Domain().bigEnd());
       inputs_domain[lev] = bx;
    }

    if (max_level >= mx_lev) {

       for (int i(0); i <= mx_lev; ++i) { is >> Geom(i);      }
       for (int i(0); i <  mx_lev; ++i) { is >> ref_ratio[i]; }
       for (int i(0); i <= mx_lev; ++i) { is >> dt_level[i];  }

       if (new_checkpoint_format)
       {
           for (int i(0); i <= mx_lev; ++i) { is >> dt_min[i]; }
       }
       else
       {
           for (int i(0); i <= mx_lev; ++i) { dt_min[i] = dt_level[i]; }
       }

       Vector<int>  n_cycle_in;
       n_cycle_in.resize(mx_lev+1);  
       for (int i(0); i <= mx_lev; ++i) { is >> n_cycle_in[i]; }
       bool any_changed = false;

       for (int i(0); i <= mx_lev; ++i) {
           if (n_cycle[i] != n_cycle_in[i]) {
               any_changed = true;
               if (verbose > 0) {
		   amrex::Print() << "Warning: n_cycle has changed at level " << i << 
		       " from " << n_cycle_in[i] << " to " << n_cycle[i] << "\n";
	       }
           }
       }

       // If we change n_cycle then force a full regrid from level 0 up
       if (max_level > 0 && any_changed)
       {
           level_count[0] = regrid_int[0];
           if (verbose > 0) {
	       amrex::Print() << "Warning: This forces a full regrid \n";
	   }
       }


       for (int i(0); i <= mx_lev; ++i) { is >> level_steps[i]; }
       for (int i(0); i <= mx_lev; ++i) { is >> level_count[i]; }

       //
       // Set bndry conditions.
       //
       if (max_level > mx_lev)
       {
           for (int i(mx_lev + 1); i <= max_level; ++i)
           {
               dt_level[i]    = dt_level[i-1]/n_cycle[i];
               level_steps[i] = n_cycle[i]*level_steps[i-1];
               level_count[i] = 0;
           }

           // This is just an error check
           if ( ! sub_cycle)
           {
               for (int i(1); i <= finest_level; ++i)
               {
                   if (dt_level[i] != dt_level[i-1]) {
                      amrex::Error("restart: must have same dt at all levels if not subcycling");
		   }
               }
           }
       }

       if (regrid_on_restart && max_level > 0)
       {
           if (regrid_int[0] > 0) {
               level_count[0] = regrid_int[0];
	   } else {
               amrex::Error("restart: can't have regrid_on_restart and regrid_int <= 0");
	   }
       }

       checkInput();
       //
       // Read levels.
       //
       for (int lev(0); lev <= finest_level; ++lev)
       {
	   amr_level[lev].reset((*levelbld)());
           amr_level[lev]->restart(*this, is);
	   this->SetBoxArray(lev, amr_level[lev]->boxArray());
	   this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
       }
       //
       // Build any additional data structures.
       //
       for (int lev = 0; lev <= finest_level; lev++) {
           amr_level[lev]->post_restart();
       }

    } else {

       if (ParallelDescriptor::IOProcessor()) {
          amrex::Warning("Amr::restart(): max_level is lower than before");
       }

       int new_finest_level = std::min(max_level,finest_level);

       finest_level = new_finest_level;
 
       // These are just used to hold the extra stuff we have to read in.
       Geometry   geom_dummy;
       Real       real_dummy;
       int         int_dummy;
       IntVect intvect_dummy;

       for (int i(0)            ; i <= max_level; ++i) { is >> Geom(i); }
       for (int i(max_level + 1); i <= mx_lev   ; ++i) { is >> geom_dummy; }

       for (int i(0)        ; i <  max_level; ++i) { is >> ref_ratio[i]; }
       for (int i(max_level); i <  mx_lev   ; ++i) { is >> intvect_dummy; }

       for (int i(0)            ; i <= max_level; ++i) { is >> dt_level[i]; }
       for (int i(max_level + 1); i <= mx_lev   ; ++i) { is >> real_dummy; }

       if (new_checkpoint_format) {
           for (int i(0)            ; i <= max_level; ++i) { is >> dt_min[i]; }
           for (int i(max_level + 1); i <= mx_lev   ; ++i) { is >> real_dummy; }
       } else {
           for (int i(0); i <= max_level; ++i) { dt_min[i] = dt_level[i]; }
       }

       for (int i(0)            ; i <= max_level; ++i) { is >> n_cycle[i]; }
       for (int i(max_level + 1); i <= mx_lev   ; ++i) { is >> int_dummy; }

       for (int i(0)            ; i <= max_level; ++i) { is >> level_steps[i]; }
       for (int i(max_level + 1); i <= mx_lev   ; ++i) { is >> int_dummy; }

       for (int i(0)            ; i <= max_level; ++i) { is >> level_count[i]; }
       for (int i(max_level + 1); i <= mx_lev   ; ++i) { is >> int_dummy; }

       if (regrid_on_restart && max_level > 0) {
           if (regrid_int[0] > 0)  {
               level_count[0] = regrid_int[0];
	   } else {
               amrex::Error("restart: can't have regrid_on_restart and regrid_int <= 0");
	   }
       }

       checkInput();

       //
       // Read levels.
       //
       for (int lev = 0; lev <= new_finest_level; lev++)
       {
	   amr_level[lev].reset((*levelbld)());
           amr_level[lev]->restart(*this, is);
	   this->SetBoxArray(lev, amr_level[lev]->boxArray());
	   this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
       }
       //
       // Build any additional data structures.
       //
       for (int lev = 0; lev <= new_finest_level; lev++) {
           amr_level[lev]->post_restart();
       }
    }

    // Save the number of steps taken so far. This mainly
    // helps in the edge case where we end up not taking
    // any timesteps before the run terminates, so that
    // we know not to unnecessarily overwrite the old file.
    last_checkpoint = level_steps[0];
    last_plotfile = level_steps[0];

    for (int lev = 0; lev <= finest_level; ++lev)
    {
	Box restart_domain(Geom(lev).Domain());
       if ( ! (inputs_domain[lev] == restart_domain) )
       {
           std::ostringstream ss;
           ss  << "Problem at level " << lev << '\n'
	       << "Domain according to     inputs file is " <<  inputs_domain[lev] << '\n'
	       << "Domain according to checkpoint file is " << restart_domain      << '\n'
	       << "Amr::restart() failed -- box from inputs file does not "
	       << "equal box from restart file. \n";
           amrex::Abort(ss.str());
       }
    }

    if (verbose > 0)
    {
        Real dRestartTime = amrex::second() - dRestartTime0;

        ParallelDescriptor::ReduceRealMax(dRestartTime,ParallelDescriptor::IOProcessorNumber());

	amrex::Print() << "Restart time = " << dRestartTime << " seconds." << '\n';
    }
    BL_PROFILE_REGION_STOP("Amr::restart()");
}

void
Amr::checkPoint ()
{
    if( ! checkpoint_files_output) {
      return;
    }

    BL_PROFILE_REGION_START("Amr::checkPoint()");
    BL_PROFILE("Amr::checkPoint()");

    VisMF::SetNOutFiles(checkpoint_nfiles);
    //
    // In checkpoint files always write out FABs in NATIVE format.
    //
    FABio::Format thePrevFormat = FArrayBox::getFormat();

    FArrayBox::setFormat(FABio::FAB_NATIVE);

    VisMF::Header::Version currentVersion(VisMF::GetHeaderVersion());
    VisMF::SetHeaderVersion(checkpoint_headerversion);

    Real dCheckPointTime0 = amrex::second();

    const std::string& ckfile = amrex::Concatenate(check_file_root,level_steps[0],file_name_digits);

    if(verbose > 0) {
	amrex::Print() << "CHECKPOINT: file = " << ckfile << "\n";
    }

    if(record_run_info && ParallelDescriptor::IOProcessor()) {
        runlog << "CHECKPOINT: file = " << ckfile << '\n';
    }


  amrex::StreamRetry sretry(ckfile, abort_on_stream_retry_failure,
                             stream_max_tries);

  const std::string ckfileTemp(ckfile + ".temp");

  while(sretry.TryFileOutput()) {

    StateData::ClearFabArrayHeaderNames();

    //
    //  if either the ckfile or ckfileTemp exists, rename them
    //  to move them out of the way.  then create ckfile
    //  with the temporary name, then rename it back when
    //  it is finished writing.  then stream retry can rename
    //  it to a bad suffix if there were stream errors.
    //

    if (precreateDirectories) {    // ---- make all directories at once
      amrex::UtilRenameDirectoryToOld(ckfile, false);      // dont call barrier
      amrex::UtilCreateCleanDirectory(ckfileTemp, false);  // dont call barrier
      for (int i(0); i <= finest_level; ++i) 
      {
        amr_level[i]->CreateLevelDirectory(ckfileTemp);
      }
      ParallelDescriptor::Barrier("Amr::precreateDirectories");
    } else {
      amrex::UtilRenameDirectoryToOld(ckfile, false);     // dont call barrier
      amrex::UtilCreateCleanDirectory(ckfileTemp, true);  // call barrier
    }

    std::string HeaderFileName = ckfileTemp + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec = 0;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out | std::ios::trunc |
	                                        std::ios::binary);

        if ( ! HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
	}

        old_prec = HeaderFile.precision(17);

        HeaderFile << CheckPointVersion << '\n'
                   << AMREX_SPACEDIM       << '\n'
                   << cumtime           << '\n'
                   << max_level         << '\n'
                   << finest_level      << '\n';
        //
        // Write out problem domain.
        //
        for (int i(0); i <= max_level; ++i) { HeaderFile << Geom(i)        << ' '; }
        HeaderFile << '\n';
        for (int i(0); i < max_level; ++i)  { HeaderFile << ref_ratio[i]   << ' '; }
        HeaderFile << '\n';
        for (int i(0); i <= max_level; ++i) { HeaderFile << dt_level[i]    << ' '; }
        HeaderFile << '\n';
        for (int i(0); i <= max_level; ++i) { HeaderFile << dt_min[i]      << ' '; }
        HeaderFile << '\n';
        for (int i(0); i <= max_level; ++i) { HeaderFile << n_cycle[i]     << ' '; }
        HeaderFile << '\n';
        for (int i(0); i <= max_level; ++i) { HeaderFile << level_steps[i] << ' '; }
        HeaderFile << '\n';
        for (int i(0); i <= max_level; ++i) { HeaderFile << level_count[i] << ' '; }
        HeaderFile << '\n';
    }

    for (int i = 0; i <= finest_level; ++i) {
        amr_level[i]->checkPointPre(ckfileTemp, HeaderFile);
    }

    for (int i = 0; i <= finest_level; ++i) {
        amr_level[i]->checkPoint(ckfileTemp, HeaderFile);
    }

    for (int i = 0; i <= finest_level; ++i) {
        amr_level[i]->checkPointPost(ckfileTemp, HeaderFile);
    }

    if (ParallelDescriptor::IOProcessor()) {
	const Vector<std::string> &FAHeaderNames = StateData::FabArrayHeaderNames();
	if(FAHeaderNames.size() > 0) {
          std::string FAHeaderFilesName = ckfileTemp + "/FabArrayHeaders.txt";
          std::ofstream FAHeaderFile(FAHeaderFilesName.c_str(),
	                             std::ios::out | std::ios::trunc |
	                             std::ios::binary);
          if ( ! FAHeaderFile.good()) {
              amrex::FileOpenFailed(FAHeaderFilesName);
	  }

	  for(int i(0); i < FAHeaderNames.size(); ++i) {
	    FAHeaderFile << FAHeaderNames[i] << '\n';
	  }
	}
    }

    if(ParallelDescriptor::IOProcessor()) {
        HeaderFile.precision(old_prec);

        if( ! HeaderFile.good()) {
            amrex::Error("Amr::checkpoint() failed");
	}
    }

    last_checkpoint = level_steps[0];

    if (verbose > 0)
    {
        Real dCheckPointTime = amrex::second() - dCheckPointTime0;

        ParallelDescriptor::ReduceRealMax(dCheckPointTime,
	                            ParallelDescriptor::IOProcessorNumber());

	amrex::Print() << "checkPoint() time = " << dCheckPointTime << " secs." << '\n';
    }
    ParallelDescriptor::Barrier("Amr::checkPoint::end");

    if(ParallelDescriptor::IOProcessor()) {
      std::rename(ckfileTemp.c_str(), ckfile.c_str());
    }
    ParallelDescriptor::Barrier("Renaming temporary checkPoint file.");

  }  // end while

  //
  // Restore the previous FAB format.
  //
  FArrayBox::setFormat(thePrevFormat);

  VisMF::SetHeaderVersion(currentVersion);

  BL_PROFILE_REGION_STOP("Amr::checkPoint()");
}

void
Amr::RegridOnly (Real time, bool do_io)
{
    BL_ASSERT(regrid_on_restart == 1);

    int lev_top = std::min(finest_level, max_level-1);

    for (int i = 0; i <= lev_top; i++)
       regrid(i,time);

    if (do_io) {

        if (plotfile_on_restart)
            writePlotFile();

        if (checkpoint_on_restart)
            checkPoint();

        if (insitu_on_restart)
            updateInSitu();

    }
}

void
Amr::timeStep (int  level,
               Real time,
               int  iteration,
               int  niter,
               Real stop_time)
{
#if defined(USE_PERILLA_PTHREADS) || defined(USE_PERILLA_OMP)
    perilla::syncAllWorkerThreads();
    if(perilla::isMasterThread())
    {
#endif
    BL_PROFILE("Amr::timeStep()");
    BL_COMM_PROFILE_NAMETAG("Amr::timeStep TOP");

    // This is used so that the AmrLevel functions can know which level is being advanced 
    //      when regridding is called with possible lbase > level.
    which_level_being_advanced = level;


    // Update so that by default, we don't force a post-step regrid.
    amr_level[level]->setPostStepRegrid(0);

    //
    // Allow regridding of level 0 calculation on restart.
    //
    if (max_level == 0 && regrid_on_restart)
    {
	regrid_level_0_on_restart();
    }
    else
    {
        int lev_top = std::min(finest_level, max_level-1);

#ifdef USE_PERILLA
        int cnt=0;
        bool *metadataChanged=new bool[finest_level+1];
        for (int l=0; l <= finest_level; l++)
            metadataChanged[l]=false;
#endif

        for (int i(level); i <= lev_top; ++i)
        {
            const int old_finest = finest_level;

            if (okToRegrid(i))
            {
#ifdef USE_PERILLA
#if defined(USE_PERILLA_PTHREADS) || defined(USE_PERILLA_OMP)
		//ask the communication thread to stop so that I can update the metadata
                Perilla::updateMetadata_request=1;
		while(!Perilla::updateMetadata_noticed){

		}
#endif
                //for (int k(i>0?i-1:0); k <= finest_level; ++k) {
                for (int k=0; k <= finest_level; ++k) {
                    if(metadataChanged[k]==false){
                        graphArray[k].clear();
                        getLevel(k).finalizePerilla(time);
                        metadataChanged[k]=true;
			cnt++;
		    }
		}
#endif
                regrid(i,time);

                //
                // Compute new dt after regrid if at level 0 and compute_new_dt_on_regrid.
                //
                if ( compute_new_dt_on_regrid && (i == 0) )
                {
                    int post_regrid_flag = 1;
                    amr_level[0]->computeNewDt(finest_level,
					       sub_cycle,
					       n_cycle,
					       ref_ratio,
					       dt_min,
					       dt_level,
					       stop_time, 
					       post_regrid_flag);
                }

                for (int k(i); k <= finest_level; ++k) {
                    level_count[k] = 0;
		}

                if (old_finest < finest_level)
                {
                    //
                    // The new levels will not have valid time steps
                    // and iteration counts.
                    //
                    for (int k(old_finest + 1); k <= finest_level; ++k)
                    {
                        dt_level[k]    = dt_level[k-1]/n_cycle[k];
                    }
                }
            }
            if (old_finest > finest_level) {
                lev_top = std::min(finest_level, max_level - 1);
	    }
        }
#ifdef USE_PERILLA
	if(cnt){
	    if(ParallelDescriptor::NProcs()>1){
	        Perilla::clearTagMap();
	        Perilla::clearMyTagMap();
	        Perilla::genTags=true;
	        Perilla::uTags=0;
	        Perilla::pTagCnt.clear();
            }
            for(int i=0; i<= finest_level; i++){
                getLevel(i).initPerilla(cumtime);
	    }
 	    Perilla::updateMetadata_done++;
	}
        delete metadataChanged;
#endif

        if (max_level == 0 && loadbalance_level0_int > 0 && loadbalance_with_workestimates)
        {
            if (level_steps[0] == 1 || level_count[0] >= loadbalance_level0_int) {
                LoadBalanceLevel0(time);
                level_count[0] = 0;
            }
        }
    }
    //
    // Check to see if should write plotfile.
    // This routine is here so it is done after the restart regrid.
    //
    if (plotfile_on_restart && ! (restart_chkfile.empty()) )
    {
	plotfile_on_restart = 0;
	writePlotFile();
    }
    //
    // Advance grids at this level.
    //
    if (verbose > 0)
    {
	amrex::Print() << "[Level " << level << " step " << level_steps[level]+1 << "] "
		       << "ADVANCE with dt = " << dt_level[level] << "\n";
    }

#if defined(USE_PERILLA_PTHREADS) || defined(USE_PERILLA_OMP)
    }
    perilla::syncAllWorkerThreads();
#endif

    BL_PROFILE_REGION_START("amr_level.advance");
    Real dt_new = amr_level[level]->advance(time,dt_level[level],iteration,niter);
    BL_PROFILE_REGION_STOP("amr_level.advance");

#if defined(USE_PERILLA_PTHREADS) || defined(USE_PERILLA_OMP)
    perilla::syncWorkerThreads();
    if(perilla::isMasterThread())
    {
#endif

    dt_min[level] = iteration == 1 ? dt_new : std::min(dt_min[level],dt_new);

    level_steps[level]++;
    level_count[level]++;

    if (verbose > 0)
    {
	amrex::Print() << "[Level " << level << " step " << level_steps[level] << "] "
		       << "Advanced " << amr_level[level]->countCells() << " cells\n";
    }

    // If the level signified that it wants a regrid after the advance has
    // occurred, do that now.
    if (amr_level[level]->postStepRegrid()) {

	int old_finest = finest_level;

	regrid(level, time);

	if (old_finest < finest_level)
	{
	    //
	    // The new levels will not have valid time steps.
	    //
	    for (int k = old_finest + 1; k <= finest_level; ++k)
	    {
		dt_level[k] = dt_level[k-1] / n_cycle[k];
	    }
	}
#ifdef USE_PERILLA
//        getLevel(level).finalizePerilla(cumtime);
//        getLevel(level).initPerilla(cumtime);
#endif
    }

#if defined(USE_PERILLA_PTHREADS) || defined(USE_PERILLA_OMP)
    }
    perilla::syncAllWorkerThreads();
#endif

    //
    // Advance grids at higher level.
    //
    if (level < finest_level)
    {
        const int lev_fine = level+1;

        if (sub_cycle)
        {
            const int ncycle = n_cycle[lev_fine];

            BL_COMM_PROFILE_NAMETAG("Amr::timeStep timeStep subcycle");
            for (int i = 1; i <= ncycle; i++)
                timeStep(lev_fine,time+(i-1)*dt_level[lev_fine],i,ncycle,stop_time);
        }
        else
        {
            BL_COMM_PROFILE_NAMETAG("Amr::timeStep timeStep nosubcycle");
            timeStep(lev_fine,time,1,1,stop_time);
        }
    }

#if defined(USE_PERILLA_PTHREADS) || defined(USE_PERILLA_OMP)
    perilla::syncAllWorkerThreads();
#endif

    amr_level[level]->post_timestep(iteration);

#if defined(USE_PERILLA_PTHREADS) || defined(USE_PERILLA_OMP)
    perilla::syncAllWorkerThreads();
    if(perilla::isMasterThread())
    {
#endif
    // Set this back to negative so we know whether we are in fact in this routine
    which_level_being_advanced = -1;
#if defined(USE_PERILLA_PTHREADS) || defined(USE_PERILLA_OMP)
    }
    perilla::syncAllWorkerThreads();
#endif
}

Real
Amr::coarseTimeStepDt (Real stop_time)
{
    coarseTimeStep(stop_time);
    return dt_level[0];
}

void
Amr::coarseTimeStep (Real stop_time)
{
    Real      run_stop;
    Real run_strt;
#ifdef USE_PERILLA_PTHREADS
    //mpi+pthreads (default) or upcxx+pthreads
    std::vector<RegionGraph*> flattenedGraphArray;
    perilla::syncAllThreads();
    if(perilla::isMasterThread())
    {
#endif
    BL_PROFILE_REGION_START("Amr::coarseTimeStep()");
    BL_PROFILE("Amr::coarseTimeStep()");
    std::stringstream stepName;
    stepName << "timeStep STEP " << level_steps[0];

    run_strt = amrex::second() ;

    //
    // Compute new dt.
    //
    if (levelSteps(0) > 0)
    {
        int post_regrid_flag = 0;
        amr_level[0]->computeNewDt(finest_level,
				   sub_cycle,
				   n_cycle,
				   ref_ratio,
				   dt_min,
				   dt_level,
				   stop_time,
				   post_regrid_flag);
    }
    else
    {
        amr_level[0]->computeInitialDt(finest_level,
				       sub_cycle,
				       n_cycle,
				       ref_ratio,
				       dt_level,
				       stop_time);
    }

    BL_PROFILE_REGION_START(stepName.str());

#ifdef USE_PERILLA
    std::vector<RegionGraph*> flattenedGraphArray;
#ifdef USE_PERILLA_PTHREADS
    //mpi+pthreads (default) or upcxx+pthreads
    }
    perilla::syncAllThreads();

    if(perilla::isMasterThread()){
        Perilla::updateMetadata_request = 0;
        Perilla::updateMetadata_noticed = 0;
        Perilla::updateMetadata_done = 0;
        Perilla::numTeamsFinished = 0;
        RegionGraph::graphCnt = 0;
        if(levelSteps(0)==0){
	    graphArray.resize(finest_level+1);
            for(int i=0; i<= finest_level; i++)
                getLevel(i).initPerilla(cumtime);
	    if(ParallelDescriptor::NProcs()>1){
  	        Perilla::syncProcesses();
                Perilla::communicateTags();
	        Perilla::syncProcesses();
	    }
        }
    }
    perilla::syncAllThreads();

    if(perilla::isCommunicationThread())
    {
        Perilla::flattenGraphHierarchy(graphArray, flattenedGraphArray);
	bool doublechecked=false;
        while(true){
   	    if(!Perilla::updateMetadata_request){
                Perilla::serviceMultipleGraphCommDynamic(flattenedGraphArray,true,perilla::tid());
                if( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
		{
                    Perilla::syncProcesses();
	            flattenedGraphArray.clear();
                    Perilla::syncProcesses();
                    break;
		}
            }else{
	        Perilla::syncProcesses();
        	for(int g=0; g<flattenedGraphArray.size(); g++)
          	{
		       //cancel messages preposted previously
		       flattenedGraphArray[g]->graphTeardown();
		}
#ifdef PERILLA_USE_UPCXX
                    pthread_mutex_lock(&(rMsgMap.lock));
                    for(int i=0; i<rMsgMap.map.size(); i++){
                        for(int j=0; j<rMsgMap.map[i].size(); j++){
                            while(rMsgMap.map[i][j].size()>0){
                               rMsgMap.map[i][j].pop_front();
                               rMsgMap.size--;
                            }
                        }
                    }
                    pthread_mutex_unlock(&(rMsgMap.lock));
                    while(sMsgMap.size>0){
                    }
#endif
	        Perilla::syncProcesses();
	        Perilla::updateMetadata_noticed=1;
	        while(Perilla::updateMetadata_done==0){//!= (max_level+1)){
		
	        }
	        Perilla::updateMetadata_request=0;
	        Perilla::updateMetadata_noticed=0;
	        Perilla::updateMetadata_done=0;
                if(ParallelDescriptor::NProcs()>1){
	            Perilla::syncProcesses();
                    Perilla::communicateTags();
	            Perilla::syncProcesses();
		}
	        flattenedGraphArray.clear();
		Perilla::flattenGraphHierarchy(graphArray, flattenedGraphArray);
	        Perilla::serviceMultipleGraphCommDynamic(flattenedGraphArray,true,perilla::tid());

                if( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
		{
	 	    Perilla::syncProcesses();
  	            flattenedGraphArray.clear();
	 	    Perilla::syncProcesses();
                    break;
		}
 	    }
        }  
    }else{
        timeStep(0,cumtime,1,1,stop_time);
        if(perilla::isMasterWorkerThread()){
            pthread_mutex_lock(&teamFinishLock);
            Perilla::numTeamsFinished++;
            pthread_mutex_unlock(&teamFinishLock);
        }
    }

    perilla::syncAllThreads();
    if(perilla::isMasterThread()){
        if(!okToContinue() || (level_steps[0] == Perilla::max_step) || (stop_time -(dt_level[0] + cumTime())<=0)){
            for(int i=0; i<= finest_level; i++){
                getLevel(i).finalizePerilla(cumtime);
            }
        }
    }
#else
    Perilla::numTeamsFinished = 0;
    RegionGraph::graphCnt = 0;
    if(levelSteps(0)==0){
	graphArray.resize(finest_level+1);
        for(int i=0; i<= finest_level; i++)
            getLevel(i).initPerilla(cumtime);
        if(ParallelDescriptor::NProcs()>1){
            Perilla::communicateTags();
        }
    }

#ifdef USE_PERILLA_OMP
#pragma omp parallel
    {
        if(perilla::isCommunicationThread())
        {
            while(true){
                Perilla::flattenGraphHierarchy(graphArray, flattenedGraphArray);
                Perilla::serviceMultipleGraphCommDynamic(flattenedGraphArray,true,perilla::tid());
                if( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
                {
                    flattenedGraphArray.clear();
                    break;
                }
            }
        }
        else{
            timeStep(0,cumtime,1,1,stop_time);
            if(perilla::isMasterWorkerThread()){
		#pragma omp atomic
                Perilla::numTeamsFinished++;
            }
        }
    }
#elif defined(USE_PERILLA_ON_DEMAND)
    //RTS on-demand
    timeStep(0,cumtime,1,1,stop_time);
#else
    cout<<"Undefined Async Mode"<<endl;
    exit(0);
#endif

#if 0
    if(!okToContinue() || (level_steps[0] == Perilla::max_step) || (stop_time -(dt_level[0] + cumTime())<=0)){
        for(int i=0; i<= finest_level; i++){
            getLevel(i).finalizePerilla(cumtime);
        }
    }
#endif
//end nonPthreads backends
#endif
//end Perilla backends
#else
    //synchronous
    timeStep(0,cumtime,1,1,stop_time);
#endif

#ifdef USE_PERILLA_PTHREADS
    perilla::syncAllThreads();
    if(perilla::isMasterThread())
    {
#endif

    BL_PROFILE_REGION_STOP(stepName.str());

    cumtime += dt_level[0];

    amr_level[0]->postCoarseTimeStep(cumtime);


    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        run_stop = amrex::second() - run_strt;
	const int istep    = level_steps[0];

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_stop,IOProc);
	amrex::Print() << "\n[STEP " << istep << "] Coarse TimeStep time: " << run_stop << '\n';
#ifdef BL_LAZY
	});
#endif

#ifndef BL_MEM_PROFILING
        long min_fab_kilobytes  = amrex::TotalBytesAllocatedInFabsHWM()/1024;
        long max_fab_kilobytes  = min_fab_kilobytes;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceLongMin(min_fab_kilobytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_kilobytes, IOProc);

	amrex::Print() << "[STEP " << istep << "] FAB kilobyte spread across MPI nodes: ["
		       << min_fab_kilobytes << " ... " << max_fab_kilobytes << "]\n";
#ifdef BL_LAZY
	amrex::Print() << "\n";
	});
#endif
#endif
    }

#ifdef BL_MEM_PROFILING
    {
	std::ostringstream ss;
	ss << "[STEP " << level_steps[0] << "]";
	MemProfiler::report(ss.str());
    }
#endif

    BL_PROFILE_ADD_STEP(level_steps[0]);
    BL_PROFILE_REGION_STOP("Amr::coarseTimeStep()");
    BL_COMM_PROFILE_NAMETAG(stepName.str());
    //BL_PROFILE_FLUSH();
    BL_TRACE_PROFILE_FLUSH();
    BL_COMM_PROFILE_FLUSH();

    if (verbose > 0)
    {
	amrex::Print()
	    << "\nSTEP = " << level_steps[0]
	    << " TIME = "  << cumtime
	    << " DT = "    << dt_level[0] << "\n\n";
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "STEP = "  << level_steps[0]
               << " TIME = " << cumtime
               << " DT = "   << dt_level[0] << '\n';
    }
    if (record_run_info_terse && ParallelDescriptor::IOProcessor())
        runlog_terse << level_steps[0] << " " << cumtime << " " << dt_level[0] << '\n';

    int check_test = 0;

    if (check_per > 0.0)
    {

        // Check to see if we've crossed a check_per interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = (cumtime-dt_level[0]) / check_per;
        int num_per_new = (cumtime            ) / check_per;

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next check_per interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * std::abs(cumtime);
        const Real next_chk_time = (num_per_old + 1) * check_per;

        if ((num_per_new == num_per_old) && std::abs(cumtime - next_chk_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && std::abs((cumtime - dt_level[0]) - next_chk_time) <= eps)
        {
            num_per_old += 1;
        }

        if (num_per_old != num_per_new)
        {
            check_test = 1;
        }

    }

    int to_stop       = 0;    
    int to_checkpoint = 0;
    int to_plot       = 0;
    int to_small_plot = 0;
    if (message_int > 0 && level_steps[0] % message_int == 0) {
	if (ParallelDescriptor::IOProcessor())
	{
	    FILE *fp;
	    if ((fp=fopen("dump_and_continue","r")) != 0)
	    {
		remove("dump_and_continue");
		to_checkpoint = 1;
		fclose(fp);
	    }
	    else if ((fp=fopen("stop_run","r")) != 0)
	    {
		remove("stop_run");
		to_stop = 1;
		fclose(fp);
	    }
	    else if ((fp=fopen("dump_and_stop","r")) != 0)
	    {
		remove("dump_and_stop");
		to_checkpoint = 1;
		to_stop = 1;
		fclose(fp);
	    }

	    if ((fp=fopen("plot_and_continue","r")) != 0)
	    {
		remove("plot_and_continue");
		to_plot = 1;
		fclose(fp);
	    }

            if ((fp=fopen("small_plot_and_continue","r")) != 0)
            {
                remove("small_plot_and_continue");
                to_small_plot = 1;
                fclose(fp);
            }
	}
        int packed_data[4];
	packed_data[0] = to_stop;
	packed_data[1] = to_checkpoint;
        packed_data[2] = to_plot;
        packed_data[3] = to_small_plot;
	ParallelDescriptor::Bcast(packed_data, 4, ParallelDescriptor::IOProcessorNumber());
	to_stop = packed_data[0];
	to_checkpoint = packed_data[1];
        to_plot = packed_data[2];
        to_small_plot = packed_data[3];

    }

    if(to_stop == 1 && to_checkpoint == 0) {  // prevent main from writing files
      last_checkpoint = level_steps[0];
      last_plotfile   = level_steps[0];
    }

    if (to_checkpoint && write_plotfile_with_checkpoint) {
      to_plot = 1;
      to_small_plot = 1;
    }

    if ((check_int > 0 && level_steps[0] % check_int == 0) || check_test == 1
	|| to_checkpoint)
    {
        checkPoint();
    }


    if (writePlotNow() || to_plot)
    {
        writePlotFile();
    }

    if (writeSmallPlotNow() || to_small_plot)
    {
        writeSmallPlotFile();
    }

    updateInSitu();

    bUserStopRequest = to_stop;
    if (to_stop)
    {
        ParallelDescriptor::Barrier("Amr::coarseTimeStep::to_stop");
        if(ParallelDescriptor::IOProcessor()) {
          if (to_checkpoint)
          {
            amrex::ErrorStream() << "Stopped by user w/ checkpoint" << std::endl;
          }
          else
          {
            amrex::ErrorStream() << "Stopped by user w/o checkpoint" << std::endl;
          }
	}
    }


#ifdef USE_PERILLA_PTHREADS
    }
#endif

}

bool
Amr::writePlotNow() noexcept
{
    int plot_test = 0;
    if (plot_per > 0.0)
    {

        // Check to see if we've crossed a plot_per interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = (cumtime-dt_level[0]) / plot_per;
        int num_per_new = (cumtime            ) / plot_per;

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next plot_per interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * std::abs(cumtime);
        const Real next_plot_time = (num_per_old + 1) * plot_per;

        if ((num_per_new == num_per_old) && std::abs(cumtime - next_plot_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && std::abs((cumtime - dt_level[0]) - next_plot_time) <= eps)
        {
            num_per_old += 1;
        }

        if (num_per_old != num_per_new)
        {
            plot_test = 1;
        }

    }

    if (plot_log_per > 0.0)
    {

        // Check to see if we've crossed a plot_log_per interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.
        // This only works when cumtime > 0.

        int num_per_old = 0;
        int num_per_new = 0;

        if (cumtime-dt_level[0] > 0.) {
            num_per_old = log10(cumtime-dt_level[0]) / plot_log_per;
        }
        if (cumtime > 0.) {
            num_per_new = log10(cumtime) / plot_log_per;
        }

        if (num_per_old != num_per_new)
        {
            plot_test = 1;
        }

    }

    return ( (plot_int > 0 && level_steps[0] % plot_int == 0) || 
              plot_test == 1 ||
              amr_level[0]->writePlotNow());
} 

bool
Amr::writeSmallPlotNow() noexcept
{
    int plot_test = 0;
    if (small_plot_per > 0.0)
    {

        // Check to see if we've crossed a small_plot_per interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = (cumtime-dt_level[0]) / small_plot_per;
        int num_per_new = (cumtime            ) / small_plot_per;

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next small_plot_per interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * std::abs(cumtime);
        const Real next_plot_time = (num_per_old + 1) * small_plot_per;

        if ((num_per_new == num_per_old) && std::abs(cumtime - next_plot_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && std::abs((cumtime - dt_level[0]) - next_plot_time) <= eps)
        {
            num_per_old += 1;
        }

        if (num_per_old != num_per_new)
	{
            plot_test = 1;
	}

    }

    if (small_plot_log_per > 0.0)
    {

        // Check to see if we've crossed a small_plot_log_per interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.
        // This only works when cumtime > 0.

        int num_per_old = 0;
        int num_per_new = 0;

        if (cumtime-dt_level[0] > 0.) {
            num_per_old = log10(cumtime-dt_level[0]) / small_plot_log_per;
        }
        if (cumtime > 0.) {
            num_per_new = log10(cumtime) / small_plot_log_per;
        }

        if (num_per_old != num_per_new)
        {
            plot_test = 1;
        }

    }

    return ( (small_plot_int > 0 && level_steps[0] % small_plot_int == 0) || 
              plot_test == 1 ||
              amr_level[0]->writeSmallPlotNow());
} 

void
Amr::defBaseLevel (Real              strt_time, 
                   const BoxArray*   lev0_grids,
                   const Vector<int>* pmap)
{
    BL_PROFILE("Amr::defBaseLevel()");
    // Just initialize this here for the heck of it
    which_level_being_advanced = -1;

    //
    // Check that base domain has even number of zones in all directions.
    //
    const Box& domain   = Geom(0).Domain();
    const IntVect& d_len = domain.size();

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
        if (d_len[idir]%2 != 0)
            amrex::Error("defBaseLevel: must have even number of cells");

    BoxArray lev0;

    if (lev0_grids != 0 && lev0_grids->size() > 0)
    {
        BL_ASSERT(pmap != 0);

        BoxArray domain_ba(domain);
        if (!domain_ba.contains(*lev0_grids))
            amrex::Error("defBaseLevel: domain does not contain lev0_grids!");
        if (!lev0_grids->contains(domain_ba))
            amrex::Error("defBaseLevel: lev0_grids does not contain domain");

        lev0 = *lev0_grids;

	if (refine_grid_layout) {
	    ChopGrids(0,lev0,ParallelDescriptor::NProcs());
	}
    }
    else
    {
	lev0 = MakeBaseGrids();
    }

    this->SetBoxArray(0, lev0);
    this->SetDistributionMap(0, DistributionMapping(lev0));

    //
    // Now build level 0 grids.
    //
    amr_level[0].reset((*levelbld)(*this,0,Geom(0),grids[0],dmap[0],strt_time));
    //
    // Now init level 0 grids with data.
    //
    amr_level[0]->initData();
}

void
Amr::regrid (int  lbase,
             Real time,
             bool initial)
{
    BL_PROFILE("Amr::regrid()");

    if (lbase > std::min(finest_level,max_level-1)) return;

    if (verbose > 0)
	amrex::Print() << "Now regridding at level lbase = " << lbase << "\n";

    //
    // Compute positions of new grids.
    //
    int             new_finest;
    Vector<BoxArray> new_grid_places(max_level+1);
    Vector<DistributionMapping> new_dmap(max_level+1);

    grid_places(lbase,time,new_finest, new_grid_places);

    bool regrid_level_zero = (!initial) && (lbase == 0)
        && ( loadbalance_with_workestimates || (new_grid_places[0] != amr_level[0]->boxArray()));

    const int start = regrid_level_zero ? 0 : lbase+1;

    bool grids_unchanged = finest_level == new_finest;
    for (int lev = start, End = std::min(finest_level,new_finest); lev <= End; lev++) {
	if (new_grid_places[lev] == amr_level[lev]->boxArray()) {
	    new_grid_places[lev] = amr_level[lev]->boxArray();  // to avoid duplicates
	    new_dmap[lev] = amr_level[lev]->DistributionMap(); 
	} else {
	    grids_unchanged = false;
	}
    }

    //
    // If use_efficient_regrid flag is set and grids are unchanged, then don't do anything more here.
    //
    if (use_efficient_regrid == 1 && grids_unchanged )
    {
	if (verbose > 0) {
	    amrex::Print() << "Regridding at level lbase = " << lbase 
			   << " but grids unchanged\n";
	}
	return;
    }

    //
    // Reclaim old-time grid space for all remain levels > lbase.
    //
    for(int lev = start; lev <= finest_level; ++lev) {
	amr_level[lev]->removeOldData();
    }
    //
    // Reclaim all remaining storage for levels > new_finest.
    //
    for(int lev = new_finest + 1; lev <= finest_level; ++lev) {
	amr_level[lev].reset();
	this->ClearBoxArray(lev);
	this->ClearDistributionMap(lev);
    }

    finest_level = new_finest;

    //
    // Define the new grids from level start up to new_finest.
    //
    for(int lev = start; lev <= new_finest; ++lev) {
        //
        // Construct skeleton of new level.
        //

        if (loadbalance_with_workestimates && !initial) {
            new_dmap[lev] = makeLoadBalanceDistributionMap(lev, time, new_grid_places[lev]);
        }
        else if (new_dmap[lev].empty()) {
	    new_dmap[lev].define(new_grid_places[lev]);
	}

        AmrLevel* a = (*levelbld)(*this,lev,Geom(lev),new_grid_places[lev],
				  new_dmap[lev],cumtime);

        if (initial)
        {
            //
            // We're being called on startup from bldFineLevels().
            // NOTE: The initData function may use a filPatch, and so needs to
            //       be officially inserted into the hierarchy prior to the call.
            //
            amr_level[lev].reset(a);
	    this->SetBoxArray(lev, amr_level[lev]->boxArray());
	    this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
            amr_level[lev]->initData();
        }
        else if (amr_level[lev])
        {
            //
            // Init with data from old structure then remove old structure.
            // NOTE: The init function may use a filPatch from the old level,
            //       which therefore needs remain in the hierarchy during the call.
            //
            a->init(*amr_level[lev]);
            amr_level[lev].reset(a);
	    this->SetBoxArray(lev, amr_level[lev]->boxArray());
	    this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
	}
        else
        {
            a->init();
            amr_level[lev].reset(a);
	    this->SetBoxArray(lev, amr_level[lev]->boxArray());
	    this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
        }

    }


    //
    // Check at *all* levels whether we need to do anything special now that the grids
    //       at levels lbase+1 and higher may have changed.  
    //
    for(int lev(0); lev <= new_finest; ++lev) {
      amr_level[lev]->post_regrid(lbase,new_finest);
    }

    //
    // Report creation of new grids.
    //

    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "REGRID: at level lbase = " << lbase << '\n';
        printGridInfo(runlog,start,finest_level);
    }

    if (record_grid_info && ParallelDescriptor::IOProcessor())
    {
        if (lbase == 0)
            gridlog << "STEP = " << level_steps[0] << ' ';

        gridlog << "TIME = "
                << time
                << " : REGRID  with lbase = "
                << lbase
                << '\n';

        printGridInfo(gridlog,start,finest_level);
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        if (lbase == 0) {
            amrex::Print() << "STEP = " << level_steps[0] << ' ';
        }

        amrex::Print() << "TIME = "
                       << time
                       << " : REGRID  with lbase = "
                       << lbase
                       << std::endl;

        if (verbose > 1)
        {
            printGridInfo(amrex::OutStream(),start,finest_level);
        }
        else
        {
            printGridSummary(amrex::OutStream(),start,finest_level);
        }
    }
}

DistributionMapping
Amr::makeLoadBalanceDistributionMap (int lev, Real time, const BoxArray& ba) const
{
    BL_PROFILE("makeLoadBalanceDistributionMap()");

    if (verbose) {
        amrex::Print() << "Load balance on level " << lev << " at t = " << time << "\n";
    }

    DistributionMapping newdm;

    const int work_est_type = amr_level[0]->WorkEstType();

    if (work_est_type < 0) {
        if (verbose) {
            amrex::Print() << "\nAMREX WARNING: work estimates type does not exist!\n\n";
        }
        newdm.define(ba);
    }
    else if (amr_level[lev])
    {
        DistributionMapping dmtmp;
        if (ba.size() == boxArray(lev).size()) {
            dmtmp = DistributionMap(lev);
        } else {
            dmtmp.define(ba);
        }

        MultiFab workest(ba, dmtmp, 1, 0, MFInfo(), FArrayBoxFactory());
        AmrLevel::FillPatch(*amr_level[lev], workest, 0, time, work_est_type, 0, 1, 0);

        Real navg = static_cast<Real>(ba.size()) / static_cast<Real>(ParallelDescriptor::NProcs());
        int nmax = std::max(std::round(loadbalance_max_fac*navg), std::ceil(navg));

        newdm = DistributionMapping::makeKnapSack(workest, nmax);
    }
    else
    {
        newdm.define(ba);
    }

    return newdm;
}

void
Amr::LoadBalanceLevel0 (Real time)
{
    BL_PROFILE("LoadBalanceLevel0()");
    const auto& dm = makeLoadBalanceDistributionMap(0, time, boxArray(0));
    InstallNewDistributionMap(0, dm);
    amr_level[0]->post_regrid(0,time);
}

void
Amr::InstallNewDistributionMap (int lev, const DistributionMapping& newdm)
{
    BL_PROFILE("InstallNewDistributionMap()");

    AmrLevel* a = (*levelbld)(*this,lev,Geom(lev),boxArray(lev),newdm,cumtime);
    a->init(*amr_level[lev]);
    amr_level[lev].reset(a);

    this->SetBoxArray(lev, amr_level[lev]->boxArray());
    this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
}

void
Amr::regrid_level_0_on_restart()
{
    regrid_on_restart = 0;
    //
    // Coarsening before we split the grids ensures that each resulting
    // grid will have an even number of cells in each direction.
    //
    BoxArray lev0(amrex::coarsen(Geom(0).Domain(),2));
    //
    // Now split up into list of grids within max_grid_size[0] limit.
    //
    lev0.maxSize(max_grid_size[0]/2);
    //
    // Now refine these boxes back to level 0.
    //
    lev0.refine(2);
    
    //
    // If use_efficient_regrid flag is set, then test to see whether we in fact 
    //    have just changed the level 0 grids. If not, then don't do anything more here.
    //
    if ( !( (use_efficient_regrid == 1) && (lev0 == amr_level[0]->boxArray()) ) ) 
    {
	//
	// Construct skeleton of new level.
	//
	DistributionMapping dm(lev0);
	AmrLevel* a = (*levelbld)(*this,0,Geom(0),lev0,dm,cumtime);
	
	a->init(*amr_level[0]);
	amr_level[0].reset(a);
	
	this->SetBoxArray(0, amr_level[0]->boxArray());
	this->SetDistributionMap(0, amr_level[0]->DistributionMap());

	amr_level[0]->post_regrid(0,0);
	
	if (ParallelDescriptor::IOProcessor())
	{
	    if (verbose > 1)
	    {
		printGridInfo(amrex::OutStream(),0,finest_level);
	    }
	    else if (verbose > 0)
	    {
		printGridSummary(amrex::OutStream(),0,finest_level);
	    }
	}
	
	if (record_grid_info && ParallelDescriptor::IOProcessor())
	    printGridInfo(gridlog,0,finest_level);
    }
    else
    {
	if (verbose > 0)
	    amrex::Print() << "Regridding at level 0 but grids unchanged \n";
    }
}

void
Amr::printGridInfo (std::ostream& os,
                    int           min_lev,
                    int           max_lev)
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray&           bs      = amr_level[lev]->boxArray();
        int                       numgrid = bs.size();
        long                      ncells  = amr_level[lev]->countCells();
        double                    ntot    = Geom(lev).Domain().d_numPts();
        Real                      frac    = 100.0*(Real(ncells) / ntot);
        const DistributionMapping& map    = amr_level[lev]->get_new_data(0).DistributionMap();

        os << "  Level "
           << lev
           << "   "
           << numgrid
           << " grids  "
           << ncells
           << " cells  "
           << frac
           << " % of domain"
           << '\n';


        for (int k = 0; k < numgrid; k++)
        {
            const Box& b = bs[k];

            os << ' ' << lev << ": " << b << "   ";
                
            for (int i = 0; i < AMREX_SPACEDIM; i++)
                os << b.length(i) << ' ';

            os << ":: " << map[k] << '\n';
        }
    }

    os << std::endl; // Make sure we flush!
}


void
Amr::grid_places (int              lbase,
                  Real             time,
                  int&             new_finest,
                  Vector<BoxArray>& new_grids)
{
    BL_PROFILE("Amr::grid_places()");

    const Real strttime = amrex::second();

    if (lbase == 0)
    {
	new_grids[0] = MakeBaseGrids();
    }

    if ( time == 0. && !initial_grids_file.empty() && !use_fixed_coarse_grids)
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min<int>(new_finest,initial_ba.size());

        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = initial_ba[lev-1].size();
            for (int i = 0; i < ngrid; i++)
            {
                Box bx(initial_ba[lev-1][i]);
                if (lev > lbase)
                    bl.push_back(bx);
            }
            if (lev > lbase)
                new_grids[lev].define(bl);
        }
        return;
    }

    // Use grids in initial_grids_file as fixed coarse grids.
    if ( ! initial_grids_file.empty() && use_fixed_coarse_grids)
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min<int>(new_finest,initial_ba.size());

        for (int lev = lbase+1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = initial_ba[lev-1].size();
            for (int i = 0; i < ngrid; i++)
            {
                Box bx(initial_ba[lev-1][i]);

                if (lev > lbase)
                    bl.push_back(bx);

            }
            if (lev > lbase)
                new_grids[lev].define(bl);
            new_grids[lev].maxSize(max_grid_size[lev]);
        }
    }
    else if ( !regrid_grids_file.empty() )     // Use grids in regrid_grids_file 
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min<int>(new_finest,regrid_ba.size());
        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = regrid_ba[lev-1].size();
            for (int i = 0; i < ngrid; i++)
            {
                Box bx(regrid_ba[lev-1][i]);
                if (lev > lbase)
                    bl.push_back(bx);
            }
            if (lev > lbase)
                new_grids[lev].define(bl);
        }
        return;
    }

    MakeNewGrids(lbase, time, new_finest, new_grids);

    if (verbose > 0)
    {
        Real stoptime = amrex::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
	amrex::Print() << "grid_places() time: " << stoptime << " new finest: " << new_finest<< '\n';
#ifdef BL_LAZY
	});
#endif
    }
}

void
Amr::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    amr_level[lev]->errorEst(tags,TagBox::CLEAR,TagBox::SET,time, n_error_buf[lev][0], ngrow);
}

BoxArray
Amr::GetAreaNotToTag (int lev)
{
    return BoxArray(amr_level[lev]->getAreaNotToTag());
}

void
Amr::ManualTagsPlacement (int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev)
{
    amr_level[lev]->manual_tags_placement(tags, bf_lev);
}

void
Amr::bldFineLevels (Real strt_time)
{
    BL_PROFILE("Amr::bldFineLevels()");
    finest_level = 0;

    Vector<BoxArray> new_grids(max_level+1);
    //
    // Get initial grid placement.
    //
    do
    {
        int new_finest;

        grid_places(finest_level,strt_time,new_finest,new_grids);

        if (new_finest <= finest_level) break;
        //
        // Create a new level and link with others.
        //
        finest_level = new_finest;

	DistributionMapping new_dm {new_grids[new_finest]};

        AmrLevel* level = (*levelbld)(*this,
                                      new_finest,
                                      Geom(new_finest),
                                      new_grids[new_finest],
				      new_dm,
                                      strt_time);

        amr_level[new_finest].reset(level);
	this->SetBoxArray(new_finest, new_grids[new_finest]);
	this->SetDistributionMap(new_finest, new_dm);

        amr_level[new_finest]->initData();
    }
    while (finest_level < max_level);
    //
    // Iterate grids to ensure fine grids encompass all interesting gunk.
    //     but only iterate if we did not provide a grids file.
    //
    if ( regrid_grids_file.empty() || (strt_time == 0.0 && !initial_grids_file.empty()) )  
    {
	bool grids_the_same;

	const int MaxCnt = 4;

	int count = 0;

	do
	{
	    for (int i = 0; i <= finest_level; i++) {
		new_grids[i] = amr_level[i]->boxArray();
	    }

	    regrid(0,strt_time,true);

	    grids_the_same = true;

	    for (int i = 0; i <= finest_level && grids_the_same; i++) {
		if (!(new_grids[i] == amr_level[i]->boxArray())) {
		    grids_the_same = false;
		}
	    }

	    count++;
	}
	while (!grids_the_same && count < MaxCnt);
    }
}

void
Amr::initSubcycle ()
{
    BL_PROFILE("Amr::initSubcycle()");
    ParmParse pp("amr");
    sub_cycle = true;
    if (pp.contains("nosub"))
    {
        if (verbose) {
            amrex::Print() << "Warning: The nosub flag has been deprecated.\n "
                           << "... please use subcycling_mode to control subcycling.\n";
        }
        int nosub;
        pp.query("nosub",nosub);
        if (nosub > 0)
            sub_cycle = false;
        else
            amrex::Error("nosub <= 0 not allowed.\n");
        subcycling_mode = "None";
    }
    else 
    {
        subcycling_mode = "Auto";
        pp.query("subcycling_mode",subcycling_mode);
    }
    
    if (subcycling_mode == "None")
    {
        sub_cycle = false;
        for (int i = 0; i <= max_level; i++)
        {
            n_cycle[i] = 1;
        }
    }
    else if (subcycling_mode == "Manual")
    {
        int cnt = pp.countval("subcycling_iterations");

        if (cnt == 1)
        {
            //
            // Set all values to the single available value.
            //
            int cycles = 0;

            pp.get("subcycling_iterations",cycles);

            n_cycle[0] = 1; // coarse level is always 1 cycle
            for (int i = 1; i <= max_level; i++)
            {
                n_cycle[i] = cycles;
            }
        }
        else if (cnt > 1)
        {
            //
            // Otherwise we expect a vector of max_grid_size values.
            //
            pp.getarr("subcycling_iterations",n_cycle,0,max_level+1);
            if (n_cycle[0] != 1)
            {
                amrex::Error("First entry of subcycling_iterations must be 1");
            }
        }
        else
        {
            amrex::Error("Must provide a valid subcycling_iterations if mode is Manual");
        }
        for (int i = 1; i <= max_level; i++)
        {
            if (n_cycle[i] > MaxRefRatio(i-1))
                amrex::Error("subcycling iterations must always be <= ref_ratio");
            if (n_cycle[i] <= 0)
                amrex::Error("subcycling iterations must always be > 0");
        }
    }
    else if (subcycling_mode == "Auto")
    {
        n_cycle[0] = 1;
        for (int i = 1; i <= max_level; i++)
        {
            n_cycle[i] = MaxRefRatio(i-1);
        } 
    }
    else if (subcycling_mode == "Optimal")
    {
        // if subcycling mode is Optimal, n_cycle is set dynamically.
        // We'll initialize it to be Auto subcycling.
        n_cycle[0] = 1;
        for (int i = 1; i <= max_level; i++)
        {
            n_cycle[i] = MaxRefRatio(i-1);
        } 
    }
    else
    {
        std::string err_message = "Unrecognzied subcycling mode: " + subcycling_mode + "\n";
        amrex::Error(err_message.c_str());
    }
}

void
Amr::initPltAndChk ()
{
    ParmParse pp("amr");

    pp.query("checkpoint_files_output", checkpoint_files_output);
    pp.query("plot_files_output", plot_files_output);

    pp.query("plot_nfiles", plot_nfiles);
    pp.query("checkpoint_nfiles", checkpoint_nfiles);
    //
    // -1 ==> use ParallelDescriptor::NProcs().
    //
    if (plot_nfiles       == -1) plot_nfiles       = ParallelDescriptor::NProcs();
    if (checkpoint_nfiles == -1) checkpoint_nfiles = ParallelDescriptor::NProcs();
    
    check_file_root = "chk";
    pp.query("check_file",check_file_root);

    check_int = -1;
    pp.query("check_int",check_int);

    check_per = -1.0;
    pp.query("check_per",check_per);

    if (check_int > 0 && check_per > 0)
    {
        if (ParallelDescriptor::IOProcessor())
	    amrex::Warning("Warning: both amr.check_int and amr.check_per are > 0.");
    }

    plot_file_root = "plt";
    pp.query("plot_file",plot_file_root);

    plot_int = -1;
    pp.query("plot_int",plot_int);

    plot_per = -1.0;
    pp.query("plot_per",plot_per);

    plot_log_per = -1.0;
    pp.query("plot_log_per",plot_log_per);

    if (plot_int > 0 && plot_per > 0)
    {
        if (ParallelDescriptor::IOProcessor())
            amrex::Warning("Warning: both amr.plot_int and amr.plot_per are > 0.");
    }

    small_plot_file_root = "smallplt";
    pp.query("small_plot_file",small_plot_file_root);

    small_plot_int = -1;
    pp.query("small_plot_int",small_plot_int);

    small_plot_per = -1.0;
    pp.query("small_plot_per",small_plot_per);

    small_plot_log_per = -1.0;
    pp.query("small_plot_log_per",small_plot_log_per);

    if (small_plot_int > 0 && small_plot_per > 0)
    {
        if (ParallelDescriptor::IOProcessor())
            amrex::Warning("Warning: both amr.small_plot_int and amr.small_plot_per are > 0.");
    }

    write_plotfile_with_checkpoint = 1;
    pp.query("write_plotfile_with_checkpoint",write_plotfile_with_checkpoint);

    stream_max_tries = 4;
    pp.query("stream_max_tries",stream_max_tries);
    stream_max_tries = std::max(stream_max_tries, 1);

    abort_on_stream_retry_failure = false;
    pp.query("abort_on_stream_retry_failure",abort_on_stream_retry_failure);

    pp.query("precreateDirectories", precreateDirectories);
    pp.query("prereadFAHeaders", prereadFAHeaders);

    int phvInt(plot_headerversion), chvInt(checkpoint_headerversion);
    pp.query("plot_headerversion", phvInt);
    if(phvInt != plot_headerversion) {
      plot_headerversion = static_cast<VisMF::Header::Version> (phvInt);
    }
    pp.query("checkpoint_headerversion", chvInt);
    if(chvInt != checkpoint_headerversion) {
      checkpoint_headerversion = static_cast<VisMF::Header::Version> (chvInt);
    }
}


bool
Amr::okToRegrid(int level) noexcept
{
    if (regrid_int[level] < 0)
        return false;
    else
        return level_count[level] >= regrid_int[level] && amr_level[level]->okToRegrid();
}

Real
Amr::computeOptimalSubcycling(int n, int* best, Real* dt_max, Real* est_work, int* cycle_max)
{
    BL_ASSERT(cycle_max[0] == 1);
    // internally these represent the total number of steps at a level, 
    // not the number of cycles
    std::vector<int> cycles(n);
    Real best_ratio = 1e200;
    Real best_dt = 0;
    Real ratio;
    Real dt;
    Real work;
    int limit = 1;
    // This provides a memory efficient way to test all candidates
    for (int i = 1; i < n; i++)
        limit *= cycle_max[i];
    for (int candidate = 0; candidate < limit; candidate++)
    {
        int temp_cand = candidate;
        cycles[0] = 1;
        dt = dt_max[0];
        work = est_work[0];
        for (int i  = 1; i < n; i++)
        {
            // grab the relevant "digit" and shift over.
            cycles[i] = (1 + temp_cand%cycle_max[i]) * cycles[i-1];
            temp_cand /= cycle_max[i];
            dt = std::min(dt, cycles[i]*dt_max[i]);
            work += cycles[i]*est_work[i];
        }
        ratio = work/dt;
        if (ratio < best_ratio) 
        {
            for (int i  = 0; i < n; i++)
                best[i] = cycles[i];
            best_ratio = ratio;
            best_dt = dt;
        }
    }
    //
    // Now we convert best back to n_cycles format
    //
    for (int i = n-1; i > 0; i--)
        best[i] /= best[i-1];
    return best_dt;
}

const Vector<BoxArray>& Amr::getInitialBA() noexcept
{
  return initial_ba;
}

#ifdef AMREX_PARTICLES
void 
Amr::RedistributeParticles ()
{
    amr_level[0]->particle_redistribute(0,true);
}
#endif

}

