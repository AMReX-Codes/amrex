
#include <winstd.H>
#include <algorithm>
#include <cstdio>
#include <list>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <iomanip>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <Geometry.H>
#include <TagBox.H>
#include <Array.H>
#include <CoordSys.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <Cluster.H>
#include <LevelBld.H>
#include <AmrLevel.H>
#include <PROB_AMR_F.H>
#include <Amr.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <DistributionMapping.H>
#include <FabSet.H>

#ifdef MG_USE_FBOXLIB
#include <mg_cpp_f.h>
#endif

#ifdef USE_PARTICLES
#include <AmrParGDB.H>
#endif

#ifdef BL_LAZY
#include <Lazy.H>
#endif

#ifdef BL_MEM_PROFILING
#include <MemProfiler.H>
#endif

#ifdef BL_USE_ARRAYVIEW
#include <DatasetClient.H>
#endif
//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif
//
// Static class members.  Set defaults in Initialize()!!!
//
std::list<std::string> Amr::state_plot_vars;
std::list<std::string> Amr::state_small_plot_vars;
std::list<std::string> Amr::derive_plot_vars;
bool                   Amr::first_plotfile;
bool                   Amr::first_smallplotfile;
Array<BoxArray>        Amr::initial_ba;
Array<BoxArray>        Amr::regrid_ba;
bool                   Amr::useFixedCoarseGrids;
int                    Amr::useFixedUpToLevel;


namespace
{
    const std::string CheckPointVersion("CheckPointVersion_1.0");

    bool initialized = false;
}

namespace
{
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
    bool refine_grid_layout;
    int  plotfile_on_restart;
    int  checkpoint_on_restart;
    bool checkpoint_files_output;
    int  compute_new_dt_on_regrid;
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
    probinit_natonce         = 32;
    plot_files_output        = true;
    checkpoint_nfiles        = 64;
    regrid_on_restart        = 0;
    use_efficient_regrid     = 0;
    refine_grid_layout       = true;
    plotfile_on_restart      = 0;
    checkpoint_on_restart    = 0;
    checkpoint_files_output  = true;
    compute_new_dt_on_regrid = 0;
    Amr::useFixedCoarseGrids = false;
    Amr::useFixedUpToLevel   = 0;

    BoxLib::ExecOnFinalize(Amr::Finalize);

    VisMF::Initialize();

    initialized = true;
}

void
Amr::Finalize ()
{
    Amr::state_plot_vars.clear();
    Amr::derive_plot_vars.clear();
    Amr::regrid_ba.clear();
    Amr::initial_ba.clear();

    initialized = false;
}

bool Amr::Plot_Files_Output () { return plot_files_output; }

std::ostream&
Amr::DataLog (int i)
{
    return datalog[i];
}

int
Amr::NumDataLogs ()
{
    return datalog.size();
}

bool
Amr::RegridOnRestart () const
{
    return regrid_on_restart;
}

const BoxArray&
Amr::boxArray (int lev) const
{
    return amr_level[lev].boxArray();
}

#ifdef USE_PARTICLES
const BoxArray&
Amr::ParticleBoxArray (int lev) const
{
    return amr_level[lev].ParticleBoxArray();
}
#endif

void
Amr::setDtMin (const Array<Real>& dt_min_in)
{
    for (int i = 0; i <= finest_level; i++)
        dt_min[i] = dt_min_in[i];
}

PArray<AmrLevel>&
Amr::getAmrLevels ()
{
    return amr_level;
}

long
Amr::cellCount (int lev)
{
    return amr_level[lev].countCells();
}

int
Amr::numGrids (int lev)
{
    return amr_level[lev].numGrids();
}

MultiFab*
Amr::derive (const std::string& name,
             Real               time,
             int                lev,
             int                ngrow)
{
    return amr_level[lev].derive(name,time,ngrow);
}

int
Amr::MaxRefRatio (int level) const
{
    int maxval = 0;
    for (int n = 0; n<BL_SPACEDIM; n++) 
        maxval = std::max(maxval,ref_ratio[level][n]);
    return maxval;
}

Amr::Amr ()
    :
    amr_level(PArrayManage),
    datalog(PArrayManage)
{
    Initialize();
    Geometry::Setup();
    int max_level_in = -1;
    Array<int> n_cell_in(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++) n_cell_in[i] = -1;
    InitAmr(max_level_in,n_cell_in);
}

Amr::Amr (const RealBox* rb, int max_level_in, Array<int> n_cell_in, int coord)
    :
    amr_level(PArrayManage),
    datalog(PArrayManage)
{
    Initialize();
    Geometry::Setup(rb,coord);
    InitAmr(max_level_in,n_cell_in);
}

void
Amr::InitAmr (int max_level_in, Array<int> n_cell_in)
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
    grid_eff               = 0.7;
    plot_int               = -1;
    small_plot_int         = -1;
    n_proper               = 1;
    max_level              = -1;
    last_plotfile          = 0;
    last_smallplotfile     = 0;
    last_checkpoint        = 0;
    record_run_info        = false;
    record_grid_info       = false;
    file_name_digits       = 5;
    record_run_info_terse  = false;
    bUserStopRequest       = false;
    message_int            = 10;
    
    int i;
    for (i = 0; i < BL_SPACEDIM; i++)
        isPeriodic[i] = false;

    ParmParse pp("amr");
    //
    // Check for command line flags.
    //
    verbose = 0;
    pp.query("v",verbose);

    pp.query("regrid_on_restart",regrid_on_restart);
    pp.query("use_efficient_regrid",use_efficient_regrid);
    pp.query("plotfile_on_restart",plotfile_on_restart);
    pp.query("checkpoint_on_restart",checkpoint_on_restart);

    pp.query("compute_new_dt_on_regrid",compute_new_dt_on_regrid);

    pp.query("refine_grid_layout", refine_grid_layout);

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
    //
    // Read max_level and alloc memory for container objects.
    //
    if (max_level_in == -1) 
       pp.get("max_level", max_level);
    else
       max_level = max_level_in;

    int nlev     = max_level+1;
    geom.resize(nlev);
    dt_level.resize(nlev);
    level_steps.resize(nlev);
    level_count.resize(nlev);
    n_cycle.resize(nlev);
    dt_min.resize(nlev);
    blocking_factor.resize(nlev);
    max_grid_size.resize(nlev);
    n_error_buf.resize(nlev);
    amr_level.resize(nlev);
    //
    // Set bogus values.
    //
    for (i = 0; i < nlev; i++)
    {
        dt_level[i]    = 1.e200; // Something nonzero so old & new will differ
        level_steps[i] = 0;
        level_count[i] = 0;
        n_cycle[i]     = 0;
        dt_min[i]      = 0.0;
        n_error_buf[i] = 1;
        blocking_factor[i] = 8;
        max_grid_size[i] = (BL_SPACEDIM == 2) ? 128 : 32;
    }

    // Make the default regrid_int = 1 for all levels.
    if (max_level > 0) 
    {
       regrid_int.resize(max_level);
       for (i = 0; i < max_level; i++)
           regrid_int[i]  = 1;
    }

    // Make the default ref_ratio = 2 for all levels.
    ref_ratio.resize(max_level);
    for (i = 0; i < max_level; i++)
        ref_ratio[i] = 2 * IntVect::TheUnitVector();

    //
    // Read other amr specific values.
    //
    pp.query("n_proper",n_proper);
    pp.query("grid_eff",grid_eff);
    pp.queryarr("n_error_buf",n_error_buf,0,max_level);
    //
    // Read in the refinement ratio IntVects as integer BL_SPACEDIM-tuples.
    //
    if (max_level > 0)
    {
        const int nratios_vect = max_level*BL_SPACEDIM;

        Array<int> ratios_vect(nratios_vect);

        int got_vect = pp.queryarr("ref_ratio_vect",ratios_vect,0,nratios_vect);

        Array<int> ratios(max_level);

        const int got_int = pp.queryarr("ref_ratio",ratios,0,max_level);
   
        if (got_int == 1 && got_vect == 1 && ParallelDescriptor::IOProcessor())
        {
            BoxLib::Warning("Only input *either* ref_ratio or ref_ratio_vect");
        }
        else if (got_vect == 1)
        {
            int k = 0;
            for (i = 0; i < max_level; i++)
            {
                for (int n = 0; n < BL_SPACEDIM; n++,k++)
                    ref_ratio[i][n] = ratios_vect[k];
            }
        }
        else if (got_int == 1)
        {
            for (i = 0; i < max_level; i++)
            {
                for (int n = 0; n < BL_SPACEDIM; n++)
                    ref_ratio[i][n] = ratios[i];
            }
        }
        else
        {
            if (ParallelDescriptor::IOProcessor())
                BoxLib::Warning("Using default ref_ratio = 2 at all levels");
        }
    }

    //
    // Are we going to keep the grids at certain levels fixed?
    //
    pp.query("useFixedCoarseGrids", useFixedCoarseGrids);
    
    //
    // Setup plot and checkpoint controls.
    //
    initPltAndChk();
    
    //
    // Setup subcycling controls.
    //
    initSubcycle();

    //
    // Read in max_grid_size.  Use defaults if not explicitly defined.
    //
    int cnt = pp.countval("max_grid_size");

    if (cnt == 1)
    {
        //
        // Set all values to the single available value.
        //
        int the_max_grid_size = 0;

        pp.get("max_grid_size",the_max_grid_size);

        for (i = 0; i <= max_level; i++)
        {
            max_grid_size[i] = the_max_grid_size;
        }
    }
    else if (cnt > 1)
    {
        //
        // Otherwise we expect a vector of max_grid_size values.
        //
        pp.getarr("max_grid_size",max_grid_size,0,max_level+1);
    }
    //
    // Read in the blocking_factors.  Use defaults if not explicitly defined.
    //
    cnt = pp.countval("blocking_factor");

    if (cnt == 1)
    {
        //
        // Set all values to the single available value.
        //
        int the_blocking_factor = 0;

        pp.get("blocking_factor",the_blocking_factor);

        for (i = 0; i <= max_level; i++)
        {
            blocking_factor[i] = the_blocking_factor;
        }
    }
    else if (cnt > 1)
    {
        //
        // Otherwise we expect a vector of blocking factors.
        //
        pp.getarr("blocking_factor",blocking_factor,0,max_level+1);
    }
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
           for (i = 0; i < max_level; i++)
           {
               regrid_int[i] = the_regrid_int;
           }
       }
       else if (numvals == 0)
       {
           if (ParallelDescriptor::IOProcessor())
               BoxLib::Warning("Using default regrid_int = 1 at all levels");
       }
       else if (numvals < max_level)
       {
           BoxLib::Error("You did not specify enough values of regrid_int");
       }
       else 
       {
           //
           // Otherwise we expect a vector of max_level values
           //
           pp.queryarr("regrid_int",regrid_int,0,max_level);
       }
    }
    //
    // Read computational domain and set geometry.
    //
    Array<int> n_cell(BL_SPACEDIM);
    if (n_cell_in[0] == -1)
    {
       pp.getarr("n_cell",n_cell,0,BL_SPACEDIM);
    }
    else
    {
       for (i = 0; i < BL_SPACEDIM; i++) n_cell[i] = n_cell_in[i];
    }
    
    IntVect lo(IntVect::TheZeroVector()), hi(n_cell);
    hi -= IntVect::TheUnitVector();
    Box index_domain(lo,hi);
    for (i = 0; i <= max_level; i++)
    {
        geom[i].define(index_domain);
        if (i < max_level)
            index_domain.refine(ref_ratio[i]);
    }
    //
    // Now define offset for CoordSys.
    //
    Real offset[BL_SPACEDIM];
    for (i = 0; i < BL_SPACEDIM; i++)
    {
        const Real delta = Geometry::ProbLength(i)/(Real)n_cell[i];
        offset[i]        = Geometry::ProbLo(i) + delta*lo[i];
    }
    CoordSys::SetOffset(offset);

    if (max_level > 0 && !initial_grids_file.empty())
    {
#define STRIP while( is.get() != '\n' )
        std::ifstream is(initial_grids_file.c_str(),std::ios::in);

        if (!is.good())
            BoxLib::FileOpenFailed(initial_grids_file);

        int in_finest,ngrid;

        is >> in_finest;
        STRIP;
        initial_ba.resize(in_finest);

        useFixedUpToLevel = in_finest;
        if (in_finest > max_level)
           BoxLib::Error("You have fewer levels in your inputs file then in your grids file!");

        for (int lev = 1; lev <= in_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (i = 0; i < ngrid; i++)
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
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Read initial_ba. Size is " << initial_ba.size() << std::endl;

#undef STRIP
    }

    if (max_level > 0 && !regrid_grids_file.empty())
    {
#define STRIP while( is.get() != '\n' )
        std::ifstream is(regrid_grids_file.c_str(),std::ios::in);

        if (!is.good())
            BoxLib::FileOpenFailed(regrid_grids_file);

        int in_finest,ngrid;

        is >> in_finest;
        STRIP;
        regrid_ba.resize(in_finest);
        for (int lev = 1; lev <= in_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (i = 0; i < ngrid; i++)
            {
                Box bx;
                is >> bx;
                STRIP;
                 bx.refine(ref_ratio[lev-1]);
                 if (bx.longside() > max_grid_size[lev])
                 {
                     std::cout << "Grid " << bx << " too large" << '\n';
                     BoxLib::Error();
                 }
                 bl.push_back(bx);
            }
            regrid_ba[lev-1].define(bl);
        }
        is.close();
#undef STRIP
    }

    rebalance_grids = 0;
    pp.query("rebalance_grids", rebalance_grids);

#ifdef USE_PARTICLES
    m_gdb = new AmrParGDB(this);
#endif
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
Amr::isDerivePlotVar (const std::string& name)
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
Amr::clearDerivePlotVarList ()
{
    derive_plot_vars.clear();
}

void
Amr::addDerivePlotVar (const std::string& name)
{
    if (!isDerivePlotVar(name))
        derive_plot_vars.push_back(name);
}

void
Amr::deleteDerivePlotVar (const std::string& name)
{
    if (isDerivePlotVar(name))
        derive_plot_vars.remove(name);
}

Amr::~Amr ()
{
    levelbld->variableCleanUp();

#ifdef USE_PARTICLES
    delete m_gdb;
#endif

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
            BoxLib::FileOpenFailed(filename);
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
            BoxLib::FileOpenFailed(filename);
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
            BoxLib::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordRunInfoTerse");
}

void
Amr::setRecordDataInfo (int i, const std::string& filename)
{
    if (ParallelDescriptor::IOProcessor())
    {
        datalog.set(i,new std::fstream);
        datalog[i].open(filename.c_str(),std::ios::out|std::ios::app);
        if (!datalog[i].good())
            BoxLib::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier("Amr::setRecordDataInfo");
}

void
Amr::setDtLevel (const Array<Real>& dt_lev)
{
    for (int i = 0; i <= finest_level; i++)
        dt_level[i] = dt_lev[i];
}

void
Amr::setDtLevel (Real dt, int lev)
{
    dt_level[lev] = dt;
}

void
Amr::setNCycle (const Array<int>& ns)
{
    for (int i = 0; i <= finest_level; i++)
        n_cycle[i] = ns[i];
}

long
Amr::cellCount ()
{
    long cnt = 0;
    for (int i = 0; i <= finest_level; i++)
        cnt += amr_level[i].countCells();
    return cnt;
}

int
Amr::numGrids ()
{
    int cnt = 0;
    for (int i = 0; i <= finest_level; i++)
        cnt += amr_level[i].numGrids();
    return cnt;
}

int
Amr::okToContinue ()
{
    int ok = true;
    for (int i = 0; ok && (i <= finest_level); i++)
        ok = ok && amr_level[i].okToContinue();
    if(bUserStopRequest) {
      ok = false;
    }
    return ok;
}

void
Amr::writePlotFile ()
{
    if ( ! Plot_Files_Output()) return;

    BL_PROFILE_REGION_START("Amr::writePlotFile()");
    BL_PROFILE("Amr::writePlotFile()");

    VisMF::SetNOutFiles(plot_nfiles);

    if (first_plotfile) 
    {
        first_plotfile = false;
        amr_level[0].setPlotVariables();
    }

    Real dPlotFileTime0 = ParallelDescriptor::second();

    const std::string& pltfile = BoxLib::Concatenate(plot_file_root,level_steps[0],file_name_digits);

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "PLOTFILE: file = " << pltfile << '\n';

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "PLOTFILE: file = " << pltfile << '\n';

  BoxLib::StreamRetry sretry(pltfile, abort_on_stream_retry_failure,
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
    BoxLib::UtilRenameDirectoryToOld(pltfile, false);     // dont call barrier
    BoxLib::UtilCreateCleanDirectory(pltfileTemp, true);  // call barrier

    std::string HeaderFileName(pltfileTemp + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec(0);

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out | std::ios::trunc |
	                                        std::ios::binary);
        if ( ! HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);
        old_prec = HeaderFile.precision(15);
    }

    for (int k(0); k <= finest_level; ++k)
        amr_level[k].writePlotFile(pltfileTemp, HeaderFile);

    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.precision(old_prec);
        if ( ! HeaderFile.good())
            BoxLib::Error("Amr::writePlotFile() failed");
    }

    last_plotfile = level_steps[0];

    if (verbose > 0)
    {
        const int IOProc        = ParallelDescriptor::IOProcessorNumber();
        Real      dPlotFileTime = ParallelDescriptor::second() - dPlotFileTime0;

        ParallelDescriptor::ReduceRealMax(dPlotFileTime,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Write plotfile time = " << dPlotFileTime << "  seconds" << "\n\n";
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
  BL_PROFILE_REGION_STOP("Amr::writePlotFile()");
}

void
Amr::writeSmallPlotFile ()
{
    if ( ! Plot_Files_Output()) return;

    BL_PROFILE_REGION_START("Amr::writeSmallPlotFile()");
    BL_PROFILE("Amr::writeSmallPlotFile()");

    VisMF::SetNOutFiles(plot_nfiles);

    if (first_smallplotfile) 
    {
        first_smallplotfile = false;
        amr_level[0].setSmallPlotVariables();
    }

    // Don't continue if we have no variables to plot.
    
    if (stateSmallPlotVars().size() == 0) return;

    Real dPlotFileTime0 = ParallelDescriptor::second();

    const std::string& pltfile = BoxLib::Concatenate(small_plot_file_root,level_steps[0],file_name_digits);

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "SMALL PLOTFILE: file = " << pltfile << '\n';

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "SMALL PLOTFILE: file = " << pltfile << '\n';

  BoxLib::StreamRetry sretry(pltfile, abort_on_stream_retry_failure,
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
    BoxLib::UtilRenameDirectoryToOld(pltfile, false);     // dont call barrier
    BoxLib::UtilCreateCleanDirectory(pltfileTemp, true);  // call barrier

    std::string HeaderFileName(pltfileTemp + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec(0);

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out | std::ios::trunc |
	                                        std::ios::binary);
        if ( ! HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);
        old_prec = HeaderFile.precision(15);
    }

    for (int k(0); k <= finest_level; ++k)
        amr_level[k].writeSmallPlotFile(pltfileTemp, HeaderFile);

    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.precision(old_prec);
        if ( ! HeaderFile.good())
            BoxLib::Error("Amr::writePlotFile() failed");
    }

    last_smallplotfile = level_steps[0];

    if (verbose > 0)
    {
        const int IOProc        = ParallelDescriptor::IOProcessorNumber();
        Real      dPlotFileTime = ParallelDescriptor::second() - dPlotFileTime0;

        ParallelDescriptor::ReduceRealMax(dPlotFileTime,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Write small plotfile time = " << dPlotFileTime << "  seconds" << "\n\n";
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
  BL_PROFILE_REGION_STOP("Amr::writeSmallPlotFile()");
}

void
Amr::checkInput ()
{
    if (max_level < 0)
        BoxLib::Error("checkInput: max_level not set");
    //
    // Check that blocking_factor is a power of 2.
    //
    for (int i = 0; i < max_level; i++)
    {
        int k = blocking_factor[i];
        while ( k > 0 && (k%2 == 0) )
            k /= 2;
        if (k != 1)
            BoxLib::Error("Amr::checkInputs: blocking_factor not power of 2");
    }
    //
    // Check level dependent values.
    //
    int i;
    for (i = 0; i < max_level; i++)
    {
        if (MaxRefRatio(i) < 2 || MaxRefRatio(i) > 12)
            BoxLib::Error("checkInput bad ref_ratios");
    }
    const Box& domain = geom[0].Domain();
    if (!domain.ok())
        BoxLib::Error("level 0 domain bad or not set");
    //
    // Check that domain size is a multiple of blocking_factor[0].
    //
    for (i = 0; i < BL_SPACEDIM; i++)
    {
        int len = domain.length(i);
        if (len%blocking_factor[0] != 0)
            BoxLib::Error("domain size not divisible by blocking_factor");
    }
    //
    // Check that max_grid_size is even.
    //
    for (i = 0; i < max_level; i++)
    {
        if (max_grid_size[i]%2 != 0)
            BoxLib::Error("max_grid_size is not even");
    }

    //
    // Check that max_grid_size is a multiple of blocking_factor at every level.
    //
    for (i = 0; i < max_level; i++)
    {
        if (max_grid_size[i]%blocking_factor[i] != 0)
            BoxLib::Error("max_grid_size not divisible by blocking_factor");
    }

    if( ! Geometry::ProbDomain().ok()) {
        BoxLib::Error("checkInput: bad physical problem size");
    }

    if(verbose > 0 && ParallelDescriptor::IOProcessor()) {
       std::cout << "Successfully read inputs file ... " << '\n';
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
        if(plot_int > 0 || plot_per > 0) {
            writePlotFile();
	}
	if (small_plot_int > 0 || small_plot_per > 0)
	    writeSmallPlotFile();
    }
#ifdef HAS_XGRAPH
    if (first_plotfile)
    {
        first_plotfile = false;
        amr_level[0].setPlotVariables();
    }
#endif

#ifdef BL_COMM_PROFILING
    Array<Box> probDomain(geom.size());
    for(int i(0); i < probDomain.size(); ++i) {
      probDomain[i] = geom[i].Domain();
    }
    BL_COMM_PROFILE_INITAMR(finest_level, max_level, ref_ratio, probDomain);
#endif
    BL_PROFILE_REGION_STOP("Amr::init()");
}

void
Amr::readProbinFile (int& init)
{
    BL_PROFILE("Amr::readProbinFile()");
    //
    // Populate integer array with name of probin file.
    //
    int probin_file_length = probin_file.length();

    Array<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
        probin_file_name[i] = probin_file[i];

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
       std::cout << "Starting to read probin ... " << std::endl;

    const int nAtOnce = probinit_natonce;
    const int MyProc  = ParallelDescriptor::MyProc();
    const int NProcs  = ParallelDescriptor::NProcs();
    const int NSets   = (NProcs + (nAtOnce - 1)) / nAtOnce;
    const int MySet   = MyProc/nAtOnce;

    Real piStart = 0, piEnd = 0, piStartAll = ParallelDescriptor::second();

    for (int iSet = 0; iSet < NSets; ++iSet)
    {
        if (MySet == iSet)
        {
            //
            // Call the pesky probin reader.
            //
            piStart = ParallelDescriptor::second();

#ifdef DIMENSION_AGNOSTIC

            FORT_PROBINIT(&init,
                          probin_file_name.dataPtr(),
                          &probin_file_length,
                          ZFILL(Geometry::ProbLo()),
                          ZFILL(Geometry::ProbHi()));

#else

            FORT_PROBINIT(&init,
                          probin_file_name.dataPtr(),
                          &probin_file_length,
                          Geometry::ProbLo(),
                          Geometry::ProbHi());
#endif

            piEnd = ParallelDescriptor::second();
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
        Real      piTotalAll = ParallelDescriptor::second() - piStartAll;

        ParallelDescriptor::ReduceRealMax(piTotal,    IOProc);
        ParallelDescriptor::ReduceRealMax(piTotalAll, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "MFRead::: PROBINIT max time   = " << piTotal    << '\n';
            std::cout << "MFRead::: PROBINIT total time = " << piTotalAll << '\n';
        }
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "Successfully read probin file: \"" << probin_file << "\"\n";
}

void
Amr::initialInit (Real              strt_time,
                  Real              stop_time,
                  const BoxArray*   lev0_grids,
                  const Array<int>* pmap)
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
                    const Array<int>* pmap)
{
    BL_PROFILE("Amr::InitializeInit()");
    BL_COMM_PROFILE_NAMETAG("Amr::InitializeInit TOP");
    checkInput();
    //
    // Generate internal values from user-supplied values.
    //
    finest_level = 0;
    //
    // Init problem dependent data.
    //
    int init = true;

    if (!probin_file.empty()) {
        readProbinFile(init);
    }

#ifdef BL_SYNC_RANTABLES
    int iGet(0), iSet(1);
    const int iTableSize(64);
    Real *RanAmpl = new Real[iTableSize];
    Real *RanPhase = new Real[iTableSize];
    FORT_SYNC_RANTABLES(RanPhase, RanAmpl, &iGet);
    ParallelDescriptor::Bcast(RanPhase, iTableSize);
    ParallelDescriptor::Bcast(RanAmpl, iTableSize);
    FORT_SYNC_RANTABLES(RanPhase, RanAmpl, &iSet);
    delete [] RanAmpl;
    delete [] RanPhase;
#endif

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
    amr_level[0].computeInitialDt(finest_level,
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

    for (int lev = 1; lev <= max_level; lev++)
    {
        dt0           /= n_cycle[lev];
        dt_level[lev]  = dt0;
        dt_min[lev]    = dt_level[lev];
    }

    if (max_level > 0)
        bldFineLevels(strt_time);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].setTimeLevel(strt_time,dt_level[lev],dt_level[lev]);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_regrid(0,finest_level);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        level_steps[lev] = 0;
        level_count[lev] = 0;
    }

    //
    // Perform any special post_initialization operations.
    //
    for(int lev(0); lev <= finest_level; ++lev) {
      amr_level[lev].post_init(stop_time);
    }

    if (ParallelDescriptor::IOProcessor())
    {
       if (verbose > 1)
       {
           std::cout << "INITIAL GRIDS \n";
           printGridInfo(std::cout,0,finest_level);
       }
       else if (verbose > 0)
       { 
           std::cout << "INITIAL GRIDS \n";
           printGridSummary(std::cout,0,finest_level);
       }
    }

    if (record_grid_info && ParallelDescriptor::IOProcessor())
    {
        gridlog << "INITIAL GRIDS \n";
        printGridInfo(gridlog,0,finest_level);
    }

#ifdef USE_STATIONDATA 
    station.init(amr_level, finestLevel());
    station.findGrid(amr_level,geom);
#endif
    BL_COMM_PROFILE_NAMETAG("Amr::initialInit BOTTOM");
}

void
Amr::restart (const std::string& filename)
{
    BL_PROFILE_REGION_START("Amr::restart()");
    BL_PROFILE("Amr::restart()");

    // Just initialize this here for the heck of it
    which_level_being_advanced = -1;

    Real dRestartTime0 = ParallelDescriptor::second();

    DistributionMapping::Initialize();

    VisMF::SetMFFileInStreams(mffile_nstreams);

    int i;

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "restarting calculation from file: " << filename << std::endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';
    //
    // Init problem dependent data.
    //
    int init = false;

    readProbinFile(init);
    //
    // Start calculation from given restart file.
    //
    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';
    //
    // Open the checkpoint header file for reading.
    //
    std::string File = filename;

    File += '/';
    File += "Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    Array<char> fileCharPtr;
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

    if (spdim != BL_SPACEDIM)
    {
        std::cerr << "Amr::restart(): bad spacedim = " << spdim << '\n';
        BoxLib::Abort();
    }

    is >> cumtime;
    int mx_lev;
    is >> mx_lev;
    is >> finest_level;

    Array<Box> inputs_domain(max_level+1);
    for (int lev = 0; lev <= max_level; lev++)
    {
       Box bx(geom[lev].Domain().smallEnd(),geom[lev].Domain().bigEnd());
       inputs_domain[lev] = bx;
    }

    if (max_level >= mx_lev) {

       for (i = 0; i <= mx_lev; i++) is >> geom[i];
       for (i = 0; i <  mx_lev; i++) is >> ref_ratio[i];
       for (i = 0; i <= mx_lev; i++) is >> dt_level[i];

       if (new_checkpoint_format)
       {
           for (i = 0; i <= mx_lev; i++) is >> dt_min[i];
       }
       else
       {
           for (i = 0; i <= mx_lev; i++) dt_min[i] = dt_level[i];
       }

       Array<int>  n_cycle_in;
       n_cycle_in.resize(mx_lev+1);  
       for (i = 0; i <= mx_lev; i++) is >> n_cycle_in[i];
       bool any_changed = false;

       for (i = 0; i <= mx_lev; i++) 
           if (n_cycle[i] != n_cycle_in[i])
           {
               any_changed = true;
               if (verbose > 0 && ParallelDescriptor::IOProcessor())
                   std::cout << "Warning: n_cycle has changed at level " << i << 
                                " from " << n_cycle_in[i] << " to " << n_cycle[i] << std::endl;;
           }

       // If we change n_cycle then force a full regrid from level 0 up
       if (max_level > 0 && any_changed)
       {
           level_count[0] = regrid_int[0];
           if ((verbose > 0) && ParallelDescriptor::IOProcessor())
               std::cout << "Warning: This forces a full regrid " << std::endl;
       }


       for (i = 0; i <= mx_lev; i++) is >> level_steps[i];
       for (i = 0; i <= mx_lev; i++) is >> level_count[i];

       //
       // Set bndry conditions.
       //
       if (max_level > mx_lev)
       {
           for (i = mx_lev+1; i <= max_level; i++)
           {
               dt_level[i]    = dt_level[i-1]/n_cycle[i];
               level_steps[i] = n_cycle[i]*level_steps[i-1];
               level_count[i] = 0;
           }

           // This is just an error check
           if (!sub_cycle)
           {
               for (i = 1; i <= finest_level; i++)
               {
                   if (dt_level[i] != dt_level[i-1])
                      BoxLib::Error("restart: must have same dt at all levels if not subcycling");
               }
           }
       }

       if (regrid_on_restart && max_level > 0)
       {
           if (regrid_int[0] > 0) 
               level_count[0] = regrid_int[0];
           else
               BoxLib::Error("restart: can't have regrid_on_restart and regrid_int <= 0");
       }

       checkInput();
       //
       // Read levels.
       //
       int lev;
       for (lev = 0; lev <= finest_level; lev++)
       {
           amr_level.set(lev,(*levelbld)());
           amr_level[lev].restart(*this, is);
       }
       //
       // Build any additional data structures.
       //
       for (lev = 0; lev <= finest_level; lev++)
           amr_level[lev].post_restart();

    } else {

       if (ParallelDescriptor::IOProcessor())
          BoxLib::Warning("Amr::restart(): max_level is lower than before");

       int new_finest_level = std::min(max_level,finest_level);

       finest_level = new_finest_level;
 
       // These are just used to hold the extra stuff we have to read in.
       Geometry   geom_dummy;
       Real       real_dummy;
       int         int_dummy;
       IntVect intvect_dummy;

       for (i = 0          ; i <= max_level; i++) is >> geom[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> geom_dummy;

       for (i = 0        ; i <  max_level; i++) is >> ref_ratio[i];
       for (i = max_level; i <  mx_lev   ; i++) is >> intvect_dummy;

       for (i = 0          ; i <= max_level; i++) is >> dt_level[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> real_dummy;

       if (new_checkpoint_format)
       {
           for (i = 0          ; i <= max_level; i++) is >> dt_min[i];
           for (i = max_level+1; i <= mx_lev   ; i++) is >> real_dummy;
       }
       else
       {
           for (i = 0; i <= max_level; i++) dt_min[i] = dt_level[i];
       }

       for (i = 0          ; i <= max_level; i++) is >> n_cycle[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> int_dummy;

       for (i = 0          ; i <= max_level; i++) is >> level_steps[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> int_dummy;

       for (i = 0          ; i <= max_level; i++) is >> level_count[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> int_dummy;

       if (regrid_on_restart && max_level > 0)
       {
           if (regrid_int[0] > 0) 
               level_count[0] = regrid_int[0];
           else
               BoxLib::Error("restart: can't have regrid_on_restart and regrid_int <= 0");
       }

       checkInput();

       //
       // Read levels.
       //
       int lev;
       for (lev = 0; lev <= new_finest_level; lev++)
       {
           amr_level.set(lev,(*levelbld)());
           amr_level[lev].restart(*this, is);
       }
       //
       // Build any additional data structures.
       //
       for (lev = 0; lev <= new_finest_level; lev++)
           amr_level[lev].post_restart();

    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
       Box restart_domain(geom[lev].Domain());
       if (! (inputs_domain[lev] == restart_domain) )
       {
          if (ParallelDescriptor::IOProcessor())
          {
             std::cout << "Problem at level " << lev << '\n';
             std::cout << "Domain according to     inputs file is " <<  inputs_domain[lev] << '\n';
             std::cout << "Domain according to checkpoint file is " << restart_domain      << '\n';
             std::cout << "Amr::restart() failed -- box from inputs file does not equal box from restart file" << std::endl;
          }
          BoxLib::Abort();
       }
    }

#ifdef USE_STATIONDATA
    station.init(amr_level, finestLevel());
    station.findGrid(amr_level,geom);
#endif

    if (verbose > 0)
    {
        Real dRestartTime = ParallelDescriptor::second() - dRestartTime0;

        ParallelDescriptor::ReduceRealMax(dRestartTime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Restart time = " << dRestartTime << " seconds." << '\n';
    }
    BL_PROFILE_REGION_STOP("Amr::restart()");
}

void
Amr::checkPoint ()
{
    if ( ! checkpoint_files_output) return;

    BL_PROFILE_REGION_START("Amr::checkPoint()");
    BL_PROFILE("Amr::checkPoint()");

    VisMF::SetNOutFiles(checkpoint_nfiles);
    //
    // In checkpoint files always write out FABs in NATIVE format.
    //
    FABio::Format thePrevFormat = FArrayBox::getFormat();

    FArrayBox::setFormat(FABio::FAB_NATIVE);

    Real dCheckPointTime0 = ParallelDescriptor::second();

    const std::string& ckfile = BoxLib::Concatenate(check_file_root,level_steps[0],file_name_digits);

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "CHECKPOINT: file = " << ckfile << std::endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "CHECKPOINT: file = " << ckfile << '\n';


  BoxLib::StreamRetry sretry(ckfile, abort_on_stream_retry_failure,
                             stream_max_tries);

  const std::string ckfileTemp(ckfile + ".temp");

  while(sretry.TryFileOutput()) {
    //
    //  if either the ckfile or ckfileTemp exists, rename them
    //  to move them out of the way.  then create ckfile
    //  with the temporary name, then rename it back when
    //  it is finished writing.  then stream retry can rename
    //  it to a bad suffix if there were stream errors.
    //
    BoxLib::UtilRenameDirectoryToOld(ckfile, false);     // dont call barrier
    BoxLib::UtilCreateCleanDirectory(ckfileTemp, true);  // call barrier

    std::string HeaderFileName = ckfileTemp + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec = 0, i;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out | std::ios::trunc |
	                                        std::ios::binary);

        if ( ! HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);

        old_prec = HeaderFile.precision(17);

        HeaderFile << CheckPointVersion << '\n'
                   << BL_SPACEDIM       << '\n'
                   << cumtime           << '\n'
                   << max_level         << '\n'
                   << finest_level      << '\n';
        //
        // Write out problem domain.
        //
        for (i = 0; i <= max_level; i++) HeaderFile << geom[i]        << ' ';
        HeaderFile << '\n';
        for (i = 0; i < max_level; i++)  HeaderFile << ref_ratio[i]   << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << dt_level[i]    << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << dt_min[i]      << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << n_cycle[i]     << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << level_steps[i] << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << level_count[i] << ' ';
        HeaderFile << '\n';
    }

    for (i = 0; i <= finest_level; ++i)
        amr_level[i].checkPoint(ckfileTemp, HeaderFile);

    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.precision(old_prec);

        if ( ! HeaderFile.good())
            BoxLib::Error("Amr::checkpoint() failed");
    }

    last_checkpoint = level_steps[0];

#ifdef USE_SLABSTAT
    //
    // Dump out any SlabStats MultiFabs.
    //
    AmrLevel::get_slabstat_lst().checkPoint(getAmrLevels(), level_steps[0]);
#endif

    if (verbose > 0)
    {
        Real dCheckPointTime = ParallelDescriptor::second() - dCheckPointTime0;

        ParallelDescriptor::ReduceRealMax(dCheckPointTime,
	                            ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "checkPoint() time = " << dCheckPointTime << " secs." << '\n';
    }
    ParallelDescriptor::Barrier("Amr::checkPoint::end");

    if(ParallelDescriptor::IOProcessor()) {
      std::rename(ckfileTemp.c_str(), ckfile.c_str());
    }
    ParallelDescriptor::Barrier("Renaming temporary checkPoint file.");

  }  // end while

  //
  // Don't forget to reset FAB format.
  //
  FArrayBox::setFormat(thePrevFormat);

  BL_PROFILE_REGION_STOP("Amr::checkPoint()");
}

void
Amr::RegridOnly (Real time)
{
    BL_ASSERT(regrid_on_restart == 1);

    int lev_top = std::min(finest_level, max_level-1);

    for (int i = 0; i <= lev_top; i++)
       regrid(i,time);

    if (plotfile_on_restart)
	writePlotFile();

    if (checkpoint_on_restart)
       checkPoint();

}

void
Amr::timeStep (int  level,
               Real time,
               int  iteration,
               int  niter,
               Real stop_time)
{
    BL_PROFILE("Amr::timeStep()");
    BL_COMM_PROFILE_NAMETAG("Amr::timeStep TOP");

    // This is used so that the AmrLevel functions can know which level is being advanced 
    //      when regridding is called with possible lbase > level.
    which_level_being_advanced = level;

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

        for (int i(level); i <= lev_top; ++i)
        {
            const int old_finest = finest_level;

            if (okToRegrid(i))
            {
                regrid(i,time);

                //
                // Compute new dt after regrid if at level 0 and compute_new_dt_on_regrid.
                //
                if ( compute_new_dt_on_regrid && (i == 0) )
                {
                    int post_regrid_flag = 1;
                    amr_level[0].computeNewDt(finest_level,
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
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
	std::cout << "[Level " << level 
		  << " step " << level_steps[level]+1 << "] ";
        std::cout << "ADVANCE with dt = "
                  << dt_level[level]
                  << std::endl;
    }
    BL_PROFILE_REGION_START("amr_level.advance");
    Real dt_new = amr_level[level].advance(time,dt_level[level],iteration,niter);
    BL_PROFILE_REGION_STOP("amr_level.advance");

    dt_min[level] = iteration == 1 ? dt_new : std::min(dt_min[level],dt_new);

    level_steps[level]++;
    level_count[level]++;

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
	std::cout << "[Level " << level
		  << " step " << level_steps[level] << "] ";
        std::cout << "Advanced "
                  << amr_level[level].countCells()
                  << " cells"
                  << std::endl;
    }

#ifdef USE_STATIONDATA
    station.report(time+dt_level[level],level,amr_level[level]);
#endif

#ifdef USE_SLABSTAT
    AmrLevel::get_slabstat_lst().update(amr_level[level],time,dt_level[level]);
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

    amr_level[level].post_timestep(iteration);

    // Set this back to negative so we know whether we are in fact in this routine
    which_level_being_advanced = -1;
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
    BL_PROFILE_REGION_START("Amr::coarseTimeStep()");
    BL_PROFILE("Amr::coarseTimeStep()");
    std::stringstream stepName;
    stepName << "timeStep STEP " << level_steps[0];

    const Real run_strt = ParallelDescriptor::second() ;
    //
    // Compute new dt.
    //
    if (levelSteps(0) > 0)
    {
        int post_regrid_flag = 0;
        amr_level[0].computeNewDt(finest_level,
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
        amr_level[0].computeInitialDt(finest_level,
                                      sub_cycle,
                                      n_cycle,
                                      ref_ratio,
                                      dt_level,
                                      stop_time);
    }

    BL_PROFILE_REGION_START(stepName.str());

    timeStep(0,cumtime,1,1,stop_time);

    BL_PROFILE_REGION_STOP(stepName.str());

    cumtime += dt_level[0];

    amr_level[0].postCoarseTimeStep(cumtime);

#ifdef BL_PROFILING
#ifdef DEBUG
    std::stringstream dfss;
    dfss << "BytesPerProc.STEP_" << std::setw(5) << std::setfill('0')
         << level_steps[0] - 1 << ".xgr";
    DistributionMapping::PrintDiagnostics(dfss.str());
#endif
#endif

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_stop = ParallelDescriptor::second() - run_strt;
	const int istep    = level_steps[0];

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_stop,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "\n[STEP " << istep << "] Coarse TimeStep time: " << run_stop << '\n' ;
#ifdef BL_LAZY
	});
#endif

#ifndef BL_MEM_PROFILING
        long min_fab_kilobytes  = BoxLib::TotalBytesAllocatedInFabsHWM()/1024;
        long max_fab_kilobytes  = min_fab_kilobytes;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceLongMin(min_fab_kilobytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_kilobytes, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "[STEP " << istep << "] FAB kilobyte spread across MPI nodes: ["
                      << min_fab_kilobytes
                      << " ... "
                      << max_fab_kilobytes
                      << "]\n";
        }
#ifdef BL_LAZY
	if (ParallelDescriptor::IOProcessor()) std::cout << "\n";
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
    BL_TRACE_PROFILE_FLUSH();
    BL_COMM_PROFILE_NAMETAG(stepName.str());
    BL_COMM_PROFILE_FLUSH();

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nSTEP = "
                 << level_steps[0]
                  << " TIME = "
                  << cumtime
                  << " DT = "
                  << dt_level[0]
                  << '\n'
                  << std::endl;
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "STEP = "
               << level_steps[0]
               << " TIME = "
               << cumtime
               << " DT = "
               << dt_level[0]
               << '\n';
    }
    if (record_run_info_terse && ParallelDescriptor::IOProcessor())
        runlog_terse << level_steps[0] << " " << cumtime << " " << dt_level[0] << '\n';

    int check_test = 0;

    if (check_per > 0.0)
    {
      const int num_per_old = (cumtime-dt_level[0]) / check_per;
      const int num_per_new = (cumtime            ) / check_per;

      if (num_per_old != num_per_new)
      {
	check_test = 1;
      }
    }

    int to_stop       = 0;    
    int to_checkpoint = 0;
    int to_plot       = 0;
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
	}
	int packed_data[2];
	packed_data[0] = to_stop;
	packed_data[1] = to_checkpoint;
	ParallelDescriptor::Bcast(packed_data, 2, ParallelDescriptor::IOProcessorNumber());
	to_stop = packed_data[0];
	to_checkpoint = packed_data[1];

    }

    if(to_stop == 1 && to_checkpoint == 0) {  // prevent main from writing files
      last_checkpoint = level_steps[0];
      last_plotfile   = level_steps[0];
    }

    if (to_checkpoint && write_plotfile_with_checkpoint)
      to_plot = 1;

    if ((check_int > 0 && level_steps[0] % check_int == 0) || check_test == 1
	|| to_checkpoint)
    {
        checkPoint();
    }


    if (writePlotNow() || to_plot)
    {
        writePlotFile();
    }

    if (writeSmallPlotNow())
    {
        writeSmallPlotFile();
    }

    bUserStopRequest = to_stop;
    if (to_stop)
    {
        ParallelDescriptor::Barrier("Amr::coarseTimeStep::to_stop");
        if(ParallelDescriptor::IOProcessor()) {
          if (to_checkpoint)
          {
            std::cerr << "Stopped by user w/ checkpoint" << std::endl;
          }
          else
          {
            std::cerr << "Stopped by user w/o checkpoint" << std::endl;
          }
	}
    }
}

bool
Amr::writePlotNow()
{
    int plot_test = 0;
    if (plot_per > 0.0)
    {
#ifdef BL_USE_NEWPLOTPER
      Real rN(0.0);
      Real rR = modf(cumtime/plot_per, &rN);
      if (rR < (dt_level[0]*0.001))
#else
      const int num_per_old = (cumtime-dt_level[0]) / plot_per;
      const int num_per_new = (cumtime            ) / plot_per;

      if (num_per_old != num_per_new)
#endif
	{
	  plot_test = 1;
	}
    }

    return ( (plot_int > 0 && level_steps[0] % plot_int == 0) || 
              plot_test == 1 ||
              amr_level[0].writePlotNow());
} 

bool
Amr::writeSmallPlotNow()
{
    int plot_test = 0;
    if (small_plot_per > 0.0)
    {
#ifdef BL_USE_NEWPLOTPER
      Real rN(0.0);
      Real rR = modf(cumtime/small_plot_per, &rN);
      if (rR < (dt_level[0]*0.001))
#else
      const int num_per_old = (cumtime-dt_level[0]) / small_plot_per;
      const int num_per_new = (cumtime            ) / small_plot_per;

      if (num_per_old != num_per_new)
#endif
	{
	  plot_test = 1;
	}
    }

    return ( (small_plot_int > 0 && level_steps[0] % small_plot_int == 0) || 
              plot_test == 1 ||
              amr_level[0].writeSmallPlotNow());
} 

void
Amr::defBaseLevel (Real              strt_time, 
                   const BoxArray*   lev0_grids,
                   const Array<int>* pmap)
{
    BL_PROFILE("Amr::defBaseLevel()");
    // Just initialize this here for the heck of it
    which_level_being_advanced = -1;

    //
    // Check that base domain has even number of zones in all directions.
    //
    const Box& domain   = geom[0].Domain();
    IntVect    d_length = domain.size();
    const int* d_len    = d_length.getVect();

    for (int idir = 0; idir < BL_SPACEDIM; idir++)
        if (d_len[idir]%2 != 0)
            BoxLib::Error("defBaseLevel: must have even number of cells");

    BoxArray lev0(1);

    if (lev0_grids != 0 && lev0_grids->size() > 0)
    {
        BL_ASSERT(pmap != 0);

        BoxArray domain_ba(domain);
        if (!domain_ba.contains(*lev0_grids))
            BoxLib::Error("defBaseLevel: domain does not contain lev0_grids!");
        if (!lev0_grids->contains(domain_ba))
            BoxLib::Error("defBaseLevel: lev0_grids does not contain domain");

        lev0 = *lev0_grids;
        //
        // Make sure that the grids all go on the processor as specified by pmap;
        //      we store this in cache so all level 0 MultiFabs will use this DistributionMap.
        //
        const bool put_in_cache = true;
        DistributionMapping dmap(*pmap,put_in_cache);
    }
    else
    {
        //
        // Coarsening before we split the grids ensures that each resulting
        // grid will have an even number of cells in each direction.
        //
        lev0.set(0,BoxLib::coarsen(domain,2));
        //
        // Now split up into list of grids within max_grid_size[0] limit.
        //
        lev0.maxSize(max_grid_size[0]/2);
        //
        // Now refine these boxes back to level 0.
        //
        lev0.refine(2);
    }

    //
    // If (refine_grid_layout == 1) and (Nprocs > ngrids) then break up the 
    //    grids into smaller chunks
    //
    Array<BoxArray> new_grids(1);
    new_grids[0] = lev0;
    impose_refine_grid_layout(0,0,new_grids);
    lev0 = new_grids[0];

    //
    // Now build level 0 grids.
    //
    amr_level.set(0,(*levelbld)(*this,0,geom[0],lev0,strt_time));

    //
    // Now init level 0 grids with data.
    //
    amr_level[0].initData();
}

void
Amr::regrid (int  lbase,
             Real time,
             bool initial)
{
    BL_PROFILE("Amr::regrid()");

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "Now regridding at level lbase = " << lbase << std::endl;

    //
    // Compute positions of new grids.
    //
    int             new_finest;
    Array<BoxArray> new_grid_places(max_level+1);

    if (lbase <= std::min(finest_level,max_level-1)) {
      grid_places(lbase,time,new_finest, new_grid_places);
    }

    bool regrid_level_zero =
        (lbase == 0 && new_grid_places[0] != amr_level[0].boxArray()) && (!initial);

    const int start = regrid_level_zero ? 0 : lbase+1;

    bool grids_unchanged = finest_level == new_finest;
    for (int lev = start, End = std::min(finest_level,new_finest); lev <= End; lev++) {
	if (new_grid_places[lev] == amr_level[lev].boxArray()) {
	    new_grid_places[lev] = amr_level[lev].boxArray();  // to avoid duplicates
	} else {
	    grids_unchanged = false;
	}
    }

    //
    // If use_efficient_regrid flag is set and grids are unchanged, then don't do anything more here.
    //
    if (use_efficient_regrid == 1 && !regrid_level_zero && grids_unchanged )
    {
	if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
	    std::cout << "Regridding at level lbase = " << lbase 
		      << " but grids unchanged " << std::endl;
	}
	return;
    }

    //
    // Reclaim old-time grid space for all remain levels > lbase.
    //
    for(int lev = start; lev <= finest_level; ++lev) {
      amr_level[lev].removeOldData();
    }
    //
    // Reclaim all remaining storage for levels > new_finest.
    //
    for(int lev = new_finest + 1; lev <= finest_level; ++lev) {
      amr_level.clear(lev);
    }

    finest_level = new_finest;
    //
    // Flush the caches.
    //
//    MultiFab::flushFBCache();  no need to flush these
//    Geometry::flushFPBCache();
//    FabArrayBase::flushCPCache();
    DistributionMapping::FlushCache();
#ifdef MG_USE_FBOXLIB
    mgt_flush_copyassoc_cache();
#endif

    //
    // Define the new grids from level start up to new_finest.
    //
    for(int lev = start; lev <= new_finest; ++lev) {
        //
        // Construct skeleton of new level.
        //

        AmrLevel* a = (*levelbld)(*this,lev,geom[lev],new_grid_places[lev],cumtime);

        if (initial)
        {
            //
            // We're being called on startup from bldFineLevels().
            // NOTE: The initData function may use a filPatch, and so needs to
            //       be officially inserted into the hierarchy prior to the call.
            //
            amr_level.clear(lev);
            amr_level.set(lev,a);
            amr_level[lev].initData();
        }
        else if (amr_level.defined(lev))
        {
            //
            // Init with data from old structure then remove old structure.
            // NOTE: The init function may use a filPatch from the old level,
            //       which therefore needs remain in the hierarchy during the call.
            //
            a->init(amr_level[lev]);
            amr_level.clear(lev);
            amr_level.set(lev,a);
       }
        else
        {
            a->init();
            amr_level.clear(lev);
            amr_level.set(lev,a);
        }
    }


    //
    // Check at *all* levels whether we need to do anything special now that the grids
    //       at levels lbase+1 and higher may have changed.  
    //
    for(int lev(0); lev <= new_finest; ++lev) {
      amr_level[lev].post_regrid(lbase,new_finest);
    }

    if(rebalance_grids > 0) {
      DistributionMapping::InitProximityMap();
      DistributionMapping::Initialize();

        Array<BoxArray> allBoxes(amr_level.size());
	for(int ilev(0); ilev < allBoxes.size(); ++ilev) {
	  allBoxes[ilev] = boxArray(ilev);
	}
        Array<Array<int> > mLDM;
	if(rebalance_grids == 1) {
          mLDM = DistributionMapping::MultiLevelMapPFC(ref_ratio, allBoxes, maxGridSize(0));
	} else if(rebalance_grids == 2) {
          mLDM = DistributionMapping::MultiLevelMapRandom(ref_ratio, allBoxes, maxGridSize(0));
	} else if(rebalance_grids == 3) {
          mLDM = DistributionMapping::MultiLevelMapKnapSack(ref_ratio, allBoxes, maxGridSize(0));
	} else if(rebalance_grids == 4) {  // ---- move all grids to proc zero
          mLDM = DistributionMapping::MultiLevelMapRandom(ref_ratio, allBoxes, maxGridSize(0), 0);
	} else {
	}

        for(int iMap(0); iMap < mLDM.size(); ++iMap) {
          MultiFab::MoveAllFabs(mLDM[iMap]);
        }
	Geometry::flushFPBCache();
    }

#ifdef USE_STATIONDATA
    station.findGrid(amr_level,geom);
#endif
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
        if (lbase == 0)
            std::cout << "STEP = " << level_steps[0] << ' ';

        std::cout << "TIME = "
                  << time
                  << " : REGRID  with lbase = "
                  << lbase
                  << std::endl;

        if (verbose > 1)
        {
           printGridInfo(std::cout,start,finest_level);
        }
        else
        {
           printGridSummary(std::cout,start,finest_level);
        }
    }
}

void
Amr::regrid_level_0_on_restart()
{
    regrid_on_restart = 0;
    //
    // Coarsening before we split the grids ensures that each resulting
    // grid will have an even number of cells in each direction.
    //
    BoxArray lev0(BoxLib::coarsen(geom[0].Domain(),2));
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
    if ( !( (use_efficient_regrid == 1) && (lev0 == amr_level[0].boxArray()) ) ) 
    {
	//
	// Construct skeleton of new level.
	//
	AmrLevel* a = (*levelbld)(*this,0,geom[0],lev0,cumtime);
	
	a->init(amr_level[0]);
	amr_level.clear(0);
	amr_level.set(0,a);
	
	amr_level[0].post_regrid(0,0);
	
	if (ParallelDescriptor::IOProcessor())
	{
	    if (verbose > 1)
	    {
		printGridInfo(std::cout,0,finest_level);
	    }
	    else if (verbose > 0)
	    {
		printGridSummary(std::cout,0,finest_level);
	    }
	}
	
	if (record_grid_info && ParallelDescriptor::IOProcessor())
	    printGridInfo(gridlog,0,finest_level);
    }
    else
    {
	if (verbose > 0 && ParallelDescriptor::IOProcessor())
	    std::cout << "Regridding at level 0 but grids unchanged " << std::endl;
    }
}

void
Amr::printGridInfo (std::ostream& os,
                    int           min_lev,
                    int           max_lev)
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray&           bs      = amr_level[lev].boxArray();
        int                       numgrid = bs.size();
        long                      ncells  = amr_level[lev].countCells();
        double                    ntot    = geom[lev].Domain().d_numPts();
        Real                      frac    = 100.0*(Real(ncells) / ntot);
        const DistributionMapping& map    = amr_level[lev].get_new_data(0).DistributionMap();

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
                
            for (int i = 0; i < BL_SPACEDIM; i++)
                os << b.length(i) << ' ';

            os << ":: " << map[k] << '\n';
        }
    }

    os << std::endl; // Make sure we flush!
}

void
Amr::printGridSummary (std::ostream& os,
                       int           min_lev,
                       int           max_lev)
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray&           bs      = amr_level[lev].boxArray();
        int                       numgrid = bs.size();
        long                      ncells  = amr_level[lev].countCells();
        double                    ntot    = geom[lev].Domain().d_numPts();
        Real                      frac    = 100.0*(Real(ncells) / ntot);

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

	if (numgrid > 1) {
	    long vmin = std::numeric_limits<long>::max();
	    long vmax = -1;
	    int lmax = -1;
	    int smin = std::numeric_limits<int>::max();
	    int imax, imin;
#ifdef _OPENMP
#pragma omp parallel
#endif	    
	    {
		long vmin_this = std::numeric_limits<long>::max();
		long vmax_this = -1;
		int lmax_this = -1;
		int smin_this = std::numeric_limits<int>::max();
		int imax_this, imin_this;
#ifdef _OPENMP
#pragma omp for
#endif	    	    
		for (int k = 0; k < numgrid; k++) {
		    const Box& bx = bs[k];
		    long v = bx.volume();
		    int ss = bx.shortside();
		    int ls = bx.longside();
		    if (v < vmin_this || (v == vmin_this && ss < smin_this)) {
			vmin_this = v;
			smin_this = ss;
			imin_this = k;
		    }
		    if (v > vmax_this || (v == vmax_this && ls > lmax_this)) {
			vmax_this = v;
			lmax_this = ls;
			imax_this = k;
		    }
		}
#ifdef _OPENMP
#pragma omp critical (amr_prtgs)
#endif	    	    
		{
		    if (vmin_this < vmin || (vmin_this == vmin && smin_this < smin)) {
			vmin = vmin_this;
			smin = smin_this;
			imin = imin_this;
		    }
		    if (vmax_this > vmax || (vmax_this == vmax && lmax_this > lmax)) {
			vmax = vmax_this;
			lmax = lmax_this;
			imax = imax_this;
		    }
		}
	    }
	    const Box& bmin = bs[imin];
	    const Box& bmax = bs[imax];
	    os << "           "
	       << " smallest grid: "
		D_TERM(<< bmin.length(0),
		       << " x " << bmin.length(1),
		       << " x " << bmin.length(2))
	       << "  biggest grid: "
		D_TERM(<< bmax.length(0),
		       << " x " << bmax.length(1),
		       << " x " << bmax.length(2))
	       << '\n';
	}
    }

    os << std::endl; // Make sure we flush!
}

void
Amr::ProjPeriodic (BoxList&        blout,
                   const Geometry& geom)
{
    //
    // Add periodic translates to blout.
    //
    Box domain = geom.Domain();

    BoxList blorig(blout);

    int nist,njst,nkst;
    int niend,njend,nkend;
    nist = njst = nkst = 0;
    niend = njend = nkend = 0;
    D_TERM( nist , =njst , =nkst ) = -1;
    D_TERM( niend , =njend , =nkend ) = +1;

    int ri,rj,rk;
    for (ri = nist; ri <= niend; ri++)
    {
        if (ri != 0 && !geom.isPeriodic(0))
            continue;
        if (ri != 0 && geom.isPeriodic(0))
            blorig.shift(0,ri*domain.length(0));
        for (rj = njst; rj <= njend; rj++)
        {
            if (rj != 0 && !geom.isPeriodic(1))
                continue;
            if (rj != 0 && geom.isPeriodic(1))
                blorig.shift(1,rj*domain.length(1));
            for (rk = nkst; rk <= nkend; rk++)
            {
                if (rk != 0 && !geom.isPeriodic(2))
                    continue;
                if (rk != 0 && geom.isPeriodic(2))
                    blorig.shift(2,rk*domain.length(2));

                BoxList tmp(blorig);
                tmp.intersect(domain);
                blout.catenate(tmp);
 
                if (rk != 0 && geom.isPeriodic(2))
                    blorig.shift(2,-rk*domain.length(2));
            }
            if (rj != 0 && geom.isPeriodic(1))
                blorig.shift(1,-rj*domain.length(1));
        }
        if (ri != 0 && geom.isPeriodic(0))
            blorig.shift(0,-ri*domain.length(0));
    }
}

void
Amr::grid_places (int              lbase,
                  Real             time,
                  int&             new_finest,
                  Array<BoxArray>& new_grids)
{
    BL_PROFILE("Amr::grid_places()");

    int i, max_crse = std::min(finest_level,max_level-1);

    const Real strttime = ParallelDescriptor::second();

    if (lbase == 0)
    {
        //
        // Recalculate level 0 BoxArray in case max_grid_size[0] has changed.
        // This is done exactly as in defBaseLev().
        //
        BoxArray lev0(1);

        lev0.set(0,BoxLib::coarsen(geom[0].Domain(),2));
        //
        // Now split up into list of grids within max_grid_size[0] limit.
        //
        lev0.maxSize(max_grid_size[0]/2);
        //
        // Now refine these boxes back to level 0.
        //
        lev0.refine(2);

        new_grids[0] = lev0;

       // If Nprocs > Ngrids and refine_grid_layout == 1 then break up the grids
       //    into smaller chunks for better load balancing
       // We need to impose this here in the event that we return with fixed_grids
       //    and never have a chance to impose it later
       impose_refine_grid_layout(lbase,lbase,new_grids);
    }

    if ( time == 0. && !initial_grids_file.empty() && !useFixedCoarseGrids)
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min(new_finest,initial_ba.size());

        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = initial_ba[lev-1].size();
            for (i = 0; i < ngrid; i++)
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
    if ( ! initial_grids_file.empty() && useFixedCoarseGrids)
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min(new_finest,initial_ba.size());

        for (int lev = lbase+1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = initial_ba[lev-1].size();
            for (i = 0; i < ngrid; i++)
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
        new_finest = std::min(new_finest,regrid_ba.size());
        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = regrid_ba[lev-1].size();
            for (i = 0; i < ngrid; i++)
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

    //
    // Construct problem domain at each level.
    //
    Array<IntVect> bf_lev(max_level); // Blocking factor at each level.
    Array<IntVect> rr_lev(max_level);
    Array<Box>     pc_domain(max_level);  // Coarsened problem domain.
    for (i = 0; i <= max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            bf_lev[i][n] = std::max(1,blocking_factor[i+1]/ref_ratio[i][n]);
    }
    for (i = lbase; i < max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            rr_lev[i][n] = (ref_ratio[i][n]*bf_lev[i][n])/bf_lev[i+1][n];
    }
    for(i = lbase; i <= max_crse; i++) {
      pc_domain[i] = BoxLib::coarsen(geom[i].Domain(),bf_lev[i]);
    }
    //
    // Construct proper nesting domains.
    //
    Array<BoxList> p_n(max_level);      // Proper nesting domain.
    Array<BoxList> p_n_comp(max_level); // Complement proper nesting domain.

    BoxList bl(amr_level[lbase].boxArray());
    bl.simplify();
    bl.coarsen(bf_lev[lbase]);
    p_n_comp[lbase].complementIn(pc_domain[lbase],bl);
    p_n_comp[lbase].simplify();
    p_n_comp[lbase].accrete(n_proper);
    Amr::ProjPeriodic(p_n_comp[lbase], Geometry(pc_domain[lbase]));
    p_n[lbase].complementIn(pc_domain[lbase],p_n_comp[lbase]);
    p_n[lbase].simplify();
    bl.clear();

    for (i = lbase+1; i <= max_crse; i++)
    {
        p_n_comp[i] = p_n_comp[i-1];

        // Need to simplify p_n_comp or the number of grids can too large for many levels.
        p_n_comp[i].simplify();

        p_n_comp[i].refine(rr_lev[i-1]);
        p_n_comp[i].accrete(n_proper);

        Amr::ProjPeriodic(p_n_comp[i], Geometry(pc_domain[i]));

        p_n[i].complementIn(pc_domain[i],p_n_comp[i]);
        p_n[i].simplify();
    }
    //
    // Now generate grids from finest level down.
    //
    new_finest = lbase;

    for (int levc = max_crse; levc >= lbase; levc--)
    {
        int levf = levc+1;
        //
        // Construct TagBoxArray with sufficient grow factor to contain
        // new levels projected down to this level.
        //
        int ngrow = 0;

        if (levf < new_finest)
        {
            BoxArray ba_proj(new_grids[levf+1]);

            ba_proj.coarsen(ref_ratio[levf]);
            ba_proj.grow(n_proper);
            ba_proj.coarsen(ref_ratio[levc]);

            BoxArray levcBA = amr_level[levc].boxArray();

            while (!levcBA.contains(ba_proj))
            {
                BoxArray tmp = levcBA;
                tmp.grow(1);
                levcBA = tmp;
                ngrow++;
            }
        }
        TagBoxArray tags(amr_level[levc].boxArray(),n_error_buf[levc]+ngrow);
    
        //
        // Only use error estimation to tag cells for the creation of new grids
        //      if the grids at that level aren't already fixed.
        //

        if ( ! (useFixedCoarseGrids && levc < useFixedUpToLevel) ) {
            amr_level[levc].errorEst(tags,
                                     TagBox::CLEAR,TagBox::SET,time,
                                     n_error_buf[levc],ngrow);
	}

        //
        // If new grids have been constructed above this level, project
        // those grids down and tag cells on intersections to ensure
        // proper nesting.
        //
        // NOTE: this loop replaces the previous code:
        //      if (levf < new_finest) 
        //          tags.setVal(ba_proj,TagBox::SET);
        // The problem with this code is that it effectively 
        // "buffered the buffer cells",  i.e., the grids at level
        // levf+1 which were created by buffering with n_error_buf[levf]
        // are then coarsened down twice to define tagging at
        // level levc, which will then also be buffered.  This can
        // create grids which are larger than necessary.
        //
        if (levf < new_finest)
        {
            int nerr = n_error_buf[levf];

            BoxList bl_tagged(new_grids[levf+1]);
            bl_tagged.simplify();
            bl_tagged.coarsen(ref_ratio[levf]);
            //
            // This grows the boxes by nerr if they touch the edge of the
            // domain in preparation for them being shrunk by nerr later.
            // We want the net effect to be that grids are NOT shrunk away
            // from the edges of the domain.
            //
            for (BoxList::iterator blt = bl_tagged.begin(), End = bl_tagged.end();
                 blt != End;
                 ++blt)
            {
                for (int idir = 0; idir < BL_SPACEDIM; idir++)
                {
                    if (blt->smallEnd(idir) == geom[levf].Domain().smallEnd(idir))
                        blt->growLo(idir,nerr);
                    if (blt->bigEnd(idir) == geom[levf].Domain().bigEnd(idir))
                        blt->growHi(idir,nerr);
                }
            }
            Box mboxF = BoxLib::grow(bl_tagged.minimalBox(),1);
            BoxList blFcomp;
            blFcomp.complementIn(mboxF,bl_tagged);
            blFcomp.simplify();
            bl_tagged.clear();

            const IntVect& iv = IntVect(D_DECL(nerr/ref_ratio[levf][0],
                                               nerr/ref_ratio[levf][1],
                                               nerr/ref_ratio[levf][2]));
            blFcomp.accrete(iv);
            BoxList blF;
            blF.complementIn(mboxF,blFcomp);
            BoxArray baF(blF);
            blF.clear();
            baF.grow(n_proper);
            //
            // We need to do this in case the error buffering at
            // levc will not be enough to cover the error buffering
            // at levf which was just subtracted off.
            //
            for (int idir = 0; idir < BL_SPACEDIM; idir++) 
            {
                if (nerr > n_error_buf[levc]*ref_ratio[levc][idir]) 
                    baF.grow(idir,nerr-n_error_buf[levc]*ref_ratio[levc][idir]);
            }

            baF.coarsen(ref_ratio[levc]);

            tags.setVal(baF,TagBox::SET);
        }
        //
        // Buffer error cells.
        //
        tags.buffer(n_error_buf[levc]+ngrow);

        if (useFixedCoarseGrids)
        {
           if (levc>=useFixedUpToLevel)
           {
               BoxArray bla(amr_level[levc].getAreaNotToTag());
               tags.setVal(bla,TagBox::CLEAR);
           }
           if (levc<useFixedUpToLevel)
              new_finest = std::max(new_finest,levf);
        }

        //
        // Coarsen the taglist by blocking_factor/ref_ratio.
        //
        int bl_max = 0;
        for (int n=0; n<BL_SPACEDIM; n++)
            bl_max = std::max(bl_max,bf_lev[levc][n]);
        if (bl_max > 1) 
            tags.coarsen(bf_lev[levc]);
        //
        // Remove or add tagged points which violate/satisfy additional 
        // user-specified criteria.
        //
        amr_level[levc].manual_tags_placement(tags, bf_lev);
        //
        // Map tagged points through periodic boundaries, if any.
        //
        tags.mapPeriodic(Geometry(pc_domain[levc]));
        //
        // Remove cells outside proper nesting domain for this level.
        //
        tags.setVal(p_n_comp[levc],TagBox::CLEAR);
        //
        // Create initial cluster containing all tagged points.
        //
	std::vector<IntVect> tagvec;
	tags.collate(tagvec);
        tags.clear();

        if (tagvec.size() > 0)
        {
            //
            // Created new level, now generate efficient grids.
            //
            if ( !(useFixedCoarseGrids && levc<useFixedUpToLevel) )
                new_finest = std::max(new_finest,levf);
            //
            // Construct initial cluster.
            //
            ClusterList clist(&tagvec[0], tagvec.size());
            clist.chop(grid_eff);
            BoxDomain bd;
            bd.add(p_n[levc]);
            clist.intersect(bd);
            bd.clear();
            //
            // Efficient properly nested Clusters have been constructed
            // now generate list of grids at level levf.
            //
            BoxList new_bx;
            clist.boxList(new_bx);
            new_bx.refine(bf_lev[levc]);
            new_bx.simplify();
            BL_ASSERT(new_bx.isDisjoint());
	    
	    if (new_bx.size()>0) {
	      if ( !(geom[levc].Domain().contains(BoxArray(new_bx).minimalBox())) ) {
		// Chop new grids outside domain, note that this is likely to result in
		//  new grids that violate blocking_factor....see warning checking below
		new_bx = BoxLib::intersect(new_bx,geom[levc].Domain());
	      }
	    }

            IntVect largest_grid_size;
            for (int n = 0; n < BL_SPACEDIM; n++)
                largest_grid_size[n] = max_grid_size[levf] / ref_ratio[levc][n];
            //
            // Ensure new grid boxes are at most max_grid_size in index dirs.
            //
            new_bx.maxSize(largest_grid_size);

#ifdef BL_FIX_GATHERV_ERROR
	      int wcount = 0, iLGS = largest_grid_size[0];

              while (new_bx.size() < 64 && wcount++ < 4)
              {
                  iLGS /= 2;
                  if (ParallelDescriptor::IOProcessor())
                  {
                      std::cout << "BL_FIX_GATHERV_ERROR:  using iLGS = " << iLGS
                                << "   largest_grid_size was:  " << largest_grid_size[0]
                                << '\n';
                      std::cout << "BL_FIX_GATHERV_ERROR:  new_bx.size() was:   "
                                << new_bx.size() << '\n';
                  }

                  new_bx.maxSize(iLGS);

                  if (ParallelDescriptor::IOProcessor())
                  {
                      std::cout << "BL_FIX_GATHERV_ERROR:  new_bx.size() now:   "
                                << new_bx.size() << '\n';
                  }
	      }
#endif
            //
            // Refine up to levf.
            //
            new_bx.refine(ref_ratio[levc]);
            BL_ASSERT(new_bx.isDisjoint());

	    if (new_bx.size()>0) {
	      if ( !(geom[levf].Domain().contains(BoxArray(new_bx).minimalBox())) ) {
		new_bx = BoxLib::intersect(new_bx,geom[levf].Domain());
	      }
	      if (ParallelDescriptor::IOProcessor()) {
		for (int d=0; d<BL_SPACEDIM; ++d) {
		  bool ok = true;
		  for (BoxList::const_iterator bli = new_bx.begin(); bli != new_bx.end(); ++bli) {
		    int len = bli->length(d);
		    int bf = blocking_factor[levf];
		    ok &= (len/bf) * bf == len;
		  }
		  if (!ok) {
		    BoxLib::Warning("WARNING: New grids violate blocking factor near upper boundary");
		  }
		}
	      }
	    }
            if(levf > useFixedUpToLevel) {
              new_grids[levf].define(new_bx);
	    }
        }
    }

    // If Nprocs > Ngrids and refine_grid_layout == 1 then break up the grids
    //    into smaller chunks for better load balancing
    impose_refine_grid_layout(lbase,new_finest,new_grids);

    if (verbose > 0)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor())
            std::cout << "grid_places() time: " << stoptime << " new finest: " << new_finest<< '\n';
#ifdef BL_LAZY
	});
#endif
    }
}

void
Amr::bldFineLevels (Real strt_time)
{
    BL_PROFILE("Amr::bldFineLevels()");
    finest_level = 0;

    Array<BoxArray> grids(max_level+1);
    //
    // Get initial grid placement.
    //
    do
    {
        int new_finest;

        grid_places(finest_level,strt_time,new_finest,grids);

        if (new_finest <= finest_level) break;
        //
        // Create a new level and link with others.
        //
        finest_level = new_finest;

        AmrLevel* level = (*levelbld)(*this,
                                      new_finest,
                                      geom[new_finest],
                                      grids[new_finest],
                                      strt_time);

        amr_level.set(new_finest,level);

        amr_level[new_finest].initData();
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
	    for (int i = 0; i <= finest_level; i++)
	      grids[i] = amr_level[i].boxArray();

	    regrid(0,strt_time,true);

	    grids_the_same = true;

	    for (int i = 0; i <= finest_level && grids_the_same; i++)
	      if (!(grids[i] == amr_level[i].boxArray()))
                grids_the_same = false;

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
    int i;
    sub_cycle = true;
    if (pp.contains("nosub"))
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "Warning: The nosub flag has been deprecated.\n ";
            std::cout << "... please use subcycling_mode to control subcycling.\n";
        }
        int nosub;
        pp.query("nosub",nosub);
        if (nosub > 0)
            sub_cycle = false;
        else
            BoxLib::Error("nosub <= 0 not allowed.\n");
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
        for (i = 0; i <= max_level; i++)
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
            for (i = 1; i <= max_level; i++)
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
                BoxLib::Error("First entry of subcycling_iterations must be 1");
            }
        }
        else
        {
            BoxLib::Error("Must provide a valid subcycling_iterations if mode is Manual");
        }
        for (i = 1; i <= max_level; i++)
        {
            if (n_cycle[i] > MaxRefRatio(i-1))
                BoxLib::Error("subcycling iterations must always be <= ref_ratio");
            if (n_cycle[i] <= 0)
                BoxLib::Error("subcycling iterations must always be > 0");
        }
    }
    else if (subcycling_mode == "Auto")
    {
        n_cycle[0] = 1;
        for (i = 1; i <= max_level; i++)
        {
            n_cycle[i] = MaxRefRatio(i-1);
        } 
    }
    else if (subcycling_mode == "Optimal")
    {
        // if subcycling mode is Optimal, n_cycle is set dynamically.
        // We'll initialize it to be Auto subcycling.
        n_cycle[0] = 1;
        for (i = 1; i <= max_level; i++)
        {
            n_cycle[i] = MaxRefRatio(i-1);
        } 
    }
    else
    {
        std::string err_message = "Unrecognzied subcycling mode: " + subcycling_mode + "\n";
        BoxLib::Error(err_message.c_str());
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
	    BoxLib::Warning("Warning: both amr.check_int and amr.check_per are > 0.");
    }

    plot_file_root = "plt";
    pp.query("plot_file",plot_file_root);

    plot_int = -1;
    pp.query("plot_int",plot_int);

    plot_per = -1.0;
    pp.query("plot_per",plot_per);

    if (plot_int > 0 && plot_per > 0)
    {
        if (ParallelDescriptor::IOProcessor())
            BoxLib::Warning("Warning: both amr.plot_int and amr.plot_per are > 0.");
    }

    small_plot_file_root = "smallplt";
    pp.query("small_plot_file",small_plot_file_root);

    small_plot_int = -1;
    pp.query("small_plot_int",small_plot_int);

    small_plot_per = -1.0;
    pp.query("small_plot_per",small_plot_per);

    if (small_plot_int > 0 && small_plot_per > 0)
    {
        if (ParallelDescriptor::IOProcessor())
            BoxLib::Warning("Warning: both amr.small_plot_int and amr.small_plot_per are > 0.");
    }

    write_plotfile_with_checkpoint = 1;
    pp.query("write_plotfile_with_checkpoint",write_plotfile_with_checkpoint);

    stream_max_tries = 4;
    pp.query("stream_max_tries",stream_max_tries);
    stream_max_tries = std::max(stream_max_tries, 1);

    abort_on_stream_retry_failure = false;
    pp.query("abort_on_stream_retry_failure",abort_on_stream_retry_failure);

}


bool
Amr::okToRegrid(int level)
{
    if (regrid_int[level] < 0)
        return false;
    else
        return level_count[level] >= regrid_int[level] && amr_level[level].okToRegrid();
}

Real
Amr::computeOptimalSubcycling(int n, int* best, Real* dt_max, Real* est_work, int* cycle_max)
{
    BL_ASSERT(cycle_max[0] == 1);
    // internally these represent the total number of steps at a level, 
    // not the number of cycles
    int cycles[n];
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

const Array<BoxArray>& Amr::getInitialBA()
{
  return initial_ba;
}

void 
Amr::impose_refine_grid_layout (int lbase, int new_finest, Array<BoxArray>& new_grids)
{
    const int NProcs = ParallelDescriptor::NProcs();
    if (NProcs > 1 && refine_grid_layout)
    {
        //
        // Chop up grids if fewer grids at level than CPUs.
        // The idea here is to make more grids on a given level
        // to spread the work around.
        //
        for (int cnt = 1; cnt <= 4; cnt *= 2)
        {
            for (int i = lbase; i <= new_finest; i++)
            {
                const int ChunkSize = max_grid_size[i]/cnt;

                IntVect chunk(D_DECL(ChunkSize,ChunkSize,ChunkSize));
                //
                // We go from Z -> Y -> X to promote cache-efficiency.
                //
                for (int j = BL_SPACEDIM-1; j >= 0 ; j--)
                {
                    chunk[j] /= 2;

                    if ( (new_grids[i].size() < NProcs) && (chunk[j]%blocking_factor[i] == 0) )
                    {
                        new_grids[i].maxSize(chunk);
                    }
                }
            }
        }
    }
}

#ifdef USE_PARTICLES
void 
Amr::addOneParticle (int id_in, int cpu_in, 
                     std::vector<double>& xloc, std::vector<double>& attributes)
{
    amr_level[0].addOneParticle(id_in,cpu_in,xloc,attributes);
}
void 
Amr::RedistributeParticles ()
{
    // Call Redistribute with where_already_called = false
    //                        full_where           = false
    //                        lev_min              = 0
    amr_level[0].particle_redistribute(0,true);
}
void
Amr::GetParticleIDs (Array<int>& part_ids)
{
    //
    // The AmrLevel class is where we have access to the particle container.
    //
    amr_level[0].GetParticleIDs(part_ids);
}
void
Amr::GetParticleCPU (Array<int>& part_cpu)
{
    //
    // The AmrLevel class is where we have access to the particle container.
    //
    amr_level[0].GetParticleCPU(part_cpu);
}
void
Amr::GetParticleLocations (Array<Real>& part_data)
{
    //
    // The AmrLevel class is where we have access to the particle container.
    //
    amr_level[0].GetParticleLocations(part_data);
}
void 
Amr::GetParticleData (Array<Real>& part_data, int start_comp, int num_comp)
{
    //
    // The AmrLevel class is where we have access to the particle container.
    //
    amr_level[0].GetParticleData(part_data,start_comp,num_comp);
}
#endif


void
Amr::AddProcsToSidecar(int nSidecarProcs, int prevSidecarProcs) {

//    MultiFab::flushFBCache();
//    Geometry::flushFPBCache();
//    FabArrayBase::flushCPCache();
    DistributionMapping::FlushCache();

    Array<BoxArray> allBoxes(finest_level + 1);

    for(int ilev(0); ilev < allBoxes.size(); ++ilev) {
      allBoxes[ilev] = boxArray(ilev);
    }

    Array<Array<int> > mLDM;
    // ---- just use the random map for now
    int maxRank(ParallelDescriptor::NProcsAll() - nSidecarProcs - 1);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "_______ maxRank = " << maxRank << std::endl;
    }

    mLDM = DistributionMapping::MultiLevelMapRandom(ref_ratio, allBoxes, maxGridSize(0), maxRank);

    for(int iMap(0); iMap < mLDM.size(); ++iMap) {
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "_in Amr::AddProcsToSidecar:  calling MoveAllFabs:" << std::endl;
      }
      MultiFab::MoveAllFabs(mLDM[iMap]);
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "_in Amr::AddProcsToSidecar:  after calling MoveAllFabs:" << std::endl;
      }
    }
    Geometry::flushFPBCache();
    VisMF::SetNOutFiles(checkpoint_nfiles);

#ifdef USE_PARTICLES
    RedistributeParticles();
#endif

    bool inSidecar(ParallelDescriptor::MyProc() > maxRank);
    if(inSidecar) {
      DistributionMapping::DeleteCache();
    }
}


void
Amr::AddProcsToComp(int nSidecarProcs, int prevSidecarProcs) {
#if BL_USE_MPI
//    MultiFab::flushFBCache();
//    Geometry::flushFPBCache();
//    FabArrayBase::CPC::flushCPCache();
    //FabArrayBase::flushTileArrayCache();
    DistributionMapping::FlushCache();

    MPI_Group scsGroup, allGroup;
    MPI_Comm  scsComm;

    BL_ASSERT(nSidecarProcs < prevSidecarProcs);
    int nProcsAll(ParallelDescriptor::NProcsAll());
    int ioProcNumAll(ParallelDescriptor::IOProcessorNumberAll());
    int ioProcNumSCS(-1);

    // ---- make a group with ioprocnum and the new comp ranks that were part of the sidecar
    // ---- then initialize all required data for the new ranks (amr, amrlevels, ...)
    Array<int> groupRanks(prevSidecarProcs - nSidecarProcs + 1, -1);  // ---- + 1 for ioprocnum
    groupRanks[0] = ioProcNumAll;
    int ngStart(nProcsAll - prevSidecarProcs);
    for(int ip(1); ip < groupRanks.size(); ++ip) {
      groupRanks[ip] = ngStart++;
    }
    if(ParallelDescriptor::IOProcessor()) {
      for(int ip(0); ip < groupRanks.size(); ++ip) {
        std::cout << "_in AddProcsToComp:  groupRanks[" << ip << "] = " << groupRanks[ip] << std::endl;
      }
    }
    BL_MPI_REQUIRE( MPI_Comm_group(ParallelDescriptor::CommunicatorAll(), &allGroup) );
    BL_MPI_REQUIRE( MPI_Group_incl(allGroup, groupRanks.size(), groupRanks.dataPtr(), &scsGroup) );
    BL_MPI_REQUIRE( MPI_Comm_create(ParallelDescriptor::CommunicatorAll(), scsGroup, &scsComm) );

    // ---- dont always assume ioprocnum == 0 everywhere
    BL_MPI_REQUIRE( MPI_Group_translate_ranks(allGroup, 1, &ioProcNumAll, scsGroup, &ioProcNumSCS) );

    int scsMyId;
    BL_MPI_REQUIRE( MPI_Group_rank(scsGroup, &scsMyId) );

    // ---- send all amr data from ioprocnum to the new comp ranks
    if(scsMyId != MPI_UNDEFINED) {
      int currentSeqNumber(-4);
      if(scsMyId == ioProcNumSCS) {
        currentSeqNumber = ParallelDescriptor::SeqNum(1);
      }
      ParallelDescriptor::Bcast(&currentSeqNumber, 1, ioProcNumAll, scsComm);
      if(scsMyId != ioProcNumSCS) {
        ParallelDescriptor::SeqNum(2, currentSeqNumber);
      }


      // ---- pack up the ints
      Array<int> allInts;
      int allIntsSize(0);
      int dt_level_Size(dt_level.size()), dt_min_Size(dt_min.size());
      int ref_ratio_Size(ref_ratio.size()), amr_level_Size(amr_level.size()), geom_Size(geom.size());
      int state_plot_vars_Size(state_plot_vars.size()), derive_plot_vars_Size(derive_plot_vars.size());
      int state_small_plot_vars_Size(state_small_plot_vars.size());
      if(scsMyId == ioProcNumSCS) {
        allInts.push_back(max_level);
        allInts.push_back(finest_level);
        allInts.push_back(n_proper);
        allInts.push_back(last_checkpoint); 
        allInts.push_back(check_int);
        allInts.push_back(last_plotfile);   
        allInts.push_back(last_smallplotfile);   
        allInts.push_back(plot_int);
        allInts.push_back(small_plot_int);
        allInts.push_back(write_plotfile_with_checkpoint);
        allInts.push_back(file_name_digits);
        allInts.push_back(message_int);
        allInts.push_back(which_level_being_advanced);
        allInts.push_back(verbose);
        allInts.push_back(record_grid_info);
        allInts.push_back(record_run_info);
        allInts.push_back(record_run_info_terse);
        allInts.push_back(sub_cycle);
        allInts.push_back(stream_max_tries);
        allInts.push_back(rebalance_grids);

	// ---- these are parmparsed in
        allInts.push_back(plot_nfiles);
        allInts.push_back(mffile_nstreams);
        allInts.push_back(probinit_natonce);
        allInts.push_back(checkpoint_nfiles);
        allInts.push_back(regrid_on_restart);
        allInts.push_back(use_efficient_regrid);
        allInts.push_back(plotfile_on_restart);
        allInts.push_back(checkpoint_on_restart);
        allInts.push_back(compute_new_dt_on_regrid);
        allInts.push_back(useFixedUpToLevel);

        allInts.push_back(level_steps.size());
        for(int i(0); i < level_steps.size(); ++i)     { allInts.push_back(level_steps[i]); }
        allInts.push_back(level_count.size());
        for(int i(0); i < level_count.size(); ++i)     { allInts.push_back(level_count[i]); }
        allInts.push_back(n_cycle.size());
        for(int i(0); i < n_cycle.size(); ++i)         { allInts.push_back(n_cycle[i]); }
        allInts.push_back(regrid_int.size());
        for(int i(0); i < regrid_int.size(); ++i)      { allInts.push_back(regrid_int[i]); }
        allInts.push_back(n_error_buf.size());
        for(int i(0); i < n_error_buf.size(); ++i)    { allInts.push_back(n_error_buf[i]); }
        allInts.push_back(blocking_factor.size());
        for(int i(0); i < blocking_factor.size(); ++i) { allInts.push_back(blocking_factor[i]); }
        allInts.push_back(max_grid_size.size());
        for(int i(0); i < max_grid_size.size(); ++i)   { allInts.push_back(max_grid_size[i]); }

	// ---- for non-int arrays
        allInts.push_back(dt_level.size());
        allInts.push_back(dt_min.size());
        allInts.push_back(ref_ratio.size());
        allInts.push_back(amr_level.size());
        allInts.push_back(geom.size());
        allInts.push_back(state_plot_vars.size());
        allInts.push_back(state_small_plot_vars.size());
        allInts.push_back(derive_plot_vars.size());

        allIntsSize = allInts.size();
      }

      BoxLib::BroadcastArray(allInts, scsMyId, ioProcNumAll, scsComm);

      // ---- unpack the ints
      if(scsMyId != ioProcNumSCS) {
	int count(0), aSize(-1);
        max_level                  = allInts[count++];
        finest_level               = allInts[count++];
        n_proper                   = allInts[count++];
        last_checkpoint            = allInts[count++]; 
        check_int                  = allInts[count++];
        last_plotfile              = allInts[count++];   
        last_smallplotfile         = allInts[count++];   
        plot_int                   = allInts[count++];
        small_plot_int             = allInts[count++];
        write_plotfile_with_checkpoint = allInts[count++];
        file_name_digits           = allInts[count++];
        message_int                = allInts[count++];
        which_level_being_advanced = allInts[count++];
        verbose                    = allInts[count++];
        record_grid_info           = allInts[count++];
        record_run_info            = allInts[count++];
        record_run_info_terse      = allInts[count++];
        sub_cycle                  = allInts[count++];
        stream_max_tries           = allInts[count++];
        rebalance_grids            = allInts[count++];

        plot_nfiles                = allInts[count++];
        mffile_nstreams            = allInts[count++];
        probinit_natonce           = allInts[count++];
        checkpoint_nfiles          = allInts[count++];
        regrid_on_restart          = allInts[count++];
        use_efficient_regrid       = allInts[count++];
        plotfile_on_restart        = allInts[count++];
        checkpoint_on_restart      = allInts[count++];
        compute_new_dt_on_regrid   = allInts[count++];
        useFixedUpToLevel          = allInts[count++];

        aSize                      = allInts[count++];
	level_steps.resize(aSize);
        for(int i(0); i < level_steps.size(); ++i)     { level_steps[i] = allInts[count++]; }
        aSize                      = allInts[count++];
        level_count.resize(aSize);
        for(int i(0); i < level_count.size(); ++i)     { level_count[i] = allInts[count++]; }
        aSize                      = allInts[count++];
        n_cycle.resize(aSize);
        for(int i(0); i < n_cycle.size(); ++i)         { n_cycle[i] = allInts[count++]; }
        aSize                      = allInts[count++];
        regrid_int.resize(aSize);
        for(int i(0); i < regrid_int.size(); ++i)      { regrid_int[i] = allInts[count++]; }
        aSize                      = allInts[count++];
        n_error_buf.resize(aSize);
        for(int i(0); i < n_error_buf.size(); ++i)     { n_error_buf[i] = allInts[count++]; }
        aSize                      = allInts[count++];
        blocking_factor.resize(aSize);
        for(int i(0); i < blocking_factor.size(); ++i) { blocking_factor[i] = allInts[count++]; }
        aSize                      = allInts[count++];
        max_grid_size.resize(aSize);
        for(int i(0); i < max_grid_size.size(); ++i)   { max_grid_size[i] = allInts[count++]; }

        dt_level_Size              = allInts[count++];
        dt_min_Size                = allInts[count++];
        ref_ratio_Size             = allInts[count++];
        amr_level_Size             = allInts[count++];
        geom_Size                  = allInts[count++];
        state_plot_vars_Size       = allInts[count++];
        state_small_plot_vars_Size = allInts[count++];
        derive_plot_vars_Size      = allInts[count++];

	BL_ASSERT(count == allInts.size());
      }


      // ---- pack up the Reals
      Array<Real> allReals;
      int allRealsSize(0);
      if(scsMyId == ioProcNumSCS) {
        allReals.push_back(cumtime);
        allReals.push_back(start_time);
        allReals.push_back(grid_eff);
        allReals.push_back(check_per);
        allReals.push_back(plot_per);
        allReals.push_back(small_plot_per);

        for(int i(0); i < dt_level.size(); ++i)   { allReals.push_back(dt_level[i]); }
        for(int i(0); i < dt_min.size(); ++i)     { allReals.push_back(dt_min[i]); }

	allRealsSize = allReals.size();
      }

      BoxLib::BroadcastArray(allReals, scsMyId, ioProcNumAll, scsComm);

      // ---- unpack the Reals
      if(scsMyId != ioProcNumSCS) {
	int count(0);
        cumtime    = allReals[count++];
        start_time = allReals[count++];
        grid_eff   = allReals[count++];
        check_per  = allReals[count++];
        plot_per   = allReals[count++];
        small_plot_per = allReals[count++];

	dt_level.resize(dt_level_Size);
        for(int i(0); i < dt_level.size(); ++i)  { dt_level[i] = allReals[count++]; }
	dt_min.resize(dt_min_Size);
        for(int i(0); i < dt_min.size(); ++i)    { dt_min[i]   = allReals[count++]; }
      }


      // ---- pack up the bools
      Array<int> allBools;  // ---- just use ints here
      int allBoolsSize(0);
      if(scsMyId == ioProcNumSCS) {
        allBools.push_back(abort_on_stream_retry_failure);
        allBools.push_back(bUserStopRequest);
        for(int i(0); i < BL_SPACEDIM; ++i)    { allBools.push_back(isPeriodic[i]); }
        allBools.push_back(first_plotfile);

        allBools.push_back(plot_files_output);
        allBools.push_back(refine_grid_layout);
        allBools.push_back(checkpoint_files_output);
        allBools.push_back(initialized);
        allBools.push_back(useFixedCoarseGrids);
        allBools.push_back(first_smallplotfile);

	allBoolsSize = allBools.size();
      }

      BoxLib::BroadcastArray(allBools, scsMyId, ioProcNumAll, scsComm);

      // ---- unpack the bools
      if(scsMyId != ioProcNumSCS) {
	int count(0);

        abort_on_stream_retry_failure = allBools[count++];
        bUserStopRequest              = allBools[count++];
        for(int i(0); i < BL_SPACEDIM; ++i)    { isPeriodic[i] = allBools[count++]; }
        first_plotfile                = allBools[count++];

        plot_files_output             = allBools[count++];
        refine_grid_layout            = allBools[count++];
        checkpoint_files_output       = allBools[count++];
        initialized                   = allBools[count++];
        useFixedCoarseGrids           = allBools[count++];
        first_smallplotfile           = allBools[count++];
      }


      // ---- pack up the strings
      Array<std::string> allStrings;
      Array<char> serialStrings;
      int serialStringsSize(0);
      if(scsMyId == ioProcNumSCS) {
        allStrings.push_back(regrid_grids_file);
        allStrings.push_back(initial_grids_file);
        allStrings.push_back(check_file_root);
        allStrings.push_back(subcycling_mode);
        allStrings.push_back(plot_file_root);
        allStrings.push_back(small_plot_file_root);
        allStrings.push_back(restart_chkfile);
        allStrings.push_back(restart_pltfile);
        allStrings.push_back(probin_file);

        std::list<std::string>::iterator lit;
	for( lit = state_plot_vars.begin(); lit != state_plot_vars.end(); ++lit) {
          allStrings.push_back(*lit);
	}
	for( lit = state_small_plot_vars.begin(); lit != state_small_plot_vars.end(); ++lit) {
          allStrings.push_back(*lit);
	}
	for( lit = derive_plot_vars.begin(); lit != derive_plot_vars.end(); ++lit) {
          allStrings.push_back(*lit);
	}

	serialStrings = BoxLib::SerializeStringArray(allStrings);
	serialStringsSize = serialStrings.size();
      }

      BoxLib::BroadcastArray(serialStrings, scsMyId, ioProcNumAll, scsComm);

      // ---- unpack the strings
      if(scsMyId != ioProcNumSCS) {
	int count(0);
        allStrings = BoxLib::UnSerializeStringArray(serialStrings);

        regrid_grids_file  = allStrings[count++];
        initial_grids_file = allStrings[count++];
        check_file_root    = allStrings[count++];
        subcycling_mode    = allStrings[count++];
        plot_file_root     = allStrings[count++];
        small_plot_file_root = allStrings[count++];
        restart_chkfile    = allStrings[count++];
        restart_pltfile    = allStrings[count++];
        probin_file        = allStrings[count++];

        for(int i(0); i < state_plot_vars_Size; ++i) {
          state_plot_vars.push_back(allStrings[count++]);
	}
        for(int i(0); i < state_small_plot_vars_Size; ++i) {
          state_small_plot_vars.push_back(allStrings[count++]);
	}
        for(int i(0); i < derive_plot_vars_Size; ++i) {
          derive_plot_vars.push_back(allStrings[count++]);
	}
      }


      // ---- pack up the IntVects
      Array<int> allIntVects;
      int allIntVectsSize(0);
      if(scsMyId == ioProcNumSCS) {
        for(int lev(0); lev < ref_ratio.size(); ++lev) {
          for(int i(0); i < BL_SPACEDIM; ++i)    { allIntVects.push_back(ref_ratio[lev][i]); }
	}

	allIntVectsSize = allIntVects.size();
	BL_ASSERT(allIntVectsSize == ref_ratio_Size * BL_SPACEDIM);
      }

      ParallelDescriptor::Bcast(&allIntVectsSize, 1, ioProcNumAll, scsComm);
      if(allIntVectsSize > 0) {
        if(scsMyId != ioProcNumSCS) {
          allIntVects.resize(allIntVectsSize);
        }
        ParallelDescriptor::Bcast(allIntVects.dataPtr(), allIntVectsSize, ioProcNumAll, scsComm);

        // ---- unpack the IntVects
        if(scsMyId != ioProcNumSCS) {
	  int count(0);
	  BL_ASSERT(allIntVectsSize == ref_ratio_Size * BL_SPACEDIM);

	  ref_ratio.resize(ref_ratio_Size);
          for(int lev(0); lev < ref_ratio.size(); ++lev) {
            for(int i(0); i < BL_SPACEDIM; ++i)    { ref_ratio[lev][i] = allIntVects[count++]; }
	  }
        }
      }



      // ---- BoxArrays
      for(int i(0); i < initial_ba.size(); ++i) {
        BoxLib::BroadcastBoxArray(initial_ba[i], scsMyId, ioProcNumAll, scsComm);
      }
      for(int i(0); i < regrid_ba.size(); ++i) {
        BoxLib::BroadcastBoxArray(regrid_ba[i], scsMyId, ioProcNumAll, scsComm);
      }


      if(scsMyId != ioProcNumSCS) {
        levelbld = getLevelBld();
        levelbld->variableSetUpForNewCompProcs();
      }

      // ---- handle amrlevels
      if(scsMyId == ioProcNumSCS) {
        MultiFab::LockAllFAPointers();
      }

      if(scsMyId != ioProcNumSCS) {
        amr_level.resize(0);
        amr_level.resize(amr_level_Size, PArrayManage);
        for(int lev(0); lev < amr_level.size(); ++lev) {
	  amr_level.set(lev,(*levelbld)());
	}
      }

      for(int lev(0); lev <= finest_level; ++lev) {
        amr_level[lev].AddProcsToComp(this, nSidecarProcs, prevSidecarProcs,
	                              ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
      }


      // ---- handle geom
      if(scsMyId != ioProcNumSCS) {
        geom.resize(geom_Size);
      }
      for(int lev(0); lev < geom.size(); ++lev) {
        Geometry::BroadcastGeometry(geom[lev], ioProcNumSCS, scsComm);
      }

      // ---- handle BoundaryPointLists
      BroadcastBoundaryPointList(intersect_lox, scsMyId, ioProcNumSCS, scsComm);
      BroadcastBoundaryPointList(intersect_loy, scsMyId, ioProcNumSCS, scsComm);
      BroadcastBoundaryPointList(intersect_loz, scsMyId, ioProcNumSCS, scsComm);
      BroadcastBoundaryPointList(intersect_hix, scsMyId, ioProcNumSCS, scsComm);
      BroadcastBoundaryPointList(intersect_hiy, scsMyId, ioProcNumSCS, scsComm);
      BroadcastBoundaryPointList(intersect_hiz, scsMyId, ioProcNumSCS, scsComm);

#ifdef USE_STATIONDATA
      BoxLib::Abort("**** Error:  USE_STATIONDATA not yet supported in sidecar resize.");
      // ---- handle station
      if(scsMyId != ioProcNumSCS) {
      }
#endif

      // ---- initialize fortran data
      if(scsMyId != ioProcNumSCS) {
        int probin_file_length(probin_file.length());
	int init(true);
        Array<int> probin_file_name(probin_file_length);
        for(int i(0); i < probin_file_length; ++i) {
          probin_file_name[i] = probin_file[i];
        }
        std::cout << "Starting to read probin ... " << std::endl;
        FORT_PROBINIT(&init, probin_file_name.dataPtr(), &probin_file_length,
                      Geometry::ProbLo(), Geometry::ProbHi());
      }

    }  // ---- end if(scsMyId != MPI_UNDEFINED)


    if(scsComm != MPI_COMM_NULL) {
      BL_MPI_REQUIRE( MPI_Comm_free(&scsComm) );
    }
    if(scsGroup != MPI_GROUP_NULL) {
      BL_MPI_REQUIRE( MPI_Group_free(&scsGroup) );
    }

    VisMF::SetNOutFiles(checkpoint_nfiles);

#ifdef USE_PARTICLES
    RedistributeParticles();
#endif

    bool abortOnError(false);
    MultiFab::CheckFAPointers(abortOnError);

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "%%%%%%%% finished AddProcsToComp." << std::endl;
    }

#endif
}


void
Amr::RedistributeGrids(int how) {
//    MultiFab::flushFBCache();
//    Geometry::flushFPBCache();
//    FabArrayBase::CPC::flushCPCache();
    DistributionMapping::FlushCache();
    if( ! ParallelDescriptor::InCompGroup()) {
      return;
    }

    if(how >= 0) {
      DistributionMapping::InitProximityMap();
      DistributionMapping::Initialize();

        Array<BoxArray> allBoxes(finest_level + 1);
        for(int ilev(0); ilev < allBoxes.size(); ++ilev) {
          allBoxes[ilev] = boxArray(ilev);
        }
        Array<Array<int> > mLDM;
        if(how == 1) {
          mLDM = DistributionMapping::MultiLevelMapPFC(ref_ratio, allBoxes, maxGridSize(0));
        } else if(how == 2) {
          mLDM = DistributionMapping::MultiLevelMapRandom(ref_ratio, allBoxes, maxGridSize(0));
        } else if(how == 3) {
          mLDM = DistributionMapping::MultiLevelMapKnapSack(ref_ratio, allBoxes, maxGridSize(0));
        } else if(how == 0) {   // ---- move all grids to proc zero
	  int minRank(0), maxRank(0);
          mLDM = DistributionMapping::MultiLevelMapRandom(ref_ratio, allBoxes, maxGridSize(0),
	                                                  maxRank, minRank);
        } else if(how == 8) {   // ---- move all grids to proc 8
	  int minRank(8), maxRank(8);
          mLDM = DistributionMapping::MultiLevelMapRandom(ref_ratio, allBoxes, maxGridSize(0),
	                                                  maxRank, minRank);
        } else if(how == 13) {  // ---- move all grids to proc 13
	  int minRank(13), maxRank(13);
          mLDM = DistributionMapping::MultiLevelMapRandom(ref_ratio, allBoxes, maxGridSize(0),
	                                                  maxRank, minRank);
        } else {
	  return;
        }

        for(int iMap(0); iMap < mLDM.size(); ++iMap) {
          MultiFab::MoveAllFabs(mLDM[iMap]);
        }
	Geometry::flushFPBCache();
    }
#ifdef USE_PARTICLES
    RedistributeParticles();
#endif
}


void
Amr::BroadcastBoundaryPointList(BoundaryPointList &bpl, int myLocalId, int rootId, MPI_Comm comm) {
  bool bcastSource(ParallelDescriptor::MyProc() == rootId);
  Array<int> pF, pS;
  Array<double> bplD;
  if(bcastSource) {  // ---- initialize the source data
    std::multimap< std::pair<int, int>, double >::iterator it;
    for(it = bpl.begin(); it != bpl.end(); ++it) {
      pF.push_back(it->first.first);
      pS.push_back(it->first.second);
      bplD.push_back(it->second);
    }
  }
  BoxLib::BroadcastArray(pF, myLocalId, rootId, comm);
  BoxLib::BroadcastArray(pS, myLocalId, rootId, comm);
  BoxLib::BroadcastArray(bplD, myLocalId, rootId, comm);

  BL_ASSERT(pF.size() == pS.size());
  BL_ASSERT(pS.size() == bplD.size());

  if( ! bcastSource) {
    for(int i(0); i < pF.size(); ++i) {
      bpl.insert(std::make_pair(std::make_pair(pF[i],pS[i]),bplD[i]));
    }
  }
}


void
Amr::PrintData(std::ostream& os) {
using std::endl;
#define SHOWVAL(val) { os << #val << " = " << val << std::endl; }
  os << "---------------------------------------------" << std::endl;
  //SHOWVAL(regrid_grids_file);
  //SHOWVAL(initial_grids_file);
  SHOWVAL(max_level);
  SHOWVAL(finest_level);
  SHOWVAL(cumtime);
  SHOWVAL(start_time);
  //SHOWVAL(verbose);
  //SHOWVAL(plot_file_root);
  //SHOWVAL(check_int);
  //SHOWVAL(ref_ratio.size());
  //for(int i(0); i < ref_ratio.size(); ++i) {
    //os << "ref_ratio[" << i << "] = " << ref_ratio[i] << endl;
  //}
  SHOWVAL(amr_level.size());
  os << endl;
  for(int i(0); i < amr_level.size(); ++i) {
    AmrLevel &amrlev = amr_level[i];
    os << "amr_level[" << i << "] = " << &(amr_level[i]) << endl;
    SHOWVAL(amrlev.numGrids());
    SHOWVAL(amrlev.nStep());
    SHOWVAL(amrlev.countCells());
    MultiFab &mf0 = amrlev.get_new_data(0);
    SHOWVAL(mf0.DistributionMap());
    SHOWVAL(mf0.boxArray());
    SHOWVAL(mf0.NFabArrays());
    SHOWVAL(mf0.AllocatedFAPtrID());
  }
  SHOWVAL(geom.size());
  for(int i(0); i < geom.size(); ++i) {
    os << "geom[" << i << "] = " << geom[i] << endl;
  }

  /*
  std::cout << "state_plot_vars = " << endl;
  for(std::list<std::string>::const_iterator li = state_plot_vars.begin(), End = state_plot_vars.end();
      li != End; ++li)
  {
    os << ":::: " << *li << endl;
  }
  */
  os << "=============================================" << endl;
}


void
Amr::BroadcastBCRec(BCRec &bcrec, int myLocalId, int rootId, MPI_Comm localComm)
{
  int bcvect[bcrec.vectSize()];
  if(myLocalId == rootId) {
    for(int i(0); i < bcrec.vectSize(); ++i) {
      bcvect[i] = bcrec.vect()[i];
    }
  }
  ParallelDescriptor::Bcast(bcvect, bcrec.vectSize(), rootId, localComm);
  if(myLocalId != rootId) {
    bcrec.setVect(bcvect);
  }
}


#if 0
void
Amr::SendDataToNewProcs() {

    if (ParallelDescriptor::IOProcessor()) {
      Array<std::string> origSA;
      origSA.push_back("string0");
      origSA.push_back("__string1");
      origSA.push_back("string222");
      std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
      for(int i(0); i < origSA.size(); ++i) {
        std::cout << "origSA[" << i << "] = " << origSA[i] << std::endl;
      }
      Array<char> charArray(BoxLib::SerializeStringArray(origSA));
      Array<std::string> unSA(BoxLib::UnSerializeStringArray(charArray));
      for(int i(0); i < unSA.size(); ++i) {
        std::cout << "unSA[" << i << "] = " << unSA[i] << std::endl;
      }
      std::cout << "~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "<<<<<<<<<<<<<<<<<<" << std::endl;
    }

}
#endif





