//BL_COPYRIGHT_NOTICE

//
// $Id: Amr.cpp,v 1.98 1999-10-21 17:27:21 lijewski Exp $
//

#include <TagBox.H>
#include <Array.H>
#include <CoordSys.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <Cluster.H>
#include <RunStats.H>
#include <LevelBld.H>
#include <AmrLevel.H>
#include <PROB_AMR_F.H>
#include <Amr.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <DistributionMapping.H>
#include <FabSet.H>

#ifdef BL_USE_ARRAYVIEW
#include <DatasetClient.H>
#endif

#ifdef BL_USE_NEW_HFILES
#include <cstdio>
using std::ifstream;
using std::ios;
#else
#include <stdio.h>
#endif

#include <sys/types.h>
#include <sys/wait.h>

//
// This is supposed to be in <sys/types.h>
//
extern "C" pid_t fork();

//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

//
// Static objects.
//
List<aString> Amr::state_plot_vars;
List<aString> Amr::derive_plot_vars;
bool          Amr::first_plotfile = true;
//
// I now want to add a version string to the checkpoint file.
//
static const aString CheckPointVersion = "CheckPointVersion_1.0";
//
// Force immediate full (level 0) regrid() on restart?
//
static int regrid_on_restart = 0;

static int plotfile_on_restart = 0;
//
// Tar up and then remove checkpoint and plot files?
//
#if defined(BL_T3E)
static int tar_and_rm_files = 1;
#else
static int tar_and_rm_files = 0;
#endif
//
// The pid's of the last `tar_and_rm_files' commands for chk and plt files.
//
static pid_t pid_chk = -1;
static pid_t pid_plt = -1;

//
// Run cmd in a child process.  Returns child pid.  Ignores errors.
//

static
pid_t
DoIt (const char* cmd)
{
    pid_t pid = fork();

    if (pid == 0)
    {
        system(cmd);

        exit(0);
    }

    return pid;
}

void
Amr::setDtMin (const Array<REAL>& dt_min_in)
{
    for (int i = 0; i <= finest_level; i++)
        dt_min[i] = dt_min_in[i];
}

AmrLevel&
Amr::getLevel (int lev)
{
    return amr_level[lev];
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

const BoxArray&
Amr::boxArray (int lev) const
{
    return amr_level[lev].boxArray();
}

MultiFab*
Amr::derive (const aString& name,
             Real           time,
             int            lev,
             int            ngrow)
{
    return amr_level[lev].derive(name,time,ngrow);
}

int
Amr::MaxRefRatio (int level) const
{
    int maxval = 0;
    for (int n = 0; n<BL_SPACEDIM; n++) 
        maxval = Max(maxval,ref_ratio[level][n]);
    return maxval;
}

Amr::Amr ()
    :
    amr_level(PArrayManage)
{
    //
    // Setup Geometry from ParmParse file.
    // May be needed for variableSetup or even getLevelBld.
    //
    Geometry::Setup();
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
    max_level        = -1;
    record_run_info  = false;
    record_grid_info = false;
    grid_eff         = 0.7;
    blocking_factor  = 1;
    last_checkpoint  = 0;
    last_plotfile    = 0;
    plot_int         = -1;
    n_proper         = 1;
    max_grid_size    = (BL_SPACEDIM == 2) ? 128 : 32;

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
    pp.query("plotfile_on_restart",plotfile_on_restart);

    if (pp.contains("tar_and_rm_files"))
        //
        // Only query() if it's in the file.
        // Otherwise, we could overwrite any non-zero default.
        //
        pp.query("tar_and_rm_files",tar_and_rm_files);

    sub_cycle = true;
    if (pp.contains("nosub"))
        sub_cycle = false;

    pp.query("regrid_file",grids_file);
    if (pp.contains("run_log"))
    {
        aString log_file_name;
        pp.get("run_log",log_file_name);
        setRecordRunInfo(log_file_name);
    }
    if (pp.contains("grid_log"))
    {
        aString grid_file_name;
        pp.get("grid_log",grid_file_name);
        setRecordGridInfo(grid_file_name);
    }
    if (pp.contains("data_log"))
    {
	aString data_file_name;
	pp.get("data_log",data_file_name);
	setRecordDataInfo(data_file_name);
    }
    if (pp.contains("probin_file"))
    {
        pp.get("probin_file",probin_file);
    }
    //
    // Restart or run from scratch?
    //
    pp.query("restart", restart_file);
    //
    // Read max_level and alloc memory for container objects.
    //
    pp.get("max_level", max_level);
    int nlev     = max_level+1;
    geom.resize(nlev);
    dt_level.resize(nlev);
    level_steps.resize(nlev);
    level_count.resize(nlev);
    regrid_int.resize(nlev);
    n_cycle.resize(nlev);
    dt_min.resize(nlev);
    n_error_buf.resize(nlev);
    amr_level.resize(nlev);
    //
    // Set bogus values.
    //
    for (i = 0; i < nlev; i++)
    {
        dt_level[i]    = 1.e20; // Something nonzero so old & new will differ
        level_steps[i] = 0;
        level_count[i] = 0;
        regrid_int[i]  = 0;
        n_cycle[i]     = 0;
        dt_min[i]      = 0.0;
        n_error_buf[i] = 1;
    }
    ref_ratio.resize(max_level);
    for (i = 0; i < max_level; i++)
        ref_ratio[i] = IntVect::TheZeroVector();
    //
    // Read other amr specific values.
    //
    check_file_root = "chk";
    pp.query("check_file",check_file_root);

    check_int = -1;
    int got_check_int = pp.query("check_int",check_int);

    check_per = -1.0;
    int got_check_per = pp.query("check_per",check_per);

    if (got_check_int == 1 && got_check_per == 1)
    {
        BoxLib::Error("Must only specify amr.check_int OR amr.check_per");
    }
    else if (got_check_per == 1 && ParallelDescriptor::IOProcessor())
    {
        BoxLib::Warning("Specifying amr.check_per will change the time step");
    }

    plot_file_root = "plt";
    pp.query("plot_file",plot_file_root);

    plot_int = -1;
    int got_plot_int = pp.query("plot_int",plot_int);

    plot_per = -1.0;
    int got_plot_per = pp.query("plot_per",plot_per);

    if (got_plot_int == 1 && got_plot_per == 1)
    {
        BoxLib::Error("Must only specify amr.plot_int OR amr.plot_per");
    }
    else if (got_plot_per == 1 && ParallelDescriptor::IOProcessor())
    {
        BoxLib::Warning("Specifying amr.plot_per will change the time step");
    }

    pp.query("max_grid_size",max_grid_size);
    pp.query("n_proper",n_proper);
    pp.query("blocking_factor",blocking_factor);
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
          BoxLib::Error("Must input *either* ref_ratio or ref_ratio_vect");
      }
    }
    //
    // Read computational domain and set geometry.
    //
    Array<int> n_cell(BL_SPACEDIM);
    pp.getarr("n_cell",n_cell,0,BL_SPACEDIM);
    BL_ASSERT(n_cell.length() == BL_SPACEDIM);
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
    //
    // Set regrid interval.
    //
    int ri;
    pp.get("regrid_int",ri);

    for (int k = 0; k <= max_level; k++)
        regrid_int[k] = ri;
}

bool
Amr::isStatePlotVar (const aString& name)
{
    for (ListIterator<aString> li(state_plot_vars); li; ++li)
        if (li() == name)
            return true;
  
    return false;
}

void
Amr::fillStatePlotVarList()
{
    state_plot_vars.clear();
    const DescriptorList& desc_lst = AmrLevel::get_desc_lst();
    for (int typ = 0; typ < desc_lst.length(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (desc_lst[typ].getType() == IndexType::TheCellType())
                state_plot_vars.append(desc_lst[typ].name(comp));
}

void
Amr::clearStatePlotVarList()
{
    state_plot_vars.clear();
}

void
Amr::addStatePlotVar (const aString& name)
{
    if (state_plot_vars.isEmpty())
    {
        state_plot_vars.append(name);
    }
    else 
    {
        if (isDerivePlotVar(name) == false)
            state_plot_vars.append(name);
    }
}

void
Amr::deleteStatePlotVar (const aString& name)
{
    if (state_plot_vars.isNotEmpty())
        state_plot_vars.remove(name);
}

bool
Amr::isDerivePlotVar (const aString& name)
{
    for (ListIterator<aString> li(derive_plot_vars); li; ++li)
        if (li() == name)
            return true;
  
    return false;
}

void 
Amr::fillDerivePlotVarList()
{
    derive_plot_vars.clear();
    DeriveList& derive_lst = AmrLevel::get_derive_lst();
    List<DeriveRec>& dlist = derive_lst.dlist();
    for (ListIterator<DeriveRec> it(dlist); it; ++it)
        if (it().deriveType() == IndexType::TheCellType())
            derive_plot_vars.append(it().name());
}

void
Amr::clearDerivePlotVarList()
{
    derive_plot_vars.clear();
}

void
Amr::addDerivePlotVar (const aString& name)
{
    if (derive_plot_vars.isEmpty())
    {
        derive_plot_vars.append(name);
    }
    else 
    {
        if (isDerivePlotVar(name) == false)
            derive_plot_vars.append(name);
    }
}

void
Amr::deleteDerivePlotVar (const aString& name)
{
    if (derive_plot_vars.isNotEmpty())
        derive_plot_vars.remove(name);
}

Amr::~Amr ()
{
    if (level_steps[0] > last_checkpoint)
        checkPoint();

    if (level_steps[0] > last_plotfile)
        writePlotFile();

    levelbld->variableCleanUp();

    if (ParallelDescriptor::IOProcessor() && tar_and_rm_files)
    {
        //
        // Wait for the last chk and plt file tar_and_rm processes.
        //
        int status;

        if (pid_chk != -1) waitpid(pid_chk, &status, 0);
        if (pid_plt != -1) waitpid(pid_plt, &status, 0);
    }
}

void
Amr::setRecordGridInfo (const aString& filename)
{
    record_grid_info= true;
    gridlog.open(filename.c_str(),ios::out|ios::app);
    if (!gridlog.good())
        Utility::FileOpenFailed(filename);
}

void
Amr::setRecordRunInfo (const aString& filename)
{
    record_run_info= true;
    runlog.open(filename.c_str(),ios::out|ios::app);
    if (!runlog.good())
        Utility::FileOpenFailed(filename);
}

void
Amr::setRecordDataInfo (const aString& filename)
{
    datalog.open(filename.c_str(),ios::out|ios::app);
    if (!datalog.good())
        Utility::FileOpenFailed(filename);
}

void
Amr::setDtLevel (const Array<Real>& dt_lev)
{
    for (int i = 0; i <= finest_level; i++)
        dt_level[i] = dt_lev[i];
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
    return ok;
}

void
Amr::writePlotFile (const aString& root,
                    int            num)
{
    if (first_plotfile) 
    {
        first_plotfile = false;
        amr_level[0].setPlotVariables();
    }

    static RunStats stats("write_pltfile");

    stats.start();

    Real dPlotFileTime0 = ParallelDescriptor::second();

    const aString pltfile = Utility::Concatenate(root,num);

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "PLOTFILE: file = " << pltfile << endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "PLOTFILE: file = " << pltfile << '\n';
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!Utility::UtilCreateDirectory(pltfile, 0755))
            Utility::CreateDirectoryFailed(pltfile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    aString HeaderFileName = pltfile + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());

    int old_prec;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), ios::out|ios::trunc);

        if (!HeaderFile.good())
            Utility::FileOpenFailed(HeaderFileName);

        old_prec = HeaderFile.precision(30);
    }

    for (int k = 0; k <= finest_level; k++)
        amr_level[k].writePlotFile(pltfile, HeaderFile);

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Accumulate # of bytes written to header file.
        //
        RunStats::addBytes(VisMF::FileOffset(HeaderFile));

        HeaderFile.precision(old_prec);

        if (!HeaderFile.good())
            BoxLib::Error("Amr::writePlotFile() failed");
    }

    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real dPlotFileTime1 = ParallelDescriptor::second();
    Real dPlotFileTime  = dPlotFileTime1 - dPlotFileTime0;

    ParallelDescriptor::ReduceRealMax(dPlotFileTime,IOProc);

    stats.end();

    if (ParallelDescriptor::IOProcessor())
    {
        cout << "Write plotfile time = " << dPlotFileTime << "  seconds." << endl;

        if (tar_and_rm_files)
        {
            int status;

            if (pid_plt != -1) waitpid(pid_plt, &status, 0);

            aString cmd;

            cmd += "tar cf ";
            cmd += pltfile;
            cmd += ".tar ";
            cmd += pltfile;
            cmd += "; rm -rf ";
            cmd += pltfile;

            pid_plt = DoIt(cmd.c_str());
        }
    }
}

void
Amr::checkInput ()
{
    if (max_level < 0)
        BoxLib::Error("checkInput: max_level not set");
    //
    // Check that multigrid factor is a power of 2.
    //
    int k = blocking_factor;
    while ( k > 0 && (k%2 == 0) )
        k = k/2;
    if (k != 1)
        BoxLib::Error("Amr::checkInputs: multiGrid factor not a power of 2");
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
    // Check that domain size has a factor of blocking_factor.
    //
    for (i = 0; i < BL_SPACEDIM; i++)
    {
        int len = domain.length(i);
        if (len%blocking_factor != 0)
            BoxLib::Error("domain size not divisible by blocking_factor");
    }
    //
    // Check that max_grid_size has a factor of blocking_factor.
    //
    for (i = 0; i < max_level; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
        {
            int lratio = blocking_factor*ref_ratio[i][n];
            if (max_grid_size%lratio != 0)
                BoxLib::Error("max_grid_size not divisible by blocking_factor*ref_ratio");
        }
    }
    if (!Geometry::ProbDomain().ok())
        BoxLib::Error("checkInput: bad physical problem size");

    if (regrid_int[0] <= 0)
        BoxLib::Error("checkinput: regrid_int not defined");
}

void
Amr::init (Real strt_time,
           Real stop_time)
{
    if (!restart_file.isNull())
    {
        restart(restart_file);
    }
    else
    {
        initialInit(strt_time,stop_time);
        checkPoint();
        if (plot_int > 0 || plot_per > 0)
            writePlotFile();
    }
}

void
Amr::initialInit (Real strt_time,
                  Real stop_time)
{
    checkInput();
    //
    // Generate internal values from user-supplied values.
    //
    finest_level = 0;
    //
    // Init problem dependent data.
    //
    int  init = true;
    //
    // Populate integer array with name of `probin' file.
    //
    int probin_file_length = probin_file.length();

    Array<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
        probin_file_name[i] = probin_file[i];

    FORT_PROBINIT(&init,
                  probin_file_name.dataPtr(),
                  &probin_file_length,
                  Geometry::ProbLo(),
                  Geometry::ProbHi());

    cumtime = strt_time;
    //
    // Define base level grids.
    //
    defBaseLevel(strt_time);
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
        const int fact = sub_cycle ? ref_ratio[lev-1][0] : 1;

        dt0           /= Real(fact);
        dt_level[lev]  = dt0;
        dt_min[lev]    = dt_level[lev];
        n_cycle[lev]   = fact;
    }

    if (max_level > 0)
        bldFineLevels(strt_time);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].setTimeLevel(strt_time,dt_level[lev],dt_level[lev]);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_regrid(0,finest_level);
    //
    // Perform any special post_initialization operations.
    //
    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_init(stop_time);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        level_count[lev] = 0;
        level_steps[lev] = 0;
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "INITIAL GRIDS \n";
        printGridInfo(cout,0,finest_level);
    }

    if (record_grid_info && ParallelDescriptor::IOProcessor())
    {
        gridlog << "INITIAL GRIDS \n";
        printGridInfo(gridlog,0,finest_level);
    }

    station.init();
    station.findGrid(amr_level,geom);
}

void
Amr::restart (const aString& filename)
{
    Real dRestartTime0 = ParallelDescriptor::second();

    int i;

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "restarting calculation from file: " << filename << endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';
    //
    // Init problem dependent data.
    //
    int init = false;
    //
    // Populate integer array with name of `probin' file.
    //
    int probin_file_length = probin_file.length();

    Array<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
        probin_file_name[i] = probin_file[i];

    FORT_PROBINIT(&init,
                  probin_file_name.dataPtr(),
                  &probin_file_length,
                  Geometry::ProbLo(),
                  Geometry::ProbHi());
    //
    // Start calculation from given restart file.
    //
    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';
    //
    // Open the checkpoint header file for reading.
    //
    aString File = filename;

    File += '/';
    File += "Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    ifstream is;

    is.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());

    is.open(File.c_str(), ios::in);

    if (!is.good())
        Utility::FileOpenFailed(File);
    //
    // Read global data.
    //
    // Attempt to differentiate between old and new CheckPointFiles.
    //
    int     spdim;
    bool    new_checkpoint_format = false;
    aString first_line;

    first_line.getline(is);

    if (first_line == CheckPointVersion)
    {
        new_checkpoint_format = true;
        is >> spdim;
    }
    else
    {
        spdim = first_line.toInteger();
    }

    if (spdim != BL_SPACEDIM)
    {
        cerr << "Amr::restart(): bad spacedim = " << spdim << '\n';
        BoxLib::Abort();
    }

    is >> cumtime;
    int mx_lev;
    is >> mx_lev;
    if (max_level < mx_lev)
        BoxLib::Error("Amr::restart(): different max_level");

    is >> finest_level;

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

    for (i = 0; i <= mx_lev; i++) is >> n_cycle[i];
    for (i = 0; i <= mx_lev; i++) is >> level_steps[i];
    for (i = 0; i <= mx_lev; i++) is >> level_count[i];
    //
    // Set bndry conditions.
    //
    if (max_level > mx_lev)
    {
        for (i = mx_lev+1; i <= max_level; i++)
        {
            const int rat  = MaxRefRatio(i-1);
            const int mult = sub_cycle ? rat : 1;

            dt_level[i]    = dt_level[i-1]/Real(rat);
            n_cycle[i]     = mult;
            level_steps[i] = mult*level_steps[i-1];
            level_count[i] = 0;
        }
        if (!sub_cycle)
        {
            for (i = 0; i <= max_level; i++)
                dt_level[i] = dt_level[max_level];
        }
    }

    if (regrid_on_restart)
        level_count[0] = regrid_int[0];

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
    // Initialize local stats indicating this is a restart.
    //
    RunStats::readStats(is,true);
    //
    // Build any additional data structures.
    //
    for (lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_restart();

    station.init();
    station.findGrid(amr_level,geom);

    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real dRestartTime1 = ParallelDescriptor::second();
    Real dRestartTime  = dRestartTime1 - dRestartTime0;

    ParallelDescriptor::ReduceRealMax(dRestartTime,IOProc);

    if (ParallelDescriptor::IOProcessor())
        cout << "Restart time = " << dRestartTime << " seconds." << endl;

}

void
Amr::checkPoint ()
{
    static RunStats stats("write_chkfile");

    stats.start();
    //
    // In checkpoint files always write out FABs in NATIVE format.
    //
    FABio::Format thePrevFormat = FArrayBox::getFormat();

    FArrayBox::setFormat(FABio::FAB_NATIVE);

    Real dCheckPointTime0 = ParallelDescriptor::second();

    const aString ckfile = Utility::Concatenate(check_file_root,level_steps[0]);

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "CHECKPOINT: file = " << ckfile << endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "CHECKPOINT: file = " << ckfile << '\n';
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!Utility::UtilCreateDirectory(ckfile, 0755))
            Utility::CreateDirectoryFailed(ckfile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    aString HeaderFileName = ckfile + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());

    int old_prec, i;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), ios::out|ios::trunc);

        if (!HeaderFile.good())
            Utility::FileOpenFailed(HeaderFileName);

        old_prec = HeaderFile.precision(30);

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

    for (i = 0; i <= finest_level; i++)
        amr_level[i].checkPoint(ckfile, HeaderFile);

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Accumulate # of bytes written to header file.
        //
        RunStats::addBytes(VisMF::FileOffset(HeaderFile));

        HeaderFile.precision(old_prec);

        if (!HeaderFile.good())
            BoxLib::Error("Amr::checkpoint() failed");
    }
    //
    // Don't forget to reset FAB format.
    //
    FArrayBox::setFormat(thePrevFormat);

    stats.end();

    RunStats::dumpStats(HeaderFile);

    const int IOProc     = ParallelDescriptor::IOProcessorNumber();
    Real dCheckPointTime = ParallelDescriptor::second() - dCheckPointTime0;

    ParallelDescriptor::ReduceRealMax(dCheckPointTime,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        cout << "checkPoint() time = " << dCheckPointTime << " secs." << endl;

        if (tar_and_rm_files)
        {
            int status;

            if (pid_chk != -1) waitpid(pid_chk, &status, 0);

            aString cmd;

            cmd += "tar cf ";
            cmd += ckfile;
            cmd += ".tar ";
            cmd += ckfile;
            cmd += "; rm -rf ";
            cmd += ckfile;

            pid_chk = DoIt(cmd.c_str());
        }
    }
}

void
Amr::timeStep (int  level,
               Real time,
               int  iteration,
               int  niter)
{
    //
    // Time to regrid?
    //
    int lev_top = Min(finest_level, max_level-1);

    for (int i = level; i <= lev_top; i++)
    {
        const int old_finest = finest_level;

        if (level_count[i] >= regrid_int[i] && amr_level[i].okToRegrid())
        {
            regrid(i,time);

            for (int k = i; k <= finest_level; k++)
                level_count[k] = 0;

            if (old_finest < finest_level)
            {
                //
                // The new levels will not have valid time steps
                // and iteration counts.
                //
                for (int k = old_finest+1; k <= finest_level; k++)
                {
                    const int fact = sub_cycle ? MaxRefRatio(k-1) : 1;
                    dt_level[k]    = dt_level[k-1]/Real(fact);
                    n_cycle[k]     = fact;
                }
            }
        }
        if (old_finest > finest_level)
          lev_top = Min(finest_level, max_level-1);
    }

    //
    // check to see if should write plotfile
    // This routine is here so it is done after the restart regrid.
    //
    if (plotfile_on_restart && !(restart_file.isNull()) )
    {
	plotfile_on_restart = 0;
	writePlotFile();
    }
    //
    // Advance grids at this level.
    //
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "ADVANCE grids at level "
             << level
             << " with dt = "
             << dt_level[level]
             << endl;
    }
    Real dt_new = amr_level[level].advance(time,dt_level[level],iteration,niter);

    dt_min[level] = iteration == 1 ? dt_new : Min(dt_min[level],dt_new);

    level_steps[level]++;
    level_count[level]++;

    RunStats::addCells(level,amr_level[level].countCells());

    station.report(time+dt_level[level],level,amr_level[level]);
    //
    // Advance grids at higher level.
    //
    if (level < finest_level)
    {
        const int lev_fine = level+1;

        if (sub_cycle)
        {
            const int ncycle = n_cycle[lev_fine];

            for (int i = 1; i <= ncycle; i++)
                timeStep(lev_fine,time+(i-1)*dt_level[lev_fine],i,ncycle);
        }
        else
        {
            timeStep(lev_fine,time,1,1);
        }
    }

    amr_level[level].post_timestep();
}

void
Amr::coarseTimeStep (Real stop_time)
{
    //
    // Compute new dt.
    //
    if (level_steps[0] > 0)
        amr_level[0].computeNewDt(finest_level,
                                  sub_cycle,
                                  n_cycle,
                                  ref_ratio,
                                  dt_min,
                                  dt_level,
                                  stop_time);
    timeStep(0,cumtime,1,1);
    cumtime += dt_level[0];

    amr_level[0].postCoarseTimeStep(cumtime);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "\nSTEP = "
             << level_steps[0]
             << " TIME = "
             << cumtime
             << " DT = "
             << dt_level[0]
             << '\n'
             << endl;
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

    int check_test = 0;
    if (check_per > 0.0)
    {
        const int num_per = int((cumtime+.001*dt_level[0]) / check_per);
        const Real resid  = cumtime - num_per * check_per;

        if (resid < .001*dt_level[0]) check_test = 1;
    }

    if ((check_int > 0 && level_steps[0] % check_int == 0) || check_test == 1)
    {
        last_checkpoint = level_steps[0];
        checkPoint();
    }

    int plot_test = 0;
    if (plot_per > 0.0)
    {
        const int num_per = int((cumtime+.001*dt_level[0]) / plot_per);
        const Real resid  = cumtime - num_per * plot_per;

        if (resid < .001*dt_level[0]) plot_test = 1;
    }

    if ((plot_int > 0 && level_steps[0] % plot_int == 0) || plot_test == 1)
    {
        last_plotfile = level_steps[0];
        writePlotFile();
    }
}

void
Amr::defBaseLevel (Real strt_time)
{
    //
    // Check that base domain has even number of zones in all directions.
    //
    const Box& domain = geom[0].Domain();
    const int* d_len  = domain.length().getVect();

    for (int idir = 0; idir < BL_SPACEDIM; idir++)
        if (d_len[idir]%2 != 0)
            BoxLib::Error("defBaseLevel: must have even number of cells");
    //
    // Coarsening before we split the grids ensures that each resulting
    // grid will have an even number of cells in each direction.
    //
    BoxArray lev0(1);

    lev0.set(0,::coarsen(domain,2));
    //
    // Now split up into list of grids within max_grid_size limit.
    //
    lev0.maxSize(max_grid_size/2);
    //
    // Now refine these boxes back to level 0.
    //
    lev0.refine(2);
    //
    // Now build level 0 grids.
    //
    amr_level.set(0,(*levelbld)(*this,0,geom[0],lev0,strt_time));

    lev0.clear();
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
    static RunStats stats("regrid");

    stats.start();

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "REGRID: at level lbase = " << lbase << endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "REGRID: at level lbase = " << lbase << endl;
    //
    // Remove old-time grid space at highest level.
    //
    if (finest_level == max_level)
        amr_level[finest_level].removeOldData();
    //
    // Compute positions of new grids.
    //
    int             new_finest;
    Array<BoxArray> new_grid_places(max_level+1);

    if (lbase <= Min(finest_level,max_level-1))
      grid_places(lbase,time,new_finest, new_grid_places);

    bool regrid_level_zero =
        lbase == 0 && new_grid_places[0] != amr_level[0].boxArray();

    const int start = regrid_level_zero ? 0 : lbase+1;
    //
    // Reclaim old-time grid space for all remain levels > lbase.
    //
    for (int lev = start; lev <= finest_level; lev++)
        amr_level[lev].removeOldData();
    //
    // Reclaim all remaining storage for levels > new_finest.
    //
    for (int lev = new_finest+1; lev <= finest_level; lev++)
        amr_level.clear(lev);

    finest_level = new_finest;

    if (lbase == 0)
    {
        FabSet::FlushCache();
        MultiFab::FlushSICache();
        Geometry::FlushPIRMCache();
        DistributionMapping::FlushCache();

        if (!regrid_level_zero)
        {
            //
            // Since we're not regridding level 0 we've got to recache its
            // procmap.  Otherwise we'd have a problem if a higher level had
            // the same number of grids as does level 0.
            //
            MultiFab dummy(amr_level[0].boxArray(),1,0,Fab_noallocate);
        }
    }
    //
    // Define the new grids from level start up to new_finest.
    //
    for (int lev = start; lev <= new_finest; lev++)
    {
        //
        // Construct skeleton of new level.
        //
        AmrLevel* a = (*levelbld)(*this,
                                  lev,
                                  geom[lev],
                                  new_grid_places[lev],
                                  cumtime);
        if (initial)
        {
            //
            // We're being called on startup from bldFineLevels().
            //
            a->initData();
        }
        else if (amr_level.defined(lev))
        {
            //
            // Init with data from old structure then remove old structure.
            //
            a->init(amr_level[lev]);
        }
        else
        {
            a->init();
        }

        amr_level.clear(lev);

        amr_level.set(lev,a);
    }
    //
    // Build any additional data structures after grid generation.
    //
    for (int lev = 0; lev <= new_finest; lev++)
        amr_level[lev].post_regrid(lbase,new_finest);

    station.findGrid(amr_level,geom);
    //
    // Report creation of new grids.
    //
    if ((verbose || record_run_info) && ParallelDescriptor::IOProcessor())
    {
        for (int lev = start; lev <= finest_level; lev++)
        {
            const int  numgrids = amr_level[lev].numGrids();
            const long ncells   = amr_level[lev].countCells();
            const long ntot     = geom[lev].Domain().numPts();
            const Real frac     = 100.0*(Real(ncells) / Real(ntot));

            if (verbose)
            {
                cout << "   level "
                     << lev
                     << ": "
                     << numgrids
                     << " grids, "
                     << ncells
                     << " cells  = "
                     << frac
                     << " % of domain"
                     << endl;
            }
            if (record_run_info)
            {
                runlog << "   level "
                       << lev
                       << ": "
                       << numgrids
                       << " grids, "
                       << ncells
                       << " cells  = "
                       << frac
                       << " % of domain"
                       << endl;
            }
        }
    }
    if (record_grid_info && ParallelDescriptor::IOProcessor())
    {
        if (lbase == 0)
            gridlog << "STEP = " << level_steps[0] << ' ';

        gridlog << "TIME = "
                << time
                << " : REGRID  with lbase = "
                << lbase
                << endl;

        printGridInfo(gridlog,start,finest_level);
    }

    stats.end();
}

void
Amr::printGridInfo (ostream& os,
                    int      min_lev,
                    int      max_lev)
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray& bs      = amr_level[lev].boxArray();
        int             numgrid = bs.length();
        long            ncells  = amr_level[lev].countCells();
        long            ntot    = geom[lev].Domain().numPts();
        Real            frac    = 100.0*(Real(ncells) / Real(ntot));

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

            os << '\n';
        }
    }
    os << '\n';
}

void
Amr::ProjPeriodic (BoxDomain&      bd,
                   const Geometry& geom)
{
    Box domain(geom.Domain());
    //
    // blout will initially contain all of bd, periodic translates
    // will be added to it.
    //
    BoxList blout;  
    for (BoxDomainIterator bdi(bd); bdi; ++bdi)
        blout.add(bdi());
    //
    // blorig will be the original bd intersected with domain.
    //
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
                blout.join(tmp);
 
                if (rk != 0 && geom.isPeriodic(2))
                    blorig.shift(2,-rk*domain.length(2));
            }
            if (rj != 0 && geom.isPeriodic(1))
                blorig.shift(1,-rj*domain.length(1));
        }
        if (ri != 0 && geom.isPeriodic(0))
            blorig.shift(0,-ri*domain.length(0));
    }
    bd.clear();
    bd.add(blout);
}

void
Amr::grid_places (int              lbase,
                  Real             time,
                  int&             new_finest,
                  Array<BoxArray>& new_grids)
{
    int i, max_crse = Min(finest_level,max_level-1);

    if (lbase == 0)
    {
        //
        // Recalculate level 0 BoxArray in case max_grid_size has changed.
        // This is done exactly as in defBaseLev().
        //
        BoxArray lev0(1);

        lev0.set(0,::coarsen(geom[0].Domain(),2));
        //
        // Now split up into list of grids within max_grid_size limit.
        //
        lev0.maxSize(max_grid_size/2);
        //
        // Now refine these boxes back to level 0.
        //
        lev0.refine(2);

        new_grids[0] = lev0;
    }

    if (!grids_file.isNull())
    {
#define STRIP while( is.get() != '\n' )

        ifstream is(grids_file.c_str(),ios::in);

        if (!is.good())
            Utility::FileOpenFailed(grids_file);

        new_finest = Min(max_level,(finest_level+1));
        int in_finest;
        is >> in_finest;
        STRIP;
        new_finest = Min(new_finest,in_finest);
        int ngrid;
        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (i = 0; i < ngrid; i++)
            {
                Box bx;
                is >> bx;
                STRIP;
                if (lev > lbase)
                {
                    bx.refine(ref_ratio[lev-1]);
                    if (bx.longside() > max_grid_size)
                    {
                        cout << "Grid "
                             << bx
                             << " too large"
                             << '\n';
                        BoxLib::Error();
                    }
                    bl.append(bx);
                }
            }
            if (lev > lbase)
                new_grids[lev].define(bl);
        }
        is.close();
        return;
#undef STRIP
    }
    //
    // Construct problem domain at each level.
    //
    Array<IntVect> bf_lev(max_level); // Blocking factor at each level.
    Array<IntVect> rr_lev(max_level);
    Array<Box> pc_domain(max_level);  // Coarsened problem domain.
    for (i = lbase; i <= max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            bf_lev[i][n] = Max(1,blocking_factor/ref_ratio[i][n]);
    }
    for (i = lbase; i < max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            rr_lev[i][n] = (ref_ratio[i][n]*bf_lev[i][n])/bf_lev[i+1][n];
    }
    for (i = lbase; i <= max_crse; i++)
        pc_domain[i] = ::coarsen(geom[i].Domain(),bf_lev[i]);
    //
    // Construct proper nesting domains.
    //
    Array<BoxDomain> p_n(max_level);      // Proper nesting domain.
    Array<BoxDomain> p_n_comp(max_level); // Complement proper nesting domain.

    BoxDomain fd1;
    const BoxArray& bbase = amr_level[lbase].boxArray();
    for (i = 0; i < bbase.length(); i++)
        fd1.add(::coarsen(bbase[i],bf_lev[lbase]));
    p_n_comp[lbase].complementIn(pc_domain[lbase],fd1);
    p_n_comp[lbase].accrete(n_proper);
    Geometry tmp2(pc_domain[lbase]);
    Amr::ProjPeriodic(p_n_comp[lbase], tmp2);
    p_n_comp[lbase].minimize();
    p_n[lbase].complementIn(pc_domain[lbase],p_n_comp[lbase]);
    p_n[lbase].minimize();
    fd1.clear();

    for (i = lbase+1; i <= max_crse;  i++)
    {
        BoxList bl;
        for (BoxDomainIterator bdi(p_n_comp[i-1]); bdi; ++bdi)
            bl.add(refine(bdi(),rr_lev[i-1]));
        p_n_comp[i].clear();
        p_n_comp[i].add(bl);
        p_n_comp[i].accrete(n_proper);
        Geometry tmp3(pc_domain[i]);
        Amr::ProjPeriodic(p_n_comp[i], tmp3);
        p_n[i].complementIn(pc_domain[i],p_n_comp[i]);
        p_n[i].minimize();
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
        const BoxArray& old_grids = amr_level[levc].boxArray();
        int ngrow = 0;
        BoxArray ba_proj;
        if (levf < new_finest)
        {
            BoxList blst(old_grids);
            ba_proj.define(new_grids[levf+1]);
            ba_proj.coarsen(ref_ratio[levf]);
            ba_proj.grow(n_proper);
            ba_proj.coarsen(ref_ratio[levc]);
            while (!blst.contains(ba_proj))
            {
                blst.accrete(1);
                ngrow++;
            }
        }
        TagBoxArray tags(old_grids,n_error_buf[levc]+ngrow);
        amr_level[levc].errorEst(tags,
                                 TagBox::CLEAR,TagBox::SET,time,
				 n_error_buf[levc],ngrow);
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

            BoxList bl_tagged;
            for (int i = 0; i < new_grids[levf+1].length(); i++)
                bl_tagged.add(coarsen(new_grids[levf+1][i],ref_ratio[levf]));
            //
            // This grows the boxes by nerr if they touch the edge of the
            // domain in preparation for them being shrunk by nerr later.
            // We want the net effect to be that grids are NOT shrunk away
            // from the edges of the domain.
            //
            for (BoxListIterator blt(bl_tagged); blt; ++blt)
            {
                for (int idir = 0; idir < BL_SPACEDIM; idir++)
                {
                    if (blt().smallEnd(idir) == geom[levf].Domain().smallEnd(idir))
                        bl_tagged[blt].growLo(idir,nerr);
                    if (blt().bigEnd(idir) == geom[levf].Domain().bigEnd(idir))
                        bl_tagged[blt].growHi(idir,nerr);
                }
            }

            Box mboxF = Box(bl_tagged.minimalBox()).grow(1);
            BoxList blFcomp = complementIn(mboxF,bl_tagged);
            blFcomp.accrete(nerr);
            BoxList blF = complementIn(mboxF,blFcomp);

            BoxArray baF(blF);
            baF.grow(n_proper);
            //
            // We need to do this in case the error buffering at
            // levc will not be enough to cover the error buffering
            // at levf which was just subtracted off.
            //
            for (int idir=0; idir < BL_SPACEDIM; idir++) 
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
        tags.buffer(n_error_buf[levc]);
        //
        // Coarsen the taglist by blocking_factor.
        //
        int bl_max = 0;
        for (int n=0; n<BL_SPACEDIM; n++)
            bl_max = Max(bl_max,bf_lev[levc][n]);
        if (bl_max > 1) 
            tags.coarsen(bf_lev[levc]);
        //
        // Remove or add tagged points which violate/satisfy additional 
        // user-specified criteria.
        //
        amr_level[levc].manual_tags_placement(tags, bf_lev[levc]);
        //
        // Map tagged points through periodic boundaries, if any.
        //
        Geometry tmpgeom(pc_domain[levc]);
        tags.mapPeriodic(tmpgeom);
        //
        // Remove cells outside proper nesting domain for this level.
        //
        tags.setVal(p_n_comp[levc],TagBox::CLEAR);
        //
        // Create initial cluster containing all tagged points.
        //
        long     len = 0;
        IntVect* pts = tags.collate(len);

        tags.clear();

        if (len > 0)
        {
            //
            // Created new level, now generate efficient grids.
            //
            new_finest = Max(new_finest,levf);
            //
            // Construct initial cluster.
            //
            ClusterList clist(pts,len);
            clist.chop(grid_eff);
            clist.intersect(p_n[levc]);
            //
            // Efficient properly nested Clusters have been constructed
            // now generate list of grids at level levf.
            //
            BoxList new_bx;
            clist.boxList(new_bx);

            new_bx.refine(bf_lev[levc]);

            new_bx.minimize();

            IntVect largest_grid_size;
            for (int n = 0; n < BL_SPACEDIM; n++)
              largest_grid_size[n] = max_grid_size / ref_ratio[levc][n];
            //
            // Ensure new grid boxes are at most max_grid_size in index dirs.
            //
            new_bx.maxSize(largest_grid_size);
            //
            // Refine up to levf.
            //
            new_bx.refine(ref_ratio[levc]);

            if (!new_bx.isDisjoint())
            {
                cout << "WARNING: new grids at level "
                     << levf
                     << " not disjoint:\n"
                     << new_bx
                     << '\n';
            }
            new_grids[levf].define(new_bx);
        }
    }
}

void
Amr::bldFineLevels (Real strt_time)
{
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
    //
    bool grids_the_same;

    do
    {
        for (int i = 0; i <= finest_level; i++)
            grids[i] = amr_level[i].boxArray();

        regrid(0,strt_time,true);

        grids_the_same = true;

        for (int i = 0; i <= finest_level && grids_the_same; i++)
            if (!(grids[i] == amr_level[i].boxArray()))
                grids_the_same = false;
    }
    while (!grids_the_same);
}
