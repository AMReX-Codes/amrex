//BL_COPYRIGHT_NOTICE

//
// $Id: Amr.cpp,v 1.22 1997-12-11 06:14:02 lijewski Exp $
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
#include <Amr_auxil.H>
#include <ParallelDescriptor.H>
#include <Utility.H>

#include <DistributionMapping.H>

#ifdef BL_USE_NEW_HFILES
#include <cstdio>
using std::ifstream;
using std::ios;
#else
#include <stdio.h>
#endif

#ifndef BL_USE_SETBUF
#define setbuf pubsetbuf
#endif

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
             int            lev)
{
    return amr_level[lev].derive(name,time);
}

FArrayBox*
Amr::derive (const Box&     b,
             const aString& name,
             Real           time,
             int            lev)
{
    return amr_level[lev].derive(b,name,time);
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
    // Setup Geometry from ParmParse file.  May be needed for variableSetup or
    // even getLevelBld.
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
    max_level     = -1;
    record_run_info   = false;
    record_grid_info  = false;

    grid_eff	      = 0.7;
    blocking_factor   = 1;
    last_checkpoint   = 0;
    last_plotfile     = 0;
    plot_int          = -1;
    n_proper          = 1;

#if (BL_SPACEDIM == 2)
    max_grid_size     = 120;
#else
    max_grid_size     = 32;
#endif

    int i;
    for (i = 0; i < BL_SPACEDIM; i++) isPeriodic[i] = false;

    ParmParse pp("amr");
    //
    // Check for command line flags.
    //
    trace = pp.contains("trace");
    debug = pp.contains("debug");
    silent = pp.contains("silent");
    verbose = pp.contains("verbose");
    sub_cycle = true;
    if (pp.contains("nosub")) sub_cycle = false;

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
	dt_level[i] = 0.0;
	level_steps[i] = 0;
	level_count[i] = 0;
	regrid_int[i] = 0;
	n_cycle[i] = 0;
	dt_min[i] = 0.0;
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
        BoxLib::Error("Must only specify amr.check_int OR amr.check_per");
    else if (got_check_per == 1)
        BoxLib::Warning("Specifying amr.check_per will change the time step");

    plot_file_root = "plt";
    pp.query("plot_file",plot_file_root);

    plot_int = -1;
    int got_plot_int = pp.query("plot_int",plot_int);

    plot_per = -1.0;
    int got_plot_per = pp.query("plot_per",plot_per);

    if (got_plot_int == 1 && got_plot_per == 1)
        BoxLib::Error("Must only specify amr.plot_int OR amr.plot_per");
    else if (got_plot_per == 1)
        BoxLib::Warning("Specifying amr.plot_per will change the time step");

    pp.query("max_grid_size",max_grid_size);
    pp.query("n_proper",n_proper);
    pp.get("blocking_factor",blocking_factor);
    pp.get("grid_eff",grid_eff);

    pp.getarr("n_error_buf",n_error_buf,0,max_level);
    //
    // Read in the refinement ratio IntVects as integer BL_SPACEDIM-tuples.
    //
    if (max_level > 0)
    {
      int nratios_vect = max_level*BL_SPACEDIM;
      Array<int> ratios_vect(nratios_vect);

      int got_vect = pp.queryarr("ref_ratio_vect",ratios_vect,0,nratios_vect);

      Array<int> ratios(max_level);

      int got_int = pp.queryarr("ref_ratio",ratios,0,max_level);
   
      if (got_int == 1 && got_vect == 1)
          BoxLib::Warning("Only input *either* ref_ratio or ref_ratio_vect");
      else if (got_vect == 1)
      {
          int k = 0;
          for( i = 0; i < max_level; i++ )
          {
              for( int n=0; n<BL_SPACEDIM; n++,k++ )
                  ref_ratio[i][n] = ratios_vect[k];
          }

      }
      else if (got_int == 1)
      {
          for( i = 0; i < max_level; i++ )
          {
              for( int n=0; n<BL_SPACEDIM; n++)
                  ref_ratio[i][n] = ratios[i];
          }
      }
      else
          BoxLib::Error("Must input *either* ref_ratio or ref_ratio_vect");
    }
    //
    // Read computational domain and set geometry.
    //
    Array<int> n_cell(BL_SPACEDIM);
    pp.getarr("n_cell",n_cell,0,BL_SPACEDIM);
    assert(n_cell.length() == BL_SPACEDIM);
    IntVect lo(IntVect::TheZeroVector()), hi(n_cell);
    hi -= IntVect::TheUnitVector();
    Box index_domain(lo,hi);
    for (i = 0; i <= max_level; i++)
    {
	geom[i].define(index_domain);
	if (i < max_level) index_domain.refine(ref_ratio[i]);
    }
    //
    // Now define offset for CoordSys.
    //
    Real offset[BL_SPACEDIM];
    for (i = 0; i < BL_SPACEDIM; i++)
    {
	Real delta = Geometry::ProbLength(i)/(Real)n_cell[i];
	offset[i] = Geometry::ProbLo(i) + delta*lo[i];
    }
    CoordSys::SetOffset(offset);
    //
    // Set regrid interval.
    //
    int ri;
    pp.get("regrid_int",ri);
    int k;
    for (k = 0; k <= max_level; k++) regrid_int[k] = ri;
    if (!sub_cycle)
    {
        //
        // Must adjust regridding trigger.
        //
	int factor = 1;
	for (k = max_level-1; k >= 0; k--)
        {
	    factor *= MaxRefRatio(k);
	    regrid_int[k] = ri*factor;
	}
    }
}

Amr::~Amr ()
{
    if (level_steps[0] > last_checkpoint) checkPoint();
    if (level_steps[0] > last_plotfile) writePlotFile();

    levelbld->variableCleanUp();
}

void
Amr::setRecordGridInfo (const aString& filename)
{
    record_grid_info= true;
    gridlog.open(filename.c_str(),ios::out);
    if (!gridlog.good())
        Utility::FileOpenFailed(filename);
}

void
Amr::setRecordRunInfo (const aString& filename)
{
    record_run_info= true;
    runlog.open(filename.c_str(),ios::out);
    if (!runlog.good())
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

static
aString
Concatenate (const aString& root,
             int            num)
{
    aString result = root;
    char buf[sizeof(int) + 1];
    sprintf(buf, "%04d", num);
    result += buf;
    return result;
}

#ifdef BL_PARALLEL_IO
void
Amr::writePlotFile (const aString& root,
                    int            num)
{
    aString pltfile = Concatenate(root, num);

    if (trace && ParallelDescriptor::IOProcessor())
    {
        cout << "PLOTFILE: file = " << pltfile << endl;
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "PLOTFILE: file = " << pltfile << '\n';
    }
    if (!Utility::CreateDirectory(pltfile, 0755))
        Utility::CreateDirectoryFailed(pltfile);

    aString HeaderFileName = pltfile + "/Header";

#ifdef BL_USE_SETBUF
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
#endif
    ofstream HeaderFile;

#ifdef BL_USE_SETBUF
    HeaderFile.rdbuf()->setbuf(io_buffer.dataPtr(), io_buffer.length());
#endif
    int old_prec = HeaderFile.precision(15);
    //
    // Only the IOProcessor() writes to the header file.
    //
    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.open(HeaderFileName.c_str(), ios::out|ios::trunc);

        if (!HeaderFile.good())
            Utility::FileOpenFailed(HeaderFileName);
    }

    for (int k = 0; k <= finest_level; k++)
    {
        RunStats write_pltfile_stats("write_pltfile", k);
        write_pltfile_stats.start();
	amr_level[k].writePlotFile(pltfile, HeaderFile);
        write_pltfile_stats.end();
    }
    //
    // Accumulate # of bytes written to header file.
    //
    RunStats::addBytes(VisMF::FileOffset(HeaderFile));

    HeaderFile.precision(old_prec);

    if (!HeaderFile.good())
        BoxLib::Error("Amr::writePlotFile() failed");
}

#else

void
Amr::writePlotFile (const aString& root,
                    int            num)
{
    aString pltfile = Concatenate(root, num);

    if (trace && ParallelDescriptor::IOProcessor())
    {
        cout << "PLOTFILE: file = " << pltfile << endl;
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "PLOTFILE: file = " << pltfile << '\n';
    }

#ifdef BL_USE_SETBUF
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
#endif
    ofstream os;

#ifdef BL_USE_SETBUF
    os.rdbuf()->setbuf(io_buffer.dataPtr(), io_buffer.length());
#endif
    os.open(pltfile.c_str(), ios::out|ios::trunc);

    if (!os.good())
        Utility::FileOpenFailed(pltfile);

    for (int k = 0; k <= finest_level; k++)
    {
        RunStats write_pltfile_stats("write_pltfile", k);
        write_pltfile_stats.start();
        amr_level[k].writePlotFile(os);
        write_pltfile_stats.end();
    }
    //
    // Accumulate total # of bytes.
    //
    RunStats::addBytes(VisMF::FileOffset(os));
}
#endif /*BL_PARALLEL_IO*/

void
Amr::checkInput ()
{
    if (max_level < 0)
	BoxLib::Error("checkInput: max_level not set");
    //
    // Check that multigrid factor is a power of 2.
    //
    int k = blocking_factor;
    while ( k > 0 && (k%2 == 0) ) k = k/2;
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
Amr::init ()
{
    if (!restart_file.isNull())
    {
	restart(restart_file);
    }
    else
    {
	initialInit();
	checkPoint();
	if (plot_int > 0)
            writePlotFile();
    }
}

void
Amr::initialInit ()
{
    checkInput();
    //
    // Generate internal values from user-supplied values.
    //
    finest_level = 0;
    //
    // Init problem dependent data.
    //
    Real strt_time = 0;
    int  init = true;

    FORT_PROBINIT (&init);
    cumtime = strt_time;
    //
    // Define base level grids.
    //
    defBaseLevel(strt_time);

    if (max_level > 0) bldFineLevels(strt_time);
    //
    // Build any additional data structures after grid generation.
    //
    int lev;
    for (lev = 0; lev <= finest_level; lev++)
	amr_level[lev].post_regrid(0,finest_level);
    //
    // Compute dt and set time levels of all grid data.
    //
    amr_level[0].computeInitialDt(finest_level,sub_cycle,
				  n_cycle,ref_ratio,dt_level);
    for (lev = 0; lev <= finest_level; lev++)
	amr_level[lev].setTimeLevel(strt_time,dt_level[lev],dt_level[lev]);
    //
    // Perform any special post_initialization operations.
    //
    for (lev = 0; lev <= finest_level; lev++)
	amr_level[lev].post_init();

    for (lev = 0; lev <= finest_level; lev++)
    {
	level_count[lev] = 0;
	level_steps[lev] = 0;
    }

    if (record_grid_info)
    {
	gridlog << "INITIAL GRIDS \n";
	printGridInfo(gridlog,0,finest_level);
    }
}

void
Amr::restart (const aString& filename)
{
    int i;

    if (trace && ParallelDescriptor::IOProcessor())
    {
	cout << "restarting calculation from file " << filename << endl;
    }
    //
    // Init problem dependent data.
    //
    int init = false;
    FORT_PROBINIT (&init);
    //
    // Start calculation from given restart file.
    //
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
	runlog << "RESTART from file = " << filename << '\n';
    }
    //
    // Open the checkpoint header file for reading.
    //
    aString File = filename;

#ifdef BL_PARALLEL_IO
    File += '/';
    File += "Header";
#endif /*BL_PARALLEL_IO*/

#ifdef BL_USE_SETBUF
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
#endif
    ifstream is;

#ifdef BL_USE_SETBUF
    is.rdbuf()->setbuf(io_buffer.dataPtr(), io_buffer.length());
#endif
    is.open(File.c_str(), ios::in);

    if (!is.good())
        Utility::FileOpenFailed(File);
    //
    // Read global data.
    //
    int spdim;
    is >> spdim;
    if (spdim != BL_SPACEDIM)
    {
	cerr << "restart: bad spacedim = " << spdim << '\n';
	BoxLib::Abort();
    }
    is >> cumtime;
    int mx_lev;
    is >> mx_lev;
    if (max_level < mx_lev)
	BoxLib::Error("restart: different max_level");

    is >> finest_level;

    for (i = 0; i <= mx_lev; i++) is >> geom[i];
    for (i = 0; i <  mx_lev; i++) is >> ref_ratio[i];

    for (i = 0; i <= mx_lev; i++) is >> dt_level[i];
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
	    int rat = MaxRefRatio(i-1);
	    dt_level[i] = dt_level[i-1]/Real(rat);
            //
            // NEED SUB_CYCLE.
            //
	    if (sub_cycle)
            {
		n_cycle[i] = rat;
		level_steps[i] = rat*level_steps[i-1];
	    }
            else
            {
		n_cycle[i] = 1;
		level_steps[i] = level_steps[i-1];
	    }
	    level_count[i] = 0;
	}
	if (!sub_cycle)
        {
	    for (i = 0; i <= max_level; i++)
		dt_level[i] = dt_level[max_level];
	}
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
    // Initialize local stats indicating this is a restart.
    //
    RunStats::readStats(is,true);
    //
    // Build any additional data structures.
    //
    for (lev = 0; lev <= finest_level; lev++)
	amr_level[lev].post_restart();
}

#ifdef BL_PARALLEL_IO
void
Amr::checkPoint ()
{
    aString ckfile = Concatenate(check_file_root, level_steps[0]);

    if (trace && ParallelDescriptor::IOProcessor())
    {
        cout << "CHECKPOINT: file = " << ckfile << endl;
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "CHECKPOINT: file = " << ckfile << '\n';
    }
    if (!Utility::CreateDirectory(ckfile, 0755))
        Utility::CreateDirectoryFailed(ckfile);

    aString HeaderFileName = ckfile + "/Header";

#ifdef BL_USE_SETBUF
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
#endif
    ofstream HeaderFile;

#ifdef BL_USE_SETBUF
    HeaderFile.rdbuf()->setbuf(io_buffer.dataPtr(), io_buffer.length());
#endif
    int old_prec = HeaderFile.precision(15), i;
    //
    // Only the IOProcessor() writes to the header file.
    //
    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.open(HeaderFileName.c_str(), ios::out|ios::trunc);

        if (!HeaderFile.good())
            Utility::FileOpenFailed(HeaderFileName);

        HeaderFile << BL_SPACEDIM  << '\n'
                   << cumtime      << '\n'
                   << max_level    << '\n'
                   << finest_level << '\n';
        //
        // Write out problem domain.
        //
        for (i = 0; i <= max_level; i++) HeaderFile << geom[i] << ' ';
        HeaderFile << '\n';
        for (i = 0; i < max_level; i++) HeaderFile << ref_ratio[i] << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << dt_level[i] << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << n_cycle[i] << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << level_steps[i] << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << level_count[i] << ' ';
        HeaderFile << '\n';
    }

    for (i = 0; i <= finest_level; i++)
    {
        RunStats write_chkfile_stats("write_chkfile", i);
        write_chkfile_stats.start();
	amr_level[i].checkPoint(ckfile, HeaderFile);
        write_chkfile_stats.end();
    }

    RunStats::dumpStats(HeaderFile);
    //
    // Accumulate # of bytes written to header file.
    //
    RunStats::addBytes(VisMF::FileOffset(HeaderFile));

    HeaderFile.precision(old_prec);

    if (!HeaderFile.good())
        BoxLib::Error("Amr::checkpoint() failed");
}

#else

void
Amr::checkPoint ()
{
    aString ckfile = Concatenate(check_file_root, level_steps[0]);

    if (trace && ParallelDescriptor::IOProcessor())
    {
        cout << "CHECKPOINT: file = " << ckfile << endl;
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "CHECKPOINT: file = " << ckfile << '\n';
    }

#ifdef BL_USE_SETBUF
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
#endif
    ofstream os;

#ifdef BL_USE_SETBUF
    os.rdbuf()->setbuf(io_buffer.dataPtr(), io_buffer.length());
#endif
    os.open(ckfile.c_str(), ios::out|ios::trunc);

    if (!os.good())
        Utility::FileOpenFailed(ckfile);

    int old_prec = os.precision(15), i;

    if (ParallelDescriptor::IOProcessor())
    {
        os << BL_SPACEDIM  << '\n'
           << cumtime      << '\n'
           << max_level    << '\n'
           << finest_level << '\n';
        //
        // Write out problem domain.
        //
        for (i = 0; i <= max_level; i++) os << geom[i] << ' ';
        os << '\n';
        for (i = 0; i < max_level; i++) os << ref_ratio[i] << ' ';
        os << '\n';
        for (i = 0; i <= max_level; i++) os << dt_level[i] << ' ';
        os << '\n';
        for (i = 0; i <= max_level; i++) os << n_cycle[i] << ' ';
        os << '\n';
        for (i = 0; i <= max_level; i++) os << level_steps[i] << ' ';
        os << '\n';
        for (i = 0; i <= max_level; i++) os << level_count[i] << ' ';
        os << '\n';
    }

    for (i = 0; i <= finest_level; i++)
    {
        RunStats write_chkfile_stats("write_chkfile", i);
        write_chkfile_stats.start();
	amr_level[i].checkPoint(os);
        write_chkfile_stats.end();
    }

    RunStats::dumpStats(os);
    //
    // Accumulate # of bytes.
    //
    RunStats::addBytes(VisMF::FileOffset(os));

    os.precision(old_prec);

    if (!os.good())
        BoxLib::Error("Amr::checkpoint() failed");
}
#endif /*BL_PARALLEL_IO*/

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
	if (level_count[i] >= regrid_int[i])
        {
	    int old_finest = finest_level;
	    regrid(i,time);
            int k;
	    for (k = i; k <= finest_level; k++)
                level_count[k] = 0;

	    if (old_finest < finest_level)
            {
                //
                // The new levels will not have valid time steps
                // and iteration counts.
                //
		for (k = old_finest+1; k <= finest_level; k++)
                {
                    if (sub_cycle)
                    {
                        dt_level[k] = dt_level[k-1]/Real(MaxRefRatio(k-1));
                        n_cycle[k] = MaxRefRatio(k-1);
                    }
                    else
                    {
                        dt_level[k] = dt_level[k-1] ;
                        n_cycle[k] = 1 ;
                    }
		}
	    }
	}
    }
    //
    // Advance grids at this level.
    //
    if (trace && ParallelDescriptor::IOProcessor())
    {
	cout << "ADVANCE grids at level "
             << level
             << " with dt = "
             << dt_level[level]
             << endl;
    }
    Real dt_new = amr_level[level].advance(time,
                                           dt_level[level],
                                           iteration,
                                           niter);
    if (iteration == 1)
    {
	dt_min[level] = dt_new;
    }
    else
    {
	dt_min[level] = Min(dt_min[level],dt_new);
    }
    level_steps[level]++;
    level_count[level]++;
    RunStats::addCells(level,amr_level[level].countCells());
    //
    // Advance grids at higher level.
    //
    if (level < finest_level)
    {
	if (sub_cycle)
        {
	    int lev_fine = level+1;
	    int ncycle = n_cycle[lev_fine];
	    for (int i = 1; i <= ncycle; i++)
		timeStep(lev_fine,time + (i-1)*dt_level[lev_fine],i,ncycle);
	}
        else
        {
	    int lev_fine = level+1;
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
    {
	amr_level[0].computeNewDt(finest_level,sub_cycle,n_cycle,ref_ratio,
				  dt_min,dt_level,stop_time);
    }

    timeStep(0,cumtime,1,1);
    cumtime += dt_level[0];

    amr_level[0].postCoarseTimeStep(cumtime);

    if (!silent && ParallelDescriptor::IOProcessor())
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
      int num_per = (cumtime+.001*dt_level[0]) / check_per;
      Real resid = cumtime - num_per * check_per;
      if (resid < .001*dt_level[0])
          check_test = 1;
    }

    if ((check_int > 0 && level_steps[0] % check_int == 0) || check_test == 1)
    {
	last_checkpoint = level_steps[0];
	checkPoint();
    }

    int plot_test = 0;
    if (plot_per > 0.0)
    {
      int num_per = (cumtime+.001*dt_level[0]) / plot_per;
      Real resid = cumtime - num_per * plot_per;
      if (resid < .001*dt_level[0])
          plot_test = 1;
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
    const int* d_len = domain.length().getVect();
    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
	if (d_len[idir]%2 != 0)
        {
	    BoxLib::Error("defBaseLevel: domain not have even number of cells");
	}
    }
    //
    // Coarsening before we split the grids ensures that
    // each resulting grid will have an even number of
    // cells in each direction.
    //
    BoxArray  lev0(1);
    lev0.set(0, ::coarsen(domain,2));
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
             Real time)
{
    int new_finest;

    if (!silent && ParallelDescriptor::IOProcessor())
    {
	cout << "REGRID: at level lbase = " << lbase << endl;
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
	runlog << "REGRID: at level lbase = " << lbase << '\n';
    }
    //
    // Remove old-time grid space at highest level.
    //
    if (finest_level == max_level)
	amr_level[finest_level].removeOldData();
    //
    // Compute positions of new grids.
    //
    Array<BoxArray> new_grid_places(max_level+1);
    assert(new_grid_places.ready());

    if (lbase<=Min(finest_level,max_level-1))
      grid_places(lbase,time,new_finest, new_grid_places);
    //
    // Reclaim old-time grid space for all remain levels > lbase.
    //
    int lev;
    for (lev = lbase+1; lev <= finest_level; lev++)
	amr_level[lev].removeOldData();
    //
    // Reclaim all remaining storage for levels > new_finest.
    //
    for (lev = new_finest+1; lev <= finest_level; lev++)
	amr_level.clear(lev);

    finest_level = new_finest;
    //
    // Flush the grid -> processor map cache. 
    //
    DistributionMapping::FlushCache();
    //
    // Define the new grids from level lbase+1 up to new_finest.
    //
    for (lev = lbase+1; lev <= new_finest; lev++)
    {
        //
        // Construct skeleton of new level.
        //
	AmrLevel *a = (*levelbld)(*this,lev,geom[lev],
				  new_grid_places[lev],cumtime);
	assert(a);
        //
        // Init with data from old structure then remove old structure.
        //
	if (!amr_level.defined(lev))
        {
	    a->init();
	}
        else
        {
	    a->init(amr_level[lev]);
	    amr_level.clear(lev);
	}
        //
        // Install new structure.
        //
	amr_level.set(lev,a);
    }
    //
    // Build any additional data structures after grid generation.
    //
    for (lev = 0; lev <= new_finest; lev++)
	amr_level[lev].post_regrid(lbase,new_finest);
    //
    // Report creation of new grids.
    //
    if (!silent || record_run_info)
    {
        int lev;
        for (lev = lbase+1; lev <= finest_level; lev++)
        {
	    int numgrids= amr_level[lev].numGrids();
	    long ncells = amr_level[lev].countCells();
	    long ntot = geom[lev].Domain().numPts();
	    Real frac = 100.0*(Real(ncells) / Real(ntot));
	    if (!silent)
            {
                if (ParallelDescriptor::IOProcessor())
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
	    }
	    if (record_run_info)
            {
                if (ParallelDescriptor::IOProcessor())
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
                           << '\n';
                }
	    }
	}
    }

    if (record_grid_info && ParallelDescriptor::IOProcessor())
    {
	gridlog << "REGRID  with lbase = " << lbase << '\n';
	printGridInfo(gridlog,lbase+1,finest_level);
    }
}

void
Amr::printGridInfo (ostream& os,
                    int      min_lev,
                    int      max_lev) {
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
	const BoxArray& bs = amr_level[lev].boxArray();
	int numgrid = bs.length();
	long ncells = amr_level[lev].countCells();
	long ntot = geom[lev].Domain().numPts();
	Real frac = 100.0*(Real(ncells) / Real(ntot));

	os << "  Level "
           << lev
           << ' '
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
	    os << ' ' << lev << ": " << b << ' ';
	    for (int i = 0; i < BL_SPACEDIM; i++)
		os << b.length(i) << ' ';
	    os << '\n';
	}
    }
    os << '\n';
}

void
proj_periodic (BoxDomain&      bd,
               const Geometry& geom)
{
    Box domain( geom.Domain() );
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
    for (ri=nist; ri<=niend; ri++)
    {
        if (ri!=0 && !geom.isPeriodic(0))
            continue;
        if (ri!=0 && geom.isPeriodic(0))
            blorig.shift(0,ri*domain.length(0));
        for (rj=njst; rj<=njend; rj++)
        {
            if (rj!=0 && !geom.isPeriodic(1))
                continue;
            if (rj!=0 && geom.isPeriodic(1))
                blorig.shift(1,rj*domain.length(1));
            for (rk=nkst; rk<=nkend; rk++)
            {
                if (rk!=0 && !geom.isPeriodic(2))
                    continue;
                if (rk!=0 && geom.isPeriodic(2))
                    blorig.shift(2,rk*domain.length(2));

                BoxList tmp(blorig);
                tmp.intersect(domain);
                blout.join(tmp);
 
                if (rk!=0 && geom.isPeriodic(2))
                    blorig.shift(2,-rk*domain.length(2));
            }
            if (rj!=0 && geom.isPeriodic(1))
                blorig.shift(1,-rj*domain.length(1));
        }
        if (ri!=0 && geom.isPeriodic(0))
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
    int i;
    int  max_crse = Min(finest_level,max_level-1);

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
    Array<BoxDomain> p_n(max_level);   // proper nesting domain
    Array<BoxDomain> p_n_comp(max_level);   // complement of proper nesting domain

    BoxDomain fd1;
    const BoxArray& bbase = amr_level[lbase].boxArray();
    for (i = 0; i < bbase.length(); i++)
    {
	Box tmp(::coarsen(bbase[i],bf_lev[lbase]));
	fd1.add(tmp);
    }

    p_n_comp[lbase].complementIn(pc_domain[lbase],fd1);
    p_n_comp[lbase].accrete(n_proper);
    Geometry tmp2(pc_domain[lbase]);
    proj_periodic( p_n_comp[lbase], tmp2 );
    p_n_comp[lbase].minimize();
    p_n[lbase].complementIn(pc_domain[lbase],p_n_comp[lbase]);
    p_n[lbase].minimize();
    fd1.clear();

    for (i = lbase+1; i <= max_crse;  i++)
    {
        BoxList bl;
	for (BoxDomainIterator bdi(p_n_comp[i-1]); bdi; ++bdi)
        {
	  bl.add(refine(bdi(),rr_lev[i-1]));
	}
	p_n_comp[i].clear();
	p_n_comp[i].add(bl);
	p_n_comp[i].accrete(n_proper);
	Geometry tmp3(pc_domain[i]);
	proj_periodic( p_n_comp[i], tmp3 );
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
	amr_level[levc].errorEst(tags,TagBox::CLEAR,TagBox::SET,time);
        //
        // If new grids have been constructed above this level, project
        // those grids down and tag cells on intersections to ensure
        // proper nesting.
        //
	if (levf < new_finest)
	    tags.setVal(ba_proj,TagBox::SET);
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
	// Map tagged points through periodic boundaries, if any.
        //
	Geometry tmpgeom(pc_domain[levc]);
	tags.mapPeriodic(tmpgeom);
        //
        // Merge tagged points on overlap, remove redundant tags.
        //
	tags.mergeUnique();
        //
        // Remove cells outside proper nesting domain for this level.
        //
	tags.setVal(p_n_comp[levc],TagBox::CLEAR);
        //
        // Create initial cluster containing all tagged points.
        //
	Array<IntVect> *pts = tags.colate();

	tags.clear();

	if (pts->length() == 0)
        {
            //
            // No tagged points, no new level -- do cleanup.
            //
	    delete pts;
	}
        else
        {
            //
            // Created new level, now generate efficient grids.
            //
	    new_finest = Max(new_finest,levf);
            //
            // Construct initial cluster.
            //
	    ClusterList clist(pts);
	    clist.chop(grid_eff);
	    clist.intersect(p_n[levc]);
            //
            // Efficient properly nested Clusters have been constructed
            // now generate list of grids at level levf.
            //
	    BoxList new_bx;
	    clist.boxList(new_bx);

	    int nmerged = new_bx.minimize();

	    IntVect lratio = ref_ratio[levc]*bf_lev[levc];

	    IntVect largest_grid_size;
	    for (int n = 0; n < BL_SPACEDIM; n++)
	      largest_grid_size[n] = max_grid_size / lratio[n];
            //
            // Ensure new grid boxes are at most max_grid_size in
            // each index direction.
            //
	    maxSizeBox(new_bx,largest_grid_size);
            //
	    // Refine up to levf.
            //
            refineBoxList(new_bx,lratio);

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
    int more_levels = true;

    int num_levs = max_level+1;
    Array<BoxArray> new_grid_places(num_levs);
    while (more_levels)
    {
        //
        // Get new grid placement.
        //
	int     new_finest;
	grid_places(finest_level,strt_time,new_finest,new_grid_places);

	if (new_finest <= finest_level)
        {
	    more_levels = false;
	}
        else
        {
            //
            // Create a new level and link with others.
            //
	    int levf = new_finest;
	    finest_level = new_finest;

	    amr_level.set(levf,(*levelbld)(*this,levf,geom[levf],
					   new_grid_places[levf],
					   strt_time));
            //
            // Init data on new level.
            //
	    amr_level[levf].initData();

	    more_levels = (finest_level < max_level);
	}
    }
}

extern "C" void maxSizeBox (BoxList& bx_list, IntVect& block_size)
{
    for (BoxListIterator bli(bx_list); bli; ++bli)
    {
        const IntVect& ivlen = bli().length();
        const int* len       = ivlen.getVect();
        for (int i = 0; i < SpaceDim; i++)
        {
            if (len[i] > block_size[i])
            {
                //
                // Reduce by powers of 2.
                //
                int ratio = 1;
                int bs    = block_size[i];
                int nlen  = len[i];
                while ((bs%2==0) && (nlen%2==0))
                {
                    ratio *= 2;
                    bs /=2;
                    nlen /=2;
                }
                //
                // Determine number and size of (coarsened) cuts.
                //
                int numblk = nlen/bs + (nlen%bs ? 1 : 0);
                int size   = nlen/numblk;
                int extra  = nlen%numblk;
                //
                // Number of cuts = number of blocks - 1.
                //
                for (int k = 0; k < numblk-1; k++)
                {
                    //
                    // Compute size of this chunk, expand by power of 2.
                    //
                    int ksize = (k < extra) ? size+1 : size;
                    ksize *= ratio;
                    //
                    // Chop from high end.
                    //
                    IntVect iv = bli().bigEnd();
                    int pos    = iv[i] - ksize + 1;

                    bx_list.append(bx_list[bli].chop(i,pos));
                }
            }
        }
        //
        // b has been chopped down to size and pieces split off
        // have been added to the end of the list so that they
        // can be checked for splitting (in other directions) later.
        //
    }
}

extern "C" void refineBoxList (BoxList& bx_list, IntVect& lratio)
{
    BoxList bl;
    for (BoxListIterator bli(bx_list); bli; ++bli)
    {
      bl.add(refine(bli(),lratio));
    }
    bx_list.clear();
    bx_list = bl;
}
