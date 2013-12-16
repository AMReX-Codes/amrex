
#include <winstd.H>
#include <algorithm>
#include <cstdio>
#include <list>
#include <iostream>
#include <sstream>

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
std::list<std::string> Amr::derive_plot_vars;
bool                   Amr::first_plotfile;
Array<BoxArray>        Amr::initial_ba;
Array<BoxArray>        Amr::regrid_ba;

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
    n_proper               = 1;
    max_level              = -1;
    last_plotfile          = 0;
    last_checkpoint        = 0;
    record_run_info        = false;
    record_grid_info       = false;
    file_name_digits       = 5;
    record_run_info_terse  = false;
    bUserStopRequest       = false;

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
      Array<std::string> data_file_names(num_datalogs);
      pp.queryarr("data_log",data_file_names,0,num_datalogs);
      for (int i = 0; i < num_datalogs; i++) 
        setRecordDataInfo(i,data_file_names[i]);
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
        blocking_factor[i] = 2;
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
            BoxLib::Warning("Using default ref_ratio = 2 at all levels");
        }
    }
    
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
}

bool
Amr::isStatePlotVar (const std::string& name)
{
    for (std::list<std::string>::const_iterator li = state_plot_vars.begin(), End = state_plot_vars.end();
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
    const DescriptorList& desc_lst = AmrLevel::get_desc_lst();
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (desc_lst[typ].getType() == IndexType::TheCellType())
                state_plot_vars.push_back(desc_lst[typ].name(comp));
}

void
Amr::clearStatePlotVarList ()
{
    state_plot_vars.clear();
}

void
Amr::addStatePlotVar (const std::string& name)
{
    if (!isStatePlotVar(name))
        state_plot_vars.push_back(name);
}

void
Amr::deleteStatePlotVar (const std::string& name)
{
    if (isStatePlotVar(name))
        state_plot_vars.remove(name);
}

bool
Amr::isDerivePlotVar (const std::string& name)
{
    for (std::list<std::string>::const_iterator li = derive_plot_vars.begin(), End = derive_plot_vars.end();
         li != End;
         ++li)
    {
        if (*li == name)
            return true;
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
        datalog.set(i,new std::ofstream);
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
    if (!Plot_Files_Output()) return;

    BL_PROFILE("Amr::writePlotFile()");

    VisMF::SetNOutFiles(plot_nfiles);

    if (first_plotfile) 
    {
        first_plotfile = false;
        amr_level[0].setPlotVariables();
    }

    Real dPlotFileTime0 = ParallelDescriptor::second();

    const std::string pltfile = BoxLib::Concatenate(plot_file_root,level_steps[0],file_name_digits);

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "PLOTFILE: file = " << pltfile << '\n';

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "PLOTFILE: file = " << pltfile << '\n';

  BoxLib::StreamRetry sretry(pltfile, abort_on_stream_retry_failure,
                             stream_max_tries);

  while(sretry.TryFileOutput()) {

    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(pltfile, 0755))
            BoxLib::CreateDirectoryFailed(pltfile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier("Amr::writePlotFile::dir");

    std::string HeaderFileName = pltfile + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec = 0;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);
        old_prec = HeaderFile.precision(15);
    }

    for (int k = 0; k <= finest_level; k++)
        amr_level[k].writePlotFile(pltfile, HeaderFile);

    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.precision(old_prec);
        if (!HeaderFile.good())
            BoxLib::Error("Amr::writePlotFile() failed");
    }

    last_plotfile = level_steps[0];

    if (verbose > 0)
    {
        const int IOProc        = ParallelDescriptor::IOProcessorNumber();
        Real      dPlotFileTime = ParallelDescriptor::second() - dPlotFileTime0;

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(dPlotFileTime,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Write plotfile time = " << dPlotFileTime << "  seconds" << '\n';
#endif
    }
    ParallelDescriptor::Barrier("Amr::writePlotFile::end");

  }  // end while

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

    if (!Geometry::ProbDomain().ok())
        BoxLib::Error("checkInput: bad physical problem size");

    if (max_level > 0) 
       if (regrid_int[0] <= 0)
          BoxLib::Error("checkinput: regrid_int not defined and max_level > 0");

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
       std::cout << "Successfully read inputs file ... " << '\n';
}

void
Amr::init (Real strt_time,
           Real stop_time)
{
    if (!restart_chkfile.empty() && restart_chkfile != "init")
    {
        restart(restart_chkfile);
    }
    else
    {
        initialInit(strt_time,stop_time);
        checkPoint();
        if (plot_int > 0 || plot_per > 0)
            writePlotFile();
    }
#ifdef HAS_XGRAPH
    if (first_plotfile)
    {
        first_plotfile = false;
        amr_level[0].setPlotVariables();
    }
#endif
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
            FORT_PROBINIT(&init,
                          probin_file_name.dataPtr(),
                          &probin_file_length,
                          Geometry::ProbLo(),
                          Geometry::ProbHi());
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

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(piTotal,    IOProc);
        ParallelDescriptor::ReduceRealMax(piTotalAll, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "MFRead::: PROBINIT max time   = " << piTotal    << '\n';
            std::cout << "MFRead::: PROBINIT total time = " << piTotalAll << '\n';
        }
#endif
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
    BL_COMM_PROFILE_NAMETAG("Amr::initialInit TOP");
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
    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_init(stop_time);

    for (int lev = 0; lev <= finest_level; lev++)

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
    BL_PROFILE("Amr::restart()");

    Real dRestartTime0 = ParallelDescriptor::second();

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
                      BoxLib::Error("defBaseLevel: must have even number of cells");
               }
           }
       }

       if (regrid_on_restart && max_level > 0)
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
           level_count[0] = regrid_int[0];

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

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(dRestartTime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Restart time = " << dRestartTime << " seconds." << '\n';
#endif
    }
}

void
Amr::checkPoint ()
{
    if (!checkpoint_files_output) return;

    BL_PROFILE("Amr::checkPoint()");

    VisMF::SetNOutFiles(checkpoint_nfiles);
    //
    // In checkpoint files always write out FABs in NATIVE format.
    //
    FABio::Format thePrevFormat = FArrayBox::getFormat();

    FArrayBox::setFormat(FABio::FAB_NATIVE);

    Real dCheckPointTime0 = ParallelDescriptor::second();

    const std::string ckfile = BoxLib::Concatenate(check_file_root,level_steps[0],file_name_digits);

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "CHECKPOINT: file = " << ckfile << std::endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "CHECKPOINT: file = " << ckfile << '\n';


  BoxLib::StreamRetry sretry(ckfile, abort_on_stream_retry_failure,
                             stream_max_tries);

  while(sretry.TryFileOutput()) {

    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(ckfile, 0755))
            BoxLib::CreateDirectoryFailed(ckfile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier("Amr::checkPoint::dir");

    std::string HeaderFileName = ckfile + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec = 0, i;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);

        if (!HeaderFile.good())
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

    for (i = 0; i <= finest_level; i++)
        amr_level[i].checkPoint(ckfile, HeaderFile);

    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.precision(old_prec);

        if (!HeaderFile.good())
            BoxLib::Error("Amr::checkpoint() failed");
    }

    last_checkpoint = level_steps[0];

#ifdef USE_SLABSTAT
    //
    // Dump out any SlabStats MultiFabs.
    //
    AmrLevel::get_slabstat_lst().checkPoint(getAmrLevels(), level_steps[0]);
#endif
    //
    // Don't forget to reset FAB format.
    //
    FArrayBox::setFormat(thePrevFormat);

    if (verbose > 0)
    {
        Real dCheckPointTime = ParallelDescriptor::second() - dCheckPointTime0;

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(dCheckPointTime,
	                            ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "checkPoint() time = " << dCheckPointTime << " secs." << '\n';
#endif
    }
    ParallelDescriptor::Barrier("Amr::checkPoint::end");

  }  // end while

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
    //
    // Allow regridding of level 0 calculation on restart.
    //
    if (max_level == 0 && regrid_on_restart)
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
    else
    {
        int lev_top = std::min(finest_level, max_level-1);

        for (int i = level; i <= lev_top; i++)
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
                        dt_level[k]    = dt_level[k-1]/n_cycle[k];
                    }
                }
            }
            if (old_finest > finest_level)
                lev_top = std::min(finest_level, max_level-1);
        }
    }
    //
    // Check to see if should write plotfile.
    // This routine is here so it is done after the restart regrid.
    //
    if (plotfile_on_restart && !(restart_chkfile.empty()) )
    {
	plotfile_on_restart = 0;
	writePlotFile();
    }
    //
    // Advance grids at this level.
    //
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "ADVANCE grids at level "
                  << level
                  << " with dt = "
                  << dt_level[level]
                  << std::endl;
    }
    Real dt_new = amr_level[level].advance(time,dt_level[level],iteration,niter);

    dt_min[level] = iteration == 1 ? dt_new : std::min(dt_min[level],dt_new);

    level_steps[level]++;
    level_count[level]++;

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Advanced "
                  << amr_level[level].countCells()
                  << " cells at level "
                  << level
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
    BL_PROFILE("Amr::coarseTimeStep()");

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
    timeStep(0,cumtime,1,1,stop_time);

    cumtime += dt_level[0];

    amr_level[0].postCoarseTimeStep(cumtime);

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_stop = ParallelDescriptor::second() - run_strt;

        ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "\nCoarse TimeStep time: " << run_stop << '\n' ;

        long min_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;
        long max_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;

        ParallelDescriptor::ReduceLongMin(min_fab_bytes,IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_bytes,IOProc);
        //
        // Reset to zero to calculate high-water-mark for next timestep.
        //
        BoxLib::total_bytes_allocated_in_fabs_hwm = 0;

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "\nFAB byte spread across MPI nodes for timestep: ["
                      << min_fab_bytes
                      << " ... "
                      << max_fab_bytes
                      << "]\n";

            std::cout << "\nHigh water mark for bytes in BoxArray hash tables: "
                      << BoxLib::total_bytes_in_hashtables_hwm
                      << '\n';
        }

        BoxLib::total_bytes_in_hashtables_hwm = 0;
    }

    BL_PROFILE_ADD_STEP(level_steps[0]);
 #ifdef BL_COMM_PROFILING
    std::stringstream stepName;
    stepName << "STEP " << level_steps[0];
    BL_COMM_PROFILE_NAMETAG(stepName.str());
    BL_COMM_PROFILE_FLUSH();
 #endif

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
      const int num_per_old = cumtime / check_per;
      const int num_per_new = (cumtime+dt_level[0]) / check_per;

      if (num_per_old != num_per_new)
	{
	 check_test = 1;
	}
    }

    int to_stop       = 0;    
    int to_checkpoint = 0;
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
    }
    int packed_data[2];
    packed_data[0] = to_stop;
    packed_data[1] = to_checkpoint;
    ParallelDescriptor::Bcast(packed_data, 2, ParallelDescriptor::IOProcessorNumber());
    to_stop = packed_data[0];
    to_checkpoint = packed_data[1];

    if(to_stop == 1 && to_checkpoint == 0) {  // prevent main from writing files
      last_checkpoint = level_steps[0];
      last_plotfile   = level_steps[0];
    }
    if ((check_int > 0 && level_steps[0] % check_int == 0) || check_test == 1
	|| to_checkpoint)
    {
        checkPoint();
    }


    if (writePlotNow() || to_checkpoint)
    {
        writePlotFile();
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

void
Amr::defBaseLevel (Real              strt_time, 
                   const BoxArray*   lev0_grids,
                   const Array<int>* pmap)
{
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
    BL_PROFILE("Amr::regrid()");

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "REGRID: at level lbase = " << lbase << std::endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "REGRID: at level lbase = " << lbase << '\n';
    //
    // Compute positions of new grids.
    //
    int             new_finest;
    Array<BoxArray> new_grid_places(max_level+1);

    if (lbase <= std::min(finest_level,max_level-1))
      grid_places(lbase,time,new_finest, new_grid_places);

    bool regrid_level_zero =
        (lbase == 0 && new_grid_places[0] != amr_level[0].boxArray()) && (!initial);

    const int start = regrid_level_zero ? 0 : lbase+1;
    //
    // If use_efficient_regrid flag is set, then test to see whether we in fact 
    // have changed the grids at any of the levels through the regridding process.
    // If not, then don't do anything more here.
    //
    if (use_efficient_regrid == 1 && !regrid_level_zero && (finest_level == new_finest) )
    {
        bool grids_unchanged = true;
        for (int lev = start; lev <= finest_level && grids_unchanged; lev++)
        {
            if (new_grid_places[lev] != amr_level[lev].boxArray()) grids_unchanged = false;
        }
        if (grids_unchanged) 
        {
            if (verbose > 0 && ParallelDescriptor::IOProcessor())
                std::cout << "Regridding at level lbase = " << lbase << " but grids unchanged " << std::endl;
            return;
        }
    }

    //
    // Reclaim old-time grid space for all remain levels > lbase.
    //
    for (int lev = start; lev <= finest_level; lev++)
    {
        amr_level[lev].removeOldData();
    }
    //
    // Reclaim all remaining storage for levels > new_finest.
    //
    for (int lev = new_finest+1; lev <= finest_level; lev++)
        amr_level.clear(lev);

    finest_level = new_finest;

    if (lbase == 0)
    {
        MultiFab::FlushSICache();
        Geometry::FlushPIRMCache();
        FabArrayBase::CPC::FlushCache();
        DistributionMapping::FlushCache();
    }
    //
    // Define the new grids from level start up to new_finest.
    //
    for (int lev = start; lev <= new_finest; lev++) 
    {
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
    // Build any additional data structures at levels start and higher after grid generation.
    //
    for (int lev = 0; lev <= new_finest; lev++)
        amr_level[lev].post_regrid(lbase,new_finest);

#ifdef USE_STATIONDATA
    station.findGrid(amr_level,geom);
#endif
    //
    // Report creation of new grids.
    //
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
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

    if ( time == 0. && !initial_grids_file.empty() )
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

    else if ( !regrid_grids_file.empty() )
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
            bf_lev[i][n] = std::max(1,blocking_factor[i]/ref_ratio[i][n]);
    }
    for (i = lbase; i < max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            rr_lev[i][n] = (ref_ratio[i][n]*bf_lev[i][n])/bf_lev[i+1][n];
    }
    for (i = lbase; i <= max_crse; i++)
        pc_domain[i] = BoxLib::coarsen(geom[i].Domain(),bf_lev[i]);
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

            const IntVect iv = IntVect(D_DECL(nerr/ref_ratio[levf][0],
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
        long     len = 0;
        IntVect* pts = tags.collate(len);

        tags.clear();

        if (len > 0)
        {
            //
            // Created new level, now generate efficient grids.
            //
            new_finest = std::max(new_finest,levf);
            //
            // Construct initial cluster.
            //
            ClusterList clist(pts,len);
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
            new_grids[levf].define(new_bx);
        }
        //
        // Don't forget to get rid of space used for collate()ing.
        //
        delete [] pts;
    }

    // If Nprocs > Ngrids and refine_grid_layout == 1 then break up the grids
    //    into smaller chunks for better load balancing
    impose_refine_grid_layout(lbase,new_finest,new_grids);

    if (verbose > 0)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "grid_places() time: " << stoptime << '\n';
#endif
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
    int got_check_int = pp.query("check_int",check_int);

    check_per = -1.0;
    int got_check_per = pp.query("check_per",check_per);

    if (got_check_int == 1 && got_check_per == 1)
    {
        BoxLib::Error("Must only specify amr.check_int OR amr.check_per");
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

    stream_max_tries = 4;
    pp.query("stream_max_tries",stream_max_tries);
    stream_max_tries = std::max(stream_max_tries, 1);

    abort_on_stream_retry_failure = false;
    pp.query("abort_on_stream_retry_failure",abort_on_stream_retry_failure);

}


bool
Amr::okToRegrid(int level)
{
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
