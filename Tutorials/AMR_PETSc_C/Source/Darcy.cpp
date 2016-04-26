#include <winstd.H>

#include <algorithm>
#include <cstdio>
#include <vector>
#include <iostream>
#include <string>
using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <Utility.H>
#include <CONSTANTS.H>
#include <Darcy.H>
#include <Darcy_F.H>

#include <VisMF.H>
#include <TagBox.H>
#include <ParmParse.H>

static Box the_same_box (const Box& b) { return b; }

typedef StateDescriptor::BndryFunc BndryFunc;

static Real fixed_dt     = -1.0;
static Real initial_dt   = -1.0;
static Real dt_cutoff    = 0.0;


static bool initialized = false;
void 
Darcy::CleanupStatics ()
{
  initialized = false;
  PetscFinalize();
}

bool         Darcy::dump_old;

int          Darcy::verbose;
Real         Darcy::init_shrink;
Real         Darcy::change_max;
ErrorList    Darcy::err_list;
int          Darcy::radius_grow;
BCRec        Darcy::phys_bc;
int          Darcy::NUM_STATE;
int          Darcy::NUM_PHASES;
int          Darcy::NUM_COMPS_PER_PHASE;

int          Darcy::RhoSat;
int          Darcy::Pressure;

Real         Darcy::dt_suggest_for_next;
Real         Darcy::dt_max;

Layout*      Darcy::layout;
DarcySNES*   Darcy::snes;
MLBoundary*  Darcy::mlb;
int          Darcy::verbose_SNES;

void
Darcy::InitializeStaticVariables ()
{
  dump_old      = false;
  
  verbose             = 0;
  init_shrink         = 1.0;
  change_max          = 1.1;
  radius_grow         = 1;
  NUM_STATE           = -1;
  NUM_PHASES          = -1;
  NUM_COMPS_PER_PHASE = -1;
  
  RhoSat   = -1;
  Pressure = -1;

  dt_suggest_for_next = -1;
  dt_max = -1;

  layout   = 0;
  mlb      = 0;
  snes     = 0;
  verbose_SNES = 0;

  int argc=0;
  char** argv;

  ParmParse pp;
  std::string petsc_help="Petsc Options File";
  std::string petsc_options_file=""; pp.query("petsc_options_file",petsc_options_file);
  if (!petsc_options_file.empty()) {
    cout << "Setting PETSc options file: \"" << petsc_options_file << "\"" << endl;
  }
  PetscInitialize(&argc,&argv,petsc_options_file.c_str(),petsc_help.c_str());
}

//
// Components are:
//  0:Interior, 1:p-specified, 2:flux-specified, 3:Symmetry, 4:SlipWall-not relevant, 5:NoSlipWall-not relevant
//
// Definitions are:
// REFLECT_ODD:-1, INT_DIR:0, REFLECT_EVEN:1, FOEXTRAP:2, EXT_DIR:3, HOEXTRAP:4
//    Will use EXT_DIR to indicate pressure values specified at the boundary
//    Will use FOEXTRAP to indicate flux given at boundary
//
static int scalar_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
void
set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,scalar_bc[lo_bc[i]]);
        bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

void
Darcy::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    read_params();
    //
    // Set number of state variables and pointers to components
    //
    NUM_PHASES = 1;
    NUM_COMPS_PER_PHASE = 1;

    int cnt = 0;
    RhoSat = cnt;
    cnt += NUM_PHASES * NUM_COMPS_PER_PHASE;
    Pressure = cnt;
    cnt += NUM_PHASES;
    NUM_STATE = cnt;

    const Real run_strt = ParallelDescriptor::second() ; 

    //BL_FORT_PROC_CALL(SET_METHOD_PARAMS, set_method_params)
    //(RhoSat, Pressure, NUM_STATE);

    Real run_stop = ParallelDescriptor::second() - run_strt;
 
    ParallelDescriptor::ReduceRealMax(run_stop,ParallelDescriptor::IOProcessorNumber());
 
    if (ParallelDescriptor::IOProcessor())
        std::cout << "\nTime in set_method_params: " << run_stop << '\n' ;

    //int coord_type = Geometry::Coord();
    //int dm = BL_SPACEDIM;
    //BL_FORT_PROC_CALL(SET_PROBLEM_PARAMS, set_problem_params)
    //     (dm,phys_bc.lo(),phys_bc.hi(),Outflow,Symmetry,coord_type);

    Interpolater* interp = &cell_cons_interp;

    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,1,NUM_STATE,interp);

    Array<BCRec>       bcs(NUM_STATE);
    Array<std::string> name(NUM_STATE);
    BCRec bc;
    
    set_scalar_bc(bc,phys_bc); 
    bcs[RhoSat] = bc; name[RhoSat] = "RhoSat";
    bcs[Pressure] = bc; name[Pressure] = "Pressure";
    desc_lst.setComponent(State_Type,0,name,bcs,
                          BndryFunc(BL_FORT_PROC_CALL(STATEFILL,statefill)),
                          interp);
                          
    //
    // DEFINE DERIVED QUANTITIES
    //
    derive_lst.add("StateErr",IndexType::TheCellType(),NUM_STATE,
                   BL_FORT_PROC_CALL(DERSTATE,derstate),the_same_box);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,RhoSat,2);

    // These must correspond to those in PorousMedia::derive
    IndexType deriveType(IndexType::TheCellType());
    int nCompDerive = 1;
    std::string amr_prefix = "amr";
    ParmParse pp(amr_prefix);
    int num_derive_plot_vars = pp.countval("derive_plot_vars");
    if (num_derive_plot_vars) {
      Array<std::string> derive_plot_vars(num_derive_plot_vars);
      pp.getarr("derive_plot_vars",derive_plot_vars);

      for (int i=0; i<num_derive_plot_vars; ++i) {
        derive_lst.add(derive_plot_vars[i], deriveType, nCompDerive);
      }
    }

    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    err_list.add("StateErr",1,ErrorRec::Special,
                 BL_FORT_PROC_CALL(STATE_ERROR,state_error));
}

void
Darcy::variableCleanUp () 
{
  desc_lst.clear();
}

void
Darcy::read_params ()
{
  static bool done = false;

  if (done) return;

  done = true;

  ParmParse pp("darcy");   

  pp.query("v",verbose);
  pp.query("v_snes",verbose_SNES);

  pp.get("init_shrink",init_shrink);
  pp.query("change_max",change_max);
  pp.query("fixed_dt",fixed_dt);
  pp.query("initial_dt",initial_dt);
  pp.get("dt_cutoff",dt_cutoff);
  pp.get("dt_max",dt_max);

  pp.query("dump_old",dump_old);

  // Get boundary conditions
  Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
  pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
  pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    phys_bc.setLo(i,lo_bc[i]);
    phys_bc.setHi(i,hi_bc[i]);
  }

  //
  // Check phys_bc against possible periodic geometry
  // if periodic, must have internal BC marked.
  //
  if (Geometry::isAnyPeriodic())
  {
    //
    // Do idiot check.  Periodic means interior in those directions.
    //
    for (int dir = 0; dir<BL_SPACEDIM; dir++)
    {
      if (Geometry::isPeriodic(dir))
      {
        if (lo_bc[dir] != Interior)
        {
          std::cerr << "Darcy::read_params:periodic in direction "
                    << dir
                    << " but low BC is not Interior\n";
          BoxLib::Error();
        }
        if (hi_bc[dir] != Interior)
        {
          std::cerr << "Darcy::read_params:periodic in direction "
                    << dir
                    << " but high BC is not Interior\n";
          BoxLib::Error();
        }
      }
    }
  }
  else
  {
    //
    // Do idiot check.  If not periodic, should be no interior.
    //
    for (int dir=0; dir<BL_SPACEDIM; dir++)
    {
      if (lo_bc[dir] == Interior)
      {
        std::cerr << "Darcy::read_params:interior bc in direction "
                  << dir
                  << " but not periodic\n";
        BoxLib::Error();
      }
      if (hi_bc[dir] == Interior)
      {
        std::cerr << "Darcy::read_params:interior bc in direction "
                  << dir
                  << " but not periodic\n";
        BoxLib::Error();
      }
    }
  }

  if ( Geometry::IsRZ() && (lo_bc[0] != Symmetry) ) {
    std::cerr << "ERROR:Darcy::read_params: must set r=0 boundary condition to Symmetry for r-z\n";
    BoxLib::Error();
  }

#if (BL_SPACEDIM == 1)
  if ( Geometry::IsSPHERICAL() )
  {
    if ( (lo_bc[0] != Symmetry) && (Geometry::ProbLo(0) == 0.0) ) 
    {
      std::cerr << "ERROR:Darcy::read_params: must set r=0 boundary condition to Symmetry for spherical\n";
      BoxLib::Error();
    }
  }
#elif (BL_SPACEDIM == 2)
  if ( Geometry::IsSPHERICAL() )
  {
    BoxLib::Abort("We don't support spherical coordinate systems in 2D");
  }
#elif (BL_SPACEDIM == 3)
  if ( Geometry::IsRZ() )
  {
    BoxLib::Abort("We don't support cylindrical coordinate systems in 3D"); 
  }
  else if ( Geometry::IsSPHERICAL() )
  {
    BoxLib::Abort("We don't support spherical coordinate systems in 3D");
  }
#endif
}

void
Darcy::PostStaticInitialize () 
{
  if (!initialized) {

    BoxLib::ExecOnFinalize(Darcy::CleanupStatics);

    initialized = true;
  }

  post_restart_flag = false;
}

Darcy::Darcy ()
{
  PostStaticInitialize();
}

Darcy::Darcy (Amr&            papa,
              int             lev,
              const Geometry& level_geom,
              const BoxArray& bl,
              Real            time)
  :
  AmrLevel(papa,lev,level_geom,bl,time) 
{
  PostStaticInitialize();

  buildMetrics();
}

Darcy::~Darcy () 
{
  if (level==0) {
    BL_ASSERT(snes!=0);
    BL_ASSERT(mlb!=0);
    BL_ASSERT(layout!=0);
    delete snes; snes=0;
    delete mlb; mlb=0;
    delete layout; layout=0;
  }
}

void
Darcy::restart (Amr&     papa,
                istream& is,
                bool     bReadSpecial)
{
  post_restart_flag = true;

  AmrLevel::restart(papa,is,bReadSpecial);

  if (level == parent->finestLevel()) {
    build_layout();
    build_snes();
  }

  buildMetrics();
}

void
Darcy::checkPoint(const std::string& dir,
                  std::ostream&  os,
                  VisMF::How     how,
                  bool dump_old_default)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
}

std::string
Darcy::thePlotFileType () const
{
  //
  // Increment this whenever the writePlotFile() format changes.
  //
  static const std::string the_plot_file_type("HyperCLaw-V1.1");

  return the_plot_file_type;
}

void
Darcy::setPlotVariables ()
{
  AmrLevel::setPlotVariables();
}

void
Darcy::writePlotFile (const std::string& dir,
                      ostream&       os,
                      VisMF::How     how)
{
  //
  // The list of indices of State to write to plotfile.
  // first component of pair is state_type,
  // second component of pair is component # within the state_type
  //
  std::vector<std::pair<int,int> > plot_var_map;
  for (int typ = 0; typ < desc_lst.size(); typ++)
    for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
      if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
          desc_lst[typ].getType() == IndexType::TheCellType())
        plot_var_map.push_back(std::pair<int,int>(typ,comp));

  int num_derive = 0;
  std::list<std::string> derive_names;
  const std::list<DeriveRec>& dlist = derive_lst.dlist();

  for (std::list<DeriveRec>::const_iterator it = dlist.begin();
       it != dlist.end();
       ++it)
  {
    if (parent->isDerivePlotVar(it->name()))
    {
      derive_names.push_back(it->name());
      num_derive++;
    }
  }

  int n_data_items = plot_var_map.size() + num_derive;

  Real cur_time = state[State_Type].curTime();

  if (level == 0 && ParallelDescriptor::IOProcessor())
  {
    //
    // The first thing we write out is the plotfile type.
    //
    os << thePlotFileType() << '\n';

    if (n_data_items == 0)
      BoxLib::Error("Must specify at least one valid data item to plot");

    os << n_data_items << '\n';

    //
    // Names of variables -- first state, then derived
    //
    for (int i=0; i < plot_var_map.size(); i++)
    {
      int typ = plot_var_map[i].first;
      int comp = plot_var_map[i].second;
      os << desc_lst[typ].name(comp) << '\n';
    }

    for ( std::list<std::string>::iterator it = derive_names.begin();
          it != derive_names.end(); ++it)
    {
      const DeriveRec* rec = derive_lst.get(*it);
      os << rec->variableName(0) << '\n';
    }

    os << BL_SPACEDIM << '\n';
    os << parent->cumTime() << '\n';
    int f_lev = parent->finestLevel();
    os << f_lev << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++)
      os << Geometry::ProbLo(i) << ' ';
    os << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++)
      os << Geometry::ProbHi(i) << ' ';
    os << '\n';
    for (int i = 0; i < f_lev; i++)
      os << parent->refRatio(i)[0] << ' ';
    os << '\n';
    for (int i = 0; i <= f_lev; i++)
      os << parent->Geom(i).Domain() << ' ';
    os << '\n';
    for (int i = 0; i <= f_lev; i++)
      os << parent->levelSteps(i) << ' ';
    os << '\n';
    for (int i = 0; i <= f_lev; i++)
    {
      for (int k = 0; k < BL_SPACEDIM; k++)
        os << parent->Geom(i).CellSize()[k] << ' ';
      os << '\n';
    }
    os << (int) Geometry::Coord() << '\n';
    os << "0\n"; // Write bndry data.

  }
  // Build the directory to hold the MultiFab at this level.
  // The name is relative to the directory containing the Header file.
  //
  static const std::string BaseName = "/Cell";
  char buf[64];
  sprintf(buf, "Level_%d", level);
  std::string Level = buf;
  //
  // Now for the full pathname of that directory.
  //
  std::string FullPath = dir;
  if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
    FullPath += '/';
  FullPath += Level;
  //
  // Only the I/O processor makes the directory if it doesn't already exist.
  //
  if (ParallelDescriptor::IOProcessor())
    if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
      BoxLib::CreateDirectoryFailed(FullPath);
  //
  // Force other processors to wait till directory is built.
  //
  ParallelDescriptor::Barrier();

  if (ParallelDescriptor::IOProcessor())
  {
    os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
    os << parent->levelSteps(level) << '\n';

    for (int i = 0; i < grids.size(); ++i)
    {
      RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
      for (int n = 0; n < BL_SPACEDIM; n++)
        os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
    }
    //
    // The full relative pathname of the MultiFabs at this level.
    // The name is relative to the Header file containing this name.
    // It's the name that gets written into the Header.
    //
    if (n_data_items > 0)
    {
      std::string PathNameInHeader = Level;
      PathNameInHeader += BaseName;
      os << PathNameInHeader << '\n';
    }
  }
  //
  // We combine all of the multifabs -- state, derived, etc -- into one
  // multifab -- plotMF.
  // NOTE: we are assuming that each state variable has one component,
  // but a derived variable is allowed to have multiple components.
  int       cnt   = 0;
  const int nGrow = 0;
  MultiFab  plotMF(grids,n_data_items,nGrow);
  MultiFab* this_dat = 0;
  //
  // Cull data from state variables -- use no ghost cells.
  //
  for (int i = 0; i < plot_var_map.size(); i++)
  {
    int typ  = plot_var_map[i].first;
    int comp = plot_var_map[i].second;
    this_dat = &state[typ].newData();
    MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
    cnt++;
  }
  //
  // Cull data from derived variables.
  // 
  if (derive_names.size() > 0)
  {
    for (std::list<std::string>::iterator it = derive_names.begin();
         it != derive_names.end(); ++it) 
    {
      MultiFab* derive_dat = derive(*it,cur_time,nGrow);
      MultiFab::Copy(plotMF,*derive_dat,0,cnt,1,nGrow);
      delete derive_dat;
      cnt++;
    }
  }

  //
  // Use the Full pathname when naming the MultiFab.
  //
  std::string TheFullPath = FullPath;
  TheFullPath += BaseName;
  VisMF::Write(plotMF,TheFullPath,how,true);
}

void
Darcy::buildMetrics ()
{
  const int ngrd = grids.size();

  radius.resize(ngrd);

  const Real* dx = geom.CellSize();

  for (int i = 0; i < ngrd; i++)
  {
    const Box& b = grids[i];
    int ilo      = b.smallEnd(0)-radius_grow;
    int ihi      = b.bigEnd(0)+radius_grow;
    int len      = ihi - ilo + 1;

    radius[i].resize(len);

    Real* rad = radius[i].dataPtr();

    if (Geometry::IsCartesian())
    {
      for (int j = 0; j < len; j++)
      {
        rad[j] = 1.0;
      }
    }
    else
    {
      RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());

      const Real xlo = gridloc.lo(0) + (0.5 - radius_grow)*dx[0];

      for (int j = 0; j < len; j++)
      {
        rad[j] = xlo + j*dx[0];
      }
    }
  }
  //
  // Build volume, face area and dLogArea arrays.
  // volume is not PArrayManaged, must manually delete.
  //
  volume.clear();
  //
  // area is not PArrayManaged, must manually delete.
  //
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    area[dir].clear();
  }
  dLogArea[0].clear();

  int nGrowGeom = 0;
  geom.GetVolume(volume,grids,nGrowGeom);

  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    geom.GetFaceArea(area[dir],grids,dir,nGrowGeom);
  }
#if (BL_SPACEDIM <= 2)
  geom.GetDLogA(dLogArea[0],grids,0,nGrowGeom);
#endif
}

void
Darcy::setTimeLevel (Real time,
                     Real dt_old,
                     Real dt_new)
{
  AmrLevel::setTimeLevel(time,dt_old,dt_new);
}

void
Darcy::initData ()
{
  if (layout==0 || layout->NumLevels()!=parent->finestLevel()+1) {
    delete snes; snes = 0;
    delete mlb; mlb = 0;
    delete layout; layout = 0;

    build_layout();
    build_snes();
  }

  MultiFab& S_new = get_new_data(State_Type);
  S_new.setVal(0,Pressure,1);

  for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    const Array<Real> dsGrav = snes->gravity;
    const Real rho = snes->density;
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Box& vbox = mfi.validbox();
    FArrayBox& p = S_new[mfi];
    for (IntVect iv=vbox.smallEnd(); iv<=vbox.bigEnd(); vbox.next(iv)) {
      Real mag = 0;
      for (int d=0; d<BL_SPACEDIM; ++d) {
        Real x = plo[d] + dx[d]*(iv[d]+0.5);
        Real vd = x * dsGrav[d] * rho;
        mag += vd * vd;
      }
      p(iv,Pressure) = - std::sqrt(mag);
    }
  }

  const Layout::MultiIntFab& nodeIds = layout->NodeIds()[level];

  for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    const FArrayBox& p = S_new[mfi];
    FArrayBox& rs = S_new[mfi];
    const Box& box = mfi.validbox();
    const Layout::IntFab& ids = nodeIds[mfi];
    snes->ReducedSaturationGivenPressure(p,Pressure,rs,RhoSat,box,ids);
    snes->RhoSatGivenReducedSaturation(rs,RhoSat,rs,RhoSat,box,ids);
  }
}

void
Darcy::init (AmrLevel &old)
{
  Darcy* oldlev = (Darcy*) &old;
  //
  // Create new grid data by fillpatching from old.
  //
  Real dt_new    = parent->dtLevel(level);
  Real cur_time  = oldlev->state[State_Type].curTime();
  Real prev_time = oldlev->state[State_Type].prevTime();
  Real dt_old    = cur_time - prev_time;
  setTimeLevel(cur_time,dt_old,dt_new);

  MultiFab& S_new = get_new_data(State_Type);
    
  for (FillPatchIterator fpi(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
       fpi.isValid();
       ++fpi)
  {
    S_new[fpi].copy(fpi());
  }
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
Darcy::init ()
{
  Real dt        = parent->dtLevel(level);
  Real cur_time  = getLevel(level-1).state[State_Type].curTime();
  Real prev_time = getLevel(level-1).state[State_Type].prevTime();

  Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

  setTimeLevel(cur_time,dt_old,dt);
  MultiFab& S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);
}

void
Darcy::setPhysBoundaryValues (FArrayBox& dest,
                              int        state_indx,
                              Real       time,
                              int        dest_comp,
                              int        src_comp,
                              int        num_comp)
{
  // Currently, do nothing here
  // At the moment, the boundary conditions are handled directly inside DarcySNES
}


Real
Darcy::initialTimeStep ()
{
  Real dummy_dt = 0.0;
  return (initial_dt > 0.0) ? initial_dt : init_shrink*estTimeStep(dummy_dt);
}

Real
Darcy::estTimeStep (Real dt_old)
{
  if (fixed_dt > 0.0)
    return fixed_dt;

  Real estdt = dt_old;

  if (post_restart_flag) {

    estdt = parent->dtMin(level);

  } else if (dt_suggest_for_next > 0) {

    estdt = dt_suggest_for_next;

  }

  if (verbose > 1 && ParallelDescriptor::IOProcessor())
    cout << "Darcy::estTimeStep at level " << level << ":  estdt = " << estdt << '\n';

  return estdt;
}
void
Darcy::computeNewDt (int                   finest_level,
                     int                   sub_cycle,
                     Array<int>&           n_cycle,
                     const Array<IntVect>& ref_ratio,
                     Array<Real>&          dt_min,
                     Array<Real>&          dt_level,
                     Real                  stop_time,
                     int                   post_regrid_flag)
{
  //
  // We are at the end of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  //
  if (level > 0)
    return;

  int i;
  n_cycle[0] = 1;
  for (i = 1; i <= finest_level; i++) {
    n_cycle[i] = sub_cycle ? parent->MaxRefRatio(i-1) : 1;
  }

  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (i = 0; i <= finest_level; i++)
  {
    dt_min[i] = getLevel(i).estTimeStep(dt_level[i]);
  }

  if (fixed_dt <= 0.0)
  {
    if (post_regrid_flag == 1) 
    {
      //
      // Limit dt's by pre-regrid dt
      //
      for (i = 0; i <= finest_level; i++)
      {
        dt_min[i] = std::min(dt_min[i],dt_level[i]);
      }
    } 
    else 
    {
      //
      // Limit dt's by change_max * old dt
      //
      for (i = 0; i <= finest_level; i++)
      {
        if (verbose && ParallelDescriptor::IOProcessor())
          if (dt_min[i] > change_max*dt_level[i])
          {
            cout << "Darcy::computeNewDt : limiting dt at level " << i << std::endl;
            cout << " ... new dt computed: " << dt_min[i] << std::endl;
            cout << " ... but limiting to: " << change_max*dt_level[i] <<
              " = " << change_max << " * " << dt_level[i] << std::endl;
          }
        dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
      }
    } 
  }

  //
  // Find the minimum over all levels
  //
  for (i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_min[i]);
  }

  //
  // Limit dt's by the value of stop_time, t_max
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }
  if (dt_max > 0.0) {
    dt_0 = std::min(dt_max,dt_0);
  }

  n_factor = 1;
  for (i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

void
Darcy::computeInitialDt (int                   finest_level,
                         int                   sub_cycle,
                         Array<int>&           n_cycle,
                         const Array<IntVect>& ref_ratio,
                         Array<Real>&          dt_level,
                         Real                  stop_time)
{
  //
  // Grids have been constructed, compute dt for all levels.
  //
  if (level > 0)
    return;

  n_cycle[0] = 1;
  for (int i = 1; i <= finest_level; i++)
  {
    n_cycle[i] = sub_cycle ? parent->MaxRefRatio(i-1) : 1;
  }

  Real dt_0 = (initial_dt>0 ? initial_dt : 1.0e+100);
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor   *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_level[i]);
  }

  //
  // Limit dt's by the value of stop_time, dt_max
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }
  if (dt_max > 0.0) {
    dt_0 = std::min(dt_max,dt_0);
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

void
Darcy::post_timestep (int iteration)
{
  post_restart_flag = false;
}

void
Darcy::post_restart ()
{
}

void
Darcy::postCoarseTimeStep (Real cumtime)
{
}

void
Darcy::build_layout ()
{
  BL_ASSERT(layout == 0);
  int nLevs = parent->finestLevel() + 1;
  Array<BoxArray> aba(nLevs);
  PArray<Geometry> ag(nLevs,PArrayManage);
  const Array<IntVect>& ar = parent->refRatio();
  for (int lev=0; lev<nLevs; ++lev) {
    BL_ASSERT(parent->getAmrLevels().defined(lev));
    aba[lev] = parent->boxArray(lev);
    ag.set(lev, new Geometry(parent->Geom(lev)));
  }
  layout = new Layout(aba,ag,ar); // Nulls the Geom ptr
}

void
Darcy::build_snes ()
{
  const BCRec& pressure_bc = get_desc_lst()[State_Type].getBC(Pressure);
  BL_ASSERT(mlb == 0);
  mlb = new MLBoundary(*layout,pressure_bc);

  DarcySNES::SetVerbose(verbose_SNES);
  BL_ASSERT(snes == 0);
  snes = new DarcySNES(*layout,*mlb);
}

void
Darcy::post_regrid (int lbase,
                    int new_finest)
{
}

void
Darcy::post_init (Real stop_time)
{
  if (level > 0) {
    return;
  }
}

int
Darcy::okToContinue ()
{
  if (level > 0)
    return 1;
  return  parent->dtLevel(0) > dt_cutoff;
}

void
Darcy::allocOldData ()
{
  for (int k = 0; k < NUM_STATE_TYPE; k++)
    state[k].allocOldData();
}

void
Darcy::removeOldData()
{
  AmrLevel::removeOldData();
}

void
Darcy::errorEst (TagBoxArray& tags,
                 int          clearval,
                 int          tagval,
                 Real         time,
                 int          n_error_buf,
                 int          ngrow)
{
  return;

  const int*  domain_lo = geom.Domain().loVect();
  const int*  domain_hi = geom.Domain().hiVect();
  const Real* dx        = geom.CellSize();
  const Real* prob_lo   = geom.ProbLo();
  Array<int>  itags;

  for (int j = 0; j < err_list.size(); j++)
  {
    MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());

    BL_ASSERT(!(mf == 0));

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
      int         idx     = mfi.index();
      RealBox     gridloc = RealBox(grids[idx],geom.CellSize(),geom.ProbLo());
      itags               = tags[idx].tags();
      int*        tptr    = itags.dataPtr();
      const int*  tlo     = tags[idx].box().loVect();
      const int*  thi     = tags[idx].box().hiVect();
      const int*  lo      = mfi.validbox().loVect();
      const int*  hi      = mfi.validbox().hiVect();
      const Real* xlo     = gridloc.lo();
      Real*       dat     = (*mf)[mfi].dataPtr();
      const int*  dlo     = (*mf)[mfi].box().loVect();
      const int*  dhi     = (*mf)[mfi].box().hiVect();
      const int   ncomp   = (*mf)[mfi].nComp();

      err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                            &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                            lo,hi, &ncomp, domain_lo, domain_hi,
                            dx, xlo, prob_lo, &time, &level);
      //
      // Don't forget to set the tags in the TagBox.
      //
      tags[idx].tags(itags);
    }

    delete mf;
  }
}

MultiFab*
Darcy::derive (const std::string& name,
               Real               time,
               int                ngrow)
{
    BL_ASSERT(ngrow >= 0);

    MultiFab* mf = 0;
    const DeriveRec* rec = derive_lst.get(name);
    if (rec)
    {
        BoxArray dstBA(grids);
        mf = new MultiFab(dstBA, rec->numDerive(), ngrow);
        int dcomp = 0;
        derive(name,time,*mf,dcomp);
    }
    else
    {
        //
        // If we got here, cannot derive given name.
        //
        std::string msg("Darcy::derive(): unknown variable: ");
        msg += name;
        BoxLib::Error(msg.c_str());
    }
    return mf;
}

void
Darcy::derive (const std::string& name,
               Real               time,
               MultiFab&          mf,
               int                dcomp)
{
  bool not_found_yet = false;
  
  const DeriveRec* rec = derive_lst.get(name);

  if (name=="Seff") {
        
    BL_ASSERT(dcomp < mf.nComp());
    BL_ASSERT(rec);
    
    BoxArray dstBA(mf.boxArray());
    BL_ASSERT(rec->deriveType() == dstBA.ixType());

    const Layout::MultiIntFab& nodeIds = layout->NodeIds()[level];

    MultiFab& S_new = get_new_data(State_Type);
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
      const FArrayBox& p = S_new[mfi];
      FArrayBox& rs = mf[mfi];
      const Box& box = mfi.validbox();
      const Layout::IntFab& ids = nodeIds[mfi];
      snes->ReducedSaturationGivenPressure(p,Pressure,rs,dcomp,box,ids);

    }
  }
  else if (name=="Cell_ID") {
    
    BL_ASSERT(dcomp < mf.nComp());
    BL_ASSERT(layout!=0);
    const int ngrow = mf.nGrow();
    
    BoxArray dstBA(mf.boxArray());
    BL_ASSERT(rec->deriveType() == dstBA.ixType());
    
    mf.setVal(-1,dcomp,1,ngrow);
    Layout::IntFab ifab;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
      Box gbox = mf[mfi].box();
      ifab.resize(gbox,1);
      layout->SetNodeIds(ifab,level,mfi.index(),gbox);
      const int* idat = ifab.dataPtr();
      Real* rdat = mf[mfi].dataPtr();
      int numpts = gbox.numPts();
      for (int i=0; i<numpts; ++i) {
        rdat[i] = Real(idat[i]);
      }
    }
  } else {
    
    not_found_yet = true;
  }
  
  if (not_found_yet)
  {
    AmrLevel::derive(name,time,mf,dcomp);
  }
}

bool
Darcy::multilevel_advance(Real  t, 
                          Real  dt, 
                          Real& dt_suggest)
{
  BL_ASSERT(level==0);

  // Set state data in PArrays to pass into solver
  int nLevs = parent->finestLevel() + 1;
  PArray<MultiFab> S_pa_new(nLevs,PArrayNoManage);
  PArray<MultiFab> S_pa_old(nLevs,PArrayNoManage);
  for (int lev=0; lev<nLevs; ++lev) {
    S_pa_new.set(lev,&(getLevel(lev).get_new_data(State_Type)));
    S_pa_old.set(lev,&(getLevel(lev).get_old_data(State_Type)));
    MultiFab::Copy(S_pa_new[lev],S_pa_old[lev],Pressure,Pressure,1,0);
  }

  int ret = snes->Solve(S_pa_new,S_pa_old,RhoSat,S_pa_new,Pressure,dt);
  bool solve_successful = ret > 0;

  // Crude time step control
  //dt_suggest = (solve_successful ? 1.5 : 0.7) * dt;
  dt_suggest = (solve_successful ? 2 : 0.7) * dt;
  return solve_successful;
}

Real
Darcy::advance (Real t,
                Real dt,
                int  iteration,
                int  ncycle)
{
  bool driver_ok = true;
  Real dt_taken = -1;

  if (level == 0) 
  {
    dt_suggest_for_next = dt;
    int max_dt_iters = 20;

    Real dt_try = dt;
    Real dt_this_attempt = dt_try;
    Real subcycle_t = t;
    int dt_iter = 0;
    bool step_ok = false;
    
    bool continue_dt_iteration =
      (dt_this_attempt >= dt_cutoff) 
      && (dt_iter < max_dt_iters) 
      && (dt_taken<0 || subcycle_t < t+dt_taken);

    while (continue_dt_iteration)
    {
      for (int lev=0; lev<=parent->finestLevel(); ++lev) {
        for (int i = 0; i < desc_lst.size(); i++){
          getLevel(lev).state[i].allocOldData();
          getLevel(lev).state[i].swapTimeLevels(dt_this_attempt);
          getLevel(lev).state[i].setTimeLevel(t,dt_this_attempt,dt_this_attempt);
        }

        // Synchronize data at old time
        MultiFab& S_old = get_old_data(State_Type);

        const Layout::MultiIntFab& nodeIds = layout->NodeIds()[level];
        for (MFIter mfi(S_old); mfi.isValid(); ++mfi) {
          const FArrayBox& p = S_old[mfi];
          FArrayBox& rs = S_old[mfi];
          const Box& box = mfi.validbox();
          const Layout::IntFab& ids = nodeIds[mfi];
          snes->ReducedSaturationGivenPressure(p,Pressure,rs,RhoSat,box,ids);
          snes->RhoSatGivenReducedSaturation(rs,RhoSat,rs,RhoSat,box,ids);
        }
      }
      
      if (verbose > 0 && ParallelDescriptor::IOProcessor())
      {
        std::cout << "ADVANCE grids at time = " 
                  << subcycle_t
                  << ", attempting with dt = "
                  << dt_this_attempt
                  << std::endl;
      }

      step_ok = multilevel_advance(subcycle_t,dt_this_attempt,dt_suggest_for_next);
      ParallelDescriptor::ReduceBoolAnd(step_ok);
      ParallelDescriptor::ReduceRealMin(dt_suggest_for_next);
      
      if (step_ok) {
        dt_taken = dt_this_attempt;
        subcycle_t = subcycle_t + dt_taken;
      }
      
      // Little hack to prevent subcycled step from falling just short of target
      Real subtime_remain = t + dt - subcycle_t;
      if (std::abs(subtime_remain) < 1.e-6*dt) {
        subtime_remain = 0;
      }
      if (subtime_remain > dt_suggest_for_next && subtime_remain < 2*dt_suggest_for_next) {
        dt_suggest_for_next = subtime_remain / 2;
      }

      dt_this_attempt = dt_suggest_for_next;
      
      if (!step_ok) {
        for (int lev=0; lev<=parent->finestLevel(); ++lev) {
          for (int i = 0; i < desc_lst.size(); i++){
            getLevel(lev).state[i].swapTimeLevels(dt_this_attempt);
            getLevel(lev).state[i].setTimeLevel(subcycle_t,dt_this_attempt,dt_this_attempt);
          }
        }
      }
      
      dt_iter++;
      
      continue_dt_iteration = 
        subtime_remain > 0
        && (dt_this_attempt >= dt_cutoff) 
        && (dt_iter < max_dt_iters) 
        && (dt_taken<0 || subcycle_t < t+dt);


      ParallelDescriptor::ReduceBoolAnd(continue_dt_iteration);
    }
    
    driver_ok = step_ok && dt_taken >= dt_cutoff &&  dt_iter < max_dt_iters;
    ParallelDescriptor::ReduceBoolAnd(driver_ok);
    
    if (driver_ok && verbose > 0 && ParallelDescriptor::IOProcessor())
    {
      std::cout << "SUCCESS: grids advanced from time: " << t
                << ", with dt: " << dt_taken << ", suggested new dt: " << dt_suggest_for_next
                << std::endl;
    }
  }
  
  return dt_suggest_for_next;
}

DarcyBld Darcy_bld;

LevelBld*
getLevelBld ()
{
  return &Darcy_bld;
}

void
DarcyBld::variableSetUp ()
{
  Darcy::InitializeStaticVariables();
  Darcy::variableSetUp();
}

void
DarcyBld::variableCleanUp ()
{
  Darcy::variableCleanUp();
}

AmrLevel*
DarcyBld::operator() ()
{
  return new Darcy;
}

AmrLevel*
DarcyBld::operator() (Amr&            papa,
                      int             lev,
                      const Geometry& level_geom,
                      const BoxArray& ba,
                      Real            time)
{
  return new Darcy(papa, lev, level_geom, ba, time);
}
