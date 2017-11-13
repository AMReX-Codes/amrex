
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>

namespace amrex
{

  int      AmrLevelAdv::verbose         = 0;
  Real     AmrLevelAdv::cfl             = 0.9;
  Real     AmrLevelAdv::diffco          = 0.0;
  int      AmrLevelAdv::do_reflux       = 1;

  int      AmrLevelAdv::NUM_STATE       = 1;  // One variable in the state
  int      AmrLevelAdv::NUM_GROW        = 5;  // number of ghost cells

  using std::string;

  /////
  void
  AmrLevelAdv::
  fillGhostCellsRK3 (MultiFab & a_phi, //phi on the this level
                     const int& a_stage) //rk4 stage to fill ghost for
  {
    if(level > 0)
    {
      BL_ASSERT(a_stage >= 0);
      BL_ASSERT(a_stage <= 2);
      //this one needs to be prevTime because we just did swapTimeLevels
      Real tf      = getLevel(level  ).state[Phi_Type].prevTime(); 
      Real dt_f    = parent->dtLevel(level);
      Real dt_c    = parent->dtLevel(level-1);
      BoxArray coarGrids         = getLevel(level-1).grids;
      DistributionMapping dmCoar = getLevel(level-1).dmap;
      MultiFab phiC(coarGrids,dmCoar,NUM_STATE,0);
      MultiFab& oldC = getLevel(level-1).get_old_data(Phi_Type);
      Real tc_new  = getLevel(level-1).state[Phi_Type].curTime();
      Real tc_old  = getLevel(level-1).state[Phi_Type].prevTime();
      BL_ASSERT(tf >= tc_old);
      BL_ASSERT(tf <  tc_new);
      Real xi;//coefficient from mccorquodale
      if(a_stage == 0)
      {
        xi = (tf - tc_old)/dt_c;
      }
      else if(a_stage == 1)
      {
        xi = (tf + dt_f - tc_old)/dt_c;
      }
      else if(a_stage == 2)
      {
        xi = (tf + 0.5*dt_f - tc_old)/dt_c; 
      }
      else
      {
        amrex::Error("bogus stage number");
      }

      MultiFab & k1  = getLevel(level-1).m_k1;
      MultiFab & k2  = getLevel(level-1).m_k2;
      BoxArray bak1 = k1.boxArray();
      BoxArray bak2 = k2.boxArray();
      for (MFIter mfi(phiC); mfi.isValid(); ++mfi)
      {
        const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

        timeinterpolaterk3(xi, ARLIM_3D(lo), ARLIM_3D(hi),
                           BL_TO_FORTRAN_3D(  phiC[mfi]),
                           BL_TO_FORTRAN_3D(  oldC[mfi]),
                           BL_TO_FORTRAN_3D(    k1[mfi]),
                           BL_TO_FORTRAN_3D(    k2[mfi]));
      }
      //debug turn off time interp
      //phiC.copy(oldC);
      //end debug

      //now we can spatially interpolate 
      //these need to be in vectors but they are the same data
      //since we have already done time interpolation
      MultiFab phiSave(grids,dmap,NUM_STATE,0);
      phiSave.copy(a_phi);

      BCRec bcs;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        bcs.setLo(idir, BCType::int_dir);  // periodic uses "internal Dirichlet"
        bcs.setHi(idir, BCType::int_dir);  // periodic uses "internal Dirichlet"
      }

      Vector<Real> timevec(2, tf);
      Vector<MultiFab*> coarmf(2, &phiC);
      Vector<MultiFab*> finemf(2, &phiSave);
      int isrc = 0; int idst = 0; int inco = NUM_STATE;
      Geometry geomC = getLevel(level-1).geom;
      //I have no idea why this needs to be here for periodic
      PhysBCFunct cphysbc(geomC,bcs,BndryFunctBase());
      PhysBCFunct fphysbc(geom ,bcs,BndryFunctBase());
      IntVect refrat  = parent->refRatio(level-1);
      Interpolater* mapper = &quartic_interp;
      FillPatchTwoLevels(a_phi, tf, 
                         coarmf, timevec,
                         finemf, timevec,
                         isrc, idst, inco,
                         geomC, geom,
                         cphysbc, fphysbc, refrat,
                         mapper, bcs);
    }

    //for level >0, this might not be necessary.   I cannot tell from the FillPatchUtil documentation.
    //since we are in periodic-land, no need to call FillPatchUtil::FillPatchSingleLevel
    a_phi.FillBoundary(geom.periodicity());
  }
  /////
  void
  AmrLevelAdv::
  fillGhostCellsRK4 (MultiFab & a_phi, //phi on the this level
                     const int& a_stage) //rk4 stage to fill ghost for
  {
    if(level > 0)
    {
      BL_ASSERT(a_stage >= 0);
      BL_ASSERT(a_stage <= 3);
      //this one needs to be prevTime because we just did swapTimeLevels
      Real tf      = getLevel(level  ).state[Phi_Type].prevTime(); 
      Real dt_f    = parent->dtLevel(level);
      Real dt_c    = parent->dtLevel(level-1);
      BoxArray coarGrids         = getLevel(level-1).grids;
      DistributionMapping dmCoar = getLevel(level-1).dmap;
      MultiFab phiC(coarGrids,dmCoar,NUM_STATE,0);
      MultiFab& oldC = getLevel(level-1).get_old_data(Phi_Type);
      Real tc_new  = getLevel(level-1).state[Phi_Type].curTime();
      Real tc_old  = getLevel(level-1).state[Phi_Type].prevTime();
      BL_ASSERT(tf >= tc_old);
      BL_ASSERT(tf <  tc_new);
      Real xi;//coefficient from mccorquodale
      if(a_stage == 0)
      {
        xi = (tf - tc_old)/dt_c;
      }
      else if(a_stage == 1)
      {
        xi = (tf + 0.5*dt_f - tc_old)/dt_c;
      }
      else if(a_stage == 2)
      {
        xi = (tf + 0.5*dt_f - tc_old)/dt_c; //why, yes, this *is* the same as stage 1
      }
      else
      {
        xi = (tf + dt_f - tc_old)/dt_c;
      }

      MultiFab & k1  = getLevel(level-1).m_k1;
      MultiFab & k2  = getLevel(level-1).m_k2;
      MultiFab & k3  = getLevel(level-1).m_k3;
      MultiFab & k4  = getLevel(level-1).m_k4;
      BoxArray bak1 = k1.boxArray();
      BoxArray bak2 = k2.boxArray();
      BoxArray bak3 = k3.boxArray();
      BoxArray bak4 = k4.boxArray();
      for (MFIter mfi(phiC); mfi.isValid(); ++mfi)
      {
        const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

        timeinterpolaterk4(xi, ARLIM_3D(lo), ARLIM_3D(hi),
                           BL_TO_FORTRAN_3D(  phiC[mfi]),
                           BL_TO_FORTRAN_3D(  oldC[mfi]),
                           BL_TO_FORTRAN_3D(    k1[mfi]),
                           BL_TO_FORTRAN_3D(    k2[mfi]),
                           BL_TO_FORTRAN_3D(    k3[mfi]),
                           BL_TO_FORTRAN_3D(    k4[mfi]));
      }
      //debug turn off time interp
      //phiC.copy(oldC);
      //end debug

      //now we can spatially interpolate 
      //these need to be in vectors but they are the same data
      //since we have already done time interpolation
      MultiFab phiSave(grids,dmap,NUM_STATE,0);
      phiSave.copy(a_phi);

      BCRec bcs;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        bcs.setLo(idir, BCType::int_dir);  // periodic uses "internal Dirichlet"
        bcs.setHi(idir, BCType::int_dir);  // periodic uses "internal Dirichlet"
      }

      Vector<Real> timevec(2, tf);
      Vector<MultiFab*> coarmf(2, &phiC);
      Vector<MultiFab*> finemf(2, &phiSave);
      int isrc = 0; int idst = 0; int inco = NUM_STATE;
      Geometry geomC = getLevel(level-1).geom;
      //I have no idea why this needs to be here for periodic
      PhysBCFunct cphysbc(geomC,bcs,BndryFunctBase());
      PhysBCFunct fphysbc(geom ,bcs,BndryFunctBase());
      IntVect refrat  = parent->refRatio(level-1);
      Interpolater* mapper = &quartic_interp;
      FillPatchTwoLevels(a_phi, tf, 
                         coarmf, timevec,
                         finemf, timevec,
                         isrc, idst, inco,
                         geomC, geom,
                         cphysbc, fphysbc, refrat,
                         mapper, bcs);
    }

    //for level >0, this might not be necessary.   I cannot tell from the FillPatchUtil documentation.
    //since we are in periodic-land, no need to call FillPatchUtil::FillPatchSingleLevel
    a_phi.FillBoundary(geom.periodicity());
  }
//
//Default constructor.  Builds invalid object.
//
  AmrLevelAdv::
  AmrLevelAdv ()
  {
    flux_reg = 0;
    initSwitches();
  }

//
//The basic constructor.
//
  AmrLevelAdv::
  AmrLevelAdv (Amr&            papa,
               int             lev,
               const Geometry& level_geom,
               const BoxArray& bl,
               const DistributionMapping& dm,
               Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time) 
  {
    initSwitches();
    flux_reg = 0;
    if (level > 0 && do_reflux)
      flux_reg = new YAFluxRegister(bl, papa.boxArray(level-1),
                                    dm, papa.DistributionMap(level-1),
                                    level_geom, papa.Geom(level-1),
                                    papa.refRatio(level-1), level, NUM_STATE);
  }

//
//The destructor.
//
  AmrLevelAdv::
  ~AmrLevelAdv () 
  {
    delete flux_reg;
  }

//
//Restart from a checkpoint file.
//
  void
  AmrLevelAdv::
  restart (Amr&          papa,
           std::istream& is,
           bool          bReadSpecial)
  {
    AmrLevel::restart(papa,is,bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
      flux_reg = new YAFluxRegister(grids, papa.boxArray(level-1),
                                    dmap, papa.DistributionMap(level-1),
                                    geom, papa.Geom(level-1),
                                    papa.refRatio(level-1), level, NUM_STATE);
  }

  void 
  AmrLevelAdv::
  checkPoint (const std::string& dir,
              std::ostream&      os,
              VisMF::How         how,
              bool               dump_old) 
  {
    AmrLevel::checkPoint(dir, os, how, dump_old);
  }

//
//Write a plotfile to specified directory.
//
  void
  AmrLevelAdv::
  writePlotFile (const std::string& dir,
                 std::ostream&      os,
                 VisMF::How         how)
  {

    AmrLevel::writePlotFile (dir,os,how);
  }

//
//Define data descriptors.
//
  void
  AmrLevelAdv::
  variableSetUp ()
  {
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    read_params();

    desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NUM_STATE,
                           &quartic_interp);

    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i) 
    {
      lo_bc[i] = hi_bc[i] = INT_DIR;   // periodic boundaries
    }
    
    BCRec bc(lo_bc, hi_bc);

    desc_lst.setComponent(Phi_Type, 0, "phi", bc, 
                          StateDescriptor::BndryFunc(nullfill));
  }

//
//Cleanup data descriptors at end of run.
//
  void
  AmrLevelAdv::
  variableCleanUp () 
  {
    desc_lst.clear();
  }

//
//Initialize grid data at problem start-up.
//
  void
  AmrLevelAdv::
  initData ()
  {
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(Phi_Type);
    Real cur_time   = state[Phi_Type].curTime();

    if (verbose) 
    {
      amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
      const Box& box     = mfi.validbox();
      const int* lo      = box.loVect();
      const int* hi      = box.hiVect();

      initdata(level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
               BL_TO_FORTRAN_3D(S_new[mfi]), ZFILL(dx),
               ZFILL(prob_lo));
    }

    if (verbose) 
    {
      amrex::Print() << "Done initializing the level " << level 
                     << " data " << std::endl;
    }
  }

//
//Initialize data on this level from another AmrLevelAdv (during regrid).
//
  void
  AmrLevelAdv::
  initSwitches()
  {
    string algorithm("rk3");
    ParmParse pp;
    pp.query("algorithm", algorithm);
    static bool printedstuff = false;
    bool use_limiting = true;
    pp.query("use_limiting",use_limiting);
    m_use_fixed_dt = false;
    pp.query("use_fixed_dt", m_use_fixed_dt);
    if(m_use_fixed_dt)
    {
      pp.get("fixed_dt", m_fixed_dt);
    }

    if(!printedstuff)
    {
      amrex::Print() << "**** using the " << algorithm << " algorithm ";
      if(use_limiting)
      {
        amrex::Print() << "with limting ON ****" << endl;
      }
      else
      {
        amrex::Print() << "with limting OFF ****" << endl;
      }
      if(m_use_fixed_dt)
      {
        amrex::Print() << "using fixed dt = " << m_fixed_dt  << endl;
      }
      else
      {
        amrex::Print() << "dt calculated dynamically" << endl;
      }
      printedstuff = true;
    }

    m_algorithm = algorithm;
    m_use_limiting = use_limiting;
    if(m_use_limiting)
    {
      m_iuselimit = 1;
    }
    else
    {
      m_iuselimit = 0;
    }

  }
  void
  AmrLevelAdv::
  init (AmrLevel &old)
  {
    AmrLevelAdv* oldlev = (AmrLevelAdv*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Phi_Type].curTime();
    Real prev_time = oldlev->state[Phi_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Phi_Type);

    FillPatch(old, S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);

  }

//
//Initialize data on this level after regridding if old level did not previously exist
//
  void
  AmrLevelAdv::
  init ()
  {
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);

  }

  Real
  AmrLevelAdv::
  advance (Real time,
           Real dt,
           int  iteration,
           int  ncycle)
  {
    BL_PROFILE("AmrLevelAdv::advance()");
    MultiFab& S_new = get_new_data(Phi_Type);
    Real maxval = S_new.max(0);
    Real minval = S_new.min(0);
    amrex::Print() << "phi max = " << maxval << ", min = " << minval  << endl;
    Real dtval;

    if(m_algorithm.compare(string("godunov")) == 0)
    {
      dtval = advanceGodunov(time, dt, iteration, ncycle);
    }
    else if(m_algorithm.compare(string("rk2")) == 0)
    {
      dtval = advanceMOLRK2(time, dt, iteration, ncycle);
    }
    else if(m_algorithm.compare(string("rk3")) == 0)
    {
      dtval = advanceMOLRK3(time, dt, iteration, ncycle);
    }
    else if(m_algorithm.compare(string("rk4")) == 0)
    {
      dtval = advanceMOLRK4(time, dt, iteration, ncycle);
    }
    else
    {
      amrex::Error("bogus algorithm parameter");
    }


    if(m_use_fixed_dt);
    {
      dtval = m_fixed_dt;
    }
    return dtval;
  }

  /////
  Real
  AmrLevelAdv::
  advanceMOLRK2 (Real time,
                 Real dt,
                 int  iteration,
                 int  ncycle)
  {
    bool truncationErrorTest = false;
    ParmParse pp;
    pp.query("truncation_error_only", truncationErrorTest);
    if(truncationErrorTest && (iteration > 1))
    {
      return dt;
    }

    for (int i = 0; i < NUM_STATE_TYPE; ++i) 
    {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);
    MultiFab& S_old = get_old_data(Phi_Type);
    MultiFab dPhiDt(grids,dmap,NUM_STATE,0);
    MultiFab Sborder(grids,dmap,NUM_STATE,NUM_GROW);
  
    YAFluxRegister* fr_as_crse = nullptr;
    if (do_reflux && level < parent->finestLevel()) 
    {
      fr_as_crse = &getFluxReg(level+1);
    }

    YAFluxRegister* fr_as_fine = nullptr;
    if (do_reflux && level > 0) 
    {
      fr_as_fine = &getFluxReg(level);
    }

    if (fr_as_crse) 
    {
      fr_as_crse->reset();
    }

    // RK2 stage 1
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    compute_dPhiDt_MOL2ndOrd(Sborder, dPhiDt, time, 0.5*dt, fr_as_crse, fr_as_fine, iteration);

    if(truncationErrorTest)
    {
      S_new.copy(dPhiDt);
      return dt;
    }
    

    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dPhiDt, 0, 0, NUM_STATE, 0);
    /**/
    // RK2 stage 2
    // After fillpatch Sborder = U^n+dt*dUdt^n
    FillPatch(*this, Sborder, NUM_GROW, time+dt, Phi_Type, 0, NUM_STATE);
    compute_dPhiDt_MOL2ndOrd(Sborder, dPhiDt, time, 0.5*dt, fr_as_crse, fr_as_fine, iteration);
    // S_new = 0.5*(Sborder+S_old) = U^n + 0.5*dt*dUdt^n
    MultiFab::LinComb(S_new, 0.5, Sborder, 0, 0.5, S_old, 0, 0, NUM_STATE, 0);
    // S_new += 0.5*dt*dPhiDt
    MultiFab::Saxpy(S_new, 0.5*dt, dPhiDt, 0, 0, NUM_STATE, 0);
    // We now have S_new = U^{n+1} = (U^n+0.5*dt*dUdt^n) + 0.5*dt*dUdt^*
    
    /**/ 
    return dt;
  }

  Real
  AmrLevelAdv::
  advanceMOLRK3 (Real time,
                 Real dt,
                 int  iteration,
                 int  ncycle)
  {
    bool truncationErrorTest = false;
    ParmParse pp;
    pp.query("truncation_error_only", truncationErrorTest);

    if(truncationErrorTest && (iteration > 1))
    {
      return dt;
    }

    for (int i = 0; i < NUM_STATE_TYPE; ++i) 
    {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);
    MultiFab& S_old = get_old_data(Phi_Type);

    MultiFab u1(grids,dmap,NUM_STATE,NUM_GROW);
    MultiFab u2(grids,dmap,NUM_STATE,NUM_GROW);
    MultiFab k1(grids,dmap,NUM_STATE,0);
    MultiFab k2(grids,dmap,NUM_STATE,0);
    MultiFab k3(grids,dmap,NUM_STATE,0);
  
    YAFluxRegister* fr_as_crse = nullptr;
    if (do_reflux && level < parent->finestLevel() ) 
    {
      fr_as_crse = &getFluxReg(level+1);
    }

    YAFluxRegister* fr_as_fine = nullptr;
    if (do_reflux && level > 0) 
    {
      fr_as_fine = &getFluxReg(level);
    }

    if (fr_as_crse) 
    {
      fr_as_crse->reset();
    }


    // RK3 stage 1
    //FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    u1.copy(S_old);
    fillGhostCellsRK3(u1, 0);
    //the dt/6 is for the flux register.
    compute_dPhiDt_MOL4thOrd(u1, k1, time, dt/6., fr_as_crse, fr_as_fine, iteration);

    if(truncationErrorTest)
    {
      S_new.copy(k1);
      //cannot just return here because we still need the k coefs for finer levels
      //return dt;
    }

    //u1 is already set to u^n
    // this sets U1 = U^n + dt*dPhiDt^n
    MultiFab::Saxpy(u1, dt, k1, 0, 0, NUM_STATE, 0);
    k1.mult(dt);

    fillGhostCellsRK3(u1, 1);
    compute_dPhiDt_MOL4thOrd(u1, k2, time, dt/6., fr_as_crse, fr_as_fine, iteration);
    //           u2 = 3/4 U^n + 1/4 U^1 + 1/4 dPhiDt^1
    //this makes u2 = 3/4 u^n + 1/4 u^1
    MultiFab::LinComb(u2, 0.75, S_old, 0, 0.25, u1, 0, 0, NUM_STATE, 0);
    //this makes u2 = 3/4 u^n + 1/4 u^1 + 1/4*dt*dPhiDt^1
    MultiFab::Saxpy(u2, 0.25*dt, k2, 0, 0, NUM_STATE, 0);
    k2.mult(dt);

    fillGhostCellsRK3(u2, 2);
    //the 2*dt/3 is for the flux register.
    compute_dPhiDt_MOL4thOrd(u2, k3, time, 2.*dt/3., fr_as_crse, fr_as_fine, iteration);

    // RK3 stage 3
    //           S_new = 1/3 u^n + 2/3 u^2 + 2/3 dtL(u^2)

    if(!truncationErrorTest)
    {
      //this makes S_new = 1/3 u^n + 2/3 u^2
      MultiFab::LinComb(S_new, 1./3., S_old, 0, 2./3., u2, 0, 0, NUM_STATE, 0);
      //this makes S_new = 1/3 u^n + 2/3  u^2 + 2/3 dtL(u^2)
      MultiFab::Saxpy(S_new, 2.*dt/3., k3, 0, 0, NUM_STATE, 0);
    }

    if(level < parent->finestLevel() ) 
    {
      m_k1.define(grids,dmap,NUM_STATE,0);
      m_k2.define(grids,dmap,NUM_STATE,0);

      m_k1.copy(k1);
      m_k2.copy(k2);
    }
    /**/ 
    return dt;
  }


  Real
  AmrLevelAdv::
  advanceMOLRK4 (Real time,
                 Real dt,
                 int  iteration,
                 int  ncycle)
  {
    bool truncationErrorTest = false;
    ParmParse pp;
    pp.query("truncation_error_only", truncationErrorTest);
    if(truncationErrorTest && (iteration > 1))
    {
      return dt;
    }

    for (int i = 0; i < NUM_STATE_TYPE; ++i) 
    {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);
    MultiFab& S_old = get_old_data(Phi_Type);

    MultiFab u1(grids,dmap,NUM_STATE,NUM_GROW);
    MultiFab u2(grids,dmap,NUM_STATE,NUM_GROW);
    MultiFab u3(grids,dmap,NUM_STATE,NUM_GROW);

    MultiFab k1(grids,dmap,NUM_STATE,0);
    MultiFab k2(grids,dmap,NUM_STATE,0);
    MultiFab k3(grids,dmap,NUM_STATE,0);
    MultiFab k4(grids,dmap,NUM_STATE,0);
  
    YAFluxRegister* fr_as_crse = nullptr;
    if (do_reflux && level < parent->finestLevel()) 
    {
      fr_as_crse = &getFluxReg(level+1);
    }

    YAFluxRegister* fr_as_fine = nullptr;
    if (do_reflux && level > 0) 
    {
      fr_as_fine = &getFluxReg(level);
    }

    if (fr_as_crse) 
    {
      fr_as_crse->reset();
    }

    //phi^1 = phi^n + dt*F(phi^n)
    // RK3 stage 1
    u1.copy(S_old);
    fillGhostCellsRK4(u1, 0);
    //the dt/6 is for the flux register.
    compute_dPhiDt_MOL4thOrd(u1, k1, time, dt/6., fr_as_crse, fr_as_fine, iteration);

    if(truncationErrorTest)
    {
      S_new.copy(k1);
      //cannot just return here because we still need the k coefs for finer levels
      //return dt;
    }

    // phi^1 = phi^n + dt/2*dPhiDt^n
    // this sets U1 = U^n + dt*dPhiDt^n 
    k1.mult(dt);
    MultiFab::LinComb(u1, 1., S_old, 0, 0.5, k1, 0, 0, NUM_STATE, 0);
    fillGhostCellsRK4(u1, 1);

    //phi^2 = phi^n + dt/2*dPhiDt(phi^1)
    //the dt/3 is for the flux register.
    compute_dPhiDt_MOL4thOrd(u1, k2, time, dt/3., fr_as_crse, fr_as_fine, iteration);
    k2.mult(dt);
    MultiFab::LinComb(u2, 1., S_old, 0, 0.5, k2, 0, 0, NUM_STATE, 0);
    fillGhostCellsRK4(u2, 2);

    //phi^3 = phi^n + dt*dPhiDt(phi^2)
    //the dt/3 is for the flux register.
    compute_dPhiDt_MOL4thOrd(u2, k3, time, dt/3., fr_as_crse, fr_as_fine, iteration);
    k3.mult(dt);
    MultiFab::LinComb(u3, 1., S_old, 0, 1.0, k3, 0, 0, NUM_STATE, 0);
    fillGhostCellsRK4(u3, 3);

    //phi^4 = phi^n + dt*F(phi^3)
    //the dt/6. is for the flux register.
    compute_dPhiDt_MOL4thOrd(u3, k4, time, dt/6., fr_as_crse, fr_as_fine, iteration);
    k4.mult(dt);

    //phi^n+1  = 1/6(phi1 +  2 phi2  + 2phi3  + phi4)
    if(!truncationErrorTest)
    {
      S_new.copy(S_old);

      MultiFab::Saxpy(  S_new, 1./6., k1, 0, 0, NUM_STATE, 0);
      MultiFab::Saxpy(  S_new, 1./3., k2, 0, 0, NUM_STATE, 0);
      MultiFab::Saxpy(  S_new, 1./3., k3, 0, 0, NUM_STATE, 0);
      MultiFab::Saxpy(  S_new, 1./6., k4, 0, 0, NUM_STATE, 0);
    }
    //if we are not on the finest level, need to save the k objects for c/f time interpolation
    if(level < parent->finestLevel() ) 
    {
      m_k1.define(grids,dmap,NUM_STATE,0);
      m_k2.define(grids,dmap,NUM_STATE,0);
      m_k3.define(grids,dmap,NUM_STATE,0);
      m_k4.define(grids,dmap,NUM_STATE,0);

      m_k1.copy(k1);
      m_k2.copy(k2);
      m_k3.copy(k3);
      m_k4.copy(k4);
    }
    return dt;
  }

  Real
  AmrLevelAdv::
  advanceGodunov (Real time,
                  Real dt,
                  int  iteration,
                  int  ncycle)
  {
    BL_PROFILE("AmrLevelAdv::advance()");
        
    bool truncationErrorTest = false;
    ParmParse pp;
    pp.query("truncation_error_only", truncationErrorTest);
    if(truncationErrorTest && (iteration > 1))
    {
      return dt;
    }

    for (int i = 0; i < NUM_STATE_TYPE; ++i) 
    {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);
    MultiFab dPhiDt(grids,dmap,NUM_STATE,0);
    MultiFab Sborder(grids,dmap,NUM_STATE,NUM_GROW);
  
    YAFluxRegister* fr_as_crse = nullptr;
    if (do_reflux && level < parent->finestLevel()) 
    {
      fr_as_crse = &getFluxReg(level+1);
    }

    YAFluxRegister* fr_as_fine = nullptr;
    if (do_reflux && level > 0) 
    {
      fr_as_fine = &getFluxReg(level);
    }

    if (fr_as_crse) 
    {
      fr_as_crse->reset();
    }

/** doing this as Godunov first to make sure I did not 
    screw anything up in this reorganization
*/
    // overly complicated godunov
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    //dt needs to be there because it weights the fluxes
    compute_dPhiDt_godunov(Sborder, dPhiDt, time, dt, fr_as_crse, fr_as_fine, iteration);

    if(truncationErrorTest)
    {
      S_new.copy(dPhiDt);
      return dt;
    }

    // U^n+1 = U^n + dt*dUdt^n
    MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dPhiDt, 0, 0, NUM_STATE, 0);
    return dt;
  }
//
//Advance grids at this level in time.
// computes dphi/dt = -div(F).   Needs dt to be sent in because YAFluxRegister
// needs the fluxes to be multiplied by dt*area
  amrex::Real
  AmrLevelAdv::
  compute_dPhiDt_godunov (const MultiFab& Sborder, MultiFab& dPhiDt, Real time, Real dt,
                          YAFluxRegister* fr_as_crse, YAFluxRegister* fr_as_fine, 
                          int iteration)
  {

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      std::array<FArrayBox,AMREX_SPACEDIM> flux;
      std::array<FArrayBox,AMREX_SPACEDIM> uface;

      for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();

        const FArrayBox& statein = Sborder[mfi];
        FArrayBox& dphidtout     =  dPhiDt[mfi];

        // Allocate fabs for fluxes and Godunov velocities.
        for (int i = 0; i < BL_SPACEDIM ; i++) 
        {
          const Box& bxtmp = amrex::surroundingNodes(bx,i);
          flux[i].resize(bxtmp,NUM_STATE);
          //I have no idea why iteration comes into this but
          //it is inherited from the advection example
          uface[i].resize(amrex::grow(bxtmp, iteration), 1);
        }

        const Real prev_time = state[Phi_Type].prevTime();
        const Real cur_time = state[Phi_Type].curTime();
        //velocity for godunov is at time n+1/2 because that is wheree flux is centered
        const Real ctr_time = 0.5*(prev_time + cur_time); 

        get_face_velocity(level, ctr_time,
                          AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                                       BL_TO_FORTRAN(uface[1]),
                                       BL_TO_FORTRAN(uface[2])),
                          dx, prob_lo);

        advectDiffGodunov(time, bx.loVect(), bx.hiVect(),
                          BL_TO_FORTRAN_3D(statein), 
                          BL_TO_FORTRAN_3D(dphidtout),
                          AMREX_D_DECL(BL_TO_FORTRAN_3D(uface[0]),
                                       BL_TO_FORTRAN_3D(uface[1]),
                                       BL_TO_FORTRAN_3D(uface[2])),
                          AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]), 
                                       BL_TO_FORTRAN_3D(flux[1]), 
                                       BL_TO_FORTRAN_3D(flux[2])), 
                          dx, dt, diffco, &m_iuselimit);

        if(do_reflux)
        {
          std::array<const FArrayBox*,AMREX_SPACEDIM> fluxPtrs;
          for(int idir = 0; idir < SpaceDim; idir++)
          {
            fluxPtrs[idir] = &flux[idir];
          }
          if (fr_as_crse) 
          {
          
            fr_as_crse->CrseAdd(mfi,fluxPtrs,dx,dt);
          }

          if (fr_as_fine) 
          {
            fr_as_fine->FineAdd(mfi,fluxPtrs,dx,dt);
          }
        }

      }
    }


    return dt;
  }


//
//Advance grids at this level in time.
// computes dphi/dt = -div(F).   Needs dt to be sent in because YAFluxRegister
// needs the fluxes to be multiplied by dt*area
  Real
  AmrLevelAdv::
  compute_dPhiDt_MOL4thOrd (const MultiFab& Sborder, MultiFab& dPhiDt, Real time, Real dt,
                            YAFluxRegister* fr_as_crse, YAFluxRegister* fr_as_fine,
                            int iteration)
  {

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    int reftocoarsest=1;
    ParmParse pp;
    pp.query("ref_to_coarsest", reftocoarsest);
    Box domain = geom.Domain();

    IntVect startpt = IntVect::Zero;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      startpt[idir] = domain.size()[idir]/4;
    }
    startpt.coarsen(reftocoarsest);
    Box debboxcc(startpt, startpt);
    debboxcc.refine(reftocoarsest);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

      for (MFIter mfi(Sborder); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.validbox();

        const FArrayBox& statein = Sborder[mfi];
        FArrayBox& dphidtout     =  dPhiDt[mfi];

        Box grownBox = grow(bx,NUM_GROW);
        // Allocate fabs for fluxes and Godunov velocities.
        for (int i = 0; i < BL_SPACEDIM ; i++) 
        {
          const Box& bxtmp = amrex::surroundingNodes(bx,i);
          flux[i].resize(bxtmp,NUM_STATE);
          //here I have to stop the madness of using iteration because
          //I have a bigger stencil.
          Box velBox = surroundingNodes(grownBox, i);
          uface[i].resize(velBox, 1);
        }

        int debdir = 0;
        Box debugboxcell = debboxcc & bx;
        int printstuff = 0;
        Box debugboxfaceHi;
        Box debugboxfaceLo; 
        if(debugboxcell.ok())
        {
          debugboxfaceHi = bdryHi(debugboxcell, debdir, 1);
          debugboxfaceLo = bdryLo(debugboxcell, debdir, 1);
          printstuff = 1;
        }

        //velocity is a time because this is MOL
        const Real ctr_time = time;

        get_face_velocity(level, ctr_time,
                          AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                                       BL_TO_FORTRAN(uface[1]),
                                       BL_TO_FORTRAN(uface[2])),
                          dx, prob_lo);

        advectDiffMOL4thOrd(time, bx.loVect(), bx.hiVect(),
                            BL_TO_FORTRAN_3D(statein), 
                            BL_TO_FORTRAN_3D(dphidtout),
                            AMREX_D_DECL(BL_TO_FORTRAN_3D(uface[0]),
                                         BL_TO_FORTRAN_3D(uface[1]),
                                         BL_TO_FORTRAN_3D(uface[2])),
                            AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]), 
                                         BL_TO_FORTRAN_3D(flux[1]), 
                                         BL_TO_FORTRAN_3D(flux[2])), 
                            dx, dt, diffco, 
                            debugboxcell.loVect(), debugboxcell.hiVect(),
                            debugboxfaceHi.loVect(), debugboxfaceHi.hiVect(),
                            debugboxfaceLo.loVect(), debugboxfaceLo.hiVect(),
                           &printstuff, &m_iuselimit
                          );

        if(do_reflux)
        {
          std::array<const FArrayBox*,AMREX_SPACEDIM> fluxPtrs;
          for(int idir = 0; idir < SpaceDim; idir++)
          {
            fluxPtrs[idir] = &flux[idir];
          }
          if (fr_as_crse) 
          {
          
            fr_as_crse->CrseAdd(mfi,fluxPtrs,dx,dt);
          }

          if (fr_as_fine) 
          {
            fr_as_fine->FineAdd(mfi,fluxPtrs,dx,dt);
          }
        }
      }
    }


    return dt;
  }



//
//Advance grids at this level in time.
// computes dphi/dt = -div(F).   Needs dt to be sent in because YAFluxRegister
// needs the fluxes to be multiplied by dt*area
  Real
  AmrLevelAdv::
  compute_dPhiDt_MOL2ndOrd (const MultiFab& Sborder, MultiFab& dPhiDt, Real time, Real dt,
                            YAFluxRegister* fr_as_crse, YAFluxRegister* fr_as_fine,
                            int iteration)
  {

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

      for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();

        const FArrayBox& statein = Sborder[mfi];
        FArrayBox& dphidtout     =  dPhiDt[mfi];

        // Allocate fabs for fluxes and Godunov velocities.
        for (int i = 0; i < BL_SPACEDIM ; i++) 
        {
          const Box& bxtmp = amrex::surroundingNodes(bx,i);
          flux[i].resize(bxtmp,NUM_STATE);
          //I have no idea why iteration comes into this but
          //it is inherited from the advection example (stencil size here is the same as Goudnov)
          uface[i].resize(amrex::grow(bxtmp, iteration), 1);
        }

        //MOL
        const Real ctr_time = time;

        get_face_velocity(level, ctr_time,
                          AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                                       BL_TO_FORTRAN(uface[1]),
                                       BL_TO_FORTRAN(uface[2])),
                          dx, prob_lo);

        advectDiffMOL2ndOrd(time, bx.loVect(), bx.hiVect(),
                            BL_TO_FORTRAN_3D(statein), 
                            BL_TO_FORTRAN_3D(dphidtout),
                            AMREX_D_DECL(BL_TO_FORTRAN_3D(uface[0]),
                                         BL_TO_FORTRAN_3D(uface[1]),
                                         BL_TO_FORTRAN_3D(uface[2])),
                            AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]), 
                                         BL_TO_FORTRAN_3D(flux[1]), 
                                         BL_TO_FORTRAN_3D(flux[2])), 
                            dx, dt, diffco,  &m_iuselimit);

        if(do_reflux)
        {
          std::array<const FArrayBox*,AMREX_SPACEDIM> fluxPtrs;
          for(int idir = 0; idir < SpaceDim; idir++)
          {
            fluxPtrs[idir] = &flux[idir];
          }
          if (fr_as_crse) 
          {
          
            fr_as_crse->CrseAdd(mfi,fluxPtrs,dx,dt);
          }

          if (fr_as_fine) 
          {
            fr_as_fine->FineAdd(mfi,fluxPtrs,dx,dt);
          }
        }
      }
    }


    return dt;
  }

//
//Estimate time step.
//
  Real
  AmrLevelAdv::
  estTimeStep (Real)
  {
    // This is just a dummy value to start with 
    Real dt_est  = 1.0e+20;

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real cur_time = state[Phi_Type].curTime();
    const MultiFab& S_new = get_new_data(Phi_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
      FArrayBox uface[BL_SPACEDIM];

      for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
      {
        for (int i = 0; i < BL_SPACEDIM ; i++) 
        {
          const Box& bx = mfi.nodaltilebox(i);
          uface[i].resize(bx,1);
        }

        get_face_velocity(level, cur_time,
                          AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                                       BL_TO_FORTRAN(uface[1]),
                                       BL_TO_FORTRAN(uface[2])),
                          dx, prob_lo);

        for (int i = 0; i < BL_SPACEDIM; ++i) 
        {
          Real umax = uface[i].norm(0);
          if (umax > 1.e-100) 
          {
            dt_est = std::min(dt_est, dx[i] / umax);
          }
        }
      }
    }

    ParallelDescriptor::ReduceRealMin(dt_est);
    dt_est *= cfl;

    ParmParse pp;
    if(pp.contains("fixed_dt"))
    {
      pp.get("fixed_dt", dt_est);
    }
    if (verbose) 
    {
      amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level 
                     << ":  dt_est = " << dt_est << std::endl;
    }
    
    return dt_est;
  }

//
//Compute initial time step.
//
  Real
  AmrLevelAdv::
  initialTimeStep ()
  {
    return estTimeStep(0.0);
  }

//
//Compute initial `dt'.
//
  void
  AmrLevelAdv::
  computeInitialDt (int                   finest_level,
                    int                   sub_cycle,
                    Vector<int>&           n_cycle,
                    const Vector<IntVect>& ref_ratio,
                    Vector<Real>&          dt_level,
                    Real                  stop_time)
  {
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
      return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
      dt_level[i] = getLevel(i).initialTimeStep();
      n_factor   *= n_cycle[i];
      dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) 
    {
      if ((cur_time + dt_0) > (stop_time - eps))
        dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
      n_factor *= n_cycle[i];
      dt_level[i] = dt_0/n_factor;
    }
  }

//
//Compute new `dt'.
//
  void
  AmrLevelAdv::
  computeNewDt (int                   finest_level,
                int                   sub_cycle,
                Vector<int>&           n_cycle,
                const Vector<IntVect>& ref_ratio,
                Vector<Real>&          dt_min,
                Vector<Real>&          dt_level,
                Real                  stop_time,
                int                   post_regrid_flag)
  {
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
      return;

    for (int i = 0; i <= finest_level; i++)
    {
      AmrLevelAdv& adv_level = getLevel(i);
      dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1) 
    {
      //
      // Limit dt's by pre-regrid dt
      //
      for (int i = 0; i <= finest_level; i++)
      {
        dt_min[i] = std::min(dt_min[i],dt_level[i]);
      }
    }
    else 
    {
      //
      // Limit dt's by change_max * old dt
      //
      static Real change_max = 1.1;
      for (int i = 0; i <= finest_level; i++)
      {
        dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
      }
    }
    
    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
      n_factor *= n_cycle[i];
      dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
      if ((cur_time + dt_0) > (stop_time - eps))
        dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
      n_factor *= n_cycle[i];
      dt_level[i] = dt_0/n_factor;
    }
  }

//
//Do work after timestep().
//
  void
  AmrLevelAdv::
  post_timestep (int iteration)
  {
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    bool truncationErrorTest = false;
    ParmParse pp;
    pp.query("truncation_error_only", truncationErrorTest);

    if(!truncationErrorTest)
    {
      if (do_reflux && level < finest_level)
        reflux();

      if (level < finest_level)
        avgDown();
    }

  }

//
//Do work after regrid().
//
  void
  AmrLevelAdv::
  post_regrid (int lbase, int new_finest) 
  {
  }

//
//Do work after a restart().
//
  void
  AmrLevelAdv::
  post_restart() 
  {
  }

//
//Do work after init().
//
  void
  AmrLevelAdv::
  post_init (Real stop_time)
  {
    if (level > 0)
      return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
      getLevel(k).avgDown();
  }

//
//Error estimation for regridding.
//
  void
  AmrLevelAdv::
  errorEst (TagBoxArray& tags,
            int          clearval,
            int          tagval,
            Real         time,
            int          n_error_buf,
            int          ngrow)
  {
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(Phi_Type);

    
    ParmParse pp;
    int tag_domain_middle = 0;
    if(pp.contains("tag_domain_middle"))
    {
      pp.get("tag_domain_middle", tag_domain_middle);
    }
    Box domain = geom.Domain();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      Vector<int>  itags;
	
      for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
      {
        const Box&  tilebx  = mfi.tilebox();

        TagBox&     tagfab  = tags[mfi];
	    
        // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
        // So we are going to get a temporary integer array.
        tagfab.get_itags(itags, tilebx);
	    
        // data pointer and index space
        int*        tptr    = itags.dataPtr();
        const int*  tlo     = tilebx.loVect();
        const int*  thi     = tilebx.hiVect();
        const int*  dlo     = domain.loVect();
        const int*  dhi     = domain.hiVect();

        state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                    BL_TO_FORTRAN_3D(S_new[mfi]),
                    &tagval, &clearval, 
                    ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                    ZFILL(dx), ZFILL(prob_lo), &time, &level,
                    ARLIM_3D(dlo), ARLIM_3D(dhi), &tag_domain_middle);
        //
        // Now update the tags in the TagBox.
        //
        tagfab.tags_and_untags(itags, tilebx);
      }
    }
  }

  void
  AmrLevelAdv::
  read_params ()
  {
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");   

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.get("diffusion_coeff",diffco);
    pp.query("do_reflux",do_reflux);

    // This tutorial code only supports Cartesian coordinates.
    if (! Geometry::
        IsCartesian()) 
    {
      amrex::Abort("Please set geom.coord_sys = 0");
    }

    // This tutorial code only supports periodic boundaries.
    if (! Geometry::isAllPeriodic()) 
    {
      amrex::Abort("Please set geom.is_periodic = 1 1 1");
    }




    //
    // read tagging parameters from probin file
    //

    std::string probin_file("probin");

    ParmParse ppa("amr");
    ppa.query("probin_file",probin_file);

    int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
      probin_file_name[i] = probin_file[i];

    // use a fortran routine to
    // read in tagging parameters from probin file
    get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

  }

  void
  AmrLevelAdv::
  reflux ()
  {
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = ParallelDescriptor::second();

    getFluxReg(level+1).Reflux(get_new_data(Phi_Type));
    
    if (verbose)
    {
      const int IOProc = ParallelDescriptor::IOProcessorNumber();
      Real      end    = ParallelDescriptor::second() - strt;
	
      ParallelDescriptor::ReduceRealMax(end,IOProc);
	
      amrex::Print() << "AmrLevelAdv::reflux() at level " << level 
                     << " : time = " << end << std::endl;
    }
  }

  void
  AmrLevelAdv::
  avgDown ()
  {
    if (level == parent->finestLevel()) return;
    avgDown(Phi_Type);
  }

  void
  AmrLevelAdv::
  avgDown (int state_indx)
  {
    if (level == parent->finestLevel()) return;

    AmrLevelAdv& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);
    
    amrex::average_down(S_fine,S_crse,
                        fine_lev.geom,geom,
                        0,S_fine.nComp(),parent->refRatio(level));
  }

}
