
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

namespace amrex
{

  int      AmrLevelAdv::verbose         = 0;
  Real     AmrLevelAdv::cfl             = 0.9;
  Real     AmrLevelAdv::diffco          = 0.0;
  int      AmrLevelAdv::do_reflux       = 1;

  int      AmrLevelAdv::NUM_STATE       = 1;  // One variable in the state
  int      AmrLevelAdv::NUM_GROW        = 3;  // number of ghost cells


//
//Default constructor.  Builds invalid object.
//
  AmrLevelAdv::
  AmrLevelAdv ()
  {
    flux_reg = 0;
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
    dtval = advanceGodunov(time, dt, iteration, ncycle);
//    dtval = advanceMOLRK2(time, dt, iteration, ncycle);
//    dtval = advanceMOLRK3(time, dt, iteration, ncycle);
//    dtval = advanceMOLRK4(time, dt, iteration, ncycle);
    ParmParse pp;
    if(pp.contains("fixed_dt"))
    {
      pp.get("fixed_dt", dtval);
    }
    return dtval;
  }


  Real
  AmrLevelAdv::
  advanceMOLRK2 (Real time,
                 Real dt,
                 int  iteration,
                 int  ncycle)
  {
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

    bool truncationErrorTest = false;
    ParmParse pp;
    pp.query("truncation_error_only", truncationErrorTest);
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
    for (int i = 0; i < NUM_STATE_TYPE; ++i) 
    {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);
    MultiFab& S_old = get_old_data(Phi_Type);
    MultiFab dPhiDt(grids,dmap,NUM_STATE,0);

    MultiFab u1(grids,dmap,NUM_STATE,NUM_GROW);
    MultiFab u2(grids,dmap,NUM_STATE,NUM_GROW);
  
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


    // RK3 stage 1
    //FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    S_old.FillBoundary();
    //the dt/6 is for the flux register.
    compute_dPhiDt_MOL4thOrd(S_old, dPhiDt, time, dt/6., fr_as_crse, fr_as_fine, iteration);

    bool truncationErrorTest = false;
    ParmParse pp;
    pp.query("truncation_error_only", truncationErrorTest);
    if(truncationErrorTest)
    {
      S_new.copy(dPhiDt);
      return dt;
    }

    // U^1 = U^n + dt*dPhiDt^n
    MultiFab::LinComb(u1, 1.0, S_old, 0, dt, dPhiDt, 0, 0, NUM_STATE, 0);
    u1.FillBoundary();

    compute_dPhiDt_MOL4thOrd(u1, dPhiDt, time, dt/6., fr_as_crse, fr_as_fine, iteration);

    //           u2 = 3/4 U^n + 1/4 U^1 + 1/4 dPhiDt^1
    //this makes u2 = 3/4 u^n + 1/4 u^1
    MultiFab::LinComb(u2, 0.75, S_old, 0, 0.25, u1, 0, 0, NUM_STATE, 0);
    //this makes u2 = 3/4 u^n + 1/4 u^1 + 1/4*dt*dPhiDt^1
    MultiFab::Saxpy(u2, 0.25*dt, dPhiDt, 0, 0, NUM_STATE, 0);

    u2.FillBoundary();
    //the 2*dt/3 is for the flux register.
    compute_dPhiDt_MOL4thOrd(u2, dPhiDt, time, 2.*dt/3., fr_as_crse, fr_as_fine, iteration);

    // RK3 stage 3
    //           S_new = 1/3 u^n + 2/3 u^2 + 2/3 dtL(u^2)

    //this makes S_new = 1/3 u^n + 2/3 u^2
    MultiFab::LinComb(S_new, 1./3., S_old, 0, 2./3., u2, 0, 0, NUM_STATE, 0);
    //this makes S_new = 1/3 u^n + 2/3  u^2 + 2/3 dtL(u^2)
    MultiFab::Saxpy(S_new, 2.*dt/3., dPhiDt, 0, 0, NUM_STATE, 0);

    /**/ 
    return dt;
  }

#if 0
  Real
  AmrLevelAdv::
  advanceMOLRK4 (Real time,
                 Real dt,
                 int  iteration,
                 int  ncycle)
  {
    for (int i = 0; i < NUM_STATE_TYPE; ++i) 
    {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);
    MultiFab& S_old = get_old_data(Phi_Type);
    MultiFab dPhiDt(grids,dmap,NUM_STATE,0);
    MultiFab Sborder(grids,dmap,NUM_STATE,NUM_GROW);
    MultiFab u1(grids,dmap,NUM_STATE,NUM_GROW);
    MultiFab u2(grids,dmap,NUM_STATE,NUM_GROW);
    MultiFab u3(grids,dmap,NUM_STATE,NUM_GROW);
    MultiFab u4(grids,dmap,NUM_STATE,NUM_GROW);
  
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

    // RK3 stage 1
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE); // 
    //the dt/6 is for the flux register.
    compute_dPhiDt_MOL4thOrd(Sborder, dPhiDt, time, dt/6., fr_as_crse, fr_as_fine, iteration);

    // U^1 = U^n + dt*dPhiDt^n
    MultiFab::LinComb(u1, 1.0, S_old, 0, dt, dPhiDt, 0, 0, NUM_STATE, 0);
    ///fillpatch works with uold and unew
    ///so need to put this into S_new
    ///need to also save u1 
    S_new.copy(u1);
    
    /**/
    // RK3 stage 2
    // After fillpatch Sborder = U^n+dt*dUdt^n = u1 (including ghost)
    FillPatch(*this, Sborder, NUM_GROW, time+dt, Phi_Type, 0, NUM_STATE);
    //the dt/6 is for the flux register.
    compute_dPhiDt_MOL4thOrd(Sborder, dPhiDt, time, dt/6., fr_as_crse, fr_as_fine, iteration);

    //           u2 = 3/4 U^n + 1/4 U^1 + 1/4 dPhiDt^1
    //this makes u2 = 3/4 u^n + 1/4 u^1
    MultiFab::LinComb(u2, 0.75, S_old, 0, 0.25, u1, 0, 0, NUM_STATE, 0);
    //this makes u2 = 3/4 u^n + 1/4 u^1 + 1/4*dt*dPhiDt^1
    MultiFab::Saxpy(u2, 0.25*dt, dPhiDt, 0, 0, NUM_STATE, 0);
    ///fillpatch works with uold and unew
    ///so need to put this into S_new so u2 is the thing in the ghost cells
    S_new.copy(u2);

    // After fillpatch Sborder = u2.   I am not so sure about what dt to use here.  
    //I am using time+dt so that, away from coarse-fine boundaries, no interpolation is done.
    FillPatch(*this, Sborder, NUM_GROW, time+dt, Phi_Type, 0, NUM_STATE);
    //the 2*dt/3 is for the flux register.
    compute_dPhiDt_MOL4thOrd(Sborder, dPhiDt, time, 2.*dt/3., fr_as_crse, fr_as_fine, iteration);

    // RK3 stage 3
    //           S_new = 1/3 u^n + 2/3 u^2 + 2/3 dtL(u^2)

    //this makes S_new = 1/3 u^n + 2/3 u^2
    MultiFab::LinComb(S_new, 1./3., S_old, 0, 2./3., u2, 0, 0, NUM_STATE, 0);
    //this makes S_new = 1/3 u^n + 2/3  u^2 + 2/3 dtL(u^2)
    MultiFab::Saxpy(S_new, 2.*dt/3., dPhiDt, 0, 0, NUM_STATE, 0);

    /**/ 
    return dt;
  }
#endif
  Real
  AmrLevelAdv::
  advanceGodunov (Real time,
                  Real dt,
                  int  iteration,
                  int  ncycle)
  {
    BL_PROFILE("AmrLevelAdv::advance()");
        
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

    bool truncationErrorTest = false;
    ParmParse pp;
    pp.query("truncation_error_only", truncationErrorTest);
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
          uface[i].resize(amrex::grow(bxtmp, iteration), 1);
        }

        const Real prev_time = state[Phi_Type].prevTime();
        const Real cur_time = state[Phi_Type].curTime();
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
                          dx, dt, diffco);

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
          uface[i].resize(amrex::grow(bxtmp, iteration), 1);
        }

        const Real prev_time = state[Phi_Type].prevTime();
        const Real cur_time = state[Phi_Type].curTime();
        const Real ctr_time = 0.5*(prev_time + cur_time);

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
                            dx, dt, diffco);

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
          uface[i].resize(amrex::grow(bxtmp, iteration), 1);
        }

        const Real prev_time = state[Phi_Type].prevTime();
        const Real cur_time = state[Phi_Type].curTime();
        const Real ctr_time = 0.5*(prev_time + cur_time);

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
                            dx, dt, diffco);

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

    if (do_reflux && level < finest_level)
      reflux();

    if (level < finest_level)
      avgDown();

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

        state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                    BL_TO_FORTRAN_3D(S_new[mfi]),
                    &tagval, &clearval, 
                    ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                    ZFILL(dx), ZFILL(prob_lo), &time, &level);
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
