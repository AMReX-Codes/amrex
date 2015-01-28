#include <winstd.H>

#include "ADR.H"
#include "ADR_F.H"

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

using std::string;

Real
ADR::advance (Real time,
              Real dt,
              int  iteration,
              int  ncycle)
{
    u_gdnv = new MultiFab[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	BoxArray edge_grids(grids);
	edge_grids.surroundingNodes(dir);
	u_gdnv[dir].define(edge_grids,1,1,Fab_allocate);
	u_gdnv[dir].setVal(1.e40);
    }
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
        //
        // Set reflux registers to zero.
        //
        getFluxReg(level+1).setVal(0.0);
    }

    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // It's possible for interpolation to create very small negative values for
    //   species so we make sure here that all species are non-negative after this point
    enforce_nonnegative_species(S_old);

    const Real prev_time = state[State_Type].prevTime();

#ifdef REACTIONS
    react_first_half_dt(S_old,time,dt);
#endif

    Real cur_time = state[State_Type].curTime();
        
        //
        // Get pointers to Flux registers, or set pointer to zero if not there.
        //
        FluxRegister *fine    = 0;
        FluxRegister *current = 0;
        
        if (do_reflux && level < finest_level)
            fine = &getFluxReg(level+1);
        if (do_reflux && level > 0)
            current = &getFluxReg(level);

        AmrLevel &levelData = *this;
        Geometry g = levelData.Geom();

        // Integrate, looping through the grids.
        FArrayBox divu;
        FArrayBox flux[BL_SPACEDIM];
        
        const Real *dx = geom.CellSize();

        MultiFab fluxes[BL_SPACEDIM];

        if (do_reflux && fine)
        {
            for (int j = 0; j < BL_SPACEDIM; j++)
            {
                BoxArray ba = S_new.boxArray();
                ba.surroundingNodes(j);
                fluxes[j].define(ba, NUM_STATE, 0, Fab_allocate);
            }
        }


       MultiFab ext_src_old(grids,NUM_STATE,1,Fab_allocate);
       ext_src_old.setVal(0.0);

       if (add_ext_src)
          getOldSource(prev_time,dt,ext_src_old);

#ifdef DIFFUSION
       MultiFab OldDiffTerm(grids,NUM_STATE,1);
       OldDiffTerm.setVal(0.);
       if (diffuse_spec == 1) {
          for (int ispec = 0; ispec < NumSpec; ispec++)
             add_diffusion_to_old_source(ext_src_old,OldDiffTerm,prev_time,FirstSpec+ispec);
       }
#endif
       ext_src_old.FillBoundary();
        
       for (FillPatchIterator fpi(*this, S_new, NUM_GROW,
                                   time, State_Type, 0, NUM_STATE);
             fpi.isValid();
             ++fpi)
        {
            int mfiindex = fpi.index();

            Box bx(fpi.UngrownBox());

            // Create FAB for extended grid values (including boundaries) and fill.
            FArrayBox &state = fpi();
            FArrayBox &stateout = S_new[fpi];
            
            // Allocate fabs for fluxes.
            for (int i = 0; i < BL_SPACEDIM ; i++)  
                flux[i].resize(BoxLib::surroundingNodes(bx,i),NUM_STATE);

            const int*  domain_lo = geom.Domain().loVect();
            const int*  domain_hi = geom.Domain().hiVect();

            Real cflLoc = -1.0e+20;
            int is_finest_level = 0;
            if (level == finest_level) is_finest_level = 1;
            BL_FORT_PROC_CALL(ADVECT,advect)
                (&time,
                 bx.loVect(), bx.hiVect(),
                 domain_lo, domain_hi,
                 BL_TO_FORTRAN(state), BL_TO_FORTRAN(stateout),
		 BL_TO_FORTRAN(u_gdnv[0][fpi]),
		 BL_TO_FORTRAN(u_gdnv[1][fpi]),
#if (BL_SPACEDIM == 3)
		 BL_TO_FORTRAN(u_gdnv[2][fpi]),
#endif
                 BL_TO_FORTRAN(ext_src_old[fpi]),
                 D_DECL(BL_TO_FORTRAN(flux[0]), 
                        BL_TO_FORTRAN(flux[1]), 
                        BL_TO_FORTRAN(flux[2])), 
                 dx, &dt, verbose);

            if (do_reflux)
            {
                if (fine)
                {
                    for (int i = 0; i < BL_SPACEDIM ; i++)
                        fluxes[i][fpi].copy(flux[i]);
                }
                if (current)
                {
                    for (int i = 0; i < BL_SPACEDIM ; i++)
                        current->FineAdd(flux[i],i,mfiindex,0,0,NUM_STATE,1);
                }
            }
        }

        if (do_reflux && fine)
        {
            for (int i = 0; i < BL_SPACEDIM ; i++)
                fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1);
        }

        if (add_ext_src)  
        {
           getOldSource(prev_time,dt,ext_src_old);
           ext_src_old.mult(-0.5*dt);

           // Compute source at new time (no ghost cells needed)
           MultiFab ext_src_new(grids,NUM_STATE,0,Fab_allocate);
           ext_src_new.setVal(0.0);

           getNewSource(prev_time,cur_time,dt,ext_src_new);
           ext_src_new.mult(0.5*dt);

           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,ext_src_old,0,0,S_new.nComp(),0);
           MultiFab::Add(S_new,ext_src_new,0,0,S_new.nComp(),0);
        }

#ifdef DIFFUSION
     if (diffuse_spec == 1) {
        for (int ispec = 0; ispec < NumSpec; ispec++)
           time_center_diffusion(S_new, OldDiffTerm, cur_time, dt, FirstSpec+ispec);
     }
#endif

#ifdef REACTIONS
    react_second_half_dt(S_new,cur_time,dt);
#endif

    // Copy old velocity into new velocity -- assumes velocity unchanging.
    MultiFab::Copy(S_new,S_old,Xvel,Xvel,BL_SPACEDIM,0);

    delete [] u_gdnv;

    return dt;
}
