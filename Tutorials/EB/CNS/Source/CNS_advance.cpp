
#include <CNS.H>

using namespace amrex;

Real
CNS::advance (Real time, Real dt, int iteration, int ncycle)
{
    for (int k = 0; k < NUM_STATEDATA_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab dSdt(grids,dmap,NUM_STATE,0);
    MultiFab Sborder(grids,dmap,NUM_STATE,NUM_GROW);
  
    // RK2 stage 1
    FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);
    compute_dSdt(Sborder, dSdt, dt);
    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt, 0, 0, NUM_STATE, 0);
    
    // RK2 stage 2
    // After fillpatch Sborder = U^n+0.5*dt*dUdt^n
    FillPatch(*this, Sborder, NUM_GROW, time+0.5*dt, State_Type, 0, NUM_STATE);
    compute_dSdt(Sborder, dSdt, 0.5*dt);
    // U^{n+1} = (U^n+0.5*dt*dUdt^n) + 0.5*dt*dUdt^*
    MultiFab::LinComb(S_new, 1.0, Sborder, 0, 0.5*dt, dSdt, 0, 0, NUM_STATE, 0);
    
    return dt;
}

void
CNS::compute_dSdt (MultiFab& S, MultiFab& dSdt, Real dt)
{
    MultiFab xflux(amrex::convert(grids,IntVect(1,0,0)), dmap, NUM_STATE, 0);
    MultiFab yflux(amrex::convert(grids,IntVect(0,1,0)), dmap, NUM_STATE, 0);
    MultiFab zflux(amrex::convert(grids,IntVect(0,0,1)), dmap, NUM_STATE, 0);
    
    {
        std::array<FArrayBox,3> flux;
        for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            for (int idim = 0; idim < 3; ++dim) {
                const Box& bxtmp = amrex::surroundingNodes(bx,i);
                flux[i].resize(bxtmp,NUM_STATE);
            }

            
        }
    }
}
