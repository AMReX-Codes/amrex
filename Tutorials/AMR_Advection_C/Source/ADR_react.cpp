#include <winstd.H>

#include "ADR.H"
#include "ADR_F.H"

using std::string;

#ifdef REACTIONS
void
ADR::react_first_half_dt(MultiFab& S_old, Real time, Real dt) 
{
    if (do_react == 1) strang_chem(S_old,time,dt);
}

void
ADR::react_second_half_dt(MultiFab& S_new, Real cur_time, Real dt) 
{
    if (do_react == 1) strang_chem(S_new,cur_time,dt);
}

void
ADR::strang_chem (MultiFab&  state,
                  Real       time,
                  Real       dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::strang_chem(MultiFab&,...");
    const Real strt_time = ParallelDescriptor::second();

    for (MFIter Smfi(state); Smfi.isValid(); ++Smfi)
    {
        FArrayBox& fb   = state[Smfi];
        const Box& bx   = Smfi.validbox();
        reactState(fb, fb, bx, time, 0.5*dt);
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

//      if (ParallelDescriptor::IOProcessor()) 
//          std::cout << "strang_chem time = " << run_time << '\n';
    }
}

void
ADR::reactState(FArrayBox&        Snew,
                FArrayBox&        Sold,
                const Box&        box,
                Real              time,
                Real              dt_react)
{
    BL_FORT_PROC_CALL(REACT_STATE,react_state)
                     (box.loVect(), box.hiVect(), 
                     BL_TO_FORTRAN(Sold),
                     BL_TO_FORTRAN(Snew),
                     time,dt_react);
}
#endif
