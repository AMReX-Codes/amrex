#ifndef AMREX_EB_STATE_REDIST_TRACKER_K_H_
#define AMREX_EB_STATE_REDIST_TRACKER_K_H_

// Map itracker neighbor indices to cell offsets.  Keep the encoding in one
// place because tracker construction and both redistribution paths consume it.
namespace amrex::state_redist_detail {

#if (AMREX_SPACEDIM == 2)
AMREX_GPU_CONSTANT int const imap[9] = {0,-1, 0, 1,-1, 1,-1, 0, 1};
AMREX_GPU_CONSTANT int const jmap[9] = {0,-1,-1,-1, 0, 0, 1, 1, 1};
AMREX_GPU_CONSTANT int const kmap[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
#elif (AMREX_SPACEDIM == 3)
AMREX_GPU_CONSTANT int const imap[27] = {0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
AMREX_GPU_CONSTANT int const jmap[27] = {0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
AMREX_GPU_CONSTANT int const kmap[27] = {0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
#endif

}

#endif
