#include <AMReX_Config.H>

#ifdef AMREX_CRAYPAT
#include <pat_api.h>
#elif defined (AMREX_VTUNE)
#include <ittnotify.h>
#elif defined (AMREX_FORGE)
#include <mapsampler_api.h>
#endif
