#ifndef AMREX_COORDSYS_F_H_
#define AMREX_COORDSYS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
     void amrex_setvol (AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
		       amrex_real* vol, AMREX_ARLIM_P(vlo), AMREX_ARLIM_P(vhi),
                       const amrex_real* offset,
                       const amrex_real* dx, const int* coord);

     void amrex_setvolpt (amrex_real* vol, AMREX_ARLIM_P(vollo), AMREX_ARLIM_P(volhi),
			 amrex_real* rolo, AMREX_ARLIM_P(rolo), AMREX_ARLIM_P(rohi),
			 amrex_real* rilo, AMREX_ARLIM_P(rilo), AMREX_ARLIM_P(rihi),
			 const amrex_real* dx, const int* coord);

     void amrex_setdloga (amrex_real* vol, AMREX_ARLIM_P(vlo), AMREX_ARLIM_P(vhi),
                         const amrex_real* offset,
                         const amrex_real* dx, const int* dir,
                         const int* coord);

     void amrex_setarea (AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
			amrex_real* area,  AMREX_ARLIM_P(alo), AMREX_ARLIM_P(ahi),
                        const amrex_real* offset,
                        const amrex_real* dx, const int* dir,
                        const int* coord);
#ifdef __cplusplus
}
#endif

#endif /*_COORDSYS_F_H_*/
