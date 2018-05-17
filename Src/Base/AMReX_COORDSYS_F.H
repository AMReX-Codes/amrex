#ifndef AMREX_COORDSYS_F_H_
#define AMREX_COORDSYS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
     void amrex_setvol (ARLIM_P(rlo), ARLIM_P(rhi),
		       amrex_real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       const amrex_real* offset,
                       const amrex_real* dx, const int* coord);

     void amrex_setvolpt (amrex_real* vol, ARLIM_P(vollo), ARLIM_P(volhi),
			 amrex_real* rolo, ARLIM_P(rolo), ARLIM_P(rohi),
			 amrex_real* rilo, ARLIM_P(rilo), ARLIM_P(rihi),
			 const amrex_real* dx, const int* coord);

     void amrex_setdloga (amrex_real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                         const amrex_real* offset,
                         const amrex_real* dx, const int* dir,
                         const int* coord);

     void amrex_setarea (ARLIM_P(rlo), ARLIM_P(rhi),
			amrex_real* area,  ARLIM_P(alo), ARLIM_P(ahi),
                        const amrex_real* offset,
                        const amrex_real* dx, const int* dir,
                        const int* coord);
#ifdef __cplusplus
}
#endif

#endif /*_COORDSYS_F_H_*/
