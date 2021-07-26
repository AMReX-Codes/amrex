#include "probdata.h"

void get_face_velocity(const int* level, const amrex::real* time,
                       AMREX_D_DECL(BL_FORT_FAB_ARG(xvel),
                                    BL_FORT_FAB_ARG(yvel),
                                    BL_FORT_FAB_ARG(zvel)),
                       const amrex::real* dx, const amrex::real* problo)
{
    vx = adv_vel[0];
    vy = adv_vel[1];

    return;
}
