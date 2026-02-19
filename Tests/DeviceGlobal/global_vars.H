#pragma once

#include <AMReX_Gpu.H>
#include <AMReX_INT.H>

extern AMREX_DEVICE_GLOBAL_VARIABLE(amrex::Long,    dg_x);
extern AMREX_DEVICE_GLOBAL_VARIABLE(amrex::Long, 4, dg_y);
