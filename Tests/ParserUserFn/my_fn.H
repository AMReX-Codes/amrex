#ifndef MY_FN_H_
#define MY_FN_H_

#include <AMReX_GpuQualifiers.H>

double f1_h (double);
AMREX_GPU_DEVICE double f1_d (double);
AMREX_GPU_HOST_DEVICE double f1_hd (double);

double f2_h (double, double);
AMREX_GPU_DEVICE double f2_d (double, double);
AMREX_GPU_HOST_DEVICE double f2_hd (double, double);

double f3_h (double, double, double);
AMREX_GPU_DEVICE double f3_d (double, double, double);
AMREX_GPU_HOST_DEVICE double f3_hd (double, double, double);

double f4_h (double, double, double, double);
AMREX_GPU_DEVICE double f4_d (double, double, double, double);
AMREX_GPU_HOST_DEVICE double f4_hd (double, double, double, double);

#endif
