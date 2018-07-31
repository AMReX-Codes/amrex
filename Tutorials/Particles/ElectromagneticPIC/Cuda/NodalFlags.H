#ifndef NODAL_FLAGS_H_
#define NODAL_FLAGS_H_

#include <AMReX_IntVect.H>

class YeeGrid {
public:
    static amrex::IntVect Bx_nodal_flag;
    static amrex::IntVect By_nodal_flag;
    static amrex::IntVect Bz_nodal_flag;

    static amrex::IntVect Ex_nodal_flag;
    static amrex::IntVect Ey_nodal_flag;
    static amrex::IntVect Ez_nodal_flag;

    static amrex::IntVect jx_nodal_flag;
    static amrex::IntVect jy_nodal_flag;
    static amrex::IntVect jz_nodal_flag;
};

#endif
