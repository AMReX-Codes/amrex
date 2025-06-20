#ifndef AMREX_ORDER_H_
#define AMREX_ORDER_H_
#include <AMReX_Config.H>

namespace amrex
{
    enum struct Order { C, F, RowMajor=C, ColumnMajor=F };
}
#endif

