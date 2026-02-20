#ifndef AMREX_ORDER_H_
#define AMREX_ORDER_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Order.H
 * \brief Memory-layout order tag for multi-dimensional arrays.
 */

namespace amrex
{
/**
 * \ingroup amrex_arithmetic_types
 * \brief Memory layout order for multi-dimensional arrays.
 *
 * \c Order::C (row-major) stores the last index fastest; \c Order::F
 * (column-major / Fortran order) stores the first index fastest.
 */
    enum struct Order { C, F, RowMajor=C, ColumnMajor=F };
}
#endif

