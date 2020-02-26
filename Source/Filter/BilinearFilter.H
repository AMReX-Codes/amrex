/* Copyright 2019 Andrew Myers, Maxence Thevenet, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_BILIN_FILTER_H_
#define WARPX_BILIN_FILTER_H_

#include "Filter.H"


class BilinearFilter : public Filter
{
public:

    BilinearFilter() = default;

    void ComputeStencils();

    amrex::IntVect npass_each_dir;

};
#endif // #ifndef WARPX_BILIN_FILTER_H_
