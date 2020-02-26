/* Copyright 2019 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_GODFREY_FILTER_H_
#define WARPX_GODFREY_FILTER_H_

#include "Filter.H"


enum class godfrey_coeff_set { Ex_Ey_Bz=0, Bx_By_Ez=1 };

/**
 * \brief Class for Godrey's filter to suppress Numerical Cherenkov Instability
 *
 * It derives from the base class Filter.
 * The filter stencil is initialized in method ComputeStencils. Computing the
 * stencil requires to read parameters from a table, where each lines stands
 * for a value of c*dt/dz.
 * The filter is applied using the base class' method ApplyStencil.
 */
class NCIGodfreyFilter : public Filter
{
public:

    NCIGodfreyFilter () = default;

    NCIGodfreyFilter(godfrey_coeff_set coeff_set, amrex::Real cdtodz,
                     bool nodal_gather);

    void ComputeStencils();

    static constexpr int m_stencil_width = 4;

private:

    // Set of coefficients (different fields require to read
    // different coefficients from the table)
    godfrey_coeff_set m_coeff_set;
    // The stencil depends on parameter c*dt/dz
    amrex::Real m_cdtodz;
    // Whether the gather is from nodal fields or staggered fields
    bool m_nodal_gather;
};

#endif // #ifndef WARPX_GODFREY_FILTER_H_
