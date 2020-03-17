/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEHISTOGRAM_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEHISTOGRAM_H_

#include "ReducedDiags.H"
#include "WarpX.H"
#include <fstream>

/**
 * Reduced diagnostics that computes a histogram over particles
 * for a quantity specified by the user in the input file using the parser.
 */
class ParticleHistogram : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    ParticleHistogram(std::string rd_name);

    /// normalization type
    int m_norm;

    /// number of bins
    int m_bin_num;

    /// selected species index
    int m_selected_species_id = -1;

    /// max and min bin values
    amrex::Real m_bin_max;
    amrex::Real m_bin_min;

    /// bin size
    amrex::Real m_bin_size;

    /// Parser to read expression for particle quantity from the input file.
    /// 7 elements are t, x, y, z, ux, uy, uz
    static constexpr int m_nvars = 7;
    std::unique_ptr<ParserWrapper<m_nvars>> m_parser;

    /** This function computes a histogram of user defined quantity.
     *  \param [in] step current time step.
     */
    virtual void ComputeDiags(int step) override final;

};

#endif
