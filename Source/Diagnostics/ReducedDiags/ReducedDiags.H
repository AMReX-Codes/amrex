/* Copyright 2019-2020 Maxence Thevenet, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_

#include "AMReX_REAL.H"
#include <string>
#include <vector>
#include <fstream>

/**
 *  Base class for reduced diagnostics. Each type of reduced diagnostics is
 *  implemented in a derived class, and must override the (pure virtual)
 *  function ComputeDiags.
 */
class ReducedDiags
{
public:

    /// output path (default)
    std::string m_path = "./diags/reducedfiles/";

    /// output extension (default)
    std::string m_extension = "txt";

    /// diags name
    std::string m_rd_name;

    /// output frequency
    int m_freq = 1;

    /// check if it is a restart run
    int m_IsNotRestart = 1;

    /// separator in the output file
    std::string m_sep = " ";

    /// output data
    std::vector<amrex::Real> m_data;

    /** constructor
     *  @param[in] rd_name reduced diags name */
    ReducedDiags(std::string rd_name);

    /// function to compute diags
    virtual void ComputeDiags(int step) = 0;

    /** write to file function
     *  @param[in] step time step */
    virtual void WriteToFile(int step) const;

};

#endif
