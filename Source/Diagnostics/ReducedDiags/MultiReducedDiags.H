/* Copyright 2019-2020 Maxence Thevenet, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_

#include "ReducedDiags.H"
#include <vector>
#include <string>
#include <memory>

/**
 * This class holds multiple instances of ReducedDiagnostics, and contains
 * general functions to initialize, compute, and write these diagnostics
 * to file.
 */
class MultiReducedDiags
{
public:

    /// Bool: whether or not reduced diagnostics are activated
    int m_plot_rd = 0;

    /// names of reduced diagnostics
    std::vector<std::string> m_rd_names;

    /// m_multi_rd stores a pointer to each reduced diagnostics
    std::vector<std::unique_ptr<ReducedDiags>> m_multi_rd;

    /// constructor
    MultiReducedDiags();

    /** Loop over all ReducedDiags and call their ComputeDiags
     *  @param[in] step current iteration time */
    void ComputeDiags(int step);

    /** Loop over all ReducedDiags and call their WriteToFile
     *  @param[in] step current iteration time */
    void WriteToFile(int step);

};

#endif
