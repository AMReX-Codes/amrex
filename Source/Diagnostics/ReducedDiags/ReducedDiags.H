#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_

#include "AMReX_REAL.H"
#include <string>
#include <vector>
#include <fstream>

/// This class contains basic information for all diagnostics.
class ReducedDiags
{
public:

    /// output frequency
    int m_freq;

    /// output data
    std::vector<amrex::Real> m_data;

    /** constructor
     *  @param[in] rd_name reduced diags names */
    ReducedDiags(std::string rd_name);

    /// destructor
    ~ReducedDiags();

    /// function to compute diags
    virtual void ComputeDiags(int step) = 0;

    /** write to file function
     *  @param[in] step time step
     *  @param[in] file to write */
    virtual void WriteToFile(int step, std::ofstream & ofs);

};
///< end class ReducedDiags

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_
