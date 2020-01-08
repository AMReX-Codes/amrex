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

    /// output path (default)
    std::string m_path = "./diags/reducedfiles/";

    /// diags name
    std::string m_rd_name;

    /// output frequency
    int m_freq = 1;

    /// check if it is a restart run
    int m_IsNotRestart = 1;

    /// separator in the output file
    std::string m_sep = ",";

    /// output data
    std::vector<amrex::Real> m_data;

    /** constructor
     *  @param[in] rd_name reduced diags names */
    ReducedDiags(std::string rd_name);

    /// function to compute diags
    virtual void ComputeDiags(int step) = 0;

    /** write to file function
     *  @param[in] step time step */
    virtual void WriteToFile(int step) const;

};
///< end class ReducedDiags

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_
