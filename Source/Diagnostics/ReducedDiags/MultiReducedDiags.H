#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_

#include "ReducedDiags.H"
#include "ParticleMeanEnergy.H"
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

    /// flag
    int m_plot_rd = 0;

    /// output path (default)
    std::string m_path = "./diags/reducedfiles/";

    /// names of reduced diagnostics
    std::vector<std::string> m_rd_names;

    /// create ReducedDiags class object
    std::vector<std::unique_ptr<ReducedDiags>> m_multi_rd;

    /// constructor
    MultiReducedDiags();

    /// destructor
    ~MultiReducedDiags();

    /** Loop over all ReducedDiags and call their ComputeDiags
     *  @param[in] step current iteration time */
    void ComputeDiags(int step);

    /** Loop over all ReducedDiags and call their WriteToFile
     *  @param[in] step current iteration time */
    void WriteToFile(int step);

};
///< end class MultiReducedDiags

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_
