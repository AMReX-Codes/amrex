#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_

#include "ReducedDiags.H"
#include "ParticleKineticEnergy.H"
#include <vector>
#include <string>
#include <memory>

/// This class addresses all diagnostics
class MultiReducedDiags
{
public:

    /// names of reduced diagnostics
    std::vector<std::string> m_rd_names;

    /// create ReducedDiags class object
    std::vector<std::unique_ptr<ReducedDiags>> m_multi_rd;

    /// constructor
    MultiReducedDiags();

    /// destructor
    ~MultiReducedDiags();

    /** function to call diags computation
     *  @param[in] step current iteration time */
    void ComputeDiags(int step);

    /** function to write diags
     *  @param[in] step current iteration time */
    void WriteToFile(int step);

};
///< end class MultiReducedDiags

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_
