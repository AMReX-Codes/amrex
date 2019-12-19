#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_

#include "ReducedDiags.H"
#include "ParticleKineticEnergy.H"
#include <vector>
#include <string>
#include <memory>

class MultiReducedDiags
{
public:

    /** names of reduced diagnostics */
    std::vector<std::string> m_rd_names;

    /** */
    std::vector<std::unique_ptr<ReducedDiags>> m_multi_rd;

    /** constructor */
    MultiReducedDiags();

    /** destructor */
    ~MultiReducedDiags();

    void ComputeDiags(int step);

    void WriteToFile(int step);
};

#endif //WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_
