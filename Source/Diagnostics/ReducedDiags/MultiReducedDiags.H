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

    /** number of reduced diagnostics */
    int m_n_diags;

    /** names of reduced diagnostics */
    std::vector<std::string> m_reduced_diags_names;

    /** */
    std::vector<std::unique_ptr<ReducedDiags>> m_multi_rd;

    /** current iteration */
    int m_step;

    /** constructor */
    MultiReducedDiags(int step);

    /** destructor */
    ~MultiReducedDiags();

    void ComputeDiags();

    void WriteToFile();
};

#endif //WARPX_DIAGNOSTICS_REDUCEDDIAGS_MULTIREDUCEDDIAGS_H_
