#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_

#include <string>

class ReducedDiags
{
public:

    /** output frequency of the reduced diagnostics */
    int m_freq;

    /** number of species */
    int m_nspecies;

    /** constructor */
    ReducedDiags (std::string rd_name);

    /** destructor */
    ~ReducedDiags ();

    virtual void ComputeDiags(int step) = 0;
};

#endif //WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_
