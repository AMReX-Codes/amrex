#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEKINETICENERGY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEKINETICENERGY_H_

#include "ReducedDiags.H"

class ParticleKineticEnergy : public ReducedDiags
{
public:

    /** constructor */
    ParticleKineticEnergy();

    /** destructor */
    ~ParticleKineticEnergy();

    void ComputeDiags();

};

#endif //WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEKINETICENERGY_H_
