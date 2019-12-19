#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEKINETICENERGY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEKINETICENERGY_H_

#include "ReducedDiags.H"

class ParticleKineticEnergy : public ReducedDiags
{
public:

    /** constructor */
    ParticleKineticEnergy(std::string rd_name);

    /** destructor */
    ~ParticleKineticEnergy();

    virtual void ComputeDiags(int step) final;

};

#endif //WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEKINETICENERGY_H_
