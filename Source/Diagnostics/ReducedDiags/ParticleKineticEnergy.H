#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEKINETICENERGY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEKINETICENERGY_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class mainly contains a function that
 *  computes the particle relativistic kinetic energy
 *  of each species.
*/
class ParticleKineticEnergy : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    ParticleKineticEnergy(std::string rd_name, std::ofstream & ofs);

    /// destructor
    ~ParticleKineticEnergy();

    /** This funciton computes the particle relativistic
     *  kinetic energy (EK).
     *  EK = sqrt( p^2 c^2 + m^2 c^4 ),
     *  where p is the relativistic momentum,
     *  m is the particle rest mass. */
    virtual void ComputeDiags(int step) final;

};
///< end class ParticleKineticEnergy

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEKINETICENERGY_H_
