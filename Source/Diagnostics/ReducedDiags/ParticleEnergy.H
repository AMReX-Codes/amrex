#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEENERGY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEENERGY_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class mainly contains a function that
 *  computes the particle relativistic kinetic energy
 *  of each species.
*/
class ParticleEnergy : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    ParticleEnergy(std::string rd_name);

    /** This funciton computes the particle relativistic
     *  kinetic energy (EP).
     *  EP = sqrt( p^2 c^2 + m^2 c^4 ) - m c^2,
     *  where p is the relativistic momentum,
     *  m is the particle rest mass. */
    virtual void ComputeDiags(int step) override final;

};
///< end class ParticleEnergy

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEENERGY_H_
