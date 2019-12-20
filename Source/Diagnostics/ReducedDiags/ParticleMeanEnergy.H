#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEMEANENERGY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEMEANENERGY_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class mainly contains a function that
 *  computes the mean particle relativistic kinetic energy
 *  of each species.
*/
class ParticleMeanEnergy : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    ParticleMeanEnergy(std::string rd_name);

    /// destructor
    ~ParticleMeanEnergy();

    /** This funciton computes the mean particle relativistic
     *  kinetic energy (EK).
     *  EK = sqrt( p^2 c^2 + m^2 c^4 ) - m c^2,
     *  where p is the relativistic momentum,
     *  m is the particle rest mass. */
    virtual void ComputeDiags(int step) override final;

};
///< end class ParticleMeanEnergy

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEMEANENERGY_H_
