#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDMEANENERGY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDMEANENERGY_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class mainly contains a function that
 *  computes the mean field energy on a grid.
*/
class FieldMeanEnergy : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    FieldMeanEnergy(std::string rd_name);

    /** This funciton computes the mean field energy on a grid (EF).
     *  EF = E eps / 2 + B / mu / 2,
     *  where E is the electric field,
     *  B is the magnetic field,
     *  eps is the vacuum permittivity,
     *  mu is the vacuum permeability. */
    virtual void ComputeDiags(int step) override final;

};
///< end class FieldMeanEnergy

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDEMEANENERGY_H_
