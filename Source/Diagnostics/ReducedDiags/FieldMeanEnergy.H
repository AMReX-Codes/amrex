#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDMEANENERGY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDMEANENERGY_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class mainly contains a function that
 *  computes the mean field (E+B) energy on a grid.
*/
class FieldMeanEnergy : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names
     *  @param[in] ofs file to write */
    FieldMeanEnergy(std::string rd_name, std::ofstream & ofs);

    /// destructor
    ~FieldMeanEnergy();

    /** This funciton computes the mean field energy (EF) on a grid.
     *  EF = eps E^2 / 2 + B^2 / mu / 2,
     *  where E is the electric field,
     *  B is the magnetic field,
     *  eps is the vacuum permittivity,
     *  mu is the vacuum permeability. */
    virtual void ComputeDiags(int step) final;

};
///< end class FieldMeanEnergy

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDMEANENERGY_H_
