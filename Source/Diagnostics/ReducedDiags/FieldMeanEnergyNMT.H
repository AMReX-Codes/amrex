#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDMEANENERGYNMT_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDMEANENERGYNMT_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class mainly contains a function that
 *  computes the mean field energy on a grid.
*/
class FieldMeanEnergyNMT : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    FieldMeanEnergyNMT(std::string rd_name);

    /// destructor
    ~FieldMeanEnergyNMT();

    /** This funciton computes the mean field energy on a grid (EF).
     *  EF = E eps / 2 + B / mu / 2,
     *  where E is the electric field,
     *  B is the magnetic field,
     *  eps is the vacuum permittivity,
     *  mu is the vacuum permeability.
     *  NMT represents No Multi-Threading. */
    virtual void ComputeDiags(int step) override final;

};
///< end class FieldMeanEnergyNMT

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDEMEANENERGYNMT_H_
