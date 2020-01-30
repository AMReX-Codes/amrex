/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDENERGY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDENERGY_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class mainly contains a function that
 *  computes the field energy.
 */
class FieldEnergy : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    FieldEnergy(std::string rd_name);

    /** This funciton computes the field energy (EF).
     *  EF = E eps / 2 + B / mu / 2,
     *  where E is the electric field,
     *  B is the magnetic field,
     *  eps is the vacuum permittivity,
     *  mu is the vacuum permeability. */
    virtual void ComputeDiags(int step) override final;

};

#endif
