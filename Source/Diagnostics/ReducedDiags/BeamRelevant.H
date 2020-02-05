/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_BEAMRELEVANT_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_BEAMRELEVANT_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class contains diagnostics that are relevant to beam.
 */
class BeamRelevant : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    BeamRelevant(std::string rd_name);

    /// name of beam species
    std::string m_beam_name;

    /** This funciton computes beam relevant quantites.
     *  \param [in] step current time step
     */
    virtual void ComputeDiags(int step) override final;

};

#endif
