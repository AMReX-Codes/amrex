/* Copyright 2020 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_SPECIESPHYSICALPROPERTIES_H_
#define WARPX_SPECIESPHYSICALPROPERTIES_H_

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_AmrCore.H>

#include <map>
#include <limits>

enum struct PhysicalSpecies{unspecified=0, electron, positron, photon, hydrogen};

namespace species
{
    AMREX_FORCE_INLINE
    PhysicalSpecies from_string(std::string species)
    {
        if( species=="unspecified" )
            return PhysicalSpecies::unspecified;
        if( species=="electron" )
            return PhysicalSpecies::electron;
        if( species=="positron" )
            return PhysicalSpecies::positron;
        if( species=="photon" )
            return PhysicalSpecies::photon;
        if( species=="hydrogen" )
            return PhysicalSpecies::hydrogen;
        amrex::Abort("unknown PhysicalSpecies");
        return PhysicalSpecies::unspecified;
    }

    AMREX_FORCE_INLINE
    amrex::Real get_charge (PhysicalSpecies ps)
    {
        switch(ps) {
        case PhysicalSpecies::unspecified:
            return std::numeric_limits<amrex::Real>::quiet_NaN();
        case PhysicalSpecies::electron:
            return -PhysConst::q_e;
        case PhysicalSpecies::positron:
            return PhysConst::q_e;
        case PhysicalSpecies::photon:
            return 0.;
        case PhysicalSpecies::hydrogen:
            return PhysConst::q_e;
        default:
            amrex::Abort("unknown PhysicalSpecies");
            return 0.;
        }
    }

    AMREX_FORCE_INLINE
    amrex::Real get_mass (PhysicalSpecies ps)
    {
        switch(ps) {
        case PhysicalSpecies::unspecified:
            return std::numeric_limits<amrex::Real>::quiet_NaN();
        case PhysicalSpecies::electron:
            return PhysConst::m_e;
        case PhysicalSpecies::positron:
            return PhysConst::m_e;
        case PhysicalSpecies::photon:
            return 0.;
        case PhysicalSpecies::hydrogen:
            return PhysConst::m_p;
        default:
            amrex::Abort("unknown PhysicalSpecies");
            return 0.;
        }
    }
}

#endif // WARPX_SPECIESPHYSICALPROPERTIES_H_
