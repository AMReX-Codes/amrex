/* Copyright 2019-2020 Axel Huebl, Ligia Diana Amorim, Maxence Thevenet
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "InjectorDensity.H"
#include "PlasmaInjector.H"


using namespace amrex;

InjectorDensity::~InjectorDensity ()
{
    switch (type)
    {
    case Type::parser:
    {
        object.parser.m_parser.clear();
        break;
    }
    case Type::custom:
    {
        object.custom.clear();
        break;
    }
    case Type::predefined:
    {
        object.predefined.clear();
        break;
    }
    default:
        return;
    }
}

InjectorDensityPredefined::InjectorDensityPredefined (
    std::string const& a_species_name) noexcept
    : profile(Profile::null)
{
    ParmParse pp(a_species_name);

    std::vector<amrex::Real> v;
    // Read parameters for the predefined plasma profile,
    // and store them in managed memory
    pp.getarr("predefined_profile_params", v);
    p = static_cast<amrex::Real*>
        (amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real)*v.size()));
    for (int i = 0; i < static_cast<int>(v.size()); ++i) {
        p[i] = v[i];
    }

    // Parse predefined profile name, and update member variable profile.
    std::string which_profile_s;
    pp.query("predefined_profile_name", which_profile_s);
    std::transform(which_profile_s.begin(), which_profile_s.end(),
                   which_profile_s.begin(), ::tolower);
    if (which_profile_s == "parabolic_channel"){
        profile = Profile::parabolic_channel;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(v.size() > 5,
            "InjectorDensityPredefined::parabolic_channel: not enough parameters");
    }
}

// Note that we are not allowed to have non-trivial destructor.
// So we rely on clear() to free memory.
void InjectorDensityPredefined::clear ()
{
    amrex::The_Managed_Arena()->free(p);
}
