#include <LaserProfiles.H>

#include <WarpX_Complex.H>

void
FieldFunctionLaserProfile::init (const amrex::ParmParse& pp)
{
    // Parse the properties of the Harris profile
    pp.get("profile_waist", profile_waist);
    pp.get("profile_duration", profile_duration);
    pp.get("profile_focal_distance", profile_focal_distance);
}
