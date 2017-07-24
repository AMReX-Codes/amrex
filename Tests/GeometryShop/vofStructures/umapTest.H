#ifndef _umapTest_H
#define _umapTest_H

#ifndef BL_LANG_FORT
typedef std::array<amrex::Real,BL_SPACEDIM>  RealDIM;

struct FaceData
{
    amrex::Real m_aperature;
    RealDIM m_centroid;
};

struct EBBndryData
{
    RealDIM m_normal;
    RealDIM m_bndry_centroid;
    amrex::Real m_value;
    amrex::Real m_bndry_area;
    amrex::Real m_vol_frac;
};

#endif

#endif
