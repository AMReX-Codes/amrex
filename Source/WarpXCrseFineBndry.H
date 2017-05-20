#ifndef WARPX_CRSEFINEBNDRY_H_
#define WARPX_CRSEFINEBNDRY_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

class WarpXBndryForFine
{
public:
    WarpXBndryForFine (const amrex::BoxArray& crse_ba,
                       const amrex::BoxArray& fine_ba,
                       const amrex::DistributionMapping& crse_dm,
                       const amrex::DistributionMapping& fine_dm,
                       const amrex::Geometry& crse_gm,
                       const amrex::Geometry& fine_gm,
                       const amrex::IntVect& ref_ratio,
                       int ngrow);

    bool intersects (const amrex::Box& bx) const;

    void addFrom (amrex::FArrayBox& jx_in, amrex::FArrayBox& jy_in, amrex::FArrayBox& jz_in);

    void addTo (amrex::MultiFab& jx_out, amrex::MultiFab& jy_out, amrex::MultiFab& jz_out) const;

    void setVal (amrex::Real v);

private:
    const amrex::Geometry m_fine_geom;
    const amrex::IntVect m_ref_ratio;
    const int m_ngrow;
    amrex::MultiFab jx;
    amrex::MultiFab jy;
    amrex::MultiFab jz;
};

class WarpXBndryForCrse
{
public:
    WarpXBndryForCrse (const amrex::BoxArray& crse_ba,
                       const amrex::BoxArray& fine_ba,
                       const amrex::DistributionMapping& crse_dm,
                       const amrex::DistributionMapping& fine_dm,
                       const amrex::Geometry& crse_gm,
                       const amrex::Geometry& fine_gm,
                       const amrex::IntVect& ref_ratio,
                       int ngrow);

    bool intersects (const amrex::Box& bx) const;

    void addFrom (amrex::FArrayBox& jx_in, amrex::FArrayBox& jy_in, amrex::FArrayBox& jz_in);

    void addTo (amrex::MultiFab& jx_out, amrex::MultiFab& jy_out, amrex::MultiFab& jz_out) const;

    void setVal (amrex::Real v);

private:
    const amrex::Geometry m_crse_geom;
    const amrex::IntVect m_ref_ratio;
    const int m_ngrow;
    amrex::MultiFab jx;
    amrex::MultiFab jy;
    amrex::MultiFab jz;
};

#endif
