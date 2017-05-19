#ifndef WARPX_CRSEFINEBNDRY_H_
#define WARPX_CRSEFINEBNDRY_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

class WarpXCrseFineBndry
{
public:
    WarpXCrseFineBndry (const amrex::BoxArray& crse_ba,
                        const amrex::BoxArray& fine_ba,
                        const amrex::DistributionMapping& crse_dm,
                        const amrex::DistributionMapping& fine_dm,
                        const amrex::Geometry& crse_gm,
                        const amrex::Geometry& fine_gm,
                        const amrex::IntVect& ref_ratio,
                        int ngrow);

    amrex::Array<std::pair<amrex::Box,std::array<amrex::FArrayBox*,3> > >
    intersections (const amrex::Box& bx);    

    void AddCurrent (amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz) const;

    void setVal (amrex::Real v);

private:
    const amrex::Geometry m_fine_geom;
    amrex::MultiFab jx_crse_layout;
    amrex::MultiFab jy_crse_layout;
    amrex::MultiFab jz_crse_layout;
};

#endif
