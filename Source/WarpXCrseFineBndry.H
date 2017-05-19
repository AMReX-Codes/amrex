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

    amrex::MultiFab& GetJx () { return jx_crse_layout; }
    amrex::MultiFab& GetJy () { return jy_crse_layout; }
    amrex::MultiFab& GetJz () { return jz_crse_layout; }

private:
    amrex::BoxArray ba_fine_layout;
    amrex::DistributionMapping dm_fine_layout;

    amrex::BoxArray ba_crse_layout;
    amrex::DistributionMapping dm_crse_layout;
    amrex::MultiFab jx_crse_layout;
    amrex::MultiFab jy_crse_layout;
    amrex::MultiFab jz_crse_layout;
};

#endif
