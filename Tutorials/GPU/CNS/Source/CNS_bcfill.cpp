
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

using namespace amrex;

// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the desciptor set up in CNS::variableSetUp.

void cns_bcfill_single (Box const& bx, FArrayBox& data,
                        const int dcomp, const int numcomp,
                        Geometry const& geom, const Real time,
                        const Vector<BCRec>& bcr, const int bcomp,
                        const int scomp)
{
    AMREX_ALWAYS_ASSERT(numcomp == 1);
    amrex::Abort("cns_bcfill_single TODO");
}

void cns_bcfill_group (Box const& bx, FArrayBox& data,
                       const int dcomp, const int numcomp,
                       Geometry const& geom, const Real time,
                       const Vector<BCRec>& bcr, const int bcomp,
                       const int scomp)
{
    FArrayBox* fab = &data;
    const auto geomdata = geom.data();
    const BCRec* bp = bcr.data();
    Gpu::AsyncArray<BCRec> bcr_aa(bp+bcomp, numcomp);
    BCRec* bcr_p = bcr_aa.data();
    AMREX_LAUNCH_DEVICE_LAMBDA (bx, tbx,
    {
//        amrex_fill_bc_cc(tbx, *fab, dcomp, numcomp, geomdata, time,
//                         bcr_aa, scomp);
    });
}
