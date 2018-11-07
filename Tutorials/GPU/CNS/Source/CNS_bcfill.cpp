
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

using namespace amrex;

void cns_bcfill_single (Box const& bx, FArrayBox& data,
                        const int dcomp, const int numcomp,
                        Geometry const& geom, const Real time,
                        const Vector<BCRec>& bcr, const int bcomp,
                        const int scomp)
{
    AMREX_ALWAYS_ASSERT(numcomp == 1);
    amrex::Abort("cns_bcfill_single");
}

void cns_bcfill_group (Box const& bx, FArrayBox& data,
                       const int dcomp, const int numcomp,
                       Geometry const& geom, const Real time,
                       const Vector<BCRec>& bcr, const int bcomp,
                       const int scomp)
{
    amrex::Abort("cns_bcfill_group");
}
