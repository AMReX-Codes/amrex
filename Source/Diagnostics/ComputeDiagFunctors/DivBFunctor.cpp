#include "DivBFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

DivBFunctor::DivBFunctor(const std::array<const amrex::MultiFab* const, 3> arr_mf_src, const int lev, const amrex::IntVect crse_ratio, const int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_arr_mf_src(arr_mf_src), m_lev(lev)
{}

void
DivBFunctor::operator()(amrex::MultiFab& mf_dst, int dcomp) const
{
    auto& warpx = WarpX::GetInstance();
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_avg, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // A cell-centered divB multifab spanning the entire domain is generated
    // and divB is computed on the cell-center, with ng=1.
    MultiFab divB( warpx.boxArray(m_lev), warpx.DistributionMap(m_lev), 1, ng );
    warpx.ComputeDivB(divB, 0, m_arr_mf_src, WarpX::CellSize(m_lev) );
    // Coarsen and Interpolate from divB to coarsened/reduced_domain mf_dst
    Average::CoarsenAndInterpolate( mf_dst, divB, dcomp, 0, nComp(), 0, m_crse_ratio);
}
