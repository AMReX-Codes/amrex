#include "DivEFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

DivEFunctor::DivEFunctor(const std::array<const amrex::MultiFab* const, 3> arr_mf_src, const int lev, const int ncomp)
    : ComputeDiagFunctor(ncomp), m_arr_mf_src(arr_mf_src), m_lev(lev)
{}

void
DivEFunctor::operator()(amrex::MultiFab& mf_dst, const int dcomp) const
{
    auto& warpx = WarpX::GetInstance();

    const BoxArray& ba = amrex::convert(warpx.boxArray(m_lev),IntVect::TheUnitVector());
    MultiFab divE( ba, warpx.DistributionMap(m_lev), 1, 0 );

    warpx.ComputeDivE(divE, m_lev);
    Average::ToCellCenter ( mf_dst, divE, dcomp, 0, 0, 1 );
}
