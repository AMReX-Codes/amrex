#include "FieldMeanEnergyNMT.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "AMReX_REAL.H"
#include "AMReX_ParticleReduce.H"
#include <iostream>
#include <cmath>

using namespace amrex;

/// constructor
FieldMeanEnergyNMT::FieldMeanEnergyNMT (std::string rd_name)
: ReducedDiags{rd_name}
{

    /// open file
    std::ofstream ofs;
    ofs.open(m_path + m_rd_name + ".txt",
        std::ofstream::out | std::ofstream::app);

    // write header row
    ofs << "#";
    ofs << "step";
    ofs << m_sep;
    ofs << "time(s)";
    ofs << m_sep;
    ofs << "total(J)";
    ofs << m_sep;
    ofs << "E(J)";
    ofs << m_sep;
    ofs << "B(J)";
    ofs << std::endl;

    /// close file
    ofs.close();

}
///< end constructor

/// destructor
FieldMeanEnergyNMT::~FieldMeanEnergyNMT ()
{}
///< end destructor

/// function that computes field energy
void FieldMeanEnergyNMT::ComputeDiags (int step)
{

    /// Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    /// get WarpX class object
    auto & warpx = WarpX::GetInstance();

    /// resize data array
    /// the extra one is for total energy
    m_data.resize(3,0.0);

    /// get MultiFab data at level 0
    const MultiFab & Ex = warpx.getEfield(0,0);
    const MultiFab & Ey = warpx.getEfield(0,1);
    const MultiFab & Ez = warpx.getEfield(0,2);
    const MultiFab & Bx = warpx.getBfield(0,0);
    const MultiFab & By = warpx.getBfield(0,1);
    const MultiFab & Bz = warpx.getBfield(0,2);

    Real Es = 0.0;
    Real Bs = 0.0;

    for (MFIter mfi(Ex); mfi.isValid(); ++mfi)
    {
        const Box & box = mfi.validbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        const FArrayBox & fab = Ex[mfi];
        auto arr = fab.array();
        for (int i = lo.x; i <= hi.x; ++i) {
        for (int j = lo.y; j <  hi.y; ++j) {
        for (int k = lo.z; k <  hi.z; ++k) {
            Real F = 0.25 * ( arr(i,j,k) + arr(i,j+1,k) + arr(i,j,k+1) + arr(i,j+1,k+1) );
            Es += F*F;
        }
        }
        }
    }

    for (MFIter mfi(Ey); mfi.isValid(); ++mfi)
    {
        const Box & box = mfi.validbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        const FArrayBox & fab = Ey[mfi];
        auto arr = fab.array();
        for (int i = lo.x; i <  hi.x; ++i) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int k = lo.z; k <  hi.z; ++k) {
            Real F = 0.25 * ( arr(i,j,k) + arr(i+1,j,k) + arr(i,j,k+1) + arr(i+1,j,k+1) );
            Es += F*F;
        }
        }
        }
    }

    for (MFIter mfi(Ez); mfi.isValid(); ++mfi)
    {
        const Box & box = mfi.validbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        const FArrayBox & fab = Ez[mfi];
        auto arr = fab.array();
        for (int i = lo.x; i <  hi.x; ++i) {
        for (int j = lo.y; j <  hi.y; ++j) {
        for (int k = lo.z; k <= hi.z; ++k) {
            Real F = 0.25 * ( arr(i,j,k) + arr(i+1,j,k) + arr(i,j+1,k) + arr(i+1,j+1,k) );
            Es += F*F;
        }
        }
        }
    }

    for (MFIter mfi(Bx); mfi.isValid(); ++mfi)
    {
        const Box & box = mfi.validbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        const FArrayBox & fab = Bx[mfi];
        auto arr = fab.array();
        for (int i = lo.x; i <  hi.x; ++i) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int k = lo.z; k <= hi.z; ++k) {
            Real F = 0.5 * ( arr(i,j,k) + arr(i+1,j,k) );
            Bs += F*F;
        }
        }
        }
    }

    for (MFIter mfi(By); mfi.isValid(); ++mfi)
    {
        const Box & box = mfi.validbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        const FArrayBox & fab = By[mfi];
        auto arr = fab.array();
        for (int i = lo.x; i <= hi.x; ++i) {
        for (int j = lo.y; j <  hi.y; ++j) {
        for (int k = lo.z; k <= hi.z; ++k) {
            Real F = 0.5 * ( arr(i,j,k) + arr(i,j+1,k) );
            Bs += F*F;
        }
        }
        }
    }

    for (MFIter mfi(Bz); mfi.isValid(); ++mfi)
    {
        const Box & box = mfi.validbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        const FArrayBox & fab = Bz[mfi];
        auto arr = fab.array();
        for (int i = lo.x; i <= hi.x; ++i) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int k = lo.z; k <  hi.z; ++k) {
            Real F = 0.5 * ( arr(i,j,k) + arr(i,j,k+1) );
            Bs += F*F;
        }
        }
        }
    }

    ParallelDescriptor::ReduceRealSum(Es);
    ParallelDescriptor::ReduceRealSum(Bs);

    m_data[1] = 0.5*Es*PhysConst::ep0;
    m_data[2] = 0.5*Bs/PhysConst::mu0;
    m_data[0] = m_data[1] + m_data[2];

}
///< end void FieldMeanEnergyNMT::ComputeDiags
