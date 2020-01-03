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
    auto & Ex = warpx.getEfield(0,0);
    auto & Ey = warpx.getEfield(0,1);
    auto & Ez = warpx.getEfield(0,2);
    auto & Bx = warpx.getBfield(0,0);
    auto & By = warpx.getBfield(0,1);
    auto & Bz = warpx.getBfield(0,2);

    /// summed E and B squared
    Real Es = 0.0;
    Real Bs = 0.0;

    // summed total number of grids
    int  ng = 0;

    /// Ex
    /// loop over boxes
    for (MFIter mfi(Ex); mfi.isValid(); ++mfi)
    {
        /// This is the valid Box of the current FArrayBox.
        /// By "valid", we mean the original ungrown Box in BoxArray.
        auto box = mfi.validbox();
        /// Obtain Array4 from FArrayBox.
        auto arr = Ex.array(mfi);
        /// get indices
        auto lo = lbound(box);
        auto hi = ubound(box);
        /// loops over indices
        for (int i = lo.x; i <= hi.x; ++i) {
        for (int j = lo.y; j <  hi.y; ++j) {
        for (int k = lo.z; k <  hi.z; ++k) {
            Real F = 0.25 * ( arr(i,j,k) + arr(i,j+1,k) +
                              arr(i,j,k+1) + arr(i,j+1,k+1) );
            Es += F*F;
            ng++;
        }
        }
        }
    }

    /// Ey
    for (MFIter mfi(Ey); mfi.isValid(); ++mfi)
    {
        auto box = mfi.validbox();
        auto arr = Ey.array(mfi);
        auto lo = lbound(box);
        auto hi = ubound(box);
        for (int i = lo.x; i <  hi.x; ++i) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int k = lo.z; k <  hi.z; ++k) {
            Real F = 0.25 * ( arr(i,j,k) + arr(i+1,j,k) +
                              arr(i,j,k+1) + arr(i+1,j,k+1) );
            Es += F*F;
        }
        }
        }
    }

    /// Ez
    for (MFIter mfi(Ez); mfi.isValid(); ++mfi)
    {
        auto box = mfi.validbox();
        auto arr = Ez.array(mfi);
        auto lo = lbound(box);
        auto hi = ubound(box);
        for (int i = lo.x; i <  hi.x; ++i) {
        for (int j = lo.y; j <  hi.y; ++j) {
        for (int k = lo.z; k <= hi.z; ++k) {
            Real F = 0.25 * ( arr(i,j,k) + arr(i+1,j,k) +
                              arr(i,j+1,k) + arr(i+1,j+1,k) );
            Es += F*F;
        }
        }
        }
    }

    /// Bx
    for (MFIter mfi(Bx); mfi.isValid(); ++mfi)
    {
        auto box = mfi.validbox();
        auto arr = Bx.array(mfi);
        auto lo = lbound(box);
        auto hi = ubound(box);
        for (int i = lo.x; i <  hi.x; ++i) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int k = lo.z; k <= hi.z; ++k) {
            Real F = 0.5 * ( arr(i,j,k) + arr(i+1,j,k) );
            Bs += F*F;
        }
        }
        }
    }

    /// By
    for (MFIter mfi(By); mfi.isValid(); ++mfi)
    {
        auto box = mfi.validbox();
        auto arr = By.array(mfi);
        auto lo = lbound(box);
        auto hi = ubound(box);
        for (int i = lo.x; i <= hi.x; ++i) {
        for (int j = lo.y; j <  hi.y; ++j) {
        for (int k = lo.z; k <= hi.z; ++k) {
            Real F = 0.5 * ( arr(i,j,k) + arr(i,j+1,k) );
            Bs += F*F;
        }
        }
        }
    }

    /// Bz
    for (MFIter mfi(Bz); mfi.isValid(); ++mfi)
    {
        auto box = mfi.validbox();
        auto arr = Bz.array(mfi);
        auto lo = lbound(box);
        auto hi = ubound(box);
        for (int i = lo.x; i <= hi.x; ++i) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int k = lo.z; k <  hi.z; ++k) {
            Real F = 0.5 * ( arr(i,j,k) + arr(i,j,k+1) );
            Bs += F*F;
        }
        }
        }
    }

    /// sum over mpi ranks
    ParallelDescriptor::ReduceRealSum(Es);
    ParallelDescriptor::ReduceRealSum(Bs);
    ParallelDescriptor::ReduceIntSum(ng);

    /// save data for output
    m_data[1] = 0.5*Es*PhysConst::ep0/Real(ng);
    m_data[2] = 0.5*Bs/PhysConst::mu0/Real(ng);
    m_data[0] = m_data[1] + m_data[2];

}
///< end void FieldMeanEnergyNMT::ComputeDiags
