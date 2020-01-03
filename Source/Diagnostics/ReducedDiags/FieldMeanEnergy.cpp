#include "FieldMeanEnergy.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "AMReX_REAL.H"
#include "AMReX_ParticleReduce.H"
#include <iostream>
#include <cmath>

using namespace amrex;

/// constructor
FieldMeanEnergy::FieldMeanEnergy (std::string rd_name)
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
FieldMeanEnergy::~FieldMeanEnergy ()
{}
///< end destructor

/// function that computes field energy
void FieldMeanEnergy::ComputeDiags (int step)
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

    Es += amrex::ReduceSum(
        Ex, 0,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab) -> Real
        { return xfab.dot(bx,0,1); }
        );
    Es += amrex::ReduceSum(
        Ey, 0,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab) -> Real
        { return xfab.dot(bx,0,1); }
        );
    Es += amrex::ReduceSum(
        Ez, 0,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab) -> Real
        { return xfab.dot(bx,0,1); }
        );
    Bs += amrex::ReduceSum(
        Bx, 0,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab) -> Real
        { return xfab.dot(bx,0,1); }
        );
    Bs += amrex::ReduceSum(
        By, 0,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab) -> Real
        { return xfab.dot(bx,0,1); }
        );
    Bs += amrex::ReduceSum(
        Bz, 0,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab) -> Real
        { return xfab.dot(bx,0,1); }
        );

    ParallelDescriptor::ReduceRealSum(Es);
    ParallelDescriptor::ReduceRealSum(Bs);

    m_data[1] = 0.5*Es*PhysConst::ep0;
    m_data[2] = 0.5*Bs/PhysConst::mu0;
    m_data[0] = m_data[1] + m_data[2];

}
///< end void FieldMeanEnergy::ComputeDiags
