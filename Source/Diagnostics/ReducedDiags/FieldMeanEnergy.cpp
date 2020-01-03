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

    //auto & Ex = warpx.getEfield(0,0);
    //auto & Ey = warpx.getEfield(0,1);
    //auto & Ez = warpx.getEfield(0,2);
    //auto & Bx = warpx.getBfield(0,0);
    //auto & By = warpx.getBfield(0,1);
    //auto & Bz = warpx.getBfield(0,2);

    //Real Exs = ReduceSum( Ex,
    //[=] AMREX_GPU_HOST_DEVICE () -> Real
    //{
    //    return Ex*Ex;
    //});

    const MultiFab & Ex = warpx.getEfield(0,0);
    const MultiFab & Ey = warpx.getEfield(0,1);
    const MultiFab & Ez = warpx.getEfield(0,2);
    const MultiFab & Bx = warpx.getBfield(0,0);
    const MultiFab & By = warpx.getBfield(0,1);
    const MultiFab & Bz = warpx.getBfield(0,2);

    const int nghost = 0;
    const int xcomp = 0;
    const int ycomp = 0;
    const int numcomp = 1;
    Real Es = 0.0;
    Real Bs = 0.0;

//    Real sm = 0;
//    sm += amrex::ReduceSum(
//        Ex, nghost,
//        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab) -> Real
//        { return xfab.dot(bx,0,0,1); }
//        );

    Es += amrex::ReduceSum(
        Ex, Ex, nghost,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab, FArrayBox const& yfab) -> Real
        { return xfab.dot(bx,xcomp,yfab,bx,ycomp,numcomp); }
        );
    Es += amrex::ReduceSum(
        Ey, Ey, nghost,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab, FArrayBox const& yfab) -> Real
        { return xfab.dot(bx,xcomp,yfab,bx,ycomp,numcomp); }
        );
    Es += amrex::ReduceSum(
        Ez, Ez, nghost,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab, FArrayBox const& yfab) -> Real
        { return xfab.dot(bx,xcomp,yfab,bx,ycomp,numcomp); }
        );
    Bs += amrex::ReduceSum(
        Bx, Bx, nghost,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab, FArrayBox const& yfab) -> Real
        { return xfab.dot(bx,xcomp,yfab,bx,ycomp,numcomp); }
        );
    Bs += amrex::ReduceSum(
        By, By, nghost,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab, FArrayBox const& yfab) -> Real
        { return xfab.dot(bx,xcomp,yfab,bx,ycomp,numcomp); }
        );
    Bs += amrex::ReduceSum(
        Bz, Bz, nghost,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& xfab, FArrayBox const& yfab) -> Real
        { return xfab.dot(bx,xcomp,yfab,bx,ycomp,numcomp); }
        );
    //if (!local) ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
//    ParallelAllReduce::Sum(Es, ParallelContext::CommunicatorSub());
//    ParallelAllReduce::Sum(Bs, ParallelContext::CommunicatorSub());
    ParallelDescriptor::ReduceRealSum(Es);
    ParallelDescriptor::ReduceRealSum(Bs);

    m_data[1] = 0.5*Es*PhysConst::ep0;
    m_data[2] = 0.5*Bs/PhysConst::mu0;
    m_data[0] = m_data[1] + m_data[2];

    /// reduced sum for mpi ranks
    //ParallelDescriptor::ReduceRealSum();

}
///< end void FieldMeanEnergy::ComputeDiags
