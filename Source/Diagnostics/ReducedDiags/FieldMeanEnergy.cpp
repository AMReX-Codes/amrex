#include "FieldMeanEnergy.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "AMReX_REAL.H"
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

    /// loop over refinement levels
    for (int lev = 0; lev <= 0; ++lev)
    {
        auto & Ex = warpx.getEfield(lev,0);
        auto & Ey = warpx.getEfield(lev,1);
        auto & Ez = warpx.getEfield(lev,2);
        auto & Bx = warpx.getBfield(lev,0);
        auto & By = warpx.getBfield(lev,1);
        auto & Bz = warpx.getBfield(lev,2);
        for ( MFIter mfi(Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            //Ex(i,j,k)
        }
    }

    /// reduced sum for mpi ranks
    //ParallelDescriptor::ReduceRealSum();

}
///< end void ParticleMeanEnergy::ComputeDiags
