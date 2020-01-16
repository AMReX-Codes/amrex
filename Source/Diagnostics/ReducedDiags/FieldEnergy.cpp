#include "FieldEnergy.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "AMReX_REAL.H"
#include "AMReX_ParticleReduce.H"
#include <iostream>
#include <cmath>

using namespace amrex;

/// constructor
FieldEnergy::FieldEnergy (std::string rd_name)
: ReducedDiags{rd_name}
{

    /// resize data array
    /// the extra one is for total energy
    m_data.resize(3,0.0);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
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
    }

}
///< end constructor

/// function that computes field energy
void FieldEnergy::ComputeDiags (int step)
{

    /// Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    /// get WarpX class object
    auto & warpx = WarpX::GetInstance();

    /// get MultiFab data at level 0
    const MultiFab & Ex = warpx.getEfield(0,0);
    const MultiFab & Ey = warpx.getEfield(0,1);
    const MultiFab & Ez = warpx.getEfield(0,2);
    const MultiFab & Bx = warpx.getBfield(0,0);
    const MultiFab & By = warpx.getBfield(0,1);
    const MultiFab & Bz = warpx.getBfield(0,2);

    /// get cell size
    Geometry const & geom = warpx.Geom(0);
    auto domain_box = geom.Domain();
#if (AMREX_SPACEDIM == 2)
    auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif (AMREX_SPACEDIM == 3)
    auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif

    /// RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "FieldEnergy reduced diagnostics does not work for RZ coordinate.");
#endif

    /// compute E squared
    Real tmpx = Ex.norm2(0,geom.periodicity());
    Real tmpy = Ey.norm2(0,geom.periodicity());
    Real tmpz = Ez.norm2(0,geom.periodicity());
    Real Es = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;

    /// compute B squared
    tmpx = Bx.norm2(0,geom.periodicity());
    tmpy = By.norm2(0,geom.periodicity());
    tmpz = Bz.norm2(0,geom.periodicity());
    Real Bs = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;

    /// save data
    m_data[1] = 0.5 * Es * PhysConst::ep0 * dV;
    m_data[2] = 0.5 * Bs / PhysConst::mu0 * dV;
    m_data[0] = m_data[1] + m_data[2];

}
///< end void FieldEnergy::ComputeDiags
