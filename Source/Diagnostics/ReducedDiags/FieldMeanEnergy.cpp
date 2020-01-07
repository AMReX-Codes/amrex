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

    /// get number of grids used
    Geometry const & geom = warpx.Geom(0);
    auto domain_box = geom.Domain();
    long nx = domain_box.length(0);
    long ny = domain_box.length(1);
    long nz = domain_box.length(2);
    long nx_extra = geom.isPeriodic(0) ? 0 : 1;
    long ny_extra = geom.isPeriodic(1) ? 0 : 1;
    long nz_extra = geom.isPeriodic(2) ? 0 : 1;
    long ex_num_points = nx*(ny+ny_extra)*(nz+nz_extra);
    long ey_num_points = (nx+nx_extra)*ny*(nz+nz_extra);
    long ez_num_points = (nx+nx_extra)*(ny+ny_extra)*nz;
    long bx_num_points = (nx+nx_extra)*ny*nz;
    long by_num_points = nx*(ny+ny_extra)*nz;
    long bz_num_points = nx*ny*(nz+nz_extra);

    /// compute E squared
    Real tmpx = Ex.norm2(0,geom.periodicity());
    Real tmpy = Ey.norm2(0,geom.periodicity());
    Real tmpz = Ez.norm2(0,geom.periodicity());
    Real Es = tmpx*tmpx/Real(ex_num_points) +
              tmpy*tmpy/Real(ey_num_points) +
              tmpz*tmpz/Real(ez_num_points);

    /// compute B squared
    tmpx = Bx.norm2(0,geom.periodicity());
    tmpy = By.norm2(0,geom.periodicity());
    tmpz = Bz.norm2(0,geom.periodicity());
    Real Bs = tmpx*tmpx/Real(bx_num_points) +
              tmpy*tmpy/Real(by_num_points) +
              tmpz*tmpz/Real(bz_num_points);

    /// save data
    m_data[1] = 0.5*Es*PhysConst::ep0;
    m_data[2] = 0.5*Bs/PhysConst::mu0;
    m_data[0] = m_data[1] + m_data[2];

}
///< end void FieldMeanEnergy::ComputeDiags
