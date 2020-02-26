/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FieldEnergy.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_ParticleReduce.H>

#include <iostream>
#include <cmath>


using namespace amrex;

// constructor
FieldEnergy::FieldEnergy (std::string rd_name)
: ReducedDiags{rd_name}
{

    // RZ coordinate is not working
    #if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "FieldEnergy reduced diagnostics does not work for RZ coordinate.");
    #endif

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // read number of levels
    int nLevel = 0;
    ParmParse pp("amr");
    pp.query("max_level", nLevel);
    nLevel += 1;

    // resize data array
    m_data.resize(3*nLevel,0.0);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs;
            ofs.open(m_path + m_rd_name + "." + m_extension,
                std::ofstream::out | std::ofstream::app);
            // write header row
            ofs << "#";
            ofs << "[1]step()";
            ofs << m_sep;
            ofs << "[2]time(s)";
            for (int lev = 0; lev < nLevel; ++lev)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(3+3*lev) + "]";
                ofs << "total_lev"+std::to_string(lev)+"(J)";
                ofs << m_sep;
                ofs << "[" + std::to_string(4+3*lev) + "]";
                ofs << "E_lev"+std::to_string(lev)+"(J)";
                ofs << m_sep;
                ofs << "[" + std::to_string(5+3*lev) + "]";
                ofs << "B_lev"+std::to_string(lev)+"(J)";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }

}
// end constructor

// function that computes field energy
void FieldEnergy::ComputeDiags (int step)
{

    // Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get number of level
    auto nLevel = warpx.finestLevel() + 1;

    // loop over refinement levels
    for (int lev = 0; lev < nLevel; ++lev)
    {

        // get MultiFab data at lev
        const MultiFab & Ex = warpx.getEfield(lev,0);
        const MultiFab & Ey = warpx.getEfield(lev,1);
        const MultiFab & Ez = warpx.getEfield(lev,2);
        const MultiFab & Bx = warpx.getBfield(lev,0);
        const MultiFab & By = warpx.getBfield(lev,1);
        const MultiFab & Bz = warpx.getBfield(lev,2);

        // get cell size
        Geometry const & geom = warpx.Geom(lev);
        auto domain_box = geom.Domain();
        #if (AMREX_SPACEDIM == 2)
        auto dV = geom.CellSize(0) * geom.CellSize(1);
        #elif (AMREX_SPACEDIM == 3)
        auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
        #endif

        // compute E squared
        Real tmpx = Ex.norm2(0,geom.periodicity());
        Real tmpy = Ey.norm2(0,geom.periodicity());
        Real tmpz = Ez.norm2(0,geom.periodicity());
        Real Es = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;

        // compute B squared
        tmpx = Bx.norm2(0,geom.periodicity());
        tmpy = By.norm2(0,geom.periodicity());
        tmpz = Bz.norm2(0,geom.periodicity());
        Real Bs = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;

        // save data
        m_data[lev*3+1] = 0.5 * Es * PhysConst::ep0 * dV;
        m_data[lev*3+2] = 0.5 * Bs / PhysConst::mu0 * dV;
        m_data[lev*3+0] = m_data[lev*3+1] + m_data[lev*3+2];

    }
    // end loop over refinement levels

    /* m_data now contains up-to-date values for:
     *  [total field energy at level 0,
     *   electric field energy at level 0,
     *   magnetic field energy at level 0,
     *   total field energy at level 1,
     *   electric field energy at level 1,
     *   magnetic field energy at level 1,
     *   ......] */

}
// end void FieldEnergy::ComputeDiags
