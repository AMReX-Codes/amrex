/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleEnergy.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_ParticleReduce.H>

#include <iostream>
#include <cmath>
#include <limits>


using namespace amrex;

// constructor
ParticleEnergy::ParticleEnergy (std::string rd_name)
: ReducedDiags{rd_name}
{
    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get MultiParticleContainer class object
    auto & mypc = warpx.GetPartContainer();

    // get number of species (int)
    auto nSpecies = mypc.nSpecies();

    // resize data array
    m_data.resize(2*nSpecies+2,0.0);

    // get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

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
            ofs << m_sep;
            ofs << "[3]total(J)";
            for (int i = 0; i < nSpecies; ++i)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(4+i) + "]";
                ofs << species_names[i]+"(J)";
            }
            ofs << m_sep;
            ofs << "[" + std::to_string(4+nSpecies) + "]";
            ofs << "total_mean(J)";
            for (int i = 0; i < nSpecies; ++i)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(5+nSpecies+i) + "]";
                ofs << species_names[i]+"_mean(J)";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }

}
// end constructor

// function that computes kinetic energy
void ParticleEnergy::ComputeDiags (int step)
{

    // Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    // get MultiParticleContainer class object
    auto & mypc = WarpX::GetInstance().GetPartContainer();

    // get number of species (int)
    auto nSpecies = mypc.nSpecies();

    // get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

    // speed of light squared
    auto c2 = PhysConst::c * PhysConst::c;

    // loop over species
    for (int i_s = 0; i_s < nSpecies; ++i_s)
    {
        // get WarpXParticleContainer class object
        auto & myspc = mypc.GetParticleContainer(i_s);

        // get mass (Real)
        auto m = myspc.getMass();

        using PType = typename WarpXParticleContainer::SuperParticleType;

        // Use amex::ReduceSum to compute the sum of energies of all particles
        // held by the current MPI rank, for this species. This involves a loop over all
        // boxes held by this MPI rank.
        auto Etot = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto w  = p.rdata(PIdx::w);
            auto ux = p.rdata(PIdx::ux);
            auto uy = p.rdata(PIdx::uy);
            auto uz = p.rdata(PIdx::uz);
            auto us = (ux*ux + uy*uy + uz*uz);
            return ( std::sqrt(us*c2 + c2*c2) - c2 ) * m * w;
        });

        // Same thing for the particles weights.
        auto Wtot = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            return p.rdata(PIdx::w);
        });

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum
            (Etot, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum
            (Wtot, ParallelDescriptor::IOProcessorNumber());

        // save results for this species i_s into m_data
        m_data[i_s+1] = Etot;
        if ( Wtot > std::numeric_limits<Real>::min() )
        { m_data[nSpecies+2+i_s] = Etot / Wtot; }
        else
        { m_data[nSpecies+2+i_s] = 0.0; }

    }
    // end loop over species

    // save total energy
    // loop over species
    m_data[0] = 0.0;          // total energy
    m_data[nSpecies+1] = 0.0; // total mean energy
    for (int i_s = 0; i_s < nSpecies; ++i_s)
    {
        m_data[0] += m_data[i_s+1];
        m_data[nSpecies+1] += m_data[nSpecies+2+i_s];
    }
    // end loop over species

    /* m_data now contains up-to-date values for:
     *  [total energy (all species),
     *   total energy (species 1),
     *   ...,
     *   total energy (species n),
     *   mean energy (all species),
     *   mean energy (species 1),
     *   ...,
     *   mean energy (species n)] */

}
// end void ParticleEnergy::ComputeDiags
