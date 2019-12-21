#include "ParticleMeanEnergy.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "AMReX_REAL.H"
#include "AMReX_ParticleReduce.H"
#include <iostream>
#include <cmath>

using namespace amrex;

/// constructor
ParticleMeanEnergy::ParticleMeanEnergy (std::string rd_name)
: ReducedDiags{rd_name}
{
    /// get WarpX class object
    auto & warpx = WarpX::GetInstance();

    /// get MultiParticleContainer class object
    auto & mypc = warpx.GetPartContainer();

    /// get number of species (int)
    auto species_number = mypc.nSpecies();

    /// get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

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
    for (int i = 0; i < species_number; ++i)
    {
        ofs << m_sep;
        ofs << species_names[i]+"(J)";
    }
    ofs << std::endl;

    /// close file
    ofs.close();

}
///< end constructor

/// destructor
ParticleMeanEnergy::~ParticleMeanEnergy ()
{}
///< end destructor

/// function that computes kinetic energy
void ParticleMeanEnergy::ComputeDiags (int step)
{

    /// Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    /// get MultiParticleContainer class object
    auto & mypc = WarpX::GetInstance().GetPartContainer();

    /// get number of species (int)
    auto species_number = mypc.nSpecies();

    /// resize data array
    /// the extra one is for total energy
    m_data.resize(species_number+1,0.0);

    /// get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

    /// speed of light squared
    auto c2 = PhysConst::c * PhysConst::c;

    /// loop over species
    for (int i_s = 0; i_s < species_number; ++i_s)
    {
        /// get WarpXParticleContainer class object
        auto & myspc = mypc.GetParticleContainer(i_s);

        /// get mass (Real)
        auto m = myspc.getMass();

        using PType = typename WarpXParticleContainer::SuperParticleType;

        auto Etot = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto w  = p.rdata(PIdx::w);
            auto px = p.rdata(PIdx::ux);
            auto py = p.rdata(PIdx::uy);
            auto pz = p.rdata(PIdx::uz);
            auto ps = (px*px + py*py + pz*pz);
            return ( std::sqrt(ps*c2 + m*m*c2*c2) - m*c2 ) * w;
        });

        auto Wtot = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            return p.rdata(PIdx::w);
        });

        /// reduced sum for mpi ranks
        ParallelDescriptor::ReduceRealSum(Etot);
        ParallelDescriptor::ReduceRealSum(Wtot);

        /// save Etot to m_data
        if ( Wtot > 0.0 )
        { m_data[i_s+1] = Etot / Wtot; }
        else
        { m_data[i_s+1] = 0.0; }

    }
    ///< end loop over species

    Real E_sum = 0.0;
    /// loop over species
    for (int i_s = 0; i_s < species_number; ++i_s)
    {
        E_sum += m_data[i_s+1];
    }
    ///< end loop over species

    /// save total energy
    if ( species_number > 0 )
    { m_data[0] = E_sum / species_number; }
    else
    { m_data[0] = 0.0; }

}
///< end void ParticleMeanEnergy::ComputeDiags
