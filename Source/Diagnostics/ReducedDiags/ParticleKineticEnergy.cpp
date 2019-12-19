#include "ParticleKineticEnergy.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "AMReX_REAL.H"
#include <iostream>
#include <cmath>

/// constructor
ParticleKineticEnergy::ParticleKineticEnergy (
std::string rd_name, std::ofstream & ofs )
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

    // write header row
    ofs << "#";
    ofs << "step";
    for (int i = 0; i < species_number; ++i)
    {
        ofs << ",";
        ofs << species_names[i];
    }
    ofs << std::endl;
}
///< end constructor

/// destructor
ParticleKineticEnergy::~ParticleKineticEnergy ()
{}
///< end destructor

/// function that computes kinetic energy
void ParticleKineticEnergy::ComputeDiags (int step)
{

    /// Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    /// get WarpX class object
    auto & warpx = WarpX::GetInstance();

    /// get MultiParticleContainer class object
    auto & mypc = warpx.GetPartContainer();

    /// get number of species (int)
    auto species_number = mypc.nSpecies();

    /// resize data array
    m_data.resize(species_number,0.0);

    /// get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

    /// get number of level (int)
    auto level_number = warpx.finestLevel();

    /// speed of light squared
    auto c2 = PhysConst::c * PhysConst::c;

    /// loop over species
    for (int i_s = 0; i_s < species_number; ++i_s)
    {
        /// get WarpXParticleContainer class object
        auto & myspc = mypc.GetParticleContainer(i_s);

        /// get mass (amrex:Real)
        auto m = myspc.getMass();

        /// declare kinetic energy variable
        amrex::Real EK = 0.0;

        /// loop over refinement levels
        for (int lev = 0; lev <= level_number; ++lev)
        {

            #pragma omp parallel reduction(+:EK)
            /// Loop over boxes
            for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti)
            {

                /// get particle momentum arrays
                auto & px = pti.GetAttribs(PIdx::ux);
                auto & py = pti.GetAttribs(PIdx::uy);
                auto & pz = pti.GetAttribs(PIdx::uz);

                /// loop over particles
                for (long i = 0; i < px.size(); i++)
                {
                    /// get momentum squared
                    auto ps = (px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
                    /// get relativistic kinetic energy
                    EK += std::sqrt(ps*c2 + m*m*c2*c2);
                }
                ///< end loop over particles

            }
            ///< end loop over boxes

        }
        ///< end loop over refinement levels

        /// reduced sum for mpi ranks
        amrex::ParallelDescriptor::ReduceRealSum(EK);

        /// save EK to m_data
        m_data[i_s] = EK;

    }
    ///< end loop over species

}
///< end void ParticleKineticEnergy::ComputeDiags
