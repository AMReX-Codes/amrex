#include "ParticleMeanEnergy.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "AMReX_REAL.H"
#include <iostream>
#include <cmath>

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

    /// get WarpX class object
    auto & warpx = WarpX::GetInstance();

    /// get MultiParticleContainer class object
    auto & mypc = warpx.GetPartContainer();

    /// get number of species (int)
    auto species_number = mypc.nSpecies();

    /// resize data array
    /// the extra one is for total energy
    m_data.resize(species_number+1,0.0);

    /// get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

    /// get number of level (int)
    auto level_number = warpx.finestLevel();

    /// speed of light squared
    auto c2 = PhysConst::c * PhysConst::c;

    /// declare total kinetic energy variable
    amrex::Real EKtot = 0.0;

    /// total number of particles;
    long nptot = 0;

    /// loop over species
    for (int i_s = 0; i_s < species_number; ++i_s)
    {
        /// get WarpXParticleContainer class object
        auto & myspc = mypc.GetParticleContainer(i_s);

        /// get mass (amrex:Real)
        auto m = myspc.getMass();

        /// declare kinetic energy variable
        amrex::Real EK = 0.0;

        /// number of particles;
        long np = 0;

        /// loop over refinement levels
        for (int lev = 0; lev <= level_number; ++lev)
        {

            #pragma omp parallel reduction(+:EK, np)
            /// Loop over boxes
            for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti)
            {

                /// get particle momentum arrays
                auto & px = pti.GetAttribs(PIdx::ux);
                auto & py = pti.GetAttribs(PIdx::uy);
                auto & pz = pti.GetAttribs(PIdx::uz);

                /// sum total number of particles
                np += pti.numParticles();

                /// loop over particles
                for (long i = 0; i < px.size(); i++)
                {
                    /// get momentum squared
                    auto ps = (px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
                    /// get relativistic kinetic energy
                    EK += std::sqrt(ps*c2 + m*m*c2*c2) - m*c2;
                }
                ///< end loop over particles

            }
            ///< end loop over boxes

        }
        ///< end loop over refinement levels

        /// reduced sum for mpi ranks
        amrex::ParallelDescriptor::ReduceRealSum(EK);
        amrex::ParallelDescriptor::ReduceLongSum(np);

        /// compute total EK
        EKtot += EK;

        /// compute total np
        nptot += np;

        /// save EK to m_data
        if ( np > 0 )
        { m_data[i_s+1] = EK / np; }
        else
        { m_data[i_s+1] = 0.0; }

    }
    ///< end loop over species

    /// save total EK
    if ( nptot > 0 )
    { m_data[0] = EKtot / nptot; }
    else
    { m_data[0] = 0.0; }

}
///< end void ParticleMeanEnergy::ComputeDiags
