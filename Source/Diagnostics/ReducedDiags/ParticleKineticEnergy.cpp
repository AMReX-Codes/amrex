#include "ParticleKineticEnergy.H"
#include "WarpX.H"

#include <iostream>

ParticleKineticEnergy::ParticleKineticEnergy (std::string rd_name)
: ReducedDiags{rd_name}
{}

ParticleKineticEnergy::~ParticleKineticEnergy ()
{}

void ParticleKineticEnergy::ComputeDiags (int step)
{

    /** Judge if a reduced diags should be done */
    if ( (step+1)%m_freq != 0 ) { return; }

    /** obtain access of all data */
    WarpX& warpx = WarpX::GetInstance();

    /** obtain particle data */
    auto& mypc = warpx.GetPartContainer();

    std::vector<std::unique_ptr<WarpXParticleContainer>> allcontainers;

    /** loop over species */
    for (int i = 0; i < m_nspecies; ++i)
    {

        /** loop over refinement levels */
        for (int lev = 0; lev <= allcontainers[i]->finestLevel(); ++lev)
        {
            /** Loop over all grids/tiles at this level */
#ifdef _OPENMP
            #pragma omp parallel reduction(+:m_data)
#endif
            for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
            }
        }
        
    }

}
