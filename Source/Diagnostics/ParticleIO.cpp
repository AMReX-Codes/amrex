
#include <MultiParticleContainer.H>
#include <WarpX.H>

using namespace amrex;

void
RigidInjectedParticleContainer::ReadHeader (std::istream& is)
{
    is >> charge >> mass;
    WarpX::GotoNextLine(is);

    int nlevs;
    is >> nlevs;
    WarpX::GotoNextLine(is);

    AMREX_ASSERT(zinject_plane_levels.size() == 0);
    AMREX_ASSERT(done_injecting.size() == 0);

    for (int i = 0; i < nlevs; ++i)
    {
        int zinject_plane_tmp;
        is >> zinject_plane_tmp;
        zinject_plane_levels.push_back(zinject_plane_tmp);
        WarpX::GotoNextLine(is);        
    } 

    for (int i = 0; i < nlevs; ++i)
    {
        int done_injecting_tmp;
        is >> done_injecting_tmp;
        done_injecting.push_back(done_injecting_tmp);
        WarpX::GotoNextLine(is);
    }     
}

void
RigidInjectedParticleContainer::WriteHeader (std::ostream& os) const
{
    // no need to write species_id
    os << charge << " " << mass << "\n";
    int nlevs = zinject_plane_levels.size();
    os << nlevs << "\n";
    for (int i = 0; i < nlevs; ++i)
    {
        os << zinject_plane_levels[i] << "\n";
    }
    for (int i = 0; i < nlevs; ++i)
    {
        os << done_injecting[i] << "\n";
    }
}

void
WarpXParticleContainer::ReadHeader (std::istream& is)
{
    is >> charge >> mass;
    WarpX::GotoNextLine(is);
}

void
WarpXParticleContainer::WriteHeader (std::ostream& os) const
{
    // no need to write species_id
    os << charge << " " << mass << "\n";
}

void
MultiParticleContainer::Checkpoint (const std::string& dir) const
{
    for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
	allcontainers[i]->Checkpoint(dir, species_names[i]);
    }
}

void
MultiParticleContainer::WritePlotFile (const std::string& dir,
                                       const Vector<std::string>& real_names) const
{
    Vector<std::string> int_names;    
    Vector<int> int_flags;
    
    for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
        auto& pc = allcontainers[i];
        if (pc->plot_species) {
            // Convert momentum to SI
            pc->ConvertUnits(ConvertDirection::WarpX_to_SI);
            // real_names contains a list of all particle attributes.
            // pc->plot_flags is 1 or 0, whether quantity is dumped or not.
            pc->WritePlotFile(dir, species_names[i],
                              pc->plot_flags, int_flags,
                              real_names, int_names);
            // Convert momentum back to WarpX units
            pc->ConvertUnits(ConvertDirection::SI_to_WarpX);
        }
    }
}

void
MultiParticleContainer::Restart (const std::string& dir)
{
    for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
	allcontainers[i]->Restart(dir, species_names[i]);
    }
}

void
MultiParticleContainer::ReadHeader (std::istream& is) 
{
    for (auto& pc : allcontainers) {
	pc->ReadHeader(is);
    }
}

void
MultiParticleContainer::WriteHeader (std::ostream& os) const
{
    for (const auto& pc : allcontainers) {
	pc->WriteHeader(os);
    }
}

// Particle momentum is defined as gamma*velocity, which is neither 
// SI mass*gamma*velocity nor normalized gamma*velocity/c.
// This converts momentum to SI units (or vice-versa) to write SI data
// to file.
void
PhysicalParticleContainer::ConvertUnits(ConvertDirection convert_direction)
{
    BL_PROFILE("PPC::ConvertUnits()");

    // Compute conversion factor
    Real factor = 1;
    if (convert_direction == ConvertDirection::WarpX_to_SI){
        factor = mass;
    } else if (convert_direction == ConvertDirection::SI_to_WarpX){
        factor = 1./mass;
    }

    for (int lev=0; lev<=finestLevel(); lev++){
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // - momenta are stored as a struct of array, in `attribs`
            auto& attribs = pti.GetAttribs();
            Real* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            Real* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            Real* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
            // Loop over the particles and convert momentum
            const long np = pti.numParticles();
            ParallelFor( np,
                [=] AMREX_GPU_DEVICE (long i) {
                    ux[i] *= factor;
                    uy[i] *= factor;
                    uz[i] *= factor;
                }
            );
        }
    }
}
