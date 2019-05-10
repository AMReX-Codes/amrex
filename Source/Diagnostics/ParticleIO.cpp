
#include <MultiParticleContainer.H>
#include <WarpX.H>

using namespace amrex;

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
            // real_names contains a list of all particle attributes.
            // pc->plot_flags is 1 or 0, whether quantity is dumped or not.
            pc->WritePlotFile(dir, species_names[i],
                              pc->plot_flags, int_flags,
                              real_names, int_names);
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

