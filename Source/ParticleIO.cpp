
#include <ParticleContainer.H>
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
MultiParticleContainer::Checkpoint (const std::string& dir, 
				    bool is_checkpoint,
                                    const Vector<std::string>& varnames) const
{
    for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
	allcontainers[i]->Checkpoint(dir, species_names[i], is_checkpoint, varnames);
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

