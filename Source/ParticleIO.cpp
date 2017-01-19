
#include <ParticleContainer.H>
#include <WarpX.H>

void
MultiSpeciesContainer::Checkpoint (const std::string& dir, const std::string& name) const
{
    for (int i = 0; i < nspecies; ++i) {
	std::string namei = name + std::to_string(i);
	species[i]->Checkpoint(dir, namei);
    }
}

void
MultiSpeciesContainer::Restart (const std::string& dir, const std::string& name)
{
    for (int i = 0; i < nspecies; ++i) {
	std::string namei = name + std::to_string(i);
	species[i]->Restart(dir, namei);
    }
}

void
MultiSpeciesContainer::ReadHeader (std::istream& is) 
{
    for (int i = 0; i < nspecies; ++i) {
	species[i]->ReadHeader(is);
    }
}

void
MultiSpeciesContainer::WriteHeader (std::ostream& os) const
{
    for (int i = 0; i < nspecies; ++i) {
	species[i]->WriteHeader(os);
    }
}

void
SingleSpeciesContainer::ReadHeader (std::istream& is)
{
    is >> charge >> mass;
    WarpX::GotoNextLine(is);
}

void
SingleSpeciesContainer::WriteHeader (std::ostream& os) const
{
    // no need to write species_id
    os << charge << " " << mass << "\n";
}


