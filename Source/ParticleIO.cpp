
#include <ParticleContainer.H>
#include <WarpX.H>

void
MultiSpeciesContainer::Checkpoint (const std::string& dir, const std::string& name) const
{
    for (unsigned i = 0, n = allspecies.size(); i < n; ++i) {
	std::string namei = name + std::to_string(i);
	allspecies[i]->Checkpoint(dir, namei);
    }
}

void
MultiSpeciesContainer::Restart (const std::string& dir, const std::string& name)
{
    for (unsigned i = 0, n = allspecies.size(); i < n; ++i) {
	std::string namei = name + std::to_string(i);
	allspecies[i]->Restart(dir, namei);
    }
}

void
MultiSpeciesContainer::ReadHeader (std::istream& is) 
{
    for (auto& spec : allspecies) {
	spec->ReadHeader(is);
    }
}

void
MultiSpeciesContainer::WriteHeader (std::ostream& os) const
{
    for (const auto& spec : allspecies) {
	spec->WriteHeader(os);
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


