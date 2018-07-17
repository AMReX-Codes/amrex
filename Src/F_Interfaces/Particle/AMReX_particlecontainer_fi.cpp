
#include <AMReX_AmrParticles.H>
#include <AMReX_AmrCore.H>

using namespace amrex;

#define AMREX_FI_NSTRUCTREAL BL_SPACEDIM
#define AMREX_FI_NSTRUCTINT 0

namespace {
    using FParticleContainer = AmrParticleContainer<AMREX_FI_NSTRUCTREAL,
                                                    AMREX_FI_NSTRUCTINT>;
}

extern "C" {

    void amrex_fi_new_particlecontainer (FParticleContainer*& particlecontainer,
                                         AmrCore* amrcore)
    {
	particlecontainer = new FParticleContainer(amrcore);
    }

    void amrex_fi_delete_particlecontainer (FParticleContainer* particlecontainer)
    {
	delete particlecontainer;
    }

    void amrex_fi_init_particles_one_per_cell (FParticleContainer* particlecontainer)
    {
        FParticleContainer::ParticleInitData pdata = {AMREX_D_DECL(0.0, 0.0, 0.0)};
	particlecontainer->InitOnePerCell(0.5, 0.5, 0.5, pdata);
    }

    void amrex_fi_write_particles(FParticleContainer* particlecontainer,
                                  const char* dirname, const char* pname, int is_checkpoint)
    {
        particlecontainer->Checkpoint(dirname, pname, is_checkpoint);
    }

    void amrex_fi_particle_redistribute (FParticleContainer* particlecontainer,
                                         int lev_min, int lev_max, int ng)
    {
	particlecontainer->Redistribute(lev_min, lev_max, ng);
    }

    void amrex_fi_get_particles(FParticleContainer* particlecontainer,
                                int lev, MFIter* mfi, Real*& dp, long& np)
    {
        const int grid = mfi->index();
        const int tile = mfi->LocalTileIndex();
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
            if (np > 0) {
                auto& aos = particle_tile.GetArrayOfStructs();
                dp = aos.data();
            } else {
                dp = nullptr;
            }            
        } else {
            np = 0;
            dp = nullptr;
        }
    }

    void amrex_fi_num_particles(FParticleContainer* particlecontainer,
                                int lev, MFIter* mfi, long& np)
    {
        const int grid = mfi->index();
        const int tile = mfi->LocalTileIndex();
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {            
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
        } else {
            np = 0;
        }
    }
}
