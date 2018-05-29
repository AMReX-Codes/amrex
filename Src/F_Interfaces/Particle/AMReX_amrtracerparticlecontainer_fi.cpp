
#include <AMReX_AmrParticles.H>
#include <AMReX_AmrCore.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_amrtracerparticlecontainer (AmrTracerParticleContainer*& amrtracerparticles,
                                                  AmrCore* amrcore)
    {
	amrtracerparticles = new AmrTracerParticleContainer(amrcore);
    }

    void amrex_fi_delete_amrtracerparticlecontainer (AmrTracerParticleContainer* amrtracerparticlecontainer)
    {
	delete amrtracerparticlecontainer;
    }

    void amrex_fi_init_particles_one_per_cell (AmrTracerParticleContainer* amrtracerparticlecontainer)
    {
        AmrTracerParticleContainer::ParticleInitData pdata = {AMREX_D_DECL(0.0, 0.0, 0.0)};
	amrtracerparticlecontainer->InitOnePerCell(0.5, 0.5, 0.5, pdata);
    }

    void amrex_fi_write_particles(AmrTracerParticleContainer* amrtracerparticlecontainer,
                                  const char* dirname, const char* pname, bool is_checkpoint)
    {
        amrtracerparticlecontainer->Checkpoint(dirname, pname, is_checkpoint);
    }

    void amrex_fi_particle_redistribute (AmrTracerParticleContainer* amrtracerparticlecontainer,
                                         int lev_min, int lev_max, int ng)
    {
	amrtracerparticlecontainer->Redistribute(lev_min, lev_max, ng);
    }

    void amrex_fi_get_particles(AmrTracerParticleContainer* amrtracerparticlecontainer,
                                int lev, MFIter* mfi, Real*& dp)
    {
        const int grid = mfi->index();
        const int tile = mfi->LocalTileIndex();
        auto& particles = amrtracerparticlecontainer->ParticlesAt(lev, grid, tile);
        auto& aos = particles.GetArrayOfStructs();
        dp = aos.data();
    }
}
