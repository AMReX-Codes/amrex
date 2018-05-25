
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
        AmrTracerParticleContainer::ParticleInitData pdata = {1.0};
	amrtracerparticlecontainer->InitOnePerCell(0.5, 0.5, 0.5, pdata);
    }
}
