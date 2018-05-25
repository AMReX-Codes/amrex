
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
}
