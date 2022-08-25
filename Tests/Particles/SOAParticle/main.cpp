#include <AMReX.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_ParticleTile.H>

using namespace amrex;

template <typename T_PC>
void addParticles ()
{
    T_PC pc;
    using ParticleType = typename T_PC::ParticleType;

    const int add_num_particles = 5;

    auto& ptile1 = pc.DefineAndReturnParticleTile(0, 0, 0);

    for (int i = 0; i < add_num_particles; ++i)
    {
        ptile1.pos(i, 0) = 12.0;
        ptile1.pos(i, 1) = 12.0;
        ptile1.pos(i, 2) = 12.0;

        // TODO
        ptile1.id(i) = 1;
        ptile1.cpu(i) = 1;
    }
    //ptile1.push_back_int(3, ...std::vector);
    //ptile1.push_back_int(4, ...std::vector);

    pc.Redistribute();
}


int main(int argc, char* argv[])
 {
    amrex::Initialize(argc,argv);
    {
        //addParticles< ParticleContainer<1,2,3,4> > ();
        addParticles< ParticleContainerPureSoA<3,4> > ();
    }
    amrex::Finalize();
 }

 


