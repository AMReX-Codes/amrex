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
    //{
    //    ParticleType p(...);
    //    p.pos(0) = 12.0;
    //    p.pos(1) = 12.0;
    //    p.pos(2) = 12.0;
    //    ptile1.push_back(p);
    //}

    //DefineAndReturnParticleTile(0,0,0);

    for (int i = 0; i < add_num_particles; ++i)
    {
        ptile1.pos(i, 0) = 12.0;
        ptile1.pos(i, 1) = 12.0;
        ptile1.pos(i, 2) = 12.0;

        // TODO
        ptile1.id(i) = 1;
        //ptile1.cpu(i) = 1;
    }
    //ptile1.push_back_int(3, ...std::vector);
    //ptile1.push_back_int(4, ...std::vector);


    pc.Redistribute();
}


int main(int argc, char* argv[])
 {
    amrex::Initialize(argc,argv);
    {
        // for (int n = 0; n < BL_SPACEDIM; n++)
        // {
        //     real_box.setLo(n, 0.0);
        //     real_box.setHi(n, params.size[n]);
        // }

        // IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
        // IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
        // const Box domain(domain_lo, domain_hi);

        // int coord = 0;
        // int is_per[BL_SPACEDIM];
        // for (int i = 0; i < BL_SPACEDIM; i++)
        //     is_per[i] = params.is_periodic;
        // Geometry geom(domain, &real_box, coord, is_per);

        // BoxArray ba(domain);
        // ba.maxSize(params.max_grid_size);
        // DistributionMapping dm(ba);

        // const int ncells = 1;
    
        //addParticles< ParticleContainer<1,2,3,4> > ();
        addParticles< ParticleContainerPureSoA<3,4> > ();
    }
    amrex::Finalize();
 }

 


