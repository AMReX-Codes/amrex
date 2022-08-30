#include <AMReX.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_ParIter.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>

#include <array>

using namespace amrex;

template <typename T_PC,template<class> class Allocator=DefaultAllocator>
void addParticles ()
{
    T_PC pc;
    int const NReal = pc.NStructReal;
    int const NInt = pc.NStructInt;
    int const NArrayReal = pc.NArrayReal;
    int const NArrayInt = pc.NArrayInt;

    using ParticleType = typename T_PC::ParticleType;
    using ParticleTileDataType = typename T_PC::ParticleTileType::ParticleTileDataType;
    using RealVector = amrex::PODVector<ParticleReal, Allocator<ParticleReal> >;
    using IntVector = amrex::PODVector<int, Allocator<int> >;

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
    
    int lev=0;
    // int numparticles=0;
    using MyParIter = ParIter_impl<ParticleType, NArrayReal, NArrayInt>;
    for (MyParIter pti(pc, lev); pti.isValid(); ++pti) {
        const int np = pti.numParticles();
        //const auto t_lev = pti.GetLevel();
        //const auto index = pti.GetPairIndex();
        // ...

        // preparing access to particle data: AoS
        //using PType = ImpactXParticleContainer::ParticleType;
        auto& aos = pti.GetArrayOfStructs();
        ParticleType* AMREX_RESTRICT aos_ptr = aos().dataPtr();

        // preparing access to particle data: SoA of Reals
        auto& soa_real = pti.GetStructOfArrays().GetRealData();
        amrex::ParticleReal* const AMREX_RESTRICT part_x = soa_real[0].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_y = soa_real[1].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_z = soa_real[2].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_aaa = soa_real[3].dataPtr();
        auto& soa_int = pti.GetStructOfArrays().GetIntData();

        ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
        {
            ParticleType& AMREX_RESTRICT p = aos_ptr[ip];
            p.pos(0) += 1;
            p.pos(1) += 1;
            p.pos(2) += 1;
            
            amrex::ParticleReal & AMREX_RESTRICT x = part_x[ip];
            amrex::ParticleReal & AMREX_RESTRICT y = part_y[ip];
            amrex::ParticleReal & AMREX_RESTRICT z = part_z[ip];

            x += 1.0;
            y += 1.0;
            z += 1.0;
        });

        // new way of creating
        ParticleTileDataType ptd = pti.GetParticleTile().getParticleTileData(); 

        ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
        {
            ParticleType p(ptd, ip);
            p.pos(0) += 1;
            p.pos(1) += 1;
            p.pos(2) += 1;
        });


    }

    for (MyParIter pti(pc, lev); pti.isValid(); ++pti) {
        auto& particle_attributes = pti.GetStructOfArrays();
        RealVector& real_comp0 = particle_attributes.GetRealData(0);
        IntVector&  int_comp1  = particle_attributes.GetIntData(1);
        for (int i = 0; i < pti.numParticles(); ++i) {
            real_comp0[i] += 1;
            int_comp1[i] += 1;
        }
    }

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

 


