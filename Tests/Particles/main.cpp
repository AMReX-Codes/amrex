#include <iostream>
#include <map>
#include <vector>
#include <type_traits>

#include <AMReX_Vector.H>
#include "AMReX_Particles.H"

using namespace amrex;

struct RealIdx {
    enum {
        id = 0,
        cpu,
        nattribs
    };
};

struct IntIdx {
    enum {
        id = 0,
        cpu,
        nattribs
    };
};

/*

  This is a particle container that demonstrates how to layout the particle
  data in struct-of-arrays form. We refer to data laid out in this manner as
  particle "attributes". In this test, in addition to being stored in the particle
  struct itself, we also store the particle id and cpu in additional arrays,
  once as Reals and once as integers. The test below then moves the particles
  randomly for 10 time steps, making sure that the particle data is redistributed
  correctly.
  
 */
class StructOfArraysParticleContainer 
    : public ParticleContainer<0, 0,
                               RealIdx::nattribs,
                               IntIdx::nattribs>
{
 public:
    
    StructOfArraysParticleContainer (const Vector<Geometry>            & geom, 
                                     const Vector<DistributionMapping> & dmap,
                                     const Vector<BoxArray>            & ba,
                                     const Vector<int>                 & rr)
        : ParticleContainer<0, 0,
                            RealIdx::nattribs,
                            IntIdx::nattribs> (geom, dmap, ba, rr)
    {
        AddRealComp(true);
        AddRealComp(true);

        AddIntComp(true);
        AddIntComp(true);
    }

    void InitParticles() {
        const int lev = 0;
        const Geometry& geom = Geom(lev);
        const Real* dx  = geom.CellSize();
        
        ParticleType p;
        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const Box& tile_box = mfi.tilebox();
            const RealBox tile_real_box { tile_box, dx, geom.ProbLo() };
            
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            auto& particle_tile = DefineAndReturnParticleTile(lev, grid_id, tile_id);
            
            const auto& boxlo = tile_box.smallEnd();
            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
                
                p.id() = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
            
                AMREX_D_TERM(p.pos(0) = tile_real_box.lo(0) + (iv[0]- boxlo[0] + 0.5)*dx[0];,
                       p.pos(1) = tile_real_box.lo(1) + (iv[1]- boxlo[1] + 0.5)*dx[1];,
                       p.pos(2) = tile_real_box.lo(2) + (iv[2]- boxlo[2] + 0.5)*dx[2];);

                // set this particle's real attributes
                std::array<double, RealIdx::nattribs> real_attribs;
                real_attribs[RealIdx::id]  = p.id();
                real_attribs[RealIdx::cpu] = p.cpu();

                // set this particle's integer attributes
                std::array<int, IntIdx::nattribs> int_attribs;
                int_attribs[IntIdx::id]  = p.id();
                int_attribs[IntIdx::cpu] = p.cpu();
                                
                particle_tile.push_back(p);
                particle_tile.push_back_real(real_attribs);
                particle_tile.push_back_int(int_attribs);

                particle_tile.push_back_real(2, p.id());
                particle_tile.push_back_real(3, p.cpu());

                particle_tile.push_back_int(2, p.id());
                particle_tile.push_back_int(3, p.cpu());
            }
        }
    }
};

int main(int argc, char* argv[])
{
    
    amrex::Initialize(argc,argv);

    int ncell = 64;
    int nlevs = 2;

    RealBox real_box, fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n,0.0);
        real_box.setHi(n,1.0);
        fine_box.setLo(n,0.4);
        fine_box.setHi(n,0.6);
    }
    
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(ncell-1, ncell-1, ncell-1); 
    
    const Box domain(domain_lo, domain_hi);
    
    Vector<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;
 
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1;

    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per);
    geom[1].define(amrex::refine(geom[0].Domain(), rr[0]),
                   &real_box, CoordSys::cartesian, is_per);

    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);  

    int nfine = ncell*rr[0];
    IntVect refined_lo(nfine/4,nfine/4,nfine/4); 
    IntVect refined_hi(3*nfine/4-1,3*nfine/4-1,3*nfine/4-1);
    
    // Build a box for the level 1 domain
    Box refined_patch(refined_lo, refined_hi);
    ba[1].define(refined_patch);
    
    int max_grid_size = 32;
    for (int lev = 0; lev < nlevs; lev++)
        ba[lev].maxSize(max_grid_size);
    
    Vector<DistributionMapping> dmap(nlevs);
    for (int lev = 0; lev < nlevs; lev++)
        dmap[lev].define(ba[lev]);
    
    StructOfArraysParticleContainer MyPC(geom, dmap, ba, rr);
    
    MyPC.InitParticles();
    MyPC.Redistribute();
    
    for (int i = 0; i < 10; i ++)
        MyPC.MoveRandom();
    
    MyPC.WriteAsciiFile("particles");
    
    amrex::Finalize();
}
