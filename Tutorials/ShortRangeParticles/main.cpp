#include <iostream>
#include <random>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

using namespace amrex;

///
/// This describes the variables stored in the array-of-structs
///
struct PIdx
{
    enum {
	vx, vy, ax, ay,
	nattribs
    };
};

class MyParIter
    : public amrex::ParIter<0,0,PIdx::nattribs>
{
public:
    using amrex::ParIter<0,0,PIdx::nattribs>::ParIter;

    ///
    /// Define some convenient wrappers for accessing particle data
    ///

    ParticleType::RealType& x(int i) { return GetArrayOfStructs()[i].pos(0); }
    ParticleType::RealType& y(int i) { return GetArrayOfStructs()[i].pos(1); }

    ParticleType::RealType& vx(int i) { return GetStructOfArrays()[PIdx::vx][i]; }
    ParticleType::RealType& vy(int i) { return GetStructOfArrays()[PIdx::vy][i]; }

    ParticleType::RealType& ax(int i) { return GetStructOfArrays()[PIdx::ax][i]; }
    ParticleType::RealType& ay(int i) { return GetStructOfArrays()[PIdx::ay][i]; }

};

class MyParticleContainer
    : public ParticleContainer<0, 0, PIdx::nattribs>
{
public:

    using PairIndex = std::pair<int, int>;
    
    ///
    /// This particle container fills a mask for quickly computing
    /// neighbor grids / tiles for a given particle
    ///
    MyParticleContainer(const Geometry            & geom, 
                        const DistributionMapping & dmap,
                        const BoxArray            & ba)
    	: ParticleContainer<0, 0, PIdx::nattribs> (geom, dmap, ba)
    {
                
        mask.define(ba, dmap, 2, 1);
        mask.setVal(-1, 1);
        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const Box& box = mfi.tilebox();
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            mask.setVal(grid_id, box, 0, 1);
            mask.setVal(tile_id, box, 1, 1);
        }
        mask.FillBoundary();
    }

    ///
    /// Init one particle per cell with random velocities
    ///
    void InitParticles() {
        const Geometry& geom = Geom(lev);
        const Real* dx  = geom.CellSize();
        
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        
        ParticleType p;
        std::array<Real,PIdx::nattribs> attribs;
        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const Box& tile_box = mfi.tilebox();
            const RealBox tile_real_box { tile_box, dx, geom.ProbLo() };

            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];

            const auto& boxlo = tile_box.smallEnd();
            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {

                p.id() = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();

                p.pos(0) = tile_real_box.lo(0) + (iv[0]- boxlo[0] + 0.5)*dx[0];
                p.pos(1) = tile_real_box.lo(1) + (iv[1]- boxlo[1] + 0.5)*dx[1];

                attribs[PIdx::vx] = dist(mt);
                attribs[PIdx::vy] = dist(mt);
                attribs[PIdx::ax] = 0;
                attribs[PIdx::ay] = 0;

                particle_tile.push_back(p);
                particle_tile.push_back(attribs);
            }
        }
    }

    ///
    /// Pack a particle's data into the proper neighbor buffer
    ///
    void packNeighborParticle(const IntVect& neighbor_cell,
                              const MyParIter& pti,
                              const ParticleType& p) {
        const int neighbor_grid = mask[pti](neighbor_cell, 0);
        if (neighbor_grid >= 0) {
            const int neighbor_tile = mask[pti](neighbor_cell, 1);
            const int grid = pti.index();
            const int tile = pti.LocalTileIndex();
            PairIndex src_index(grid, tile);
            PairIndex dst_index(neighbor_grid, neighbor_tile);
            Array<char>& buffer = neighbors[src_index][dst_index];
            size_t old_size = buffer.size();
            size_t new_size = buffer.size() + pdata_size;
            buffer.resize(new_size);
            std::memcpy(&buffer[old_size], &p, pdata_size);
        }
    }

    ///
    /// This fills the neighbor buffers for each tile with the proper data
    ///
    void computeNeighbors() {
	for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            const Box& tile_box = pti.tilebox();
            const IntVect& lo = tile_box.smallEnd();
            const IntVect& hi = tile_box.bigEnd();

	    Box shrink_box = pti.tilebox();
            shrink_box.grow(-1);

            auto& particles = pti.GetArrayOfStructs();
            for (unsigned i = 0; i < pti.numParticles(); ++i) {
                const ParticleType& p = particles[i];
                const IntVect& iv = Index(p, lev);

                // if the particle is more than one cell away from 
                // the tile boundary, it's not anybody's neighbor
                if (shrink_box.contains(iv)) continue;

                // shift stores whether we are near the tile boundary in each direction.
                // -1 means lo, 1 means hi, 0 means not near the boundary
                IntVect shift = IntVect::TheZeroVector();
                for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                    if (iv[idim] == lo[idim])
                        shift[idim] = -1;
                    else if (iv[idim] == hi[idim])
                        shift[idim] = 1;
                }

                // Based on the value of shift, we add the particle to a map to be sent
                // to the neighbors. A particle can be sent to up to 3 neighbors in 2D
                // and up to 7 in 3D, depending on whether is near the tile box corners,
                // edges, or just the faces. First, add the particle for the "face" neighbors
                for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                    if (shift[idim] == 0) continue;
                    IntVect neighbor_cell = iv;
                    neighbor_cell.shift(idim, shift[idim]);
                    packNeighborParticle(neighbor_cell, pti, p);
                }

                // Now add the particle to the "edge" neighbors
                for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                    for (int jdim = 0; jdim < idim; ++jdim) {
                        if (shift[idim] != 0 and shift[jdim] != 0) {
                            IntVect neighbor_cell = iv;
                            neighbor_cell.shift(idim, shift[idim]);
                            neighbor_cell.shift(jdim, shift[jdim]);
                            packNeighborParticle(neighbor_cell, pti, p);
                        }
                    }
                }

#if (BL_SPACEDIM == 3)

                // Finally, add the particle for the "vertex" neighbors (only relevant in 3D)
                if (shift[0] != 0 and shift[1] != 0 and shift[2] != 0) {
                    IntVect neighbor_cell = iv;
                    neighbor_cell.shift(shift);
                    packNeighborParticle(neighbor_cell, pti, p);
                }
#endif
            }
        }
    }

    ///
    /// Each tile prints out the number of particles it will send to each other tile.
    ///
    void printNeighbors() {
	for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            const int grid_id = pti.index();
            const int tile_id = pti.LocalTileIndex();
            std::cout << grid_id << " " << tile_id << "\n";
            const auto& buffers = neighbors[std::make_pair(grid_id, tile_id)];
            for (auto it = buffers.begin(); it != buffers.end(); ++it) {
                const int neighbor_grid = it->first.first;
                const int neighbor_tile = it->first.second;
                const int num_particles = it->second.size() / pdata_size;  // divide by bytes per particle
                std::cout << neighbor_grid << " " << neighbor_tile << ", " << num_particles << '\n';
            }            
        }
    }

    ///
    /// Move the particles according to their forces, reflecting at domain boundaries
    ///
    void moveParticles(const Real dt) {
        const RealBox& prob_domain = Geom(lev).ProbDomain();

        const Real xlo = prob_domain.lo(0);
        const Real xhi = prob_domain.hi(0);
        const Real ylo = prob_domain.lo(1);
        const Real yhi = prob_domain.hi(1);

	for (MyParIter parts(*this, lev); parts.isValid(); ++parts) {
            for (unsigned i = 0; i < parts.numParticles(); ++i) {
                
                parts.vx(i) += parts.ax(i) * dt;
                parts.vy(i) += parts.ay(i) * dt;

                parts.x(i) += parts.vx(i) * dt;
                parts.y(i) += parts.vy(i) * dt;

                // bounce off the walls
                while( parts.x(i) < xlo || parts.x(i) > xhi) {
                    parts.x(i) = parts.x(i) < xlo ? 2*xlo - parts.x(i) : 2*xhi - parts.x(i);
                    parts.vx(i) = -parts.vx(i);
                }
                while( parts.y(i) < xlo || parts.y(i) > xhi) {
                    parts.y(i) = parts.y(i) < ylo ? 2*ylo - parts.y(i) : 2*yhi - parts.y(i);
                    parts.vy(i) = -parts.vy(i);
                }
            }
        }
    }

    void writeParticles(int n) {
        const std::string& pltfile = amrex::Concatenate("particles", n, 5);
        WriteAsciiFile(pltfile);
    }

private:

    const int lev = 0;
    const size_t pdata_size = 2*sizeof(RealType);
    FabArray<BaseFab<int> > mask;

    std::map<PairIndex, std::map<PairIndex, Array<char> > > neighbors;

};

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
  
    int nx = 32;
    int ny = 32;
    int max_step = 1000;
    Real dt = 0.0005;
    bool verbose = false;    
    int max_grid_size = 16;

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(0, 0);
    IntVect domain_hi(nx - 1, ny - 1);
    const Box domain(domain_lo, domain_hi);
    
    int coord = 0;
    
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 0; 
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    
    DistributionMapping dmap(ba);
    
    MyParticleContainer myPC(geom, dmap, ba);

    myPC.InitParticles();

    myPC.computeNeighbors();
    if (verbose) myPC.printNeighbors();

    for (int i = 0; i < max_step; i++) {
        if (verbose) myPC.writeParticles(i);
        myPC.moveParticles(dt);
        myPC.Redistribute();
    }

    if (verbose) myPC.writeParticles(max_step);
    
    amrex::Finalize();
}
