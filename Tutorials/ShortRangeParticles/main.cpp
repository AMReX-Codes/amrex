#include <iostream>
#include <random>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

#include "compute_force_F.H"

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

struct GhostCommTag {

    GhostCommTag(int pid, int gid, int tid)
        : proc_id(pid), grid_id(gid), tile_id(tid)
    {}

    int proc_id;
    int grid_id;
    int tile_id;
};

bool operator<(const GhostCommTag& l, const GhostCommTag& r) {
    return (l.proc_id < r.proc_id || 
           (l.proc_id == r.proc_id && l.grid_id < r.grid_id) ||
           (l.proc_id == r.proc_id && l.grid_id == r.grid_id && l.tile_id < r.tile_id ));
}

class MyParticleContainer
    : public ParticleContainer<0, 0, PIdx::nattribs>
{
public:

    using PairIndex = std::pair<int, int>;
    using GhostCommMap = std::map<GhostCommTag, Array<char> >;
    
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
    /// This fills the ghost buffers for each tile with the proper data
    ///
    void fillGhosts() {
        GhostCommMap ghosts_to_comm;
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
                    packGhostParticle(neighbor_cell, mask[pti], p, ghosts_to_comm);
                }

                // Now add the particle to the "edge" neighbors
                for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                    for (int jdim = 0; jdim < idim; ++jdim) {
                        if (shift[idim] != 0 and shift[jdim] != 0) {
                            IntVect neighbor_cell = iv;
                            neighbor_cell.shift(idim, shift[idim]);
                            neighbor_cell.shift(jdim, shift[jdim]);
                            packGhostParticle(neighbor_cell, mask[pti], p, ghosts_to_comm);
                        }
                    }
                }

#if (BL_SPACEDIM == 3)

                // Finally, add the particle for the "vertex" neighbors (only relevant in 3D)
                if (shift[0] != 0 and shift[1] != 0 and shift[2] != 0) {
                    IntVect neighbor_cell = iv;
                    neighbor_cell.shift(shift);
                    packGhostParticle(neighbor_cell, mask[pti], p, ghosts_to_comm);
                }
#endif
            }
        }

        fillGhostsMPI(ghosts_to_comm);
    }

    ///
    /// Each tile clears its ghosts, freeing the memory
    ///
    void clearGhosts() {
	for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            const int grid_id = pti.index();
            const int tile_id = pti.LocalTileIndex();
            auto& ghost_particles = ghosts[std::make_pair(grid_id, tile_id)];
            Array<char>().swap(ghost_particles);
        }
    }

    ///
    /// Compute the short range forces on a tile's worth of particles.
    /// fillGhosts must have already been called.
    ///
    void computeForces() {
	for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            AoS& particles = pti.GetArrayOfStructs();
            size_t Np = particles.size();
            int nstride = particles.dataShape().first;

            Array<Real>& ax = pti.GetStructOfArrays()[PIdx::ax];
            Array<Real>& ay = pti.GetStructOfArrays()[PIdx::ay];

            PairIndex index(pti.index(), pti.LocalTileIndex());            
            int nghosts = ghosts[index].size() / pdata_size;

            Array<Real> vals(nghosts, 0);
            std::memcpy(&vals[0], ghosts[index].dataPtr(), sizeof(vals));

            amrex_compute_forces(particles.data(), nstride, Np, 
                                 (RealType*) ghosts[index].dataPtr(), nghosts,
                                 ax.dataPtr(), ay.dataPtr());
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

    ///
    /// Save the particle data in an ASCII format
    ///
    void writeParticles(int n) {
        const std::string& pltfile = amrex::Concatenate("particles", n, 5);
        WriteAsciiFile(pltfile);
    }

private:

    ///
    /// Pack a particle's data into the proper neighbor buffer, or put it into
    /// the structure to be sent to the other processes
    ///
    void packGhostParticle(const IntVect& neighbor_cell,
                           const BaseFab<int>& mask,
                           const ParticleType& p,
                           GhostCommMap& ghosts_to_comm) {
        const int neighbor_grid = mask(neighbor_cell, 0);
        if (neighbor_grid >= 0) {
            const int who = ParticleDistributionMap(lev)[neighbor_grid];
            const int MyProc = ParallelDescriptor::MyProc();
            const int neighbor_tile = mask(neighbor_cell, 1);
            PairIndex dst_index(neighbor_grid, neighbor_tile);
            if (who == MyProc) {
                size_t old_size = ghosts[dst_index].size();
                size_t new_size = ghosts[dst_index].size() + pdata_size;
                ghosts[dst_index].resize(new_size);
                std::memcpy(&ghosts[dst_index][old_size], &p, pdata_size);
            } else {
                GhostCommTag tag(who, neighbor_grid, neighbor_tile);
                Array<char>& buffer = ghosts_to_comm[tag];
                size_t old_size = buffer.size();
                size_t new_size = buffer.size() + pdata_size;
                buffer.resize(new_size);
                std::memcpy(&buffer[old_size], &p, pdata_size);
            }
        }
    }

    void fillGhostsMPI(GhostCommMap& ghosts_to_comm) {

#ifdef BL_USE_MPI
        const int MyProc = ParallelDescriptor::MyProc();
        const int NProcs = ParallelDescriptor::NProcs();

        // count the number of tiles to be sent to each proc
        std::map<int, int> tile_counts;
        for (const auto& kv: ghosts_to_comm) {
            tile_counts[kv.first.proc_id] += 1;
        }

        // flatten all the data for each proc into a single buffer
        // once this is done, each dst proc will have an Array<char>
        // the buffer will be packed like:
        // ntiles, gid1, tid1, size1, data1....  gid2, tid2, size2, data2... etc. 
        std::map<int, Array<char> > send_data;
        for (const auto& kv: ghosts_to_comm) {
            Array<char>& buffer = send_data[kv.first.proc_id];
            buffer.resize(sizeof(int));
            std::memcpy(&buffer[0], &tile_counts[kv.first.proc_id], sizeof(int));
        }

        for (auto& kv : ghosts_to_comm) {
            int data_size = kv.second.size();
            Array<char>& buffer = send_data[kv.first.proc_id];
            size_t old_size = buffer.size();
            size_t new_size = buffer.size() + 2*sizeof(int) + sizeof(int) + data_size;
            buffer.resize(new_size);
            char* dst = &buffer[old_size];
            std::memcpy(dst, &(kv.first.grid_id), sizeof(int)); dst += sizeof(int);
            std::memcpy(dst, &(kv.first.tile_id), sizeof(int)); dst += sizeof(int);
            std::memcpy(dst, &data_size,          sizeof(int)); dst += sizeof(int);
            if (data_size == 0) continue;
            std::memcpy(dst, &kv.second[0], data_size);
            Array<char>().swap(kv.second);
        }

        // each proc figures out how many bytes it will send, and how
        // many it will receive
        Array<long> snds(NProcs, 0), rcvs(NProcs, 0);
        long num_snds = 0;
        for (const auto& kv : send_data) {
            num_snds      += kv.second.size();
            snds[kv.first] = kv.second.size();
        }
        ParallelDescriptor::ReduceLongMax(num_snds);
        if (num_snds == 0) return;

        // communicate that information
        BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                        ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

        BL_MPI_REQUIRE( MPI_Alltoall(snds.dataPtr(),
                                     1,
                                     ParallelDescriptor::Mpi_typemap<long>::type(),
                                     rcvs.dataPtr(),
                                     1,
                                     ParallelDescriptor::Mpi_typemap<long>::type(),
                                     ParallelDescriptor::Communicator()) );
        BL_ASSERT(rcvs[MyProc] == 0);
    
        BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                        ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

        Array<int> RcvProc;
        Array<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    
        std::size_t TotRcvBytes = 0;
        for (int i = 0; i < NProcs; ++i) {
            if (rcvs[i] > 0) {
                RcvProc.push_back(i);
                rOffset.push_back(TotRcvBytes);
                TotRcvBytes += rcvs[i];
            }
        }
    
        const int nrcvs = RcvProc.size();
        Array<MPI_Status>  stats(nrcvs);
        Array<MPI_Request> rreqs(nrcvs);
        
        const int SeqNum = ParallelDescriptor::SeqNum();
        
        // Allocate data for rcvs as one big chunk.
        Array<char> recvdata(TotRcvBytes);

        // Post receives.
        for (int i = 0; i < nrcvs; ++i) {
            const auto Who    = RcvProc[i];
            const auto offset = rOffset[i];
            const auto Cnt    = rcvs[Who];
      
            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());
            BL_ASSERT(Who >= 0 && Who < NProcs);
      
            rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum).req();
        }
    
        // Send.
        for (const auto& kv : send_data) {
            const auto Who = kv.first;
            const auto Cnt = kv.second.size();
            
            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Who >= 0 && Who < NProcs);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());
            
            ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum);
        }
    
        // unpack the received data and put them into the proper ghost buffers
        if (nrcvs > 0) {
            BL_MPI_REQUIRE( MPI_Waitall(nrcvs, rreqs.data(), stats.data()) );
            for (int i = 0; i < nrcvs; ++i) {
                const int offset = rOffset[i];
                char* buffer = &recvdata[offset];
                int num_tiles, gid, tid, size;
                std::memcpy(&num_tiles, buffer, sizeof(int)); buffer += sizeof(int);
                for (int j = 0; j < num_tiles; ++j) {
                    std::memcpy(&gid,  buffer, sizeof(int)); buffer += sizeof(int);
                    std::memcpy(&tid,  buffer, sizeof(int)); buffer += sizeof(int);
                    std::memcpy(&size, buffer, sizeof(int)); buffer += sizeof(int);

                    if (size == 0) continue;

                    PairIndex dst_index(gid, tid);
                    size_t old_size = ghosts[dst_index].size();
                    size_t new_size = ghosts[dst_index].size() + size;
                    ghosts[dst_index].resize(new_size);
                    std::memcpy(&ghosts[dst_index][old_size], buffer, size); buffer += size;
                }
            }
        }
#endif
    }

    const int lev = 0;
    const size_t pdata_size = 2*sizeof(RealType); // we communicate 2 reals (x, y) per ghost
    FabArray<BaseFab<int> > mask;
    std::map<PairIndex, Array<char> > ghosts;
};

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
  
    int nx = 32;
    int ny = 32;
    int max_step = 1000;
    Real dt = 0.0005;
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

    for (int i = 0; i < max_step; i++) {
        myPC.writeParticles(i);
        
        myPC.fillGhosts();
        myPC.computeForces();
        myPC.clearGhosts();

        myPC.moveParticles(dt);

        myPC.Redistribute();
    }

    myPC.writeParticles(max_step);
    
    amrex::Finalize();
}
