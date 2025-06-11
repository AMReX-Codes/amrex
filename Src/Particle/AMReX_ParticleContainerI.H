#include <AMReX_MakeParticle.H>

#include <algorithm>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace amrex {

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::SetParticleSize ()
{
    num_real_comm_comps  = 0;
    int comm_comps_start = 0;
    if constexpr (!ParticleType::is_soa_particle) {
        comm_comps_start += AMREX_SPACEDIM + NStructReal;
    }
    for (int i = comm_comps_start; i < comm_comps_start + NumRealComps(); ++i) {
        if (h_redistribute_real_comp[i]) {++num_real_comm_comps;}
    }

    num_int_comm_comps = 0;
    comm_comps_start   = 2 + NStructInt;
    for (int i = comm_comps_start; i < comm_comps_start + NumIntComps(); ++i) {
        if (h_redistribute_int_comp[i]) {++num_int_comm_comps;}
    }

    if constexpr (ParticleType::is_soa_particle) {
        particle_size = sizeof(uint64_t);  // idcpu
    } else {
        particle_size = sizeof(ParticleType);
    }
    superparticle_size = particle_size +
        num_real_comm_comps*sizeof(ParticleReal) + num_int_comm_comps*sizeof(int);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor> :: Initialize ()
{
    levelDirectoriesCreated = false;
    usePrePost = false;
    doUnlink = true;

    // by default communicate all components
   if constexpr (ParticleType::is_soa_particle)
   {
       h_redistribute_real_comp.resize(NArrayReal, true);
   } else {
       h_redistribute_real_comp.resize(AMREX_SPACEDIM + NStructReal + NArrayReal, true);
   }
   h_redistribute_int_comp.resize(2 + NStructInt + NArrayInt, true);

    SetParticleSize();

    // add default names for SoA Real and Int compile-time arguments
    m_soa_rdata_names.clear();
    for (int i=0; i<NArrayReal; ++i)
    {
        m_soa_rdata_names.push_back(getDefaultCompNameReal<ParticleType>(i));
    }
    m_soa_idata_names.clear();
    for (int i=0; i<NArrayInt; ++i)
    {
        m_soa_idata_names.push_back(getDefaultCompNameInt<ParticleType>(i));
    }

    static bool initialized = false;
    if ( ! initialized)
    {
        static_assert(sizeof(ParticleType)%sizeof(RealType) == 0,
                      "sizeof ParticleType is not a multiple of sizeof RealType");

        ParmParse pp("particles");
        pp.queryAdd("do_tiling", do_tiling);
        Vector<int> tilesize(AMREX_SPACEDIM);
        if (pp.queryarr("tile_size", tilesize, 0, AMREX_SPACEDIM)) {
            for (int i=0; i<AMREX_SPACEDIM; ++i) { tile_size[i] = tilesize[i]; }
        }

        static_assert(std::is_standard_layout<ParticleType>::value,
                      "Particle type must be standard layout");
        //                   && std::is_trivial<ParticleType>::value,
        //                      "Particle type must be standard layout and trivial.");

        pp.query("use_prepost", usePrePost);
        pp.query("do_unlink", doUnlink);
        pp.queryAdd("do_mem_efficient_sort", memEfficientSort);
        pp.queryAdd("use_comms_arena", use_comms_arena);

        initialized = true;
    }
}

template<
    typename ParticleType,
    int NArrayReal,
    int NArrayInt,
    template<class> class Allocator,
    class CellAssignor
>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::SetSoACompileTimeNames (
    std::vector<std::string> const & rdata_name, std::vector<std::string> const & idata_name
)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rdata_name.size() == NArrayReal, "rdata_name must be equal to NArrayReal");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(idata_name.size() == NArrayInt, "idata_name must be equal to NArrayInt");

    // ensure names for components are unique
    std::set<std::string> const unique_r_names(rdata_name.begin(), rdata_name.end());
    std::set<std::string> const unique_i_names(idata_name.begin(), idata_name.end());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rdata_name.size() == unique_r_names.size(), "SetSoACompileTimeNames: Provided names in rdata_name are not unique!");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(idata_name.size() == unique_i_names.size(), "SetSoACompileTimeNames: Provided names in idata_name are not unique!");

    for (int i=0; i<NArrayReal; ++i)
    {
        m_soa_rdata_names.at(i) = rdata_name.at(i);
    }
    for (int i=0; i<NArrayInt; ++i)
    {
        m_soa_idata_names.at(i) = idata_name.at(i);
    }
}

template<
    typename ParticleType,
    int NArrayReal,
    int NArrayInt,
    template<class> class Allocator,
    class CellAssignor
>
bool
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::HasRealComp (std::string const & name)
{
    return std::find(m_soa_rdata_names.begin(), m_soa_rdata_names.end(), name) != std::end(m_soa_rdata_names);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
        template<class> class Allocator, class CellAssignor>
bool
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::HasIntComp (std::string const & name)
{
    return std::find(m_soa_idata_names.begin(), m_soa_idata_names.end(), name) != std::end(m_soa_idata_names);
}

template<
    typename ParticleType,
    int NArrayReal,
    int NArrayInt,
    template<class> class Allocator,
    class CellAssignor
>
int
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::GetRealCompIndex (std::string const & name)
{
    const auto it = std::find(m_soa_rdata_names.begin(), m_soa_rdata_names.end(), name);

    if (it == m_soa_rdata_names.end())
    {
        throw std::runtime_error("GetRealCompIndex: Component " + name + " does not exist!");
    }
    else
    {
        return std::distance(m_soa_rdata_names.begin(), it);
    }
}

template<
    typename ParticleType,
    int NArrayReal,
    int NArrayInt,
    template<class> class Allocator,
    class CellAssignor
>
int
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::GetIntCompIndex (std::string const & name)
{
    const auto it = std::find(m_soa_idata_names.begin(), m_soa_idata_names.end(), name);

    if (it == m_soa_idata_names.end())
    {
        throw std::runtime_error("GetIntCompIndex: Component " + name + " does not exist!");
    }
    else
    {
        return std::distance(m_soa_idata_names.begin(), it);
    }
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <typename P, typename Assignor>
IntVect
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::Index (const P& p, int lev) const
{
    const Geometry& geom = Geom(lev);
    const auto& domain = geom.Domain();
    const auto& plo = geom.ProbLoArray();
    const auto& dxi = geom.InvCellSizeArray();

    return Assignor{}(p, plo, dxi, domain);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <typename P>
bool
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::Where (const P& p,
         ParticleLocData&    pld,
         int                 lev_min,
         int                 lev_max,
         int                 nGrow,
         int                 local_grid) const
{

  AMREX_ASSERT(m_gdb != nullptr);

  if (lev_max == -1) {
      lev_max = finestLevel();
  }

  AMREX_ASSERT(lev_max <= finestLevel());

  AMREX_ASSERT(nGrow == 0 || (nGrow >= 0 && lev_min == lev_max));

  std::vector< std::pair<int, Box> > isects;

  for (int lev = lev_max; lev >= lev_min; lev--) {
      const IntVect& iv = Index(p, lev);
      if (lev == pld.m_lev) {
          // The fact that we are here means this particle does not belong to any finer grids.
          if (pld.m_grid >= 0) {
              if (pld.m_grown_gridbox.contains(iv)) {
                  pld.m_cell = iv;
                  if (!pld.m_tilebox.contains(iv)) {
                      pld.m_tile = getTileIndex(iv, pld.m_gridbox, do_tiling, tile_size, pld.m_tilebox);
                  }
                  return true;
              }
          }
      }

      int grid;
      const BoxArray& ba = ParticleBoxArray(lev);
      AMREX_ASSERT(ba.ixType().cellCentered());

      if (local_grid < 0) {
          bool findfirst = (nGrow == 0) ? true : false;
          ba.intersections(Box(iv, iv), isects, findfirst, nGrow);
          grid = isects.empty() ? -1 : isects[0].first;
          if (nGrow > 0 && isects.size() > 1) {
             for (auto & isect : isects) {
                 Box bx = ba.getCellCenteredBox(isect.first);
                 for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                     Box gbx = bx;
                     IntVect gr(IntVect::TheZeroVector());
                     gr[dir] = nGrow;
                     gbx.grow(gr);
                     if (gbx.contains(iv)) {
                        grid = isect.first;
                     }
                 }
             }
          }
      } else {
          grid = (*redistribute_mask_ptr)[local_grid](iv, 0);
      }

      if (grid >= 0) {
          const Box& bx = ba.getCellCenteredBox(grid);
          pld.m_lev  = lev;
          pld.m_grid = grid;
          pld.m_tile = getTileIndex(iv, bx, do_tiling, tile_size, pld.m_tilebox);
          pld.m_cell = iv;
          pld.m_gridbox = bx;
          pld.m_grown_gridbox = amrex::grow(bx, nGrow);
          return true;
      }
  }

  return false;
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <typename P>
bool
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::EnforcePeriodicWhere (P&               p,
                        ParticleLocData& pld,
                        int              lev_min,
                        int              lev_max,
                        int              local_grid) const
{

    AMREX_ASSERT(m_gdb != nullptr);

    if (!Geom(0).isAnyPeriodic()) { return false; }

    if (lev_max == -1) {
        lev_max = finestLevel();
    }

    AMREX_ASSERT(lev_max <= finestLevel());

    // Create a copy "dummy" particle to check for periodic outs.
    Particle<0, 0> p_prime;
    AMREX_D_TERM(p_prime.pos(0) = p.pos(0);,
                 p_prime.pos(1) = p.pos(1);,
                 p_prime.pos(2) = p.pos(2));
    if (PeriodicShift(p_prime)) {
        std::vector< std::pair<int,Box> > isects;
        for (int lev = lev_max; lev >= lev_min; lev--) {

            int grid;
            IntVect iv;
            const BoxArray& ba = ParticleBoxArray(lev);
            AMREX_ASSERT(ba.ixType().cellCentered());

            if (local_grid < 0) {
                iv = Index<amrex::Particle<0, 0>, DefaultAssignor>(p_prime, lev);
                ba.intersections(Box(iv, iv), isects, true, 0);
                grid = isects.empty() ? -1 : isects[0].first;
            } else {
                iv = Index<amrex::Particle<0, 0>, DefaultAssignor>(p_prime, lev);
                if (ba[local_grid].contains(iv))
                {
                    grid = local_grid;
                }
                else
                {
                    ba.intersections(Box(iv, iv), isects, true, 0);
                    grid = isects.empty() ? -1 : isects[0].first;
                    if(grid == -1)
                    {
                        grid = (*redistribute_mask_ptr)[local_grid](Index(p, lev), 0);
                    }
                }
            }

            if (grid >= 0) {
                AMREX_D_TERM(p.pos(0) = p_prime.pos(0);,
                             p.pos(1) = p_prime.pos(1);,
                             p.pos(2) = p_prime.pos(2););

                const Box& bx = ba.getCellCenteredBox(grid);

                pld.m_lev  = lev;
                pld.m_grid = grid;
                pld.m_tile = getTileIndex(iv, bx, do_tiling, tile_size, pld.m_tilebox);
                pld.m_cell = iv;
                pld.m_gridbox = bx;
                pld.m_grown_gridbox = bx;
                return true;
            }
        }
    }

    return false;
}


template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <typename P>
bool
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::PeriodicShift (P& p) const
{
    const auto& geom = Geom(0);
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto rlo = geom.ProbLoArrayInParticleReal();
    const auto rhi = geom.ProbHiArrayInParticleReal();
    const auto is_per = geom.isPeriodicArray();

    return enforcePeriodic(p, plo, phi, rlo, rhi, is_per);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
ParticleLocData
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
Reset (ParticleType& p,
       bool          /*update*/,
       bool          verbose,
       ParticleLocData pld) const
{
    AMREX_ASSERT(m_gdb != nullptr);

    bool ok = Where(p, pld);

    if (!ok && Geom(0).isAnyPeriodic())
    {
        // Attempt to shift the particle back into the domain if it
        // crossed a periodic boundary.
      PeriodicShift(p);
      ok = Where(p, pld);
    }

    if (!ok) {
        // invalidate the particle.
        if (verbose) {
            amrex::AllPrint()<< "Invalidating out-of-domain particle: " << p << '\n';
        }

        AMREX_ASSERT(p.id() > 0);

        p.id() = -p.id();
    }

    return pld;
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::reserveData ()
{
    this->ParticleContainerBase::reserveData();
    m_particles.reserve(maxLevel()+1);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::resizeData ()
{
    this->ParticleContainerBase::resizeData();
    int nlevs = std::max(0, finestLevel()+1);
    m_particles.resize(nlevs);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <typename P>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::locateParticle (P& p, ParticleLocData& pld,
                                                                                        int lev_min, int lev_max, int nGrow, int local_grid) const
{
    bool success;
    if (Geom(0).outsideRoundoffDomain(AMREX_D_DECL(p.pos(0), p.pos(1), p.pos(2))))
    {
        // Note that EnforcePeriodicWhere may shift the particle if it is successful.
        success = EnforcePeriodicWhere(p, pld, lev_min, lev_max, local_grid);
        if (!success && lev_min == 0)
        {
            // The particle has left the domain; invalidate it.
            p.id() = -p.id();
            success = true;
        }
    }
    else
    {
        success = Where(p, pld, lev_min, lev_max, 0, local_grid);
    }

    if (!success)
    {
        success = (nGrow > 0) && Where(p, pld, lev_min, lev_min, nGrow);
        pld.m_grown_gridbox = pld.m_gridbox; // reset grown box for subsequent calls.
    }

    if (!success)
    {
        amrex::Abort("ParticleContainer::locateParticle(): invalid particle.");
    }
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
Long
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::TotalNumberOfParticles (bool only_valid, bool only_local) const
{
    Long nparticles = 0;
    for (int lev = 0; lev <= finestLevel(); lev++) {
        nparticles += NumberOfParticlesAtLevel(lev,only_valid,true);
    }
    if (!only_local) {
        ParallelAllReduce::Sum(nparticles, ParallelContext::CommunicatorSub());
    }
    return nparticles;
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
Vector<Long>
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::NumberOfParticlesInGrid (int lev, bool only_valid, bool only_local) const
{
    AMREX_ASSERT(lev >= 0 && lev < int(m_particles.size()));

    LayoutData<Long> np_per_grid_local(ParticleBoxArray(lev),
                                       ParticleDistributionMap(lev));

    for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        int gid = pti.index();
        if (only_valid)
        {
            const auto& ptile = ParticlesAt(lev, pti);
            const int np = ptile.numParticles();
            auto const ptd = ptile.getConstParticleTileData();

            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            reduce_op.eval(np, reduce_data,
                           [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                           {
                               return (ptd.id(i) > 0) ? 1 : 0;
                           });

            int np_valid = amrex::get<0>(reduce_data.value(reduce_op));
            np_per_grid_local[gid] += np_valid;
        } else
        {
            np_per_grid_local[gid] += pti.numParticles();
        }
    }

    Vector<Long> nparticles(np_per_grid_local.size(), 0);
    if (only_local)
    {
        for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            nparticles[pti.index()] = np_per_grid_local[pti.index()];
        }
    }
    else
    {
        ParallelDescriptor::GatherLayoutDataToVector(np_per_grid_local, nparticles,
                                                     ParallelContext::IOProcessorNumberSub());
        ParallelDescriptor::Bcast(nparticles.data(), nparticles.size(),
                                  ParallelContext::IOProcessorNumberSub());
    }

    return nparticles;
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
Long ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::NumberOfParticlesAtLevel (int level, bool only_valid, bool only_local) const
{
    Long nparticles = 0;

    if (level < 0 || level >= int(m_particles.size())) { return nparticles; }

    if (only_valid) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<unsigned long long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (const auto& kv : GetParticles(level)) {
            const auto& ptile = kv.second;
            auto const ptd = ptile.getConstParticleTileData();

            reduce_op.eval(ptile.numParticles(), reduce_data,
                           [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                           {
                               return (ptd.id(i) > 0) ? 1 : 0;
                           });
        }

        nparticles = static_cast<Long>(amrex::get<0>(reduce_data.value(reduce_op)));
    }
    else {
        for (const auto& kv : GetParticles(level)) {
            const auto& ptile = kv.second;
            nparticles += ptile.numParticles();
        }
    }

    if (!only_local) {
        ParallelAllReduce::Sum(nparticles, ParallelContext::CommunicatorSub());
    }

    return nparticles;
}

//
// This includes both valid and invalid particles since invalid particles still take up space.
//

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
std::array<Long, 3>
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::ByteSpread () const
{
    Long cnt = 0;

    for (unsigned lev = 0; lev < m_particles.size(); lev++) {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& ptile = kv.second;
            cnt += ptile.numParticles();
        }
    }

    Long mn = cnt, mx = mn;

    const int IOProc = ParallelContext::IOProcessorNumberSub();
    const Long sz = sizeof(ParticleType)+NumRealComps()*sizeof(ParticleReal)+NumIntComps()*sizeof(int);

#ifdef AMREX_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelReduce::Min(mn,  IOProc, ParallelContext::CommunicatorSub());
            ParallelReduce::Max(mx,  IOProc, ParallelContext::CommunicatorSub());
            ParallelReduce::Sum(cnt, IOProc, ParallelContext::CommunicatorSub());

            amrex::Print() << "ParticleContainer spread across MPI nodes - bytes (num particles): [Min: "
                           << mn*sz
                           << " (" << mn << ")"
                           << ", Max: "
                           << mx*sz
                           << " (" << mx << ")"
                           << ", Total: "
                           << cnt*sz
                           << " (" << cnt << ")]\n";
#ifdef AMREX_LAZY
        });
#endif

    return {mn*sz, mx*sz, cnt*sz};
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
std::array<Long, 3>
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::PrintCapacity () const
{
    Long cnt = 0;

    for (unsigned lev = 0; lev < m_particles.size(); lev++) {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& ptile = kv.second;
            cnt += ptile.capacity();
        }
    }

    Long mn = cnt, mx = mn;

    const int IOProc = ParallelContext::IOProcessorNumberSub();

#ifdef AMREX_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelReduce::Min(mn,  IOProc, ParallelContext::CommunicatorSub());
            ParallelReduce::Max(mx,  IOProc, ParallelContext::CommunicatorSub());
            ParallelReduce::Sum(cnt, IOProc, ParallelContext::CommunicatorSub());

            amrex::Print() << "ParticleContainer spread across MPI nodes - bytes: [Min: "
                           << mn
                           << ", Max: "
                           << mx
                           << ", Total: "
                           << cnt
                           << "]\n";
#ifdef AMREX_LAZY
        });
#endif

    return {mn, mx, cnt};
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::ShrinkToFit ()
{
    for (unsigned lev = 0; lev < m_particles.size(); lev++) {
        auto& pmap = m_particles[lev];
        for (auto& kv : pmap) {
            auto& ptile = kv.second;
            ptile.shrink_to_fit();
        }
    }
}

/**
 * Adds the number of particles in each cell to the values currently located in
 * the input MultiFab.
 */


template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::Increment (MultiFab& mf, int lev)
{
    BL_PROFILE("ParticleContainer::Increment");

    AMREX_ASSERT(OK());
    if (m_particles.empty()) { return; }
    AMREX_ASSERT(lev >= 0 && lev < int(m_particles.size()));
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto domain = geom.Domain();
    amrex::ParticleToMesh(*this, mf, lev,
        [=] AMREX_GPU_DEVICE (const typename ParticleTileType::ConstParticleTileDataType& ptd, int ip,
                              amrex::Array4<amrex::Real> const& count)
        {
            const auto p = make_particle<ConstParticleType>{}(ptd, ip);
            CellAssignor assignor;
            IntVect iv = assignor(p, plo, dxi, domain);
            amrex::Gpu::Atomic::AddNoRet(&count(iv), 1.0_rt);
        }, false);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
Long
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::IncrementWithTotal (MultiFab& mf, int lev, bool local)
{
    BL_PROFILE("ParticleContainer::IncrementWithTotal(lev)");
    Increment(mf, lev);
    return TotalNumberOfParticles(true, local);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::RemoveParticlesAtLevel (int level)
{
    BL_PROFILE("ParticleContainer::RemoveParticlesAtLevel()");
    if (level >= int(this->m_particles.size())) { return; }

    if (!this->m_particles[level].empty())
    {
        ParticleLevel().swap(this->m_particles[level]);
    }
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::RemoveParticlesNotAtFinestLevel ()
{
  BL_PROFILE("ParticleContainer::RemoveParticlesNotAtFinestLevel()");
  AMREX_ASSERT(this->finestLevel()+1 == int(this->m_particles.size()));

  Long cnt = 0;

  for (unsigned lev = 0; lev < m_particles.size() - 1; ++lev) {
      auto& pmap = m_particles[lev];
      if (!pmap.empty()) {
          for (auto& kv : pmap) {
              const auto& pbx = kv.second;
              cnt += pbx.numParticles();
          }
          ParticleLevel().swap(pmap);
      }
  }

  //
  // Print how many particles removed on each processor if any were removed.
  //
  if (this->m_verbose > 1 && cnt > 0) {
      amrex::AllPrint() << "Processor " << ParallelContext::MyProcSub() << " removed " << cnt
                        << " particles not in finest level\n";
  }
}

struct FilterVirt
{

    amrex::AssignGrid<amrex::DenseBinIteratorFactory<amrex::Box>> m_assign_buffer_grid;
    GpuArray<Real, AMREX_SPACEDIM> m_plo, m_dxi;
    Box m_domain;

    FilterVirt (const amrex::AssignGrid<amrex::DenseBinIteratorFactory<amrex::Box>>& assign_buffer_grid, const GpuArray<Real, AMREX_SPACEDIM>& plo,
                const GpuArray<Real, AMREX_SPACEDIM> & dxi, const Box& domain)
      : m_assign_buffer_grid(assign_buffer_grid), m_plo(plo), m_dxi(dxi), m_domain(domain)
    {}

    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int src_i) const noexcept
    {
        auto iv = getParticleCell(src.m_aos[src_i], m_plo, m_dxi, m_domain);
        return (m_assign_buffer_grid(iv)!=-1);
    }
};

struct TransformerVirt
{
    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst, const SrcData& src,
                     int src_i, int dst_i) const noexcept
    {
        copyParticle(dst, src, src_i, dst_i);

        (dst.m_aos[dst_i]).id() = LongParticleIds::VirtualParticleID;
        (dst.m_aos[dst_i]).cpu() = 0;
    }
};


template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::CreateVirtualParticles (int level, AoS& virts) const
{
    ParticleTileType ptile;
    CreateVirtualParticles(level, ptile);
    ptile.GetArrayOfStructs().swap(virts);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::CreateVirtualParticles (int level, ParticleTileType& virts) const
{
    BL_PROFILE("ParticleContainer::CreateVirtualParticles()");
    AMREX_ASSERT(level > 0);
    AMREX_ASSERT(virts.empty());

    if (level >= static_cast<int>(m_particles.size())) {
        return;
    }

    std::string aggregation_type = AggregationType();
    int aggregation_buffer = AggregationBuffer();

    if (aggregation_type == "None")
    {
        auto virts_offset = virts.numParticles();
        for(ParConstIterType pti(*this, level); pti.isValid(); ++pti)
        {
            const auto& src_tile = ParticlesAt(level, pti);

            auto np = src_tile.numParticles();
            virts.resize(virts_offset+np);
            transformParticles(virts, src_tile, 0, virts_offset, np, TransformerVirt());
            virts_offset += np;
        }
    }
    if (aggregation_type == "Cell")
    {
        //Components would be based on
        int nComp = AMREX_SPACEDIM + NStructReal + NArrayReal;
        // NArrayReal, NStructInt, NArrayInt behavior as before
        int nGhost = 0;
        MultiFab mf(ParticleBoxArray(level), ParticleDistributionMap(level), nComp, nGhost);

        nComp = 1 + NStructInt + NArrayInt;
        iMultiFab imf(ParticleBoxArray(level), ParticleDistributionMap(level), nComp, nGhost);

        const auto& geom = Geom(level);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto domain = geom.Domain();

        BoxList bl_buffer;
        bl_buffer.complementIn(Geom(level).Domain(), ParticleBoxArray(level));
        BoxArray buffer(std::move(bl_buffer));
        buffer.grow(aggregation_buffer);

        ParticleLocator<DenseBins<Box>> locator;
        locator.build(buffer, geom);
        AssignGrid<DenseBinIteratorFactory<Box>> assign_buffer_grid = locator.getGridAssignor();

        amrex::ParticleToMesh(*this, mf, level,
        [=] AMREX_GPU_DEVICE (const ParticleType& p,
                              amrex::Array4<amrex::Real> const& partData,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo_loc,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi_loc)
        {
            auto iv = getParticleCell(p, plo_loc, dxi_loc, domain);
            if(assign_buffer_grid(iv)==-1)
            {
               //Ordering may make this not unique
               for (int i = 0; i < NArrayReal; ++i)
               {
                   amrex::Gpu::Atomic::AddNoRet(&partData(iv,AMREX_SPACEDIM+NStructReal+i), partData(iv,AMREX_SPACEDIM)!=0.0 ? static_cast<Real>(0) : static_cast<Real>(p.rdata(NStructReal+i)));
               }
               //Add the rdata(0)-weighted sum of position
               for (int i = 0; i < AMREX_SPACEDIM; ++i)
               {
                   amrex::Gpu::Atomic::AddNoRet(&partData(iv,i), static_cast<Real>((p.rdata(0)*p.pos(i))));
               }
               //Add the rdata(0)-weighted sum of other rdata fields
               for (int i = 1; i < NStructReal; ++i)
               {
                   amrex::Gpu::Atomic::AddNoRet(&partData(iv,AMREX_SPACEDIM+i), static_cast<Real>((p.rdata(0)*p.rdata(i))));
               }
               //Add the rdata(0) sum
               for (int i = 0; i < 1; ++i)
               {
                   amrex::Gpu::Atomic::AddNoRet(&partData(iv,AMREX_SPACEDIM+i), static_cast<Real>(p.rdata(0)));
               }
            }

        }); //skipping extra false argument, doing mf.setVal(0) at beginning

        amrex::ParticleToMesh(*this, imf, level,
        [=] AMREX_GPU_DEVICE (const ParticleType& p,
                              amrex::Array4<int> const& partData,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo_loc,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi_loc)
        {

            auto iv = getParticleCell(p, plo_loc, dxi_loc, domain);
            if(assign_buffer_grid(iv)==-1)
            {
                //if this cell has no particle id info, do a straight copy to store idata
                if(partData(iv,0)==0)
                {
                    //Add 1 to indicate at least 1 particle at cell iv
                    amrex::Gpu::Atomic::AddNoRet(&partData(iv,0), 1);
                    for (int i = 0; i < NStructInt; ++i)
                    {
                        amrex::Gpu::Atomic::AddNoRet(&partData(iv,1+i), p.idata(i));
                    }
                    for (int i = 0; i < NArrayInt; ++i)
                    {
                        amrex::Gpu::Atomic::AddNoRet(&partData(iv,1+NStructInt+i), p.idata(NStructInt+i));
                    }
                }
            }
        });

    //There may be a better way to ensure virts is the right length
    virts.resize(imf.sum(0));

    int last_offset = 0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const auto bx = mfi.tilebox();
        const auto partData = mf.array(mfi);
        const auto imf_arr = imf.array(mfi);

        Gpu::DeviceVector<int> offsets(bx.numPts());
        auto *offsets_ptr = offsets.dataPtr();
        int next_offset = Scan::ExclusiveSum((int) bx.numPts(),(imf_arr.ptr(bx.smallEnd(),0)),(offsets.dataPtr()),Scan::retSum);
        auto dst = virts.getParticleTileData();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
             if(imf_arr(i,j,k,0)!=0)
             {
                 SuperParticleType p;
                 p.cpu() = 0;
                 p.id() = LongParticleIds::VirtualParticleID;

                 //Set rdata(0) first so we can normalize the weighted fields
                 p.rdata(0) = static_cast<ParticleReal>(partData(i,j,k,AMREX_SPACEDIM));
                 //Set pos with the normalized weighted field
                 for (int n = 0; n < AMREX_SPACEDIM; ++n)
                 {
                     p.pos(n) = static_cast<ParticleReal>(partData(i,j,k,n) / p.rdata(0));
                 }
                 //Set rdata(n>0) with the normalized weighted field for NStructReal
                 for (int n = 1; n < NStructReal; ++n)
                 {
                     p.rdata(n) = static_cast<ParticleReal>(partData(i,j,k,AMREX_SPACEDIM+n) / p.rdata(0));
                 }
                 //Set rdata(n>0) with the normalized weighted field for NArrayReal
                 for (int n = 0; n < NArrayReal; ++n)
                 {
                     p.rdata(NStructReal+n) = static_cast<ParticleReal>(partData(i,j,k,AMREX_SPACEDIM+NStructReal+n));
                 }
                 //Set idata with the "first" particles idata field for NStructInt
                 for (int n = 0; n < NStructInt; ++n)
                 {
                     p.idata(n) = imf_arr(i,j,k,1+n);
                 }
                 //Set idata with the "first" particles idata field for NArrayInt
                 for (int n = 0; n < NArrayInt; ++n)
                 {
                     p.idata(NStructInt+n) = imf_arr(i,j,k,1+NStructInt+n);
                 }
                 //Set the new particle in dst tile
                 dst.setSuperParticle(p, last_offset+offsets_ptr[((i-imf_arr.begin.x)+(j-imf_arr.begin.y)*imf_arr.jstride+(k-imf_arr.begin.z)*imf_arr.kstride)]);
             }

        });
        last_offset+=next_offset;
        Gpu::streamSynchronize();
    }

    // last_offset should equal virts.numParticles()
    auto virts_offset = last_offset;
    for(ParConstIterType pti(*this, level); pti.isValid(); ++pti)
    {
        const auto& src_tile = ParticlesAt(level, pti);

        auto np = src_tile.numParticles();
        virts.resize(virts_offset+np);
        virts_offset += filterAndTransformParticles(virts, src_tile, FilterVirt(assign_buffer_grid,plo,dxi,domain), TransformerVirt(),0,virts_offset);
        Gpu::streamSynchronize();
    }
    virts.resize(virts_offset);
    Gpu::streamSynchronize();
    }
}

struct AssignGridFilter
{

    int m_lev_min, m_lev_max, m_nGrow, m_gid;
    amrex::AmrAssignGrid<amrex::DenseBinIteratorFactory<amrex::Box>> m_assign_grid;

    /**
       \brief This filters based on matching grids
       *
     **/
    AssignGridFilter (amrex::AmrAssignGrid<amrex::DenseBinIteratorFactory<amrex::Box>> assign_grid, int gid, int level, int nGrow)
      : m_lev_min(level), m_lev_max(level+1), m_nGrow(nGrow), m_gid(gid), m_assign_grid(assign_grid)
    {}

    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int src_i) const noexcept
    {
        const auto tup_min = (m_assign_grid)(src.m_aos[src_i], m_lev_min, m_lev_max, m_nGrow, DefaultAssignor{});
        const auto tup_max = (m_assign_grid)(src.m_aos[src_i], m_lev_max, m_lev_max, m_nGrow, DefaultAssignor{});
        const auto p_boxes = amrex::get<0>(tup_min);
        const auto p_boxes_max = amrex::get<0>(tup_max);
        const auto p_levs_max  = amrex::get<1>(tup_max);
        return p_boxes_max >=0 && p_boxes == m_gid && p_levs_max == m_lev_max;
    }
};

struct TransformerGhost
{

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator() (DstData& dst, const SrcData& src,
                     int src_i, int dst_i) const noexcept
    {
        copyParticle(dst, src, src_i, dst_i);

        (dst.m_aos[dst_i]).id() = LongParticleIds::GhostParticleID;
        (dst.m_aos[dst_i]).cpu() = 0;
    }
};

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::CreateGhostParticles (int level, int nGrow, AoS& ghosts) const
{
    ParticleTileType ptile;
    CreateGhostParticles(level, nGrow, ptile);
    ptile.GetArrayOfStructs().swap(ghosts);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::CreateGhostParticles (int level, int nGrow, ParticleTileType& ghosts) const
{
    BL_PROFILE("ParticleContainer::CreateGhostParticles()");
    AMREX_ASSERT(ghosts.empty());
    AMREX_ASSERT(level < finestLevel());

    if (level >= static_cast<int>(m_particles.size())) {
        return;
    }

    if (! m_particle_locator.isValid(GetParGDB())) {
        m_particle_locator.build(GetParGDB());
    }

    m_particle_locator.setGeometry(GetParGDB());
    AmrAssignGrid<DenseBinIteratorFactory<Box>> assign_grid = m_particle_locator.getGridAssignor();
    auto ghost_offset = ghosts.numParticles();
    for(ParConstIterType pti(*this, level); pti.isValid(); ++pti)
    {
        const auto& src_tile = ParticlesAt(level, pti);
        int gid = pti.index();

        auto np = src_tile.numParticles();
        ghosts.resize(ghost_offset+np);
        ghost_offset += filterAndTransformParticles(ghosts, src_tile, AssignGridFilter(assign_grid,gid,level,nGrow), TransformerGhost(),0,ghost_offset);
    }
    ghosts.resize(ghost_offset);
    Gpu::streamSynchronize();
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
clearParticles ()
{
    BL_PROFILE("ParticleContainer::clearParticles()");

    for (int lev = 0; lev < static_cast<int>(m_particles.size()); ++lev)
    {
        for (auto& kv : m_particles[lev]) { kv.second.resize(0); }
        particle_detail::clearEmptyEntries(m_particles[lev]);
    }
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <class PCType, std::enable_if_t<IsParticleContainer<PCType>::value, int> foo>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
copyParticles (const PCType& other, bool local)
{
    using PData = ConstParticleTileData<typename ParticleType::StorageParticleType, NArrayReal, NArrayInt>;
    copyParticles(other, [] AMREX_GPU_HOST_DEVICE (const PData& /*data*/, int /*i*/) { return 1; }, local);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <class PCType, std::enable_if_t<IsParticleContainer<PCType>::value, int> foo>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
addParticles (const PCType& other, bool local)
{
    using PData = ConstParticleTileData<typename ParticleType::StorageParticleType, NArrayReal, NArrayInt>;
    addParticles(other, [] AMREX_GPU_HOST_DEVICE (const PData& /*data*/, int /*i*/) { return 1; }, local);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <class F, class PCType,
          std::enable_if_t<IsParticleContainer<PCType>::value, int> foo,
          std::enable_if_t<! std::is_integral_v<F>, int> bar>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
copyParticles (const PCType& other, F&& f, bool local)
{
    BL_PROFILE("ParticleContainer::copyParticles");
    clearParticles();
    addParticles(other, std::forward<F>(f), local);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <class F, class PCType,
          std::enable_if_t<IsParticleContainer<PCType>::value, int> foo,
          std::enable_if_t<! std::is_integral_v<F>, int> bar>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
addParticles (const PCType& other, F const& f, bool local)
{
    BL_PROFILE("ParticleContainer::addParticles");

    for (int lev = 0; lev < other.numLevels(); ++lev)
    {
        const auto& plevel_other = other.GetParticles(lev);
        for(MFIter mfi = other.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
            if(plevel_other.find(index) == plevel_other.end()) { continue; }

            auto& ptile = DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
            const auto& ptile_other = plevel_other.at(index);
            auto np = ptile_other.numParticles();
            if (np == 0) { continue; }

            auto dst_index = ptile.numParticles();
            ptile.resize(dst_index + np);

            auto count = amrex::filterParticles(ptile, ptile_other, f, 0, dst_index, np);

            ptile.resize(dst_index + count);
        }
    }

    if (! local) { Redistribute(); }
}

//
// This redistributes valid particles and discards invalid ones.
//
template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::Redistribute (int lev_min, int lev_max, int nGrow, int local, bool remove_negative)
{
    BL_PROFILE_SYNC_START_TIMED("SyncBeforeComms: Redist");

#ifdef AMREX_USE_GPU
    if ( Gpu::inLaunchRegion() )
    {
        RedistributeGPU(lev_min, lev_max, nGrow, local, remove_negative);
    }
    else
    {
        RedistributeCPU(lev_min, lev_max, nGrow, local, remove_negative);
    }
#else
    RedistributeCPU(lev_min, lev_max, nGrow, local, remove_negative);
#endif

    BL_PROFILE_SYNC_STOP();
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
template <class index_type>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::ReorderParticles (int lev, const MFIter& mfi, const index_type* permutations)
{
    auto& ptile           = ParticlesAt(lev, mfi);
    const size_t np       = ptile.numParticles();
    const size_t np_total = np + ptile.numNeighborParticles();

    if (memEfficientSort) {
#if defined(AMREX_USE_CUDA) && defined(_WIN32)
        if           (!ParticleType::is_soa_particle) {
#else
        if constexpr (!ParticleType::is_soa_particle) {
#endif
            static_assert(sizeof(ParticleType)%4 == 0 && sizeof(uint32_t) == 4);
            using tmp_t = std::conditional_t<sizeof(ParticleType)%8 == 0,
                                             uint64_t, uint32_t>;
            constexpr std::size_t nchunks = sizeof(ParticleType) / sizeof(tmp_t);
            Gpu::DeviceVector<tmp_t> tmp(np);
            auto* ptmp = tmp.data();
            auto* paos = (tmp_t*)(ptile.getParticleTileData().m_aos);
            for (std::size_t ichunk = 0; ichunk < nchunks; ++ichunk) {
                // Do not need to reorder neighbor particles
                AMREX_HOST_DEVICE_FOR_1D(np, i,
                {
                    ptmp[i] = paos[permutations[i]*nchunks+ichunk];
                });
                AMREX_HOST_DEVICE_FOR_1D(np, i,
                {
                    paos[i*nchunks+ichunk] = ptmp[i];
                });
            }
            Gpu::streamSynchronize();
        } else {
            typename SoA::IdCPU tmp_idcpu(np_total);

            auto src = ptile.GetStructOfArrays().GetIdCPUData().data();
            uint64_t* dst = tmp_idcpu.data();
            AMREX_HOST_DEVICE_FOR_1D( np_total, i,
            {
                dst[i] = i < np ? src[permutations[i]] : src[i];
            });

            Gpu::streamSynchronize();

            ptile.GetStructOfArrays().GetIdCPUData().swap(tmp_idcpu);
        }

        { // Create a scope for the temporary vector below
            RealVector tmp_real(np_total);
            for (int comp = 0; comp < NArrayReal + m_num_runtime_real; ++comp) {
                auto src = ptile.GetStructOfArrays().GetRealData(comp).data();
                ParticleReal* dst = tmp_real.data();
                AMREX_HOST_DEVICE_FOR_1D( np_total, i,
                {
                    dst[i] = i < np ? src[permutations[i]] : src[i];
                });

                Gpu::streamSynchronize();

                ptile.GetStructOfArrays().GetRealData(comp).swap(tmp_real);
            }
        }

        IntVector tmp_int(np_total);
        for (int comp = 0; comp < NArrayInt + m_num_runtime_int; ++comp) {
            auto src = ptile.GetStructOfArrays().GetIntData(comp).data();
            int* dst = tmp_int.data();
            AMREX_HOST_DEVICE_FOR_1D( np_total , i,
            {
                dst[i] = i < np ? src[permutations[i]] : src[i];
            });

            Gpu::streamSynchronize();

            ptile.GetStructOfArrays().GetIntData(comp).swap(tmp_int);
        }
    } else {
        ParticleTileType ptile_tmp;
        ptile_tmp.define(m_num_runtime_real, m_num_runtime_int, &m_soa_rdata_names, &m_soa_idata_names);
        ptile_tmp.resize(np_total);
        // copy re-ordered particles
        gatherParticles(ptile_tmp, ptile, np, permutations);
        // copy neighbor particles
        amrex::copyParticles(ptile_tmp, ptile, np, np, np_total-np);
        ptile.swap(ptile_tmp);
    }
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::SortParticlesByCell ()
{
    SortParticlesByBin(IntVect(AMREX_D_DECL(1, 1, 1)));
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::SortParticlesByBin (IntVect bin_size)
{
    BL_PROFILE("ParticleContainer::SortParticlesByBin()");

    if (bin_size == IntVect::TheZeroVector()) { return; }

    for (int lev = 0; lev < numLevels(); ++lev)
    {
        const Geometry& geom = Geom(lev);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();
        const auto domain = geom.Domain();

        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto& ptile           = ParticlesAt(lev, mfi);
            const size_t np       = ptile.numParticles();

            const Box& box = mfi.validbox();

            int ntiles = numTilesInBox(box, true, bin_size);

            m_bins.build(np, ptile.getParticleTileData(), ntiles,
                         GetParticleBin{plo, dxi, domain, bin_size, box});
            ReorderParticles(lev, mfi, m_bins.permutationPtr());
        }
    }
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::SortParticlesForDeposition (IntVect idx_type)
{
    BL_PROFILE("ParticleContainer::SortParticlesForDeposition()");

    for (int lev = 0; lev < numLevels(); ++lev)
    {
        const Geometry& geom = Geom(lev);

        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            const auto& ptile = ParticlesAt(lev, mfi);
            const size_t np   = ptile.numParticles();

            const Box& box = mfi.validbox();

            using index_type = typename decltype(m_bins)::index_type;
            Gpu::DeviceVector<index_type> perm;
            PermutationForDeposition<index_type>(perm, np, ptile, box, geom, idx_type);
            ReorderParticles(lev, mfi, perm.dataPtr());
        }
    }
}

//
// The GPU implementation of Redistribute
//
template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::RedistributeGPU (int lev_min, int lev_max, int nGrow, int local, bool remove_negative)
{
#ifdef AMREX_USE_GPU

    if (local) { AMREX_ASSERT(numParticlesOutOfRange(*this, lev_min, lev_max, local) == 0); }

    // sanity check
    AMREX_ALWAYS_ASSERT(do_tiling == false);

    BL_PROFILE("ParticleContainer::RedistributeGPU()");
    BL_PROFILE_VAR_NS("Redistribute_partition", blp_partition);

    int theEffectiveFinestLevel = m_gdb->finestLevel();
    while (!m_gdb->LevelDefined(theEffectiveFinestLevel)) { theEffectiveFinestLevel--; }

    if (int(m_particles.size()) < theEffectiveFinestLevel+1) {
        if (Verbose()) {
            amrex::Print() << "ParticleContainer::Redistribute() resizing containers from "
                           << m_particles.size() << " to "
                           << theEffectiveFinestLevel + 1 << '\n';
        }
        m_particles.resize(theEffectiveFinestLevel+1);
        m_dummy_mf.resize(theEffectiveFinestLevel+1);
    }

    for (int lev = 0; lev < theEffectiveFinestLevel+1; ++lev) { RedefineDummyMF(lev); }

    int finest_lev_particles;
    if (lev_max == -1) {
        lev_max = theEffectiveFinestLevel;
        finest_lev_particles = m_particles.size() - 1;
    } else {
        finest_lev_particles = lev_max;
    }
    AMREX_ASSERT(lev_max <= finestLevel());

    this->defineBufferMap();

    if (! m_particle_locator.isValid(GetParGDB())) { m_particle_locator.build(GetParGDB()); }
    m_particle_locator.setGeometry(GetParGDB());
    auto assign_grid = m_particle_locator.getGridAssignor();

    BL_PROFILE_VAR_START(blp_partition);
    ParticleCopyOp op;
    int num_levels = finest_lev_particles + 1;
    op.setNumLevels(num_levels);
    Vector<std::map<int, int> > new_sizes(num_levels);
    const auto plo    = Geom(0).ProbLoArray();
    const auto phi    = Geom(0).ProbHiArray();
    const auto rlo    = Geom(0).ProbLoArrayInParticleReal();
    const auto rhi    = Geom(0).ProbHiArrayInParticleReal();
    const auto is_per = Geom(0).isPeriodicArray();
    for (int lev = lev_min; lev <= finest_lev_particles; ++lev)
    {
        auto& plev = m_particles[lev];
        for (auto& kv : plev)
        {
            int gid = kv.first.first;
            int tid = kv.first.second;
            auto index = std::make_pair(gid, tid);

            auto& src_tile = plev[index];
            const size_t np = src_tile.numParticles();

            int num_stay = partitionParticlesByDest(src_tile, assign_grid,
                                                    std::forward<CellAssignor>(CellAssignor{}),
                                                    BufferMap(),
                                                    plo, phi, rlo, rhi, is_per, lev, gid, tid,
                                                    lev_min, lev_max, nGrow, remove_negative);

            int num_move = np - num_stay;
            new_sizes[lev][gid] = num_stay;
            op.resize(gid, lev, num_move);

            auto p_boxes = op.m_boxes[lev][gid].dataPtr();
            auto p_levs = op.m_levels[lev][gid].dataPtr();
            auto p_src_indices = op.m_src_indices[lev][gid].dataPtr();
            auto p_periodic_shift = op.m_periodic_shift[lev][gid].dataPtr();
            auto ptd = src_tile.getParticleTileData();

            AMREX_FOR_1D ( num_move, i,
            {
                const auto p = make_particle<ParticleType>{}(ptd,i + num_stay);

                if (p.id() < 0)
                {
                    p_boxes[i] = -1;
                    p_levs[i]  = -1;
                }
                else
                {
                    const auto tup = assign_grid(p, lev_min, lev_max, nGrow,
                                                 std::forward<CellAssignor>(CellAssignor{}));
                    p_boxes[i] = amrex::get<0>(tup);
                    p_levs[i]  = amrex::get<1>(tup);
                }
                p_periodic_shift[i] = IntVect(AMREX_D_DECL(0,0,0));
                p_src_indices[i] = i+num_stay;
            });
        }
    }
    BL_PROFILE_VAR_STOP(blp_partition);

    ParticleCopyPlan plan;

    plan.build(*this, op, h_redistribute_int_comp,
               h_redistribute_real_comp, local);

    // by default, this uses The_Arena();
    amrex::PODVector<char, PolymorphicArenaAllocator<char> > snd_buffer;
    amrex::PODVector<char, PolymorphicArenaAllocator<char> > rcv_buffer;

    if (use_comms_arena) {
        snd_buffer.setArena(The_Comms_Arena());
        rcv_buffer.setArena(The_Comms_Arena());
    }

    packBuffer(*this, op, plan, snd_buffer);

    // clear particles from container
    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        auto& plev = m_particles[lev];
        for (auto& kv : plev)
        {
            int gid = kv.first.first;
            int tid = kv.first.second;
            auto index = std::make_pair(gid, tid);
            auto& tile = plev[index];
            tile.resize(new_sizes[lev][gid]);
        }
    }

    for (int lev = lev_min; lev <= lev_max; lev++)
    {
        particle_detail::clearEmptyEntries(m_particles[lev]);
    }

    if (int(m_particles.size()) > theEffectiveFinestLevel+1) {
        if (m_verbose > 0) {
            amrex::Print() << "ParticleContainer::Redistribute() resizing m_particles from "
                        << m_particles.size() << " to " << theEffectiveFinestLevel+1 << '\n';
        }
        AMREX_ASSERT(int(m_particles.size()) >= 2);

        m_particles.resize(theEffectiveFinestLevel + 1);
        m_dummy_mf.resize(theEffectiveFinestLevel + 1);
    }

    if (ParallelDescriptor::UseGpuAwareMpi())
    {
        plan.buildMPIFinish(BufferMap());
        communicateParticlesStart(*this, plan, snd_buffer, rcv_buffer);
        unpackBuffer(*this, plan, snd_buffer, RedistributeUnpackPolicy());
        communicateParticlesFinish(plan);
        unpackRemotes(*this, plan, rcv_buffer, RedistributeUnpackPolicy());
    }
    else
    {
        Gpu::Device::streamSynchronize();
        Gpu::PinnedVector<char> pinned_snd_buffer;
        Gpu::PinnedVector<char> pinned_rcv_buffer;

        if (snd_buffer.arena()->isPinned()) {
            plan.buildMPIFinish(BufferMap());
            Gpu::Device::streamSynchronize();
            communicateParticlesStart(*this, plan, snd_buffer, pinned_rcv_buffer);
        } else {
            pinned_snd_buffer.resize(snd_buffer.size());
            Gpu::dtoh_memcpy_async(pinned_snd_buffer.dataPtr(), snd_buffer.dataPtr(), snd_buffer.size());
            plan.buildMPIFinish(BufferMap());
            Gpu::Device::streamSynchronize();
            communicateParticlesStart(*this, plan, pinned_snd_buffer, pinned_rcv_buffer);
        }

        rcv_buffer.resize(pinned_rcv_buffer.size());
        unpackBuffer(*this, plan, snd_buffer, RedistributeUnpackPolicy());
        communicateParticlesFinish(plan);
        Gpu::htod_memcpy_async(rcv_buffer.dataPtr(), pinned_rcv_buffer.dataPtr(), pinned_rcv_buffer.size());
        unpackRemotes(*this, plan, rcv_buffer, RedistributeUnpackPolicy());
    }

    Gpu::Device::streamSynchronize();
    AMREX_ASSERT(numParticlesOutOfRange(*this, lev_min, lev_max, nGrow) == 0);
#else
    amrex::ignore_unused(lev_min,lev_max,nGrow,local,remove_negative);
#endif
}

//
// The CPU implementation of Redistribute
//
template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::RedistributeCPU (int lev_min, int lev_max, int nGrow, int local, bool remove_negative)
{
    BL_PROFILE("ParticleContainer::RedistributeCPU()");

    const int MyProc    = ParallelContext::MyProcSub();
    auto      strttime  = amrex::second();

    if (local > 0) { BuildRedistributeMask(0, local); }

    // On startup there are cases where Redistribute() could be called
    // with a given finestLevel() where that AmrLevel has yet to be defined.
    int theEffectiveFinestLevel = m_gdb->finestLevel();

    while (!m_gdb->LevelDefined(theEffectiveFinestLevel)) {
        theEffectiveFinestLevel--;
    }

    if (int(m_particles.size()) < theEffectiveFinestLevel+1) {
        if (Verbose()) {
            amrex::Print() << "ParticleContainer::Redistribute() resizing containers from "
                           << m_particles.size() << " to "
                           << theEffectiveFinestLevel + 1 << '\n';
        }
        m_particles.resize(theEffectiveFinestLevel+1);
        m_dummy_mf.resize(theEffectiveFinestLevel+1);
    }

    // It is important to do this even if we don't have more levels because we may have changed the
    // grids at this level in a regrid.
    for (int lev = 0; lev < theEffectiveFinestLevel+1; ++lev) {
        RedefineDummyMF(lev);
    }

    int finest_lev_particles;
    if (lev_max == -1) {
        lev_max = theEffectiveFinestLevel;
        finest_lev_particles = m_particles.size() - 1;
    } else {
        finest_lev_particles = lev_max;
    }
    AMREX_ASSERT(lev_max <= finestLevel());

    // This will hold the valid particles that go to another process
    std::map<int, Vector<char> > not_ours;

    int num_threads = OpenMP::get_max_threads();

    // these are temporary buffers for each thread
    std::map<int, Vector<Vector<char> > > tmp_remote;
    Vector<std::map<std::pair<int, int>, Vector<ParticleVector> > > tmp_local;
    Vector<std::map<std::pair<int, int>, Vector<StructOfArrays<NArrayReal, NArrayInt, Allocator> > > > soa_local;
    tmp_local.resize(theEffectiveFinestLevel+1);
    soa_local.resize(theEffectiveFinestLevel+1);

    // we resize these buffers outside the parallel region
    for (int lev = lev_min; lev <= lev_max; lev++) {
        for (MFIter mfi(*m_dummy_mf[lev], this->do_tiling ? this->tile_size : IntVect::TheZeroVector());
             mfi.isValid(); ++mfi) {
            auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
            tmp_local[lev][index].resize(num_threads);
            soa_local[lev][index].resize(num_threads);
            for (int t = 0; t < num_threads; ++t) {
                soa_local[lev][index][t].define(m_num_runtime_real, m_num_runtime_int, &m_soa_rdata_names, &m_soa_idata_names);
            }
        }
    }
    if (local) {
        for (int i = 0; i < neighbor_procs.size(); ++i) {
            tmp_remote[neighbor_procs[i]].resize(num_threads);
        }
    } else {
        for (int i = 0; i < ParallelContext::NProcsSub(); ++i) {
            tmp_remote[i].resize(num_threads);
        }
    }

    // first pass: for each tile in parallel, in each thread copies the particles that
    // need to be moved into it's own, temporary buffer.
    for (int lev = lev_min; lev <= finest_lev_particles; lev++) {
        auto& pmap = m_particles[lev];

        Vector<std::pair<int, int> > grid_tile_ids;
        Vector<ParticleTileType*> ptile_ptrs;
        for (auto& kv : pmap)
        {
            grid_tile_ids.push_back(kv.first);
            ptile_ptrs.push_back(&(kv.second));
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
        for (int pmap_it = 0; pmap_it < static_cast<int>(ptile_ptrs.size()); ++pmap_it)
        {
            int thread_num = OpenMP::get_thread_num();
            int grid = grid_tile_ids[pmap_it].first;
            int tile = grid_tile_ids[pmap_it].second;
            auto& soa = ptile_ptrs[pmap_it]->GetStructOfArrays();
            auto& aos = ptile_ptrs[pmap_it]->GetArrayOfStructs();

            // AMREX_ASSERT_WITH_MESSAGE((NumRealComps() == 0 && NumIntComps() == 0)
            //                           || aos.size() == soa.size(),
            //     "The AoS and SoA data on this tile are different sizes - "
            //     "perhaps particles have not been initialized correctly?");
            unsigned npart = ptile_ptrs[pmap_it]->numParticles();
            ParticleLocData pld;

            if constexpr (!ParticleType::is_soa_particle){

                if (npart != 0) {
                    Long last = npart - 1;
                    Long pindex = 0;
                    while (pindex <= last) {
                        ParticleType& p = aos[pindex];

                        if ((remove_negative == false) && (p.id() < 0)) {
                            ++pindex;
                            continue;
                        }

                        if (p.id() < 0)
                        {
                            aos[pindex] = aos[last];
                            for (int comp = 0; comp < NumRealComps(); comp++) {
                                soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                            }
                            for (int comp = 0; comp < NumIntComps(); comp++) {
                                soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                            }
                            correctCellVectors(last, pindex, grid, aos[pindex]);
                            --last;
                            continue;
                        }

                        locateParticle(p, pld, lev_min, lev_max, nGrow, local ? grid : -1);

                        particlePostLocate(p, pld, lev);

                        if (p.id() < 0)
                        {
                            aos[pindex] = aos[last];
                            for (int comp = 0; comp < NumRealComps(); comp++) {
                                soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                            }
                            for (int comp = 0; comp < NumIntComps(); comp++) {
                                soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                            }
                            correctCellVectors(last, pindex, grid, aos[pindex]);
                            --last;
                            continue;
                        }

                        const int who = ParallelContext::global_to_local_rank(ParticleDistributionMap(pld.m_lev)[pld.m_grid]);
                        if (who == MyProc) {
                            if (pld.m_lev != lev || pld.m_grid != grid || pld.m_tile != tile) {
                                // We own it but must shift it to another place.
                                auto index = std::make_pair(pld.m_grid, pld.m_tile);
                                AMREX_ASSERT(tmp_local[pld.m_lev][index].size() == num_threads);
                                tmp_local[pld.m_lev][index][thread_num].push_back(p);
                                for (int comp = 0; comp < NumRealComps(); ++comp) {
                                    RealVector& arr = soa_local[pld.m_lev][index][thread_num].GetRealData(comp);
                                    arr.push_back(soa.GetRealData(comp)[pindex]);
                                }
                                for (int comp = 0; comp < NumIntComps(); ++comp) {
                                    IntVector& arr = soa_local[pld.m_lev][index][thread_num].GetIntData(comp);
                                    arr.push_back(soa.GetIntData(comp)[pindex]);
                                }

                                p.id() = -p.id(); // Invalidate the particle
                            }
                        }
                        else {
                            auto& particles_to_send = tmp_remote[who][thread_num];
                            auto old_size = particles_to_send.size();
                            auto new_size = old_size + superparticle_size;
                            particles_to_send.resize(new_size);
                            std::memcpy(&particles_to_send[old_size], &p, particle_size);
                            char* dst = &particles_to_send[old_size] + particle_size;
                            int array_comp_start = 0;
                            if constexpr (!ParticleType::is_soa_particle) {
                                array_comp_start = AMREX_SPACEDIM + NStructReal;
                            }
                            for (int comp = 0; comp < NumRealComps(); comp++) {
                                if (h_redistribute_real_comp[array_comp_start + comp]) {
                                    std::memcpy(dst, &soa.GetRealData(comp)[pindex], sizeof(ParticleReal));
                                    dst += sizeof(ParticleReal);
                                }
                            }
                            array_comp_start = 2 + NStructInt;
                            for (int comp = 0; comp < NumIntComps(); comp++) {
                                if (h_redistribute_int_comp[array_comp_start + comp]) {
                                    std::memcpy(dst, &soa.GetIntData(comp)[pindex], sizeof(int));
                                    dst += sizeof(int);
                                }
                            }

                            p.id() = -p.id(); // Invalidate the particle
                        }

                        if (p.id() < 0)
                        {
                            aos[pindex] = aos[last];
                            for (int comp = 0; comp < NumRealComps(); comp++) {
                                soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                            }
                            for (int comp = 0; comp < NumIntComps(); comp++) {
                                soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                            }
                            correctCellVectors(last, pindex, grid, aos[pindex]);
                            --last;
                            continue;
                        }

                        ++pindex;
                    }

                    aos().erase(aos().begin() + last + 1, aos().begin() + npart);
                    for (int comp = 0; comp < NumRealComps(); comp++) {
                        RealVector& rdata = soa.GetRealData(comp);
                        rdata.erase(rdata.begin() + last + 1, rdata.begin() + npart);
                    }
                    for (int comp = 0; comp < NumIntComps(); comp++) {
                        IntVector& idata = soa.GetIntData(comp);
                        idata.erase(idata.begin() + last + 1, idata.begin() + npart);
                    }
                }

            } else { // soa particle

                auto particle_tile = ptile_ptrs[pmap_it];
                if (npart != 0) {
                    Long last = npart - 1;
                    Long pindex = 0;
                    auto ptd = particle_tile->getParticleTileData();
                    while (pindex <= last) {
                        ParticleType p(ptd,pindex);

                        if ((remove_negative == false) && (p.id() < 0)) {
                            ++pindex;
                            continue;
                        }

                        if (p.id() < 0){
                            soa.GetIdCPUData()[pindex] = soa.GetIdCPUData()[last];
                            for (int comp = 0; comp < NumRealComps(); comp++) {
                                soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                            }
                            for (int comp = 0; comp < NumIntComps(); comp++) {
                                soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                            }
                            correctCellVectors(last, pindex, grid, ptd[pindex]);
                            --last;
                            continue;
                        }

                        locateParticle(p, pld, lev_min, lev_max, nGrow, local ? grid : -1);

                        particlePostLocate(p, pld, lev);

                        if (p.id() < 0) {
                            soa.GetIdCPUData()[pindex] = soa.GetIdCPUData()[last];
                            for (int comp = 0; comp < NumRealComps(); comp++) {
                                soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                            }
                            for (int comp = 0; comp < NumIntComps(); comp++) {
                                soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                            }
                            correctCellVectors(last, pindex, grid, ptd[pindex]);
                            --last;
                            continue;
                        }

                        const int who = ParallelContext::global_to_local_rank(ParticleDistributionMap(pld.m_lev)[pld.m_grid]);
                        if (who == MyProc) {
                            if (pld.m_lev != lev || pld.m_grid != grid || pld.m_tile != tile) {
                                // We own it but must shift it to another place.
                                auto index = std::make_pair(pld.m_grid, pld.m_tile);
                                AMREX_ASSERT(soa_local[pld.m_lev][index].size() == num_threads);
                                {
                                    auto& arr = soa_local[pld.m_lev][index][thread_num].GetIdCPUData();
                                    arr.push_back(soa.GetIdCPUData()[pindex]);
                                }
                                for (int comp = 0; comp < NumRealComps(); ++comp) {
                                    RealVector& arr = soa_local[pld.m_lev][index][thread_num].GetRealData(comp);
                                    arr.push_back(soa.GetRealData(comp)[pindex]);
                                }
                                for (int comp = 0; comp < NumIntComps(); ++comp) {
                                    IntVector& arr = soa_local[pld.m_lev][index][thread_num].GetIntData(comp);
                                    arr.push_back(soa.GetIntData(comp)[pindex]);
                                }

                                p.id() = -p.id(); // Invalidate the particle
                            }
                        }
                        else {
                            auto& particles_to_send = tmp_remote[who][thread_num];
                            auto old_size = particles_to_send.size();
                            auto new_size = old_size + superparticle_size;
                            particles_to_send.resize(new_size);

                            char* dst = &particles_to_send[old_size];
                            {
                                std::memcpy(dst, &soa.GetIdCPUData()[pindex], sizeof(uint64_t));
                                dst += sizeof(uint64_t);
                            }
                            int array_comp_start = 0;
                            if constexpr (!ParticleType::is_soa_particle) {
                                array_comp_start = AMREX_SPACEDIM + NStructReal;
                            }
                            for (int comp = 0; comp < NumRealComps(); comp++) {
                                if (h_redistribute_real_comp[array_comp_start + comp]) {
                                    std::memcpy(dst, &soa.GetRealData(comp)[pindex], sizeof(ParticleReal));
                                    dst += sizeof(ParticleReal);
                                }
                            }
                            array_comp_start = 2 + NStructInt;
                            for (int comp = 0; comp < NumIntComps(); comp++) {
                                if (h_redistribute_int_comp[array_comp_start + comp]) {
                                    std::memcpy(dst, &soa.GetIntData(comp)[pindex], sizeof(int));
                                    dst += sizeof(int);
                                }
                            }
                            p.id() = -p.id(); // Invalidate the particle
                        }

                        if (p.id() < 0){
                            soa.GetIdCPUData()[pindex] = soa.GetIdCPUData()[last];
                            for (int comp = 0; comp < NumRealComps(); comp++) {
                                soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                            }
                            for (int comp = 0; comp < NumIntComps(); comp++) {
                                soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                            }
                            correctCellVectors(last, pindex, grid, ptd[pindex]);
                            --last;
                            continue;
                        }

                        ++pindex;
                    }

                    {
                        auto& iddata = soa.GetIdCPUData();
                        iddata.erase(iddata.begin() + last + 1, iddata.begin() + npart);
                    }
                    for (int comp = 0; comp < NumRealComps(); comp++) {
                        RealVector& rdata = soa.GetRealData(comp);
                        rdata.erase(rdata.begin() + last + 1, rdata.begin() + npart);
                    }
                    for (int comp = 0; comp < NumIntComps(); comp++) {
                        IntVector& idata = soa.GetIntData(comp);
                        idata.erase(idata.begin() + last + 1, idata.begin() + npart);
                    }
                }
            }
        }
    }

    for (int lev = lev_min; lev <= lev_max; lev++) {
        particle_detail::clearEmptyEntries(m_particles[lev]);
    }

    // Second pass - for each tile in parallel, collect the particles we are owed from all thread's buffers.
    for (int lev = lev_min; lev <= lev_max; lev++) {
        typename std::map<std::pair<int, int>, Vector<ParticleVector > >::iterator pmap_it;

        if constexpr(!ParticleType::is_soa_particle) {
            Vector<std::pair<int, int> > grid_tile_ids;
            Vector<Vector<ParticleVector>* > pvec_ptrs;

            // we need to create any missing map entries in serial here
            for (pmap_it=tmp_local[lev].begin(); pmap_it != tmp_local[lev].end(); pmap_it++)
            {
                DefineAndReturnParticleTile(lev, pmap_it->first.first, pmap_it->first.second);
                grid_tile_ids.push_back(pmap_it->first);
                pvec_ptrs.push_back(&(pmap_it->second));
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
            for (int pit = 0; pit < static_cast<int>(pvec_ptrs.size()); ++pit)
            {
                auto index = grid_tile_ids[pit];
                auto& ptile = ParticlesAt(lev, index.first, index.second);
                auto& aos = ptile.GetArrayOfStructs();
                auto& soa = ptile.GetStructOfArrays();
                auto& aos_tmp = *(pvec_ptrs[pit]);
                auto& soa_tmp = soa_local[lev][index];
                for (int i = 0; i < num_threads; ++i) {
                    aos.insert(aos.end(), aos_tmp[i].begin(), aos_tmp[i].end());
                    aos_tmp[i].erase(aos_tmp[i].begin(), aos_tmp[i].end());
                    for (int comp = 0; comp < NumRealComps(); ++comp) {
                        RealVector& arr = soa.GetRealData(comp);
                        RealVector& tmp = soa_tmp[i].GetRealData(comp);
                        arr.insert(arr.end(), tmp.begin(), tmp.end());
                        tmp.erase(tmp.begin(), tmp.end());
                    }
                    for (int comp = 0; comp < NumIntComps(); ++comp) {
                        IntVector& arr = soa.GetIntData(comp);
                        IntVector& tmp = soa_tmp[i].GetIntData(comp);
                        arr.insert(arr.end(), tmp.begin(), tmp.end());
                        tmp.erase(tmp.begin(), tmp.end());
                    }
                }
            }
        } else { // soa particle
            Vector<std::pair<int, int> > grid_tile_ids;

            // we need to create any missing map entries in serial here
            for (auto soa_map_it=soa_local[lev].begin(); soa_map_it != soa_local[lev].end(); soa_map_it++)
            {
                DefineAndReturnParticleTile(lev, soa_map_it->first.first, soa_map_it->first.second);
                grid_tile_ids.push_back(soa_map_it->first);
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
            for (int pit = 0; pit < static_cast<int>(grid_tile_ids.size()); ++pit) // NOLINT(modernize-loop-convert)
            {
                auto index = grid_tile_ids[pit];
                auto& ptile = ParticlesAt(lev, index.first, index.second);
                auto& soa = ptile.GetStructOfArrays();
                auto& soa_tmp = soa_local[lev][index];
                for (int i = 0; i < num_threads; ++i) {
                    {
                        auto& arr = soa.GetIdCPUData();
                        auto& tmp = soa_tmp[i].GetIdCPUData();
                        arr.insert(arr.end(), tmp.begin(), tmp.end());
                        tmp.erase(tmp.begin(), tmp.end());
                    }
                    for (int comp = 0; comp < NumRealComps(); ++comp) {
                        RealVector& arr = soa.GetRealData(comp);
                        RealVector& tmp = soa_tmp[i].GetRealData(comp);
                        arr.insert(arr.end(), tmp.begin(), tmp.end());
                        tmp.erase(tmp.begin(), tmp.end());
                    }
                    for (int comp = 0; comp < NumIntComps(); ++comp) {
                        IntVector& arr = soa.GetIntData(comp);
                        IntVector& tmp = soa_tmp[i].GetIntData(comp);
                        arr.insert(arr.end(), tmp.begin(), tmp.end());
                        tmp.erase(tmp.begin(), tmp.end());
                    }
                }
            }
        }
    }

    for (auto& map_it : tmp_remote) {
        int who = map_it.first;
        not_ours[who];
    }

    Vector<int> dest_proc_ids;
    Vector<Vector<Vector<char> >* > pbuff_ptrs;
    for (auto& kv : tmp_remote)
    {
        dest_proc_ids.push_back(kv.first);
        pbuff_ptrs.push_back(&(kv.second));
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int pmap_it = 0; pmap_it < static_cast<int>(pbuff_ptrs.size()); ++pmap_it)
    {
        int who = dest_proc_ids[pmap_it];
        Vector<Vector<char> >& tmp = *(pbuff_ptrs[pmap_it]);
        for (int i = 0; i < num_threads; ++i) {
            not_ours[who].insert(not_ours[who].end(), tmp[i].begin(), tmp[i].end());
            tmp[i].erase(tmp[i].begin(), tmp[i].end());
        }
    }

    particle_detail::clearEmptyEntries(not_ours);

    if (int(m_particles.size()) > theEffectiveFinestLevel+1) {
        // Looks like we lost an AmrLevel on a regrid.
        if (m_verbose > 0) {
            amrex::Print() << "ParticleContainer::Redistribute() resizing m_particles from "
                           << m_particles.size() << " to " << theEffectiveFinestLevel+1 << '\n';
        }
        AMREX_ASSERT(int(m_particles.size()) >= 2);

        m_particles.resize(theEffectiveFinestLevel + 1);
        m_dummy_mf.resize(theEffectiveFinestLevel + 1);
    }

    if (ParallelContext::NProcsSub() == 1) {
        AMREX_ASSERT(not_ours.empty());
    }
    else {
        RedistributeMPI(not_ours, lev_min, lev_max, nGrow, local);
    }

    AMREX_ASSERT(OK(lev_min, lev_max, nGrow));

    if (m_verbose > 0) {
        auto stoptime = amrex::second() - strttime;

        ByteSpread();

#ifdef AMREX_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                                ParallelContext::CommunicatorSub());

            amrex::Print() << "ParticleContainer::Redistribute() time: " << stoptime << "\n\n";
#ifdef AMREX_LAZY
        });
#endif
    }
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
RedistributeMPI (std::map<int, Vector<char> >& not_ours,
                 int lev_min, int lev_max, int nGrow, int local)
{
    BL_PROFILE("ParticleContainer::RedistributeMPI()");
    BL_PROFILE_VAR_NS("RedistributeMPI_locate", blp_locate);
    BL_PROFILE_VAR_NS("RedistributeMPI_copy", blp_copy);

#ifdef AMREX_USE_MPI

    using buffer_type = unsigned long long;

    std::map<int, Vector<buffer_type> > mpi_snd_data;
    for (const auto& kv : not_ours)
    {
        auto nbt = (kv.second.size() + sizeof(buffer_type)-1)/sizeof(buffer_type);
        mpi_snd_data[kv.first].resize(nbt);
        std::memcpy((char*) mpi_snd_data[kv.first].data(), kv.second.data(), kv.second.size());
    }

    const int NProcs = ParallelContext::NProcsSub();
    const int NNeighborProcs = neighbor_procs.size();

    // We may now have particles that are rightfully owned by another CPU.
    Vector<Long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    Long NumSnds = 0;
    if (local > 0)
    {
        AMREX_ALWAYS_ASSERT(lev_min == 0);
        AMREX_ALWAYS_ASSERT(lev_max == 0);
        BuildRedistributeMask(0, local);
        NumSnds = doHandShakeLocal(not_ours, neighbor_procs, Snds, Rcvs);
    }
    else
    {
        NumSnds = doHandShake(not_ours, Snds, Rcvs);
    }

    const int SeqNum = ParallelDescriptor::SeqNum();

    if ((! local) && NumSnds == 0) {
        return;  // There's no parallel work to do.
    }

    if (local)
    {
        Long tot_snds_this_proc = 0;
        Long tot_rcvs_this_proc = 0;
        for (int i = 0; i < NNeighborProcs; ++i) {
            tot_snds_this_proc += Snds[neighbor_procs[i]];
            tot_rcvs_this_proc += Rcvs[neighbor_procs[i]];
        }
        if ( (tot_snds_this_proc == 0) && (tot_rcvs_this_proc == 0) ) {
            return; // There's no parallel work to do.
        }
    }

    Vector<int> RcvProc;
    Vector<std::size_t> rOffset; // Offset (in bytes) in the receive buffer

    std::size_t TotRcvInts = 0;
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (Rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvInts);
            TotRcvBytes += Rcvs[i];
            auto nbt = (Rcvs[i] + sizeof(buffer_type)-1)/sizeof(buffer_type);
            TotRcvInts += nbt;
        }
    }

    const auto nrcvs = static_cast<int>(RcvProc.size());
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);

    // Allocate data for rcvs as one big chunk.
    Vector<unsigned long long> recvdata(TotRcvInts);

    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt = (Rcvs[Who] + sizeof(buffer_type)-1)/sizeof(buffer_type);
        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Cnt < size_t(std::numeric_limits<int>::max()));
        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum,
                                             ParallelContext::CommunicatorSub()).req();
    }

    // Send.
    for (const auto& kv : mpi_snd_data) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum,
                                 ParallelContext::CommunicatorSub());
    }

    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);

        BL_PROFILE_VAR_START(blp_locate);

        int npart = TotRcvBytes / superparticle_size;

        Vector<int> rcv_levs(npart);
        Vector<int> rcv_grid(npart);
        Vector<int> rcv_tile(npart);

        int ipart = 0;
        ParticleLocData pld;
        for (int j = 0; j < nrcvs; ++j)
        {
            const auto offset = rOffset[j];
            const auto Who    = RcvProc[j];
            const auto Cnt    = Rcvs[Who] / superparticle_size;
            for (int i = 0; i < int(Cnt); ++i)
            {
                char* pbuf = ((char*) &recvdata[offset]) + i*superparticle_size;

                Particle<NStructReal, NStructInt> p;

                if constexpr (ParticleType::is_soa_particle) {
                    std::memcpy(&p.m_idcpu, pbuf, sizeof(uint64_t));

                    ParticleReal pos[AMREX_SPACEDIM];
                    std::memcpy(&pos[0], pbuf + sizeof(uint64_t), AMREX_SPACEDIM*sizeof(ParticleReal));
                    AMREX_D_TERM(p.pos(0) = pos[0];,
                                 p.pos(1) = pos[1];,
                                 p.pos(2) = pos[2]);
                } else {
                   std::memcpy(&p, pbuf, sizeof(ParticleType));
                }

                bool success = Where(p, pld, lev_min, lev_max, 0);
                if (!success)
                {
                    success = (nGrow > 0) && Where(p, pld, lev_min, lev_min, nGrow);
                    pld.m_grown_gridbox = pld.m_gridbox; // reset grown box for subsequent calls.
                }
                if (!success)
                {
                    amrex::Abort("RedistributeMPI_locate:: invalid particle.");
                }

                rcv_levs[ipart] = pld.m_lev;
                rcv_grid[ipart] = pld.m_grid;
                rcv_tile[ipart] = pld.m_tile;

                ++ipart;
            }
        }

        BL_PROFILE_VAR_STOP(blp_locate);

        BL_PROFILE_VAR_START(blp_copy);

#ifndef AMREX_USE_GPU
        ipart = 0;
        for (int i = 0; i < nrcvs; ++i)
        {
            const auto offset = rOffset[i];
            const auto Who    = RcvProc[i];
            const auto Cnt = Rcvs[Who] / superparticle_size;
            for (int j = 0; j < int(Cnt); ++j)
            {
                auto& ptile = m_particles[rcv_levs[ipart]][std::make_pair(rcv_grid[ipart],
                                                                          rcv_tile[ipart])];
                char* pbuf = ((char*) &recvdata[offset]) + j*superparticle_size;

                if constexpr (ParticleType::is_soa_particle) {
                    uint64_t idcpudata;
                    std::memcpy(&idcpudata, pbuf, sizeof(uint64_t));
                    pbuf += sizeof(uint64_t);
                    ptile.GetStructOfArrays().GetIdCPUData().push_back(idcpudata);
                } else {
                    ParticleType p;
                    std::memcpy(&p, pbuf, sizeof(ParticleType));
                    pbuf += sizeof(ParticleType);
                    ptile.push_back(p);
                }

                int array_comp_start = 0;
                if constexpr (!ParticleType::is_soa_particle) {
                    array_comp_start = AMREX_SPACEDIM + NStructReal;
                }
                for (int comp = 0; comp < NumRealComps(); ++comp) {
                    if (h_redistribute_real_comp[array_comp_start + comp]) {
                        ParticleReal rdata;
                        std::memcpy(&rdata, pbuf, sizeof(ParticleReal));
                        pbuf += sizeof(ParticleReal);
                        ptile.push_back_real(comp, rdata);
                    } else {
                        ptile.push_back_real(comp, 0.0);
                    }
                }

                array_comp_start = 2 + NStructInt;
                for (int comp = 0; comp < NumIntComps(); ++comp) {
                    if (h_redistribute_int_comp[array_comp_start + comp]) {
                        int idata;
                        std::memcpy(&idata, pbuf, sizeof(int));
                        pbuf += sizeof(int);
                        ptile.push_back_int(comp, idata);
                    } else {
                        ptile.push_back_int(comp, 0);
                    }
                }
                ++ipart;
            }
        }

#else
        Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
        host_particles.reserve(15);
        host_particles.resize(finestLevel()+1);

        Vector<std::map<std::pair<int, int>,
                        std::vector<Gpu::HostVector<ParticleReal> > > > host_real_attribs;
        host_real_attribs.reserve(15);
        host_real_attribs.resize(finestLevel()+1);

        Vector<std::map<std::pair<int, int>,
                        std::vector<Gpu::HostVector<int> > > > host_int_attribs;
        host_int_attribs.reserve(15);
        host_int_attribs.resize(finestLevel()+1);

        Vector<std::map<std::pair<int, int>, Gpu::HostVector<uint64_t> > > host_idcpu;
        host_idcpu.reserve(15);
        host_idcpu.resize(finestLevel()+1);

        ipart = 0;
        for (int i = 0; i < nrcvs; ++i)
        {
            const auto offset = rOffset[i];
            const auto Who    = RcvProc[i];
            const auto Cnt = Rcvs[Who] / superparticle_size;
            for (auto j = decltype(Cnt)(0); j < Cnt; ++j)
            {
                int lev = rcv_levs[ipart];
                std::pair<int, int> ind(std::make_pair(rcv_grid[ipart], rcv_tile[ipart]));

                char* pbuf = ((char*) &recvdata[offset]) + j*superparticle_size;

                host_real_attribs[lev][ind].resize(NumRealComps());
                host_int_attribs[lev][ind].resize(NumIntComps());

                if constexpr (ParticleType::is_soa_particle) {
                    uint64_t idcpudata;
                    std::memcpy(&idcpudata, pbuf, sizeof(uint64_t));
                    pbuf += sizeof(uint64_t);
                    host_idcpu[lev][ind].push_back(idcpudata);
                } else {
                    ParticleType p;
                    std::memcpy(&p, pbuf, sizeof(ParticleType));
                    pbuf += sizeof(ParticleType);
                    host_particles[lev][ind].push_back(p);
                }

                host_real_attribs[lev][ind].resize(NumRealComps());
                host_int_attribs[lev][ind].resize(NumIntComps());

                // add the real...
                int array_comp_start = 0;
                if constexpr (!ParticleType::is_soa_particle) {
                    array_comp_start = AMREX_SPACEDIM + NStructReal;
                }
                for (int comp = 0; comp < NumRealComps(); ++comp) {
                    if (h_redistribute_real_comp[array_comp_start + comp]) {
                        Real rdata;
                        std::memcpy(&rdata, pbuf, sizeof(Real));
                        pbuf += sizeof(Real);
                        host_real_attribs[lev][ind][comp].push_back(rdata);
                    } else {
                        host_real_attribs[lev][ind][comp].push_back(0.0);
                    }
                }

                // ... and int array data
                array_comp_start = 2 + NStructInt;
                for (int comp = 0; comp < NumIntComps(); ++comp) {
                    if (h_redistribute_int_comp[array_comp_start + comp]) {
                        int idata;
                        std::memcpy(&idata, pbuf, sizeof(int));
                        pbuf += sizeof(int);
                        host_int_attribs[lev][ind][comp].push_back(idata);
                    } else {
                        host_int_attribs[lev][ind][comp].push_back(0);
                    }
                }
                ++ipart;
            }
        }

        for (int host_lev = 0; host_lev < static_cast<int>(host_particles.size()); ++host_lev)
        {
            for (auto& kv : host_particles[host_lev]) {
              auto grid = kv.first.first;
              auto tile = kv.first.second;
              const auto& src_tile = kv.second;

              auto& dst_tile = GetParticles(host_lev)[std::make_pair(grid,tile)];
              auto old_size = dst_tile.size();
              auto new_size = old_size + src_tile.size();
              dst_tile.resize(new_size);

              if constexpr (ParticleType::is_soa_particle) {
                  Gpu::copyAsync(Gpu::hostToDevice,
                                 host_idcpu[host_lev][std::make_pair(grid,tile)].begin(),
                                 host_idcpu[host_lev][std::make_pair(grid,tile)].end(),
                                 dst_tile.GetStructOfArrays().GetIdCPUData().begin() + old_size);
              } else {
                  Gpu::copyAsync(Gpu::hostToDevice,
                                 src_tile.begin(), src_tile.end(),
                                 dst_tile.GetArrayOfStructs().begin() + old_size);
              }

              for (int i = 0; i < NumRealComps(); ++i) {
                  Gpu::copyAsync(Gpu::hostToDevice,
                                 host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                                 host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                                 dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
              }

              for (int i = 0; i < NumIntComps(); ++i) {
                  Gpu::copyAsync(Gpu::hostToDevice,
                                 host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                                 host_int_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                                 dst_tile.GetStructOfArrays().GetIntData(i).begin() + old_size);
              }
            }
          }

        Gpu::Device::streamSynchronize();
#endif

        BL_PROFILE_VAR_STOP(blp_copy);
    }
#else
    amrex::ignore_unused(not_ours,lev_min,lev_max,nGrow,local);
#endif
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
bool
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::OK (int lev_min, int lev_max, int nGrow) const
{
    BL_PROFILE("ParticleContainer::OK()");

    if (lev_max == -1) {
        lev_max = finestLevel();
}

    return (numParticlesOutOfRange(*this, lev_min, lev_max, nGrow) == 0);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::AddParticlesAtLevel (AoS& particles, int level, int nGrow)
{
    ParticleTileType ptile;
    ptile.GetArrayOfStructs().swap(particles);
    AddParticlesAtLevel(ptile, level, nGrow);
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::AddParticlesAtLevel (ParticleTileType& particles, int level, int nGrow)
{
    BL_PROFILE("ParticleContainer::AddParticlesAtLevel()");

    if (int(m_particles.size()) < level+1)
    {
        if (Verbose())
        {
            amrex::Print() << "ParticleContainer::AddParticlesAtLevel resizing m_particles from "
                           << m_particles.size()
                           << " to "
                           << level+1 << '\n';
        }
        m_particles.resize(level+1);
        m_dummy_mf.resize(level+1);
        for (int lev = 0; lev < level+1; ++lev) {
            RedefineDummyMF(lev);
        }
    }

    auto& ptile = DefineAndReturnParticleTile(level, 0, 0);
    int old_np = ptile.size();
    int num_to_add = particles.size();
    int new_np = old_np + num_to_add;
    ptile.resize(new_np);
    amrex::copyParticles(ptile, particles, 0, old_np, num_to_add);
    Redistribute(level, level, nGrow);
    particles.resize(0);
}

// This is the single-level version for cell-centered density
template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
AssignCellDensitySingleLevel (int rho_index,
                              MultiFab& mf_to_be_filled,
                              int       lev,
                              int       ncomp,
                              int       particle_lvl_offset) const
{
    BL_PROFILE("ParticleContainer::AssignCellDensitySingleLevel()");

    if (rho_index != 0) { amrex::Abort("AssignCellDensitySingleLevel only works if rho_index = 0"); }

    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled)) {
      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    }
    else {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new MultiFab(ParticleBoxArray(lev),
                                ParticleDistributionMap(lev),
                                ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->SumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (mf_pointer->nGrow() < 1) {
       amrex::Error("Must have at least one ghost cell when in AssignCellDensitySingleLevel");
    }

    const auto strttime = amrex::second();

    const auto dxi              = Geom(lev).InvCellSizeArray();
    const auto plo              = Geom(lev).ProbLoArray();
    const auto pdxi             = Geom(lev + particle_lvl_offset).InvCellSizeArray();

    if (Geom(lev).isAnyPeriodic() && ! Geom(lev).isAllPeriodic())
    {
        amrex::Error("AssignCellDensitySingleLevel: problem must be periodic in no or all directions");
    }

    mf_pointer->setVal(0);

    using ParConstIter = ParConstIter_impl<ParticleType, NArrayReal, NArrayInt, Allocator>;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox local_rho;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {
            const Long np = pti.numParticles();
            auto ptd = pti.GetParticleTile().getConstParticleTileData();
            FArrayBox& fab = (*mf_pointer)[pti];
            auto rhoarr = fab.array();
#ifdef AMREX_USE_OMP
            Box tile_box;
            if (Gpu::notInLaunchRegion())
            {
                tile_box = pti.tilebox();
                tile_box.grow(mf_pointer->nGrow());
                local_rho.resize(tile_box,ncomp);
                local_rho.setVal<RunOn::Host>(0.0);
                rhoarr = local_rho.array();
            }
#endif

            if (particle_lvl_offset == 0)
            {
                AMREX_HOST_DEVICE_FOR_1D( np, i,
                {
                    auto p = make_particle<ParticleType>{}(ptd,i);
                    amrex_deposit_cic(p, ncomp, rhoarr, plo, dxi);
                });
            }
            else
            {
                AMREX_HOST_DEVICE_FOR_1D( np, i,
                {
                    auto p = make_particle<ParticleType>{}(ptd,i);
                    amrex_deposit_particle_dx_cic(p, ncomp, rhoarr, plo, dxi, pdxi);
                });
            }

#ifdef AMREX_USE_OMP
            if (Gpu::notInLaunchRegion())
            {
                fab.atomicAdd<RunOn::Host>(local_rho, tile_box, tile_box, 0, 0, ncomp);
            }
#endif
        }
    }

    mf_pointer->SumBoundary(Geom(lev).periodicity());

    // If ncomp > 1, first divide the momenta (component n)
    // by the mass (component 0) in order to get velocities.
    // Be careful not to divide by zero.
    for (int n = 1; n < ncomp; n++)
    {
        for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)
        {
            (*mf_pointer)[mfi].protected_divide<RunOn::Device>((*mf_pointer)[mfi],0,n,1);
        }
    }

    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    const Real* dx = Geom(lev).CellSize();
    const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(Real(1.0)/vol, 0, 1, mf_pointer->nGrow());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.ParallelCopy(*mf_pointer,0,0,ncomp,0,0);
        delete mf_pointer;
    }

    if (m_verbose > 1)
    {
        auto stoptime = amrex::second() - strttime;

        ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                            ParallelContext::CommunicatorSub());

        amrex::Print() << "ParticleContainer::AssignCellDensitySingleLevel) time: "
                       << stoptime << '\n';
    }
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
ResizeRuntimeRealComp (int new_size, bool communicate)
{
    int old_size = m_num_runtime_real;

    m_runtime_comps_defined = (new_size > 0);
    m_num_runtime_real = new_size;
    int cur_size = h_redistribute_real_comp.size();
    h_redistribute_real_comp.resize(cur_size-old_size+new_size, communicate);
    SetParticleSize();

    for (int lev = 0; lev < numLevels(); ++lev) {
        for (ParIterType pti(*this,lev); pti.isValid(); ++pti) {
            auto& tile = DefineAndReturnParticleTile(lev, pti);
            auto np = tile.numParticles();
            if (np > 0 && new_size > old_size) {
                auto& soa = tile.GetStructOfArrays();
                soa.resize(np);
            }
        }
    }
}

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
ResizeRuntimeIntComp (int new_size, bool communicate)
{
    int old_size = m_num_runtime_int;

    m_runtime_comps_defined = (new_size > 0);
    m_num_runtime_int = new_size;
    int cur_size = h_redistribute_int_comp.size();
    h_redistribute_int_comp.resize(cur_size-old_size+new_size, communicate);
    SetParticleSize();

    for (int lev = 0; lev < numLevels(); ++lev) {
        for (ParIterType pti(*this,lev); pti.isValid(); ++pti) {
            auto& tile = DefineAndReturnParticleTile(lev, pti);
            auto np = tile.numParticles();
            if (np > 0 && new_size > old_size) {
                auto& soa = tile.GetStructOfArrays();
                soa.resize(np);
            }
        }
    }
}

}
