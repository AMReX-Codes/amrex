#include <AMReX_ParticleContainerBase.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_iMultiFab.H>

namespace amrex {

bool    ParticleContainerBase::do_tiling = false;
IntVect ParticleContainerBase::tile_size { AMREX_D_DECL(1024000,8,8) };
bool    ParticleContainerBase::memEfficientSort = true;
bool    ParticleContainerBase::use_comms_arena = false;

ParticleContainerBase::~ParticleContainerBase ()
{
#if defined(AMREX_USE_MPI) && defined(BL_USE_MPI3)
    releaseParticleHandshakeWindow();
#endif
}

ParticleContainerBase::ParticleContainerBase (ParticleContainerBase&& other) noexcept
    : m_particle_locator(std::move(other.m_particle_locator)),
      m_verbose(other.m_verbose),
      m_stable_redistribute(other.m_stable_redistribute),
      m_gdb_object(std::move(other.m_gdb_object)),
      m_gdb(other.m_gdb),
      m_dummy_mf(std::move(other.m_dummy_mf)),
      m_arena(other.m_arena),
      redistribute_mask_ptr(std::move(other.redistribute_mask_ptr)),
      redistribute_mask_nghost(other.redistribute_mask_nghost),
      neighbor_procs(std::move(other.neighbor_procs)),
      m_buffer_map(std::move(other.m_buffer_map))
#if defined(AMREX_USE_MPI) && defined(BL_USE_MPI3)
    , m_particle_handshake_win(other.m_particle_handshake_win),
      m_particle_handshake_ptr(other.m_particle_handshake_ptr),
      m_particle_handshake_nprocs(other.m_particle_handshake_nprocs),
      m_particle_handshake_comm(other.m_particle_handshake_comm)
#endif
{
    other.m_gdb = nullptr;
#if defined(AMREX_USE_MPI) && defined(BL_USE_MPI3)
    other.m_particle_handshake_win = MPI_WIN_NULL;
    other.m_particle_handshake_ptr = nullptr;
    other.m_particle_handshake_nprocs = 0;
    other.m_particle_handshake_comm = MPI_COMM_NULL;
#endif
}

ParticleContainerBase&
ParticleContainerBase::operator= (ParticleContainerBase&& other) noexcept
{
    if (this != &other)
    {
#if defined(AMREX_USE_MPI) && defined(BL_USE_MPI3)
        releaseParticleHandshakeWindow();
#endif

        m_particle_locator = std::move(other.m_particle_locator);
        m_verbose = other.m_verbose;
        m_stable_redistribute = other.m_stable_redistribute;
        m_gdb_object = std::move(other.m_gdb_object);
        m_gdb = other.m_gdb;
        m_dummy_mf = std::move(other.m_dummy_mf);
        m_arena = other.m_arena;
        redistribute_mask_ptr = std::move(other.redistribute_mask_ptr);
        redistribute_mask_nghost = other.redistribute_mask_nghost;
        neighbor_procs = std::move(other.neighbor_procs);
        m_buffer_map = std::move(other.m_buffer_map);
#if defined(AMREX_USE_MPI) && defined(BL_USE_MPI3)
        m_particle_handshake_win = other.m_particle_handshake_win;
        m_particle_handshake_ptr = other.m_particle_handshake_ptr;
        m_particle_handshake_nprocs = other.m_particle_handshake_nprocs;
        m_particle_handshake_comm = other.m_particle_handshake_comm;

        other.m_particle_handshake_win = MPI_WIN_NULL;
        other.m_particle_handshake_ptr = nullptr;
        other.m_particle_handshake_nprocs = 0;
        other.m_particle_handshake_comm = MPI_COMM_NULL;
#endif
        other.m_gdb = nullptr;
    }

    return *this;
}

void ParticleContainerBase::Define (const Geometry            & geom,
                                    const DistributionMapping & dmap,
                                    const BoxArray            & ba)
{
    *m_gdb_object = ParGDB(geom, dmap, ba);
    m_gdb = static_cast<ParGDBBase*>(m_gdb_object.get());
}

void ParticleContainerBase::Define (const Vector<Geometry>            & geom,
                                    const Vector<DistributionMapping> & dmap,
                                    const Vector<BoxArray>            & ba,
                                    const Vector<int>                 & rr)
{
    *m_gdb_object = ParGDB(geom, dmap, ba, rr);
    m_gdb = static_cast<ParGDBBase*>(m_gdb_object.get());
}

void ParticleContainerBase::Define (const Vector<Geometry>            & geom,
                                    const Vector<DistributionMapping> & dmap,
                                    const Vector<BoxArray>            & ba,
                                    const Vector<IntVect>             & rr)
{
    *m_gdb_object = ParGDB(geom, dmap, ba, rr);
    m_gdb = static_cast<ParGDBBase*>(m_gdb_object.get());
}

void ParticleContainerBase::reserveData ()
{
    m_dummy_mf.reserve(maxLevel()+1);
}

void ParticleContainerBase::resizeData ()
{
    int nlevs = std::max(0, finestLevel()+1);
    m_dummy_mf.resize(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        RedefineDummyMF(lev);
    }
}

void ParticleContainerBase::RedefineDummyMF (int lev)
{
    if (lev > m_dummy_mf.size()-1) { m_dummy_mf.resize(lev+1); }

    if (m_dummy_mf[lev] == nullptr ||
        ! BoxArray::SameRefs(m_dummy_mf[lev]->boxArray(),
                             ParticleBoxArray(lev))          ||
        ! DistributionMapping::SameRefs(m_dummy_mf[lev]->DistributionMap(),
                                        ParticleDistributionMap(lev)))
    {
        auto dm = (ParticleBoxArray(lev).size() == ParticleDistributionMap(lev).size()) ?
            ParticleDistributionMap(lev) : DistributionMapping(ParticleBoxArray(lev));
        m_dummy_mf[lev] = std::make_unique<MultiFab>(ParticleBoxArray(lev),
                                                     dm, 1,0,MFInfo().SetAlloc(false));
    };
}

void
ParticleContainerBase::defineBufferMap () const
{
    BL_PROFILE("ParticleContainer::defineBufferMap");

    if (! m_buffer_map.isValid(GetParGDB()))
    {
        m_buffer_map.define(GetParGDB());
    }
}

#if defined(AMREX_USE_MPI) && defined(BL_USE_MPI3)
void ParticleContainerBase::releaseParticleHandshakeWindow ()
{
    if (m_particle_handshake_win != MPI_WIN_NULL) {
        BL_MPI_REQUIRE(MPI_Win_free(&m_particle_handshake_win));
    }
    if (m_particle_handshake_comm != MPI_COMM_NULL) {
        BL_MPI_REQUIRE(MPI_Comm_free(&m_particle_handshake_comm));
    }
    m_particle_handshake_ptr = nullptr;
    m_particle_handshake_nprocs = 0;
}

void ParticleContainerBase::ensureParticleHandshakeWindow () const
{
    const int nprocs = ParallelContext::NProcsSub();
    MPI_Comm comm = ParallelContext::CommunicatorSub();

    bool needs_rebuild = (m_particle_handshake_win == MPI_WIN_NULL)
        || (m_particle_handshake_nprocs != nprocs)
        || (m_particle_handshake_comm == MPI_COMM_NULL);

    if (!needs_rebuild)
    {
        int cmp = MPI_UNEQUAL;
        BL_MPI_REQUIRE(MPI_Comm_compare(comm, m_particle_handshake_comm, &cmp));
        needs_rebuild = (cmp != MPI_IDENT && cmp != MPI_CONGRUENT);
    }

    if (needs_rebuild)
    {
        const_cast<ParticleContainerBase*>(this)->releaseParticleHandshakeWindow();

        Long* baseptr = nullptr;
        MPI_Win win = MPI_WIN_NULL;
        BL_MPI_REQUIRE(MPI_Win_allocate(static_cast<MPI_Aint>(nprocs*sizeof(Long)),
                                        sizeof(Long),
                                        MPI_INFO_NULL,
                                        comm,
                                        &baseptr,
                                        &win));

        MPI_Comm dup_comm = MPI_COMM_NULL;
        BL_MPI_REQUIRE(MPI_Comm_dup(comm, &dup_comm));

        m_particle_handshake_ptr = baseptr;
        m_particle_handshake_win = win;
        m_particle_handshake_nprocs = nprocs;
        m_particle_handshake_comm = dup_comm;
    }
}
#endif

void ParticleContainerBase::SetParGDB (const Geometry            & geom,
                                       const DistributionMapping & dmap,
                                       const BoxArray            & ba)
{
    *m_gdb_object = ParGDB(geom, dmap, ba);
    m_gdb = static_cast<ParGDBBase*>(m_gdb_object.get());
    resizeData();
}

void ParticleContainerBase::SetParGDB (const Vector<Geometry>            & geom,
                                       const Vector<DistributionMapping> & dmap,
                                       const Vector<BoxArray>            & ba,
                                       const Vector<int>                 & rr)
{
    *m_gdb_object = ParGDB(geom, dmap, ba, rr);
    m_gdb = static_cast<ParGDBBase*>(m_gdb_object.get());
    resizeData();
}

void ParticleContainerBase::SetParGDB (const Vector<Geometry>            & geom,
                                       const Vector<DistributionMapping> & dmap,
                                       const Vector<BoxArray>            & ba,
                                       const Vector<IntVect>             & rr)
{
    *m_gdb_object = ParGDB(geom, dmap, ba, rr);
    m_gdb = static_cast<ParGDBBase*>(m_gdb_object.get());
    resizeData();
}

void ParticleContainerBase::SetParticleBoxArray (int lev, BoxArray new_ba) // NOLINT(performance-unnecessary-value-param)
{
    // Must take the new BoxArray by value to avoid aliasing with what's
    // inside m_gdb_object
    *m_gdb_object = ParGDB(m_gdb->ParticleGeom(),
                           m_gdb->ParticleDistributionMap(),
                           m_gdb->ParticleBoxArray(),
                           m_gdb->refRatio());
    m_gdb = static_cast<ParGDBBase*>(m_gdb_object.get());
    m_gdb->SetParticleBoxArray(lev, new_ba);
    RedefineDummyMF(lev);
}

void ParticleContainerBase::SetParticleDistributionMap (int lev, DistributionMapping new_dmap) // NOLINT(performance-unnecessary-value-param)
{
    // Must take the new DistributionMapping by value to avoid aliasing with
    // what's inside m_gdb_object
    *m_gdb_object = ParGDB(m_gdb->ParticleGeom(),
                           m_gdb->ParticleDistributionMap(),
                           m_gdb->ParticleBoxArray(),
                           m_gdb->refRatio());
    m_gdb = static_cast<ParGDBBase*>(m_gdb_object.get());
    m_gdb->SetParticleDistributionMap(lev, new_dmap);
    RedefineDummyMF(lev);
}

void ParticleContainerBase::SetParticleGeometry (int lev, Geometry new_geom) // NOLINT(performance-unnecessary-value-param)
{
    // Must take the new Geometry by value to avoid aliasing with what's
    // inside m_gdb_object
    *m_gdb_object = ParGDB(m_gdb->ParticleGeom(),
                           m_gdb->ParticleDistributionMap(),
                           m_gdb->ParticleBoxArray(),
                           m_gdb->refRatio());
    m_gdb = static_cast<ParGDBBase*>(m_gdb_object.get());
    m_gdb->SetParticleGeometry(lev, new_geom);
}

const std::string& ParticleContainerBase::CheckpointVersion ()
{
    //
    // If we change the Checkpoint/Restart format we should increment this.
    //
    // Previous version strings:
    //
    //    "Version_One_Dot_Zero"
    //    "Version_One_Dot_One"
    //    "Version_Two_Dot_Zero" (before checkpoints had expanded particle ids)
    //
    static const std::string checkpoint_version("Version_Two_Dot_One");

    return checkpoint_version;
}

const std::string& ParticleContainerBase::PlotfileVersion ()
{
    //
    // If we change the plotfile format we should increment this.
    //
    // Previous version strings:
    //
    //    "Version_One_Dot_Zero"
    //    "Version_One_Dot_One"
    //
    static const std::string plotfile_version("Version_Two_Dot_Zero");

    return plotfile_version;
}

const std::string& ParticleContainerBase::DataPrefix ()
{
    //
    // The actual particle data is stored in files of the form: DATA_nnnn.
    //
    static const std::string data("DATA_");

    return data;
}

int ParticleContainerBase::MaxReaders ()
{
    const int Max_Readers_def = 64;
    static int Max_Readers;
    static bool first = true;

    if (first)
    {
        first = false;
        ParmParse pp("particles");
        Max_Readers = Max_Readers_def;
        pp.query("nreaders", Max_Readers);
        Max_Readers = std::min(ParallelDescriptor::NProcs(),Max_Readers);
        if (Max_Readers <= 0)
        {
            amrex::Abort("particles.nreaders must be positive");
        }
    }

    return Max_Readers;
}

Long ParticleContainerBase::MaxParticlesPerRead ()
{
    //
    // This is the maximum particles that "each" reader will attempt to read
    // before doing a Redistribute().
    //
    const Long Max_Particles_Per_Read_def = 100000;
    static Long Max_Particles_Per_Read;
    static bool first = true;

    if (first)
    {
        first = false;
        ParmParse pp("particles");
        Max_Particles_Per_Read = Max_Particles_Per_Read_def;
        pp.query("nparts_per_read", Max_Particles_Per_Read);
        if (Max_Particles_Per_Read <= 0)
        {
            amrex::Abort("particles.nparts_per_read must be positive");
        }
    }

    return Max_Particles_Per_Read;
}

const std::string& ParticleContainerBase::AggregationType ()
{
    static std::string aggregation_type;
    static bool first = true;

    if (first)
    {
        first = false;
        aggregation_type = "None";
        ParmParse pp("particles");
        pp.query("aggregation_type", aggregation_type);
        if (!(aggregation_type == "None" || aggregation_type == "Cell"))
        {
            amrex::Abort("particles.aggregation_type not implemented.");
        }
    }

    return aggregation_type;
}

int ParticleContainerBase::AggregationBuffer ()
{
    static int aggregation_buffer;
    static bool first = true;

    if (first)
    {
        first = false;
        aggregation_buffer = 2;
        ParmParse pp("particles");
        pp.query("aggregation_buffer", aggregation_buffer);
        if (aggregation_buffer <= 0)
        {
            amrex::Abort("particles.aggregation_buffer must be positive");
        }
    }

    return aggregation_buffer;
}

void ParticleContainerBase::BuildRedistributeMask (int lev, int nghost) const
{
    BL_PROFILE("ParticleContainer::BuildRedistributeMask");
    AMREX_ASSERT(lev == 0);

    if (redistribute_mask_ptr == nullptr ||
        redistribute_mask_nghost < nghost ||
        ! BoxArray::SameRefs(redistribute_mask_ptr->boxArray(), this->ParticleBoxArray(lev)) ||
        ! DistributionMapping::SameRefs(redistribute_mask_ptr->DistributionMap(), this->ParticleDistributionMap(lev)))
    {
        const Geometry& geom = this->Geom(lev);
        const BoxArray& ba = this->ParticleBoxArray(lev);
        const DistributionMapping& dmap = this->ParticleDistributionMap(lev);

        redistribute_mask_nghost = nghost;
        redistribute_mask_ptr = std::make_unique<iMultiFab>(ba, dmap, 2, nghost);
        redistribute_mask_ptr->setVal(-1, nghost);

        const auto tile_size_do = amrex::ParticleContainerBase::do_tiling ? amrex::ParticleContainerBase::tile_size : IntVect::TheZeroVector();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(*redistribute_mask_ptr, tile_size_do); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            (*redistribute_mask_ptr)[mfi].template setVal<RunOn::Host>(grid_id, box, 0, 1);
            (*redistribute_mask_ptr)[mfi].template setVal<RunOn::Host>(tile_id, box, 1, 1);
        }

        redistribute_mask_ptr->FillBoundary(geom.periodicity());

        neighbor_procs.clear();
        for (MFIter mfi(*redistribute_mask_ptr, tile_size_do); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.growntilebox();
            for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv))
            {
                const int grid = (*redistribute_mask_ptr)[mfi](iv, 0);
                if (grid >= 0)
                {
                    const int global_rank = this->ParticleDistributionMap(lev)[grid];
                    const int rank = ParallelContext::global_to_local_rank(global_rank);
                    if (rank != ParallelContext::MyProcSub()) {
                        neighbor_procs.push_back(rank);
                    }
                }
            }
        }
        RemoveDuplicates(neighbor_procs);
    }
}

}
