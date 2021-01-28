#include <AMReX_ParticleContainerBase.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

bool    ParticleContainerBase::do_tiling = false;
IntVect ParticleContainerBase::tile_size { AMREX_D_DECL(1024000,8,8) };

void ParticleContainerBase::Define (const Geometry            & geom,
                                    const DistributionMapping & dmap,
                                    const BoxArray            & ba)
{
    m_gdb_object = ParGDB(geom, dmap, ba);
    m_gdb = &m_gdb_object;
}

void ParticleContainerBase::Define (const Vector<Geometry>            & geom,
                                    const Vector<DistributionMapping> & dmap,
                                    const Vector<BoxArray>            & ba,
                                    const Vector<int>                 & rr)
{
    m_gdb_object = ParGDB(geom, dmap, ba, rr);
    m_gdb = &m_gdb_object;
}

void ParticleContainerBase::SetParticleBoxArray (int lev, const BoxArray& new_ba)
{
    m_gdb_object = ParGDB(m_gdb->ParticleGeom(), m_gdb->ParticleDistributionMap(),
                          m_gdb->ParticleBoxArray(), m_gdb->refRatio());
    m_gdb = &m_gdb_object;
    m_gdb->SetParticleBoxArray(lev, new_ba);
}

void ParticleContainerBase::SetParticleDistributionMap (int lev, const DistributionMapping& new_dmap)
{
    m_gdb_object = ParGDB(m_gdb->ParticleGeom(), m_gdb->ParticleDistributionMap(),
                          m_gdb->ParticleBoxArray(), m_gdb->refRatio());
    m_gdb = &m_gdb_object;
    m_gdb->SetParticleDistributionMap(lev, new_dmap);
}

void ParticleContainerBase::SetParticleGeometry (int lev, const Geometry& new_geom)
{
    m_gdb_object = ParGDB(m_gdb->ParticleGeom(), m_gdb->ParticleDistributionMap(),
                          m_gdb->ParticleBoxArray(), m_gdb->refRatio());
    m_gdb = &m_gdb_object;
    m_gdb->SetParticleGeometry(lev, new_geom);
}

const std::string& ParticleContainerBase::Version ()
{
    //
    // If we change the Checkpoint/Restart format we should increment this.
    //
    // Previous version strings:
    //
    //    "Version_One_Dot_Zero"
    //    "Version_One_Dot_One"
    //
    static const std::string version("Version_Two_Dot_Zero");

    return version;
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
        if (aggregation_type != "None" || aggregation_type != "Cell")
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
