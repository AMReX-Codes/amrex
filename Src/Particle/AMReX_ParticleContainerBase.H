#ifndef AMREX_PARTICLECONTAINERBASE_H_
#define AMREX_PARTICLECONTAINERBASE_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_INT.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParGDB.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Vector.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParticleLocator.H>
#include <AMReX_DenseBins.H>

#include <string>

namespace amrex {

class ParticleContainerBase
{
public:

    ParticleContainerBase () = default;

    ParticleContainerBase (ParGDBBase* gdb)
        :
        m_gdb(gdb)
    {}

    ParticleContainerBase (const Geometry            & geom,
                           const DistributionMapping & dmap,
                           const BoxArray            & ba)
        :
        m_gdb_object(std::make_unique<ParGDB>(geom,dmap,ba)),
        m_gdb(static_cast<ParGDBBase*>(m_gdb_object.get()))
    {
    }

    ParticleContainerBase (const Vector<Geometry>            & geom,
                           const Vector<DistributionMapping> & dmap,
                           const Vector<BoxArray>            & ba,
                           const Vector<int>                 & rr)
        :
        m_gdb_object(std::make_unique<ParGDB>(geom,dmap,ba,rr)),
        m_gdb(static_cast<ParGDBBase*>(m_gdb_object.get()))
    {
    }

    ParticleContainerBase (const Vector<Geometry>            & geom,
                           const Vector<DistributionMapping> & dmap,
                           const Vector<BoxArray>            & ba,
                           const Vector<IntVect>             & rr)
        :
        m_gdb_object
        (std::make_unique<ParGDB>(geom,dmap,ba, [&]() -> Vector<int> {
                Vector<int> ref_ratio;
                for (int i = 0; i < static_cast<int>(rr.size()); ++i)
                {
#if AMREX_SPACEDIM > 1
                    AMREX_ASSERT(rr[i][0] == rr[i][1]);
#endif
#if AMREX_SPACEDIM > 2
                    AMREX_ASSERT(rr[i][0] == rr[i][2]);
#endif
                    ref_ratio.push_back(rr[i][0]);
                }
                return ref_ratio;  }() )),
        m_gdb(static_cast<ParGDBBase*>(m_gdb_object.get()))
    {
    }

    virtual ~ParticleContainerBase () = default;

    ParticleContainerBase ( const ParticleContainerBase &) = delete;
    ParticleContainerBase& operator= ( const ParticleContainerBase & ) = delete;

    ParticleContainerBase ( ParticleContainerBase && ) = default;
    ParticleContainerBase& operator= ( ParticleContainerBase && ) = default;

    void Define (ParGDBBase* gdb) { m_gdb = gdb;}

    void Define (const Geometry            & geom,
                 const DistributionMapping & dmap,
                 const BoxArray            & ba);

    void Define (const Vector<Geometry>            & geom,
                 const Vector<DistributionMapping> & dmap,
                 const Vector<BoxArray>            & ba,
                 const Vector<int>                 & rr);

    void Define (const Vector<Geometry>            & geom,
                 const Vector<DistributionMapping> & dmap,
                 const Vector<BoxArray>            & ba,
                 const Vector<IntVect>             & rr);

    bool isDefined () const { return m_gdb != nullptr; }

    virtual void reserveData ();
    virtual void resizeData ();
    void RedefineDummyMF (int lev);

    MFIter MakeMFIter (int lev, const MFItInfo& info) const {
        AMREX_ASSERT(m_dummy_mf[lev] != nullptr);
        return MFIter(*m_dummy_mf[lev], info);
    }

    MFIter MakeMFIter (int lev) const {
        AMREX_ASSERT(m_dummy_mf[lev] != nullptr);
        return MFIter(*m_dummy_mf[lev], do_tiling ? tile_size : IntVect::TheZeroVector());
    }

    MFIter MakeMFIter (int lev, bool tile) const {
        AMREX_ASSERT(m_dummy_mf[lev] != nullptr);
        return MFIter(*m_dummy_mf[lev], tile ? tile_size : IntVect::TheZeroVector());
    }

    //! \brief Set the particle Geometry, DistributionMapping, and BoxArray. If the container
    //! was previously set to to track the AMR hierarchy of an AmrCore or AmrLevel object, that
    //! correspondence will be broken here. This is the single-level version.
    //!
    //! \param geom The new Geometry to use.
    //! \param dmap The new DistributionMapping to use.
    //! \param ba The new BoxArray to use.
    //!
    void SetParGDB (const Geometry            & geom,
                    const DistributionMapping & dmap,
                    const BoxArray            & ba);

    //! \brief Set the particle Geometry, DistributionMapping, ref ratios, and BoxArray. If the container
    //! was previously set to to track the AMR hierarchy of an AmrCore or AmrLevel object, that
    //! correspondence will be broken here. This is the multi-level version.
    //!
    //! \param geom The new vector of Geometry objects to use.
    //! \param dmap The new vector of DistributionMapping objects to use.
    //! \param ba The new vector of BoxArray objects to use.
    //! \param rr The new vector of refinement ratios to use.
    //!
    void SetParGDB (const Vector<Geometry>            & geom,
                    const Vector<DistributionMapping> & dmap,
                    const Vector<BoxArray>            & ba,
                    const Vector<IntVect>             & rr);

    //! \brief Set the particle Geometry, DistributionMapping, ref ratios, and BoxArray. If the container
    //! was previously set to to track the AMR hierarchy of an AmrCore or AmrLevel object, that
    //! correspondence will be broken here. This is the multi-level version.
    //!
    //! Like the above, except the refinement ratios are expressed as ints
    //!
    //! \param geom The new vector of Geometry objects to use.
    //! \param dmap The new vector of DistributionMapping objects to use.
    //! \param ba The new vector of BoxArray objects to use.
    //! \param rr The new vector of refinement ratios to use.
    //!
    void SetParGDB (const Vector<Geometry>            & geom,
                    const Vector<DistributionMapping> & dmap,
                    const Vector<BoxArray>            & ba,
                    const Vector<int>                 & rr);

    //! \brief Set the particle BoxArray. If the container was previously set to
    //! to track the AMR hierarchy of an AmrCore or AmrLevel object, that correspondence
    //! will be broken here.
    //!
    //! \param lev The level on which to set the BoxArray.
    //! \param new_ba The new BoxArray to use.
    //!
    void SetParticleBoxArray (int lev, BoxArray new_ba);

    //! \brief Set the particle DistributionMapping. If the container was previously set to
    //! to track the AMR hierarchy of an AmrCore or AmrLevel object, that correspondence
    //! will be broken here.
    //!
    //! \param lev The level on which to set the dm.
    //! \param new_dmap The new DistributionMapping to use.
    //!
    void SetParticleDistributionMap (int lev, DistributionMapping new_dmap);

    //! \brief Set the particle Geometry. If the container was previously set to
    //! to track the AMR hierarchy of an AmrCore or AmrLevel object, that correspondence
    //! will be broken here.
    //!
    //! \param lev The level on which to set the Geometry.
    //! \param new_geom The new Geometry to use.
    //!
    void SetParticleGeometry (int lev, Geometry new_geom);

    //! \brief Get the BoxArray for a given level
    //!
    //! \param lev The level.
    //!
    const BoxArray& ParticleBoxArray (int lev) const
        { return m_gdb->ParticleBoxArray(lev); }

    //! \brief Get the DistributionMapping for a given level
    //!
    //! \param lev The level.
    //!
    const DistributionMapping& ParticleDistributionMap (int lev) const
        { return m_gdb->ParticleDistributionMap(lev); }

    //! \brief Get the Geometry for a given level
    //!
    //! \param lev The level.
    //!
    const Geometry& Geom (int lev) const { return m_gdb->ParticleGeom(lev); }

    //! \brief Get the particle Geometry for a given level
    //!
    //! \param lev The level.
    //!
    const Geometry& ParticleGeom (int lev) const { return m_gdb->ParticleGeom(lev); }

    //! \brief the finest level actually defined for the ParticleContainer
    int finestLevel () const { return m_gdb->finestLevel(); }

    //! \brief the finest allowed level in the ParticleContainer, whether it is defined or not.
    int maxLevel ()    const { return m_gdb->maxLevel(); }

    //! \brief the number of defined levels in the ParticleContainer
    int numLevels()    const { return finestLevel() + 1; }

    //! \brief Get the ParGDB object used to define this container (const version)
    const ParGDBBase* GetParGDB () const { return m_gdb; }

    //! \brief Get the ParGDB object used to define this container
    ParGDBBase* GetParGDB ()       { return m_gdb; }

    int Verbose () const { return m_verbose; }

    void SetVerbose (int verbose) { m_verbose = verbose; }

    const ParticleBufferMap& BufferMap () const {return m_buffer_map;}

    Vector<int> NeighborProcs(int ngrow) const
    {
        return computeNeighborProcs(this->GetParGDB(), ngrow);
    }

    template <class MF>
    bool OnSameGrids (int level, const MF& mf) const { return m_gdb->OnSameGrids(level, mf); }

    static const std::string& CheckpointVersion ();
    static const std::string& PlotfileVersion ();
    static const std::string& DataPrefix ();
    static int MaxReaders ();
    static Long MaxParticlesPerRead ();
    static const std::string& AggregationType ();
    static int AggregationBuffer ();

    static AMREX_EXPORT bool do_tiling;
    static AMREX_EXPORT IntVect tile_size;
    static AMREX_EXPORT bool memEfficientSort;
    mutable AmrParticleLocator<DenseBins<Box> > m_particle_locator;

protected:

    void BuildRedistributeMask (int lev, int nghost=1) const;
    void defineBufferMap () const;

    int         m_verbose{0};
    std::unique_ptr<ParGDB> m_gdb_object = std::make_unique<ParGDB>();
    ParGDBBase* m_gdb{nullptr};
    Vector<std::unique_ptr<MultiFab> > m_dummy_mf;

    mutable std::unique_ptr<iMultiFab> redistribute_mask_ptr;
    mutable int redistribute_mask_nghost = std::numeric_limits<int>::min();
    mutable amrex::Vector<int> neighbor_procs;
    mutable ParticleBufferMap m_buffer_map;

};

} // namespace amrex

#endif // AMREX_PARTICLECONTAINERBASE_H_
