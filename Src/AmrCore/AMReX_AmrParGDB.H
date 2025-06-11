#ifndef AMREX_AmrParGDB_H
#define AMREX_AmrParGDB_H
#include <AMReX_Config.H>

#include <AMReX_ParGDB.H>
#include <AMReX_AmrCore.H>

namespace amrex {

class AmrParGDB
    : public ParGDBBase
{
    friend AmrCore;

public:

    explicit AmrParGDB (AmrCore* amr) noexcept
        : m_amrcore(amr),
          m_geom(amr->maxLevel()+1),
          m_has_geom(amr->maxLevel()+1, 0),
          m_dmap(amr->maxLevel()+1),
          m_ba(amr->maxLevel()+1)
    { }

    [[nodiscard]] const Geometry& ParticleGeom (int level) const override;
    [[nodiscard]] const Geometry&         Geom (int level) const override;

    [[nodiscard]] const Vector<Geometry>& ParticleGeom () const override;
    [[nodiscard]] const Vector<Geometry>&         Geom () const override;

    [[nodiscard]] const DistributionMapping& ParticleDistributionMap (int level) const override;
    [[nodiscard]] const DistributionMapping&         DistributionMap (int level) const override;

    [[nodiscard]] const Vector<DistributionMapping>& ParticleDistributionMap () const override;
    [[nodiscard]] const Vector<DistributionMapping>&         DistributionMap () const override;

    [[nodiscard]] const BoxArray& ParticleBoxArray (int level) const override;
    [[nodiscard]] const BoxArray&         boxArray (int level) const override;

    [[nodiscard]] const Vector<BoxArray>& ParticleBoxArray () const override;
    [[nodiscard]] const Vector<BoxArray>&         boxArray () const override;

    void SetParticleBoxArray (int level, const BoxArray& new_ba) override;
    void SetParticleDistributionMap (int level, const DistributionMapping& new_dm) override;
    void SetParticleGeometry (int level, const Geometry& new_geom) override;

    void ClearParticleBoxArray (int level) override;
    void ClearParticleDistributionMap (int level) override;
    void ClearParticleGeometry (int level) override;

    [[nodiscard]] bool LevelDefined (int level) const override;
    [[nodiscard]] int finestLevel () const override;
    [[nodiscard]] int maxLevel () const override;

    [[nodiscard]] IntVect refRatio (int level) const override;
    [[nodiscard]] int MaxRefRatio (int level) const override;

    [[nodiscard]] Vector<IntVect> refRatio () const override;

protected:

    AmrCore*                    m_amrcore;
    Vector<Geometry>            m_geom;
    Vector<int>                 m_has_geom;
    Vector<DistributionMapping> m_dmap;
    Vector<BoxArray>            m_ba;
};

inline
const Geometry&
AmrParGDB::ParticleGeom (int level) const
{
    if (! m_has_geom[level]) {
        return m_amrcore->Geom(level);
    } else {
        return m_geom[level];
    }
}

inline
const Geometry&
AmrParGDB::Geom (int level) const
{
    return m_amrcore->Geom(level);
}

inline
const Vector<Geometry>&
AmrParGDB::ParticleGeom () const
{
    if (! m_has_geom[0]) {
        return m_amrcore->Geom();
    } else {
        return m_geom;
    }
}

inline
const Vector<Geometry>&
AmrParGDB::Geom () const
{
    return m_amrcore->Geom();
}

inline
const DistributionMapping&
AmrParGDB::ParticleDistributionMap (int level) const
{
    if (m_dmap[level].empty()) {
        return m_amrcore->DistributionMap(level);
    } else {
        return m_dmap[level];
    }
}

inline
const DistributionMapping&
AmrParGDB::DistributionMap (int level) const
{
    return m_amrcore->DistributionMap(level);
}

inline
const Vector<DistributionMapping>&
AmrParGDB::ParticleDistributionMap () const
{
    if (m_dmap[0].empty()) {
        return m_amrcore->DistributionMap();
    } else {
        return m_dmap;
    }
}

inline
const Vector<DistributionMapping>&
AmrParGDB::DistributionMap () const
{
    return m_amrcore->DistributionMap();
}

inline
const BoxArray&
AmrParGDB::ParticleBoxArray (int level) const
{
    if (m_ba[level].empty()) {
        return m_amrcore->boxArray(level);
    } else {
        return m_ba[level];
    }
}

inline
const BoxArray&
AmrParGDB::boxArray (int level) const
{
    return m_amrcore->boxArray(level);
}

inline
const Vector<BoxArray>&
AmrParGDB::ParticleBoxArray () const
{
    if (m_ba[0].empty()) {
        return m_amrcore->boxArray();
    } else {
        return m_ba;
    }
}

inline
const Vector<BoxArray>&
AmrParGDB::boxArray () const
{
    return m_amrcore->boxArray();
}

inline
void AmrParGDB::SetParticleBoxArray (int level, const BoxArray& new_ba)
{
    m_ba[level] = new_ba;
}

inline
void AmrParGDB::SetParticleDistributionMap (int level, const DistributionMapping& new_dmap)
{
    m_dmap[level] = new_dmap;
}

inline
void AmrParGDB::SetParticleGeometry (int level, const Geometry& new_geom)
{
    m_has_geom[level] = 1;
    m_geom[level] = new_geom;
}

inline
void AmrParGDB::ClearParticleBoxArray (int level)
{
    m_ba[level] = BoxArray();
}

inline
void AmrParGDB::ClearParticleDistributionMap (int level)
{
    m_dmap[level] = DistributionMapping();
}

inline
void AmrParGDB::ClearParticleGeometry (int level)
{
    m_geom[level] = Geometry();
    m_has_geom[level] = 0;
}

inline
bool
AmrParGDB::LevelDefined (int level) const
{
    return m_amrcore->LevelDefined(level);
}

inline
int
AmrParGDB::finestLevel () const
{
    return m_amrcore->finestLevel();
}

inline
int
AmrParGDB::maxLevel () const
{
    return m_amrcore->maxLevel();
}

inline
IntVect
AmrParGDB::refRatio (int level) const
{
    return m_amrcore->refRatio(level);
}

inline
Vector<IntVect>
AmrParGDB::refRatio () const
{
    return m_amrcore->refRatio();
}

inline
int
AmrParGDB::MaxRefRatio (int level) const
{
    return m_amrcore->MaxRefRatio(level);
}

}

#endif
