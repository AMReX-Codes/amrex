#ifndef BL_AMRCORE_H_
#define BL_AMRCORE_H_
#include <AMReX_Config.H>

#include <AMReX_AmrMesh.H>

#include <iosfwd>
#include <memory>

namespace amrex {

#ifdef AMREX_PARTICLES
class AmrParGDB;
#endif

/**
 * \brief Provide basic functionalities to set up an AMR hierarchy
 *
 * Because it does not own any floating point data, AmrCore used
 * virtual functions to allocate, initialize and delete data.  It also
 * requires the derived class to tag cells for refinement.
 */
class AmrCore
    : public AmrMesh
{
public:

    AmrCore ();

    AmrCore (const RealBox* rb, int max_level_in,
             const Vector<int>& n_cell_in, int coord=-1,
             Vector<IntVect> ref_ratios = Vector<IntVect>(),
             const int* is_per = nullptr);

    AmrCore (const RealBox& rb, int max_level_in,
             const Vector<int>& n_cell_in, int coord,
             Vector<IntVect> const& ref_ratios,
             Array<int,AMREX_SPACEDIM> const& is_per);

    AmrCore (Geometry const& level_0_geom, AmrInfo const& amr_info);

    AmrCore (AmrCore&& rhs) noexcept;
    AmrCore& operator= (AmrCore&& rhs) noexcept;

    AmrCore (const AmrCore& rhs) = delete;
    AmrCore& operator= (const AmrCore& rhs) = delete;

    ~AmrCore () override;

#ifdef AMREX_PARTICLES
    [[nodiscard]] AmrParGDB* GetParGDB () const noexcept { return m_gdb.get(); }
#endif

    /**
     * \brief Initialize BoxArray, DistributionMapping and data from scratch.
     * Calling this function requires the derive class implement its own MakeNewLevelFromScratch
     * to allocate and initialize data.
     * Also note usually one needs to average the fine data down to coarse level after this.
     */
    void InitFromScratch (Real time);

    //! Rebuild levels finer than lbase
    virtual void regrid (int lbase, Real time, bool initial=false);

    void printGridSummary (std::ostream& os, int min_lev, int max_lev) const noexcept;

protected:

    //! Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
    void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override = 0;

    //! Make a new level from scratch using provided BoxArray and DistributionMapping.
    //! Only used during initialization.
    void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override = 0;

    //! Make a new level using provided BoxArray and DistributionMapping and fill with interpolated coarse level data.
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    //! Remake an existing level using provided BoxArray and DistributionMapping and fill with existing fine and coarse data.
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    //! Delete level data
    virtual void ClearLevel (int lev) = 0;

#ifdef AMREX_PARTICLES
    std::unique_ptr<AmrParGDB> m_gdb;
#endif

private:
    void InitAmrCore ();
};

}

#endif
