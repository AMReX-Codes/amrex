#ifndef BL_AMRCORE_H_
#define BL_AMRCORE_H_

#include <ostream>
#include <memory>

#include <AMReX_AmrMesh.H>

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
    AmrCore (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord=-1, Vector<IntVect> ref_ratios = Vector<IntVect>());

    AmrCore (const AmrCore& rhs) = delete;
    AmrCore& operator= (const AmrCore& rhs) = delete;

    virtual ~AmrCore ();

#ifdef AMREX_PARTICLES
    AmrParGDB* GetParGDB () const noexcept { return m_gdb.get(); }
#endif

    /**
     * \brief Initialize BoxArray, DistributionMapping and data from scratch.
     * Calling this function requires the derive class implement its own MakeNewLevelFromScratch
     * to allocate and intialize data.
     * Also note usually one needs to average the fine data down to coarse level after this.
     */
    void InitFromScratch (Real time);

    //! Rebuild levels finer than lbase
    virtual void regrid (int lbase, Real time, bool initial=false);

    static void Initialize ();
    static void Finalize ();

    void printGridSummary (std::ostream& os, int min_lev, int max_lev) const noexcept;

    int Verbose () const noexcept { return verbose; }

protected:

    //! Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override = 0;

    //! Make a new level from scratch using provided BoxArray and DistributionMapping.
    //! Only used during initialization.
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override = 0;

    //! Make a new level using provided BoxArray and DistributionMapping and fill with interpolated coarse level data.
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    //! Remake an existing level using provided BoxArray and DistributionMapping and fill with existing fine and coarse data.
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    //! Delete level data
    virtual void ClearLevel (int lev) = 0;

    int              verbose;

#ifdef AMREX_PARTICLES
    std::unique_ptr<AmrParGDB> m_gdb;
#endif

private:
    void InitAmrCore ();
};

}

#endif
