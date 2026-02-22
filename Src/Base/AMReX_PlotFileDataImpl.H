#ifndef AMREX_PLOT_FILE_DATA_IMPL_H_
#define AMREX_PLOT_FILE_DATA_IMPL_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <string>

/**
 * \file AMReX_PlotFileDataImpl.H
 * \brief Internal reader backing PlotFileData convenience wrapper.
 */

namespace amrex {

class PlotFileDataImpl
{
public:
    /**
     * \brief Construct a reader for an AMReX plotfile on disk.
     * \param plotfile_name Path to the pltXXXX directory or file.
     */
    PlotFileDataImpl (std::string const& plotfile_name);

    //! Spatial dimensionality encoded in the plotfile metadata.
    [[nodiscard]] int spaceDim () const noexcept { return m_spacedim; }

    //! Simulation time stored in the header.
    [[nodiscard]] Real time () const noexcept { return m_time; }

    //! Finest AMR level present in the plotfile.
    [[nodiscard]] int finestLevel () const noexcept { return m_finest_level; }

    //! Refinement ratio between level and level+1.
    [[nodiscard]] int refRatio (int level) const noexcept { return m_ref_ratio[level]; }

    //! AMR level step count written to the file.
    [[nodiscard]] int levelStep (int level) const noexcept { return m_level_steps[level]; }

    //! BoxArray describing grid decomposition at a given level.
    [[nodiscard]] const BoxArray& boxArray (int level) const noexcept { return m_ba[level]; }

    //! Distribution map used when the plotfile was written.
    [[nodiscard]] const DistributionMapping& DistributionMap (int level) const noexcept { return m_dmap[level]; }

    /**
     * \brief Copy distribution maps from \p src for all levels.
     * \param src PlotFileDataImpl supplying the new maps.
     */
    void syncDistributionMap (PlotFileDataImpl const& src) noexcept;

    /**
     * \brief Copy the distribution map for a single level.
     * \param level AMR level to update.
     * \param src   PlotFileDataImpl supplying the new map.
     */
    void syncDistributionMap (int level, PlotFileDataImpl const& src) noexcept;

    //! Coordinate-system flag stored in the header.
    [[nodiscard]] int coordSys () const noexcept { return m_coordsys; }

    //! Problem domain at the given level.
    [[nodiscard]] Box probDomain (int level) const noexcept { return m_prob_domain[level]; }

    //! Return the physical extents of the domain (probHi - probLo).
    [[nodiscard]] Array<Real,AMREX_SPACEDIM> probSize () const noexcept { return m_prob_size; }
    //! Return the physical lower corner of the domain.
    [[nodiscard]] Array<Real,AMREX_SPACEDIM> probLo () const noexcept { return m_prob_lo; }
    //! Return the physical upper corner of the domain.
    [[nodiscard]] Array<Real,AMREX_SPACEDIM> probHi () const noexcept { return m_prob_hi; }
    //! Return the cell spacing at the requested level.
    [[nodiscard]] Array<Real,AMREX_SPACEDIM> cellSize (int level) const noexcept { return m_cell_size[level]; }

    //! Variable names stored in the plotfile.
    [[nodiscard]] const Vector<std::string>& varNames () const noexcept { return m_var_names; }

    //! Number of components per MultiFab.
    [[nodiscard]] int nComp () const noexcept { return m_ncomp; }
    //! Number of guard cells stored with each level's data.
    [[nodiscard]] IntVect nGrowVect (int level) const noexcept { return m_ngrow[level]; }

    /**
     * \brief Load all components for the requested level into a new MultiFab.
     * \param level AMR level index.
     * \return MultiFab with the layout of the plotfile grids.
     */
    MultiFab get (int level) noexcept;
    /**
     * \brief Load a single component by name for the requested level.
     * \param level AMR level index.
     * \param varname Name listed in varNames().
     * \return MultiFab containing only the requested component.
     */
    MultiFab get (int level, std::string const& varname) noexcept;

private:
    std::string m_plotfile_name;
    std::string m_file_version;
    int m_ncomp;
    Vector<std::string> m_var_names;
    int m_spacedim;
    Real m_time;
    int m_finest_level, m_nlevels;
    Array<Real,AMREX_SPACEDIM> m_prob_lo {{AMREX_D_DECL(0.,0.,0.)}};
    Array<Real,AMREX_SPACEDIM> m_prob_hi {{AMREX_D_DECL(1.,1.,1.)}};
    Array<Real,AMREX_SPACEDIM> m_prob_size {{AMREX_D_DECL(1.,1.,1.)}};
    Vector<int> m_ref_ratio;
    Vector<Box> m_prob_domain;
    Vector<int> m_level_steps;
    Vector<Array<Real,AMREX_SPACEDIM> > m_cell_size;
    int m_coordsys;
    Vector<std::string> m_mf_name;
    Vector<std::unique_ptr<VisMF> > m_vismf;
    Vector<BoxArray> m_ba;
    Vector<DistributionMapping> m_dmap;
    Vector<IntVect> m_ngrow;
};

}

#endif
