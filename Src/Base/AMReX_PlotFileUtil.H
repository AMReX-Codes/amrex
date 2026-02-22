#ifndef AMREX_PLOTFILEUTIL_H_
#define AMREX_PLOTFILEUTIL_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_PlotFileUtil.H
 * \brief Utilities for writing AMReX-format plot files.
 */

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileDataImpl.H>
#include <AMReX_Vector.H>

#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif

#include <ostream>
#include <string>
#include <memory>

namespace amrex
{
    /**
     * \brief Return the directory name for a level (e.g., "Level_5").
     * \param level Zero-based AMR level index.
     * \param levelPrefix Prefix prepended to the integer index.
     */
    std::string LevelPath (int level, const std::string &levelPrefix = "Level_");

    /**
     * \brief Return the relative MultiFab path recorded in the header (e.g., "Level_5/Cell").
     * \param level Zero-based AMR level index.
     * \param levelPrefix Prefix used for directories (defaults to "Level_").
     * \param mfPrefix    Name of the data prefix ("Cell" for cell-centered data).
     */
    std::string MultiFabHeaderPath (int level,
                                    const std::string &levelPrefix = "Level_",
                                    const std::string &mfPrefix = "Cell");

    /**
     * \brief Return the absolute level directory path (e.g., "plt00005/Level_5").
     * \param level Zero-based AMR level index.
     * \param plotfilename Name of the plotfile root directory.
     * \param levelPrefix Prefix used for directories (defaults to "Level_").
     */
    std::string LevelFullPath (int level,
                               const std::string &plotfilename,
                               const std::string &levelPrefix = "Level_");

    /**
     * \brief Return the absolute path prefix for MultiFab data (e.g., "plt00005/Level_5/Cell").
     * \param level Zero-based AMR level index.
     * \param plotfilename Name of the plotfile root directory.
     * \param levelPrefix Prefix used for directories (defaults to "Level_").
     * \param mfPrefix    Component prefix ("Cell", "Node", etc.).
     */
    std::string MultiFabFileFullPrefix (int level,
                                        const std::string &plotfilename,
                                        const std::string &levelPrefix = "Level_",
                                        const std::string &mfPrefix = "Cell");


    /**
     * \brief Create a hierarchy of plotfile directories on disk.
     *
     * Builds \p dirName and renames it if it already exists, then creates the
     * numbered subdirectories `dirName/subDirPrefix_i` for `i` in `[0,nSubDirs)`.
     * MPI rank 0 (IOProcessor) performs the filesystem work; optionally, all
     * ranks can synchronize via ParallelDescriptor::Barrier().
     *
     * \param dirName     Root directory to create.
     * \param subDirPrefix Prefix for per-level subdirectories (e.g., "Level_").
     * \param nSubDirs    Number of subdirectories to create.
     * \param callBarrier Whether to synchronize all ranks after creation.
     */
    void PreBuildDirectorHierarchy (const std::string &dirName,
                                    const std::string &subDirPrefix,
                                    int nSubDirs,
                                    bool callBarrier);

    /**
     * \brief Write the ascii header that describes a multi-level plot file.
     *
     * The target plotfile directory must already exist, and the caller passes
     * header metadata describing the AMR hierarchy, component list, timesteps,
     * and refinement ratios.
     *
     * \param HeaderFile Stream already opened on `<plotfilename>/Header`.
     * \param nlevels Number of AMR levels present.
     * \param bArray Per-level BoxArrays.
     * \param varnames Component names.
     * \param geom Per-level geometries.
     * \param time Simulation time.
     * \param level_steps Per-level step counters.
     * \param ref_ratio Refinement ratios between levels.
     * \param versionName ASCII version string stored in Header.
     * \param levelPrefix Directory prefix for each level.
     * \param mfPrefix    MultiFab component prefix.
     */
    void WriteGenericPlotfileHeader (std::ostream &HeaderFile,
                                     int nlevels,
                                     const Vector<BoxArray> &bArray,
                                     const Vector<std::string> &varnames,
                                     const Vector<Geometry> &geom,
                                     Real time,
                                     const Vector<int> &level_steps,
                                     const Vector<IntVect> &ref_ratio,
                                     const std::string &versionName = "HyperCLaw-V1.1",
                                     const std::string &levelPrefix = "Level_",
                                     const std::string &mfPrefix = "Cell");

    /// \ingroup amrex_io
    /**
     * \brief Write a single-level AMReX plot file for \p mf.
     *
     * \param plotfilename Destination directory (created if needed).
     * \param mf           MultiFab holding the cell data to dump.
     * \param varnames     Component names written to the header.
     * \param geom         Geometry describing the domain.
     * \param t            Simulation time stamp stored in the header.
     * \param level_step   Step number recorded for the level.
     * \param versionName  ASCII version string stored in Header.
     * \param levelPrefix  Directory prefix for this level (defaults to "Level_").
     * \param mfPrefix     Component prefix ("Cell", "Node", etc.).
     * \param extra_dirs   Optional directories created alongside the plotfile.
     */
    void WriteSingleLevelPlotfile (const std::string &plotfilename,
                                   const MultiFab &mf,
                                   const Vector<std::string> &varnames,
                                   const Geometry &geom,
                                   Real t,
                                   int level_step,
                                   const std::string &versionName = "HyperCLaw-V1.1",
                                   const std::string &levelPrefix = "Level_",
                                   const std::string &mfPrefix = "Cell",
                                   const Vector<std::string>& extra_dirs = Vector<std::string>());

    /// \ingroup amrex_io
    /**
     * \brief Write a multi-level AMReX plot file.
     *
     * \param plotfilename Destination directory (created if needed).
     * \param nlevels      Number of AMR levels written.
     * \param mf           Vector of MultiFab pointers, one per level.
     * \param varnames     Component names recorded in the header.
     * \param geom         Per-level geometries.
     * \param time         Simulation time stamp stored in the header.
     * \param level_steps  Per-level step counters.
     * \param ref_ratio    Refinement ratios between levels.
     * \param versionName  ASCII version string stored in Header.
     * \param levelPrefix  Directory prefix for level subdirectories.
     * \param mfPrefix     Component prefix (e.g., "Cell").
     * \param extra_dirs   Optional directories created alongside the plotfile.
     */
    void WriteMultiLevelPlotfile (const std::string &plotfilename,
                                  int nlevels,
                                  const Vector<const MultiFab*> &mf,
                                  const Vector<std::string> &varnames,
                                  const Vector<Geometry> &geom,
                                  Real time,
                                  const Vector<int> &level_steps,
                                  const Vector<IntVect> &ref_ratio,
                                  const std::string &versionName = "HyperCLaw-V1.1",
                                  const std::string &levelPrefix = "Level_",
                                  const std::string &mfPrefix = "Cell",
                                  const Vector<std::string>& extra_dirs = Vector<std::string>());

    /**
     * \brief Convenience wrapper that writes a plotfile with generic variable names.
     *
     * Each component is named "Var<i>", refinement ratios are inferred from the
     * Geometry vector, and both time and level_steps are written as zero.
     *
     * \param plotfilename Destination directory.
     * \param mf           MultiFabs (one per level) to dump.
     * \param geom         Per-level geometries.
     */
    void WriteMLMF (const std::string &plotfilename,
                    const Vector<const MultiFab*>& mf,
                    const Vector<Geometry> &geom);

    /// \ingroup amrex_io
    /**
     * \brief Write only the ASCII headers for a multi-level plot file.
     *
     * The MultiFabs referenced by \p mf are used for metadata (BoxArrays, NGrow),
     * but the function assumes the raw data already exists on disk.
     *
     * \param plotfilename Destination directory.
     * \param nlevels      Number of AMR levels present.
     * \param mf           MultiFabs used to derive headers (one per level).
     * \param varnames     Component names.
     * \param geom         Per-level geometries.
     * \param time         Simulation time stamp.
     * \param level_steps  Per-level step counters.
     * \param ref_ratio    Refinement ratios between levels.
     * \param versionName  ASCII version string stored in Header.
     * \param levelPrefix  Directory prefix for levels.
     * \param mfPrefix     Component prefix.
     * \param extra_dirs   Optional directories created alongside the plotfile.
     */
    void WriteMultiLevelPlotfileHeaders (const std::string& plotfilename, int nlevels,
                                         const Vector<const MultiFab*>& mf,
                                         const Vector<std::string>& varnames,
                                         const Vector<Geometry>& geom, Real time, const Vector<int>& level_steps,
                                         const Vector<IntVect>& ref_ratio,
                                         const std::string &versionName = "HyperCLaw-V1.1",
                                         const std::string &levelPrefix = "Level_",
                                         const std::string &mfPrefix = "Cell",
                                         const Vector<std::string>& extra_dirs = Vector<std::string>());


#ifdef AMREX_USE_EB
    /// \ingroup amrex_io
    /**
     * \brief Write a single-level plotfile augmented with embedded-boundary data.
     *
     * \param plotfilename Destination directory (created if needed).
     * \param mf           Cell-centered MultiFab containing solution data.
     * \param varnames     Component names.
     * \param geom         Geometry describing the domain.
     * \param t            Simulation time stamp.
     * \param level_step   Step number recorded for the level.
     * \param versionName  ASCII version string stored in Header.
     * \param levelPrefix  Directory prefix used for levels.
     * \param mfPrefix     Component prefix written to disk.
     * \param extra_dirs   Optional sibling directories to create.
     */
    void EB_WriteSingleLevelPlotfile (const std::string &plotfilename,
                                      const MultiFab &mf,
                                      const Vector<std::string> &varnames,
                                      const Geometry &geom,
                                      Real t,
                                      int level_step,
                                      const std::string &versionName = "HyperCLaw-V1.1",
                                      const std::string &levelPrefix = "Level_",
                                      const std::string &mfPrefix = "Cell",
                                      const Vector<std::string>& extra_dirs = Vector<std::string>());

    /// \ingroup amrex_io
    /**
     * \brief Write a multi-level EB-aware plotfile.
     *
     * \param plotfilename Destination directory.
     * \param nlevels      Number of AMR levels written.
     * \param mf           Per-level MultiFabs containing solution data.
     * \param varnames     Component names.
     * \param geom         Per-level geometries.
     * \param time         Simulation time stamp.
     * \param level_steps  Per-level step counters.
     * \param ref_ratio    Refinement ratios between levels.
     * \param versionName  ASCII version string stored in Header.
     * \param levelPrefix  Directory prefix used for levels.
     * \param mfPrefix     Component prefix written to disk.
     * \param extra_dirs   Optional sibling directories to create.
     */
    void EB_WriteMultiLevelPlotfile (const std::string &plotfilename,
                                     int nlevels,
                                     const Vector<const MultiFab*> &mf,
                                     const Vector<std::string> &varnames,
                                     const Vector<Geometry> &geom,
                                     Real time,
                                     const Vector<int> &level_steps,
                                     const Vector<IntVect> &ref_ratio,
                                     const std::string &versionName = "HyperCLaw-V1.1",
                                     const std::string &levelPrefix = "Level_",
                                     const std::string &mfPrefix = "Cell",
                                     const Vector<std::string>& extra_dirs = Vector<std::string>());

#endif

    /**
     * \class amrex::PlotFileData
     * \ingroup amrex_io
     * \brief Read-only handle to an AMReX-format plot file on disk.
     *
     * Wraps PlotFileDataImpl to provide metadata queries and MultiFab extraction
     * without exposing the implementation details.
     */
    class PlotFileData
    {
    public:
        //! Open the plot file at \p plotfile_name for reading.
        PlotFileData (std::string const& plotfile_name) : m_impl(new PlotFileDataImpl(plotfile_name)) {}

        //! Return the spatial dimension of the data.
        [[nodiscard]] int spaceDim () const noexcept { return m_impl->spaceDim(); }
        //! Return the simulation time stored in the plot file.
        [[nodiscard]] Real time () const noexcept { return m_impl->time(); }
        //! Return the finest AMR level index.
        [[nodiscard]] int finestLevel () const noexcept { return m_impl->finestLevel(); }
        //! Return the refinement ratio between \p level and level+1.
        [[nodiscard]] int refRatio (int level) const noexcept { return m_impl->refRatio(level); }
        //! Return the step number recorded at \p level.
        [[nodiscard]] int levelStep (int level) const noexcept { return m_impl->levelStep(level); }
        //! Return the BoxArray at \p level.
        [[nodiscard]] const BoxArray& boxArray (int level) const noexcept { return m_impl->boxArray(level); }
        //! Return the DistributionMapping at \p level.
        [[nodiscard]] const DistributionMapping& DistributionMap (int level) const noexcept { return m_impl->DistributionMap(level); }
        //! Copy the DistributionMapping from \p src to all levels of this object.
        void syncDistributionMap (PlotFileData const& src) noexcept { m_impl->syncDistributionMap(*src.m_impl); }
        //! Copy the DistributionMapping for \p level from \p src to this object.
        void syncDistributionMap (int level, PlotFileData const& src) noexcept { m_impl->syncDistributionMap(level, *src.m_impl); }
        //! Return the coordinate system type.
        [[nodiscard]] int coordSys () const noexcept { return m_impl->coordSys(); }
        //! Return the problem domain Box at \p level.
        [[nodiscard]] Box probDomain (int level) const noexcept { return m_impl->probDomain(level); }
        //! Return the physical domain size along each dimension.
        [[nodiscard]] Array<Real,AMREX_SPACEDIM> probSize () const noexcept { return m_impl->probSize(); }
        //! Return the lower corner of the physical domain.
        [[nodiscard]] Array<Real,AMREX_SPACEDIM> probLo () const noexcept { return m_impl->probLo(); }
        //! Return the upper corner of the physical domain.
        [[nodiscard]] Array<Real,AMREX_SPACEDIM> probHi () const noexcept { return m_impl->probHi(); }
        //! Return the cell size at \p level.
        [[nodiscard]] Array<Real,AMREX_SPACEDIM> cellSize (int level) const noexcept { return m_impl->cellSize(level); }
        //! Return the list of variable names.
        [[nodiscard]] const Vector<std::string>& varNames () const noexcept { return m_impl->varNames(); }
        //! Return the number of components.
        [[nodiscard]] int nComp () const noexcept { return m_impl->nComp(); }
        //! Return the ghost cell vector at \p level.
        [[nodiscard]] IntVect nGrowVect (int level) const noexcept { return m_impl->nGrowVect(level); }
        //! Read and return all components of the MultiFab at \p level.
        MultiFab get (int level) noexcept { return m_impl->get(level); }
        //! Read and return the single component named \p varname at \p level.
        MultiFab get (int level, std::string const& varname) noexcept { return m_impl->get(level, varname); }

    private:
        std::unique_ptr<PlotFileDataImpl> m_impl;
    };
}

#endif
