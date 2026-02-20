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
    //!  return the name of the level directory,  e.g., Level_5
    std::string LevelPath (int level, const std::string &levelPrefix = "Level_");

    //!  return the path of the multifab to write to the header, e.g., Level_5/Cell
    std::string MultiFabHeaderPath (int level,
                                    const std::string &levelPrefix = "Level_",
                                    const std::string &mfPrefix = "Cell");

    //!  return the full path of the level directory, e.g., plt00005/Level_5
    std::string LevelFullPath (int level,
                               const std::string &plotfilename,
                               const std::string &levelPrefix = "Level_");

    //!  return the full path multifab prefix, e.g., plt00005/Level_5/Cell
    std::string MultiFabFileFullPrefix (int level,
                                        const std::string &plotfilename,
                                        const std::string &levelPrefix = "Level_",
                                        const std::string &mfPrefix = "Cell");


    /**
    * \brief  prebuild a hierarchy of directories
    *  dirName is built first.  if dirName exists, it is renamed.  then build
    *  dirName/subDirPrefix_0 .. dirName/subDirPrefix_nSubDirs-1
    *  if callBarrier is true, call ParallelDescriptor::Barrier()
    *  after all directories are built
    *  ParallelDescriptor::IOProcessor() creates the directories
    *
    * \param &dirName
    * \param &subDirPrefix
    * \param nSubDirs
    * \param callBarrier
    */
    void PreBuildDirectorHierarchy (const std::string &dirName,
                                    const std::string &subDirPrefix,
                                    int nSubDirs,
                                    bool callBarrier);

    /**
    *  write a generic plot file header to the file plotfilename/Header
    *  the plotfilename directory must already exist
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
    //! Write a single-level AMReX plot file for \p mf with the given variable names, geometry, and time.
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
    //! Write a multi-level AMReX plot file for the \p nlevels MultiFabs with the given variable names, geometry, and time.
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
    * \brief write a plotfile to disk given:
    * -plotfile name
    * -vector of MultiFabs
    * -vector of Geometrys
    * variable names are written as "Var0", "Var1", etc.
    * refinement ratio is computed from the Geometry vector
    * "time" and "level_steps" are set to zero
    *
    * \param &plotfilename
    * \param mf
    * \param &geom
    */
    void WriteMLMF (const std::string &plotfilename,
                    const Vector<const MultiFab*>& mf,
                    const Vector<Geometry> &geom);

    /// \ingroup amrex_io
    //! Write only the headers for a multi-level plot file (data files assumed already written).
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
    //! Write a single-level EB-aware plot file (includes volume fractions and centroid data).
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
    //! Write a multi-level EB-aware plot file (includes volume fractions and centroid data).
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
