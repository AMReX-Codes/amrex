#ifndef AMREX_PLOTFILEUTILHDF5_H_
#define AMREX_PLOTFILEUTILHDF5_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileDataImpl.H>

#include <hdf5.h>

#include <string>
#include <memory>

/**
 * \file AMReX_PlotFileUtilHDF5.H
 *
 * HDF5-backed plotfile writers mirroring the ASCII/Binary interfaces.
 */

namespace amrex
{

    /**
     * \brief Write a single-level plotfile (one dataset per level).
     *
     * \param plotfilename  Output directory (e.g., `"plt00000"`).
     * \param mf            MultiFab whose data are stored.
     * \param varnames      Names assigned to each component of \p mf.
     * \param geom          Geometry describing the level.
     * \param t             Simulation time stored in the header.
     * \param level_step    Level timestep counter.
     * \param compression   Compression descriptor (e.g., `"ZFP_RATE@4"`).
     * \param versionName   Metadata version string written to the header.
     * \param levelPrefix   Directory prefix used for each level (default `Level_`).
     * \param mfPrefix      Dataset prefix describing the MultiFab type.
     * \param extra_dirs    Optional sibling directories to create alongside the plotfile.
     */
    void WriteSingleLevelPlotfileHDF5 (const std::string &plotfilename,
                                       const MultiFab &mf,
                                       const Vector<std::string> &varnames,
                                       const Geometry &geom,
                                       Real t,
                                       int level_step,
                                       const std::string &compression = "None@0",
                                       const std::string &versionName = "HyperCLaw-V1.1",
                                       const std::string &levelPrefix = "Level_",
                                       const std::string &mfPrefix = "Cell",
                                       const Vector<std::string>& extra_dirs = Vector<std::string>());

    /**
     * \brief Write a multi-level plotfile (one dataset per level).
     *
     * \param plotfilename  Output directory (e.g., `"plt00000"`).
     * \param nlevels       Number of AMR levels provided.
     * \param mf            Vector of levelwise MultiFabs (component counts must match).
     * \param varnames      Names assigned to each component.
     * \param geom          Vector of per-level geometries.
     * \param time          Simulation time stored in the header.
     * \param level_steps   Per-level timestep counters.
     * \param ref_ratio     Refinement ratios between levels (fine/coarse).
     * \param compression   Compression descriptor (e.g., `"None@0"`).
     * \param versionName   Metadata version string written to the header.
     * \param levelPrefix   Directory prefix used for each level.
     * \param mfPrefix      Dataset prefix describing the MultiFab type.
     * \param extra_dirs    Optional sibling directories to create alongside the plotfile.
     */
    void WriteMultiLevelPlotfileHDF5 (const std::string &plotfilename,
                                      int nlevels,
                                      const Vector<const MultiFab*> &mf,
                                      const Vector<std::string> &varnames,
                                      const Vector<Geometry> &geom,
                                      Real time,
                                      const Vector<int> &level_steps,
                                      const Vector<IntVect> &ref_ratio,
                                      const std::string &compression = "None@0",
                                      const std::string &versionName = "HyperCLaw-V1.1",
                                      const std::string &levelPrefix = "Level_",
                                      const std::string &mfPrefix = "Cell",
                                      const Vector<std::string>& extra_dirs = Vector<std::string>());

    /**
     * \brief Write a single-level plotfile where all state lives in one dataset.
     *
     * \param plotfilename  Output directory (e.g., `"plt00000"`).
     * \param mf            MultiFab whose data are stored.
     * \param varnames      Names assigned to each component of \p mf.
     * \param geom          Geometry describing the level.
     * \param t             Simulation time stored in the header.
     * \param level_step    Level timestep counter.
     * \param compression   Compression descriptor (e.g., `"ZFP_RATE@4"`).
     * \param versionName   Metadata version string written to the header.
     * \param levelPrefix   Directory prefix used for each level (default `Level_`).
     * \param mfPrefix      Dataset prefix describing the MultiFab type.
     * \param extra_dirs    Optional sibling directories to create alongside the plotfile.
     */
    void WriteSingleLevelPlotfileHDF5SingleDset (const std::string &plotfilename,
                                                 const MultiFab &mf,
                                                 const Vector<std::string> &varnames,
                                                 const Geometry &geom,
                                                 Real t,
                                                 int level_step,
                                                 const std::string &compression = "None@0",
                                                 const std::string &versionName = "HyperCLaw-V1.1",
                                                 const std::string &levelPrefix = "Level_",
                                                 const std::string &mfPrefix = "Cell",
                                                 const Vector<std::string>& extra_dirs = Vector<std::string>());

    /**
     * \brief Write a multi-level plotfile using a single dataset per level.
     *
     * \param plotfilename  Output directory (e.g., `"plt00000"`).
     * \param nlevels       Number of AMR levels provided.
     * \param mf            Vector of levelwise MultiFabs (component counts must match).
     * \param varnames      Names assigned to each component.
     * \param geom          Vector of per-level geometries.
     * \param time          Simulation time stored in the header.
     * \param level_steps   Per-level timestep counters.
     * \param ref_ratio     Refinement ratios between levels (fine/coarse).
     * \param compression   Compression descriptor (e.g., `"None@0"`).
     * \param versionName   Metadata version string written to the header.
     * \param levelPrefix   Directory prefix used for each level.
     * \param mfPrefix      Dataset prefix describing the MultiFab type.
     * \param extra_dirs    Optional sibling directories to create alongside the plotfile.
     */
    void WriteMultiLevelPlotfileHDF5SingleDset (const std::string &plotfilename,
                                                int nlevels,
                                                const Vector<const MultiFab*> &mf,
                                                const Vector<std::string> &varnames,
                                                const Vector<Geometry> &geom,
                                                Real time,
                                                const Vector<int> &level_steps,
                                                const Vector<IntVect> &ref_ratio,
                                                const std::string &compression = "None@0",
                                                const std::string &versionName = "HyperCLaw-V1.1",
                                                const std::string &levelPrefix = "Level_",
                                                const std::string &mfPrefix = "Cell",
                                                const Vector<std::string>& extra_dirs = Vector<std::string>());

    /**
     * \brief Write a single-level plotfile with one HDF5 dataset per component.
     *
     * \param plotfilename  Output directory (e.g., `"plt00000"`).
     * \param mf            MultiFab whose data are stored.
     * \param varnames      Names assigned to each component of \p mf.
     * \param geom          Geometry describing the level.
     * \param t             Simulation time stored in the header.
     * \param level_step    Level timestep counter.
     * \param compression   Compression descriptor (e.g., `"ZFP_RATE@4"`).
     * \param versionName   Metadata version string written to the header.
     * \param levelPrefix   Directory prefix used for each level (default `Level_`).
     * \param mfPrefix      Dataset prefix describing the MultiFab type.
     * \param extra_dirs    Optional sibling directories to create alongside the plotfile.
     */
    void WriteSingleLevelPlotfileHDF5MultiDset (const std::string &plotfilename,
                                                const MultiFab &mf,
                                                const Vector<std::string> &varnames,
                                                const Geometry &geom,
                                                Real t,
                                                int level_step,
                                                const std::string &compression = "None@0",
                                                const std::string &versionName = "HyperCLaw-V1.1",
                                                const std::string &levelPrefix = "Level_",
                                                const std::string &mfPrefix = "Cell",
                                                const Vector<std::string>& extra_dirs = Vector<std::string>());

    /**
     * \brief Write a multi-level plotfile with one dataset per component.
     *
     * \param plotfilename  Output directory (e.g., `"plt00000"`).
     * \param nlevels       Number of AMR levels provided.
     * \param mf            Vector of levelwise MultiFabs (component counts must match).
     * \param varnames      Names assigned to each component.
     * \param geom          Vector of per-level geometries.
     * \param time          Simulation time stored in the header.
     * \param level_steps   Per-level timestep counters.
     * \param ref_ratio     Refinement ratios between levels (fine/coarse).
     * \param compression   Compression descriptor (e.g., `"None@0"`).
     * \param versionName   Metadata version string written to the header.
     * \param levelPrefix   Directory prefix used for each level.
     * \param mfPrefix      Dataset prefix describing the MultiFab type.
     * \param extra_dirs    Optional sibling directories to create alongside the plotfile.
     */
    void WriteMultiLevelPlotfileHDF5MultiDset (const std::string &plotfilename,
                                               int nlevels,
                                               const Vector<const MultiFab*> &mf,
                                               const Vector<std::string> &varnames,
                                               const Vector<Geometry> &geom,
                                               Real time,
                                               const Vector<int> &level_steps,
                                               const Vector<IntVect> &ref_ratio,
                                               const std::string &compression = "None@0",
                                               const std::string &versionName = "HyperCLaw-V1.1",
                                               const std::string &levelPrefix = "Level_",
                                               const std::string &mfPrefix = "Cell",
                                               const Vector<std::string>& extra_dirs = Vector<std::string>());
}

#endif
