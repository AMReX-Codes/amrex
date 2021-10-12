#ifndef AMREX_PLOTFILEUTIL_H_
#define AMREX_PLOTFILEUTIL_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileDataImpl.H>

#include <hdf5.h>

#include <string>
#include <memory>

namespace amrex
{

    void WriteGenericPlotfileHeaderHDF5 (hid_t fid,
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

    void WriteSingleLevelPlotfileHDF5 (const std::string &plotfilename,
                                       const MultiFab &mf,
                                       const Vector<std::string> &varnames,
                                       const Geometry &geom,
                                       Real t,
                                       int level_step,
                                       const std::string &versionName = "HyperCLaw-V1.1",
                                       const std::string &levelPrefix = "Level_",
                                       const std::string &mfPrefix = "Cell",
                                       const Vector<std::string>& extra_dirs = Vector<std::string>());

    void WriteMultiLevelPlotfileHDF5 (const std::string &plotfilename,
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

}

#endif
