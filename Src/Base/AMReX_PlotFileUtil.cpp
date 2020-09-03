
#include <fstream>
#include <iomanip>

#include <AMReX_VisMF.H>
#include <AMReX_AsyncOut.H>
#include <AMReX_PlotFileUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#ifdef AMREX_USE_HDF5
#include "hdf5.h"

#ifdef AMREX_USE_HDF5_ASYNC
#include "h5_vol_external_async_native.h"
#endif

#endif

namespace amrex {

std::string LevelPath (int level, const std::string &levelPrefix)
{
    return Concatenate(levelPrefix, level, 1);  // e.g., Level_5
}

std::string MultiFabHeaderPath (int level,
                                const std::string &levelPrefix,
                                const std::string &mfPrefix)
{
    return LevelPath(level, levelPrefix) + '/' + mfPrefix;  // e.g., Level_4/Cell
}

std::string LevelFullPath (int level,
                           const std::string &plotfilename,
                           const std::string &levelPrefix)
{
    std::string r(plotfilename);
    if ( ! r.empty() && r.back() != '/') {
	r += '/';
    }
    r += LevelPath(level, levelPrefix);  // e.g., plt00005/Level_5
    return r;
}

std::string MultiFabFileFullPrefix (int level,
                                    const std::string& plotfilename,
                                    const std::string &levelPrefix,
                                    const std::string &mfPrefix)
{
    std::string r(plotfilename);
    if ( ! r.empty() && r.back() != '/') {
	r += '/';
    }
    r += MultiFabHeaderPath(level, levelPrefix, mfPrefix);
    return r;
}


void
PreBuildDirectorHierarchy (const std::string &dirName,
                           const std::string &/*subDirPrefix*/,
                           int nSubDirs, bool callBarrier)
{
  UtilCreateCleanDirectory(dirName, false);  // ---- dont call barrier
  for(int i(0); i < nSubDirs; ++i) {
    const std::string &fullpath = LevelFullPath(i, dirName);
    UtilCreateCleanDirectory(fullpath, false);  // ---- dont call barrier
  }

  if(callBarrier) {
    ParallelDescriptor::Barrier();
  }
}


void
WriteGenericPlotfileHeader (std::ostream &HeaderFile,
                            int nlevels,
                            const Vector<BoxArray> &bArray,
                            const Vector<std::string> &varnames,
                            const Vector<Geometry> &geom,
                            Real time,
                            const Vector<int> &level_steps,
                            const Vector<IntVect> &ref_ratio,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix)
{
//        BL_PROFILE("WriteGenericPlotfileHeader()");

        BL_ASSERT(nlevels <= bArray.size());
        BL_ASSERT(nlevels <= geom.size());
        BL_ASSERT(nlevels <= ref_ratio.size()+1);
        BL_ASSERT(nlevels <= level_steps.size());

        int finest_level(nlevels - 1);

	HeaderFile.precision(17);

	// ---- this is the generic plot file type name
        HeaderFile << versionName << '\n';

        HeaderFile << varnames.size() << '\n';

        for (int ivar = 0; ivar < varnames.size(); ++ivar) {
	    HeaderFile << varnames[ivar] << "\n";
        }
        HeaderFile << AMREX_SPACEDIM << '\n';
        HeaderFile << time << '\n';
        HeaderFile << finest_level << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbLo(i) << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbHi(i) << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i < finest_level; ++i) {
            HeaderFile << ref_ratio[i][0] << ' ';
	}
        HeaderFile << '\n';
	for (int i = 0; i <= finest_level; ++i) {
	    HeaderFile << geom[i].Domain() << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << level_steps[i] << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                HeaderFile << geom[i].CellSize()[k] << ' ';
	    }
            HeaderFile << '\n';
        }
        HeaderFile << (int) geom[0].Coord() << '\n';
        HeaderFile << "0\n";

	for (int level = 0; level <= finest_level; ++level) {
	    HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
	    HeaderFile << level_steps[level] << '\n';

            const IntVect& domain_lo = geom[level].Domain().smallEnd();
            for (int i = 0; i < bArray[level].size(); ++i)
            {
                // Need to shift because the RealBox ctor we call takes the
                // physical location of index (0,0,0).  This does not affect
                // the usual cases where the domain index starts with 0.
                const Box& b = amrex::shift(bArray[level][i], -domain_lo);
                RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
                }
            }

	    HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
	}
}


void
WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
                         const Vector<const MultiFab*>& mf,
                         const Vector<std::string>& varnames,
                         const Vector<Geometry>& geom, Real time,
                         const Vector<int>& level_steps,
                         const Vector<IntVect>& ref_ratio,
                         const std::string &versionName,
                         const std::string &levelPrefix,
                         const std::string &mfPrefix,
                         const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int finest_level = nlevels-1;

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::MyProc() == ParallelDescriptor::NProcs()-1) {
        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        auto f = [=]() {
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
            std::string HeaderFileName(plotfilename + "/Header");
            std::ofstream HeaderFile;
            HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) FileOpenFailed(HeaderFileName);
            WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                       geom, time, level_steps, ref_ratio, versionName,
                                       levelPrefix, mfPrefix);
        };

        if (AsyncOut::UseAsyncOut()) {
            AsyncOut::Submit(std::move(f));
        } else {
            f();
        }
    }

    for (int level = 0; level <= finest_level; ++level)
    {
        if (AsyncOut::UseAsyncOut()) {
            VisMF::AsyncWrite(*mf[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix),
                              true);
        } else {
            const MultiFab* data;
            std::unique_ptr<MultiFab> mf_tmp;
            if (mf[level]->nGrowVect() != 0) {
                mf_tmp.reset(new MultiFab(mf[level]->boxArray(),
                                          mf[level]->DistributionMap(),
                                          mf[level]->nComp(), 0, MFInfo(),
                                          mf[level]->Factory()));
                MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
                data = mf_tmp.get();
            } else {
                data = mf[level];
            }
            VisMF::Write(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
        }
    }
}

// write a plotfile to disk given:
// -plotfile name
// -vector of MultiFabs
// -vector of Geometrys
// variable names are written as "Var0", "Var1", etc.    
// refinement ratio is computed from the Geometry vector
// "time" and "level_steps" are set to zero
void WriteMLMF (const std::string &plotfilename,
                const Vector<const MultiFab*>& mf,
                const Vector<Geometry> &geom)
{
    int nlevs = mf.size();
    int ncomp = mf[0]->nComp();

    // variables names are "Var0", "Var1", etc.
    Vector<std::string> varnames(ncomp);
    for (int i=0; i<ncomp; ++i) {
        varnames[i] = "Var" + std::to_string(i);
    }

    // compute refinement ratio by looking at hi coordinates of domain at each level from
    // the geometry object
    Vector<IntVect> ref_ratio(nlevs-1);
    for (int i = 0; i < nlevs-1; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            int rr = (geom[i+1].Domain()).bigEnd(d)/(geom[i].Domain()).bigEnd(d);
            ref_ratio[i][d] = rr;
        }
    }

    // set step_array to zero
    Vector<int> step_array(nlevs,0);

    // set time to zero
    Real time = 0.;
    
    WriteMultiLevelPlotfile(plotfilename, nlevs, mf, varnames,
                            geom, time, step_array, ref_ratio);   
    
}    


void
WriteMultiLevelPlotfileHeaders (const std::string & plotfilename, int nlevels,
                                const Vector<const MultiFab*> & mf,
                                const Vector<std::string>     & varnames,
                                const Vector<Geometry>        & geom,
                                Real time, const Vector<int>  & level_steps,
                                const Vector<IntVect>     & ref_ratio,
                                const std::string         & versionName,
                                const std::string         & levelPrefix,
                                const std::string         & mfPrefix,
                                const Vector<std::string> & extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    int finest_level = nlevels-1;

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::string HeaderFileName(plotfilename + "/Header");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                        std::ofstream::trunc |
                        std::ofstream::binary);
        if( ! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                   geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix);
    }

    for (int level = 0; level <= finest_level; ++level) {
        const MultiFab * data;
        data = mf[level];
        VisMF::WriteOnlyHeader(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

}

void
WriteSingleLevelPlotfile (const std::string& plotfilename,
                          const MultiFab& mf, const Vector<std::string>& varnames,
                          const Geometry& geom, Real time, int level_step,
                          const std::string &versionName,
                          const std::string &levelPrefix,
                          const std::string &mfPrefix,
                          const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    WriteMultiLevelPlotfile(plotfilename, 1, mfarr, varnames, geomarr, time,
                            level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}


#ifdef AMREX_USE_EB
void
EB_WriteSingleLevelPlotfile (const std::string& plotfilename,
                             const MultiFab& mf, const Vector<std::string>& varnames,
                             const Geometry& geom, Real time, int level_step,
                             const std::string &versionName,
                             const std::string &levelPrefix,
                             const std::string &mfPrefix,
                             const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    EB_WriteMultiLevelPlotfile(plotfilename, 1, mfarr, varnames, geomarr, time,
                               level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}

void
EB_WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
                            const Vector<const MultiFab*>& mf,
                            const Vector<std::string>& varnames,
                            const Vector<Geometry>& geom, Real time, const Vector<int>& level_steps,
                            const Vector<IntVect>& ref_ratio,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix,
                            const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf[0]->hasEBFabFactory(),
                                     "EB_WriteMultiLevelPlotfile: does not have EB Factory");

    int finest_level = nlevels-1;

//    int saveNFiles(VisMF::GetNOutFiles());
//    VisMF::SetNOutFiles(std::max(1024,saveNFiles));

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::string HeaderFileName(plotfilename + "/Header");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
	                                        std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        Vector<std::string> vn = varnames;
        vn.push_back("vfrac");
        WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, vn,
                                   geom, time, level_steps, ref_ratio, versionName,
                                   levelPrefix, mfPrefix);

        for (int lev = 0; lev < nlevels; ++lev) {
            HeaderFile << "1.0e-6\n";
        }
    }


    for (int level = 0; level <= finest_level; ++level)
    {
        const int nc = mf[level]->nComp();
        MultiFab mf_tmp(mf[level]->boxArray(),
                        mf[level]->DistributionMap(),
                        nc+1, 0);
        MultiFab::Copy(mf_tmp, *mf[level], 0, 0, nc, 0);
        auto const& factory = dynamic_cast<EBFArrayBoxFactory const&>(mf[level]->Factory());
        MultiFab::Copy(mf_tmp, factory.getVolFrac(), 0, nc, 1, 0);
	VisMF::Write(mf_tmp, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

//    VisMF::SetNOutFiles(saveNFiles);
}

#endif

#ifdef AMREX_USE_HDF5
static int CreateWriteHDF5AttrDouble(hid_t loc, const char *name, hsize_t n, const double *data)
{
    herr_t ret;
    hid_t attr, attr_space;
    hsize_t dims = n;

    attr_space = H5Screate_simple(1, &dims, NULL);

    attr = H5Acreate(loc, name, H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        printf("%s: Error with H5Acreate [%s]\n", __func__, name);
        return -1;
    }

    ret  = H5Awrite(attr, H5T_NATIVE_DOUBLE, (void*)data);
    if (ret < 0) {
        printf("%s: Error with H5Awrite [%s]\n", __func__, name);
        return -1;
    }
    H5Sclose(attr_space);
    H5Aclose(attr);
    return 1;
}

static int CreateWriteHDF5AttrInt(hid_t loc, const char *name, hsize_t n, const int *data)
{
    herr_t ret;
    hid_t attr, attr_space;
    hsize_t dims = n;

    attr_space = H5Screate_simple(1, &dims, NULL);

    attr = H5Acreate(loc, name, H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        printf("%s: Error with H5Acreate [%s]\n", __func__, name);
        return -1;
    }

    ret  = H5Awrite(attr, H5T_NATIVE_INT, (void*)data);
    if (ret < 0) {
        printf("%s: Error with H5Awrite [%s]\n", __func__, name);
        return -1;
    }
    H5Sclose(attr_space);
    H5Aclose(attr);
    return 1;
}

static int CreateWriteHDF5AttrString(hid_t loc, const char *name, const char* str)
{
    hid_t attr, atype, space;
    herr_t ret;

    BL_ASSERT(name);
    BL_ASSERT(str);

    space = H5Screate(H5S_SCALAR);
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(str)+1);
    H5Tset_strpad(atype,H5T_STR_NULLTERM);
    attr = H5Acreate(loc, name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        printf("%s: Error with H5Acreate [%s]\n", __func__, name);
        return -1;
    }

    ret = H5Awrite(attr, atype, str);
    if (ret < 0) {
        printf("%s: Error with H5Awrite[%s]\n", __func__, name);
        return -1;
    }

    H5Tclose(atype);
    H5Sclose(space);
    H5Aclose(attr);

    return 1;
}

static int CreateWriteDsetDouble(hid_t loc, const char *name, hsize_t n, const double *data)
{
    herr_t ret;
    hid_t dset, dset_space;
    hsize_t dims = n;

    dset_space = H5Screate_simple(1, &dims, NULL);

    dset = H5Dcreate(loc, name, H5T_NATIVE_DOUBLE, dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        printf("%s: Error with H5Dcreate [%s]\n", __func__, name);
        return -1;
    }

    ret  = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)data);
    if (ret < 0) {
        printf("%s: Error with H5Dwrite [%s]\n", __func__, name);
        return -1;
    }
    H5Sclose(dset_space);
    H5Aclose(dset);
    return 1;
}


void
WriteGenericPlotfileHeaderHDF5 (hid_t fid,
                            int nlevels,
                            const Vector<const MultiFab*>& mf,
                            const Vector<BoxArray> &bArray,
                            const Vector<std::string> &varnames,
                            const Vector<Geometry> &geom,
                            Real time,
                            const Vector<int> &level_steps,
                            const Vector<IntVect> &ref_ratio,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix, 
                            const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteGenericPlotfileHeaderHDF5()");

    BL_ASSERT(nlevels <= bArray.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());

    int finest_level(nlevels - 1);

    CreateWriteHDF5AttrString(fid, "version_name", versionName.c_str());
    CreateWriteHDF5AttrString(fid, "plotfile_type", "VanillaHDF5");

    int ncomp = varnames.size();
    CreateWriteHDF5AttrInt(fid, "num_components", 1, &ncomp);

    char comp_name[32];
    for (int ivar = 0; ivar < varnames.size(); ++ivar) {
        sprintf(comp_name, "component_%d", ivar);
        CreateWriteHDF5AttrString(fid, comp_name, varnames[ivar].c_str());
    }

    int ndim = AMREX_SPACEDIM;
    CreateWriteHDF5AttrInt(fid, "dim", 1, &ndim);
    double cur_time = (double)time;
    CreateWriteHDF5AttrDouble(fid, "time", 1, &cur_time);
    CreateWriteHDF5AttrInt(fid, "finest_level", 1, &finest_level);


    int coord = (int) geom[0].Coord();
    CreateWriteHDF5AttrInt(fid, "coordinate_system", 1, &coord);

    hid_t grp;
    char level_name[128];
    double lo[AMREX_SPACEDIM], hi[AMREX_SPACEDIM], cellsizes[AMREX_SPACEDIM];

    // For VisIt Chombo plot
    CreateWriteHDF5AttrInt(fid, "num_levels", 1, &nlevels);
    grp = H5Gcreate(fid, "Chombo_global", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    CreateWriteHDF5AttrInt(grp, "SpaceDim", 1, &ndim);
    H5Gclose(grp);

    hid_t comp_dtype;

    comp_dtype = H5Tcreate (H5T_COMPOUND, 2 * AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_j", 3 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_k", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 3 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_j", 4 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_k", 5 * sizeof(int), H5T_NATIVE_INT);
    }

    for (int level = 0; level <= finest_level; ++level) {
        sprintf(level_name, "level_%d", level);
        /* sprintf(level_name, "%s%d", levelPrefix.c_str(), level); */
        grp = H5Gcreate(fid, level_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (grp < 0) {
            std::cout << "H5Gcreate [" << level_name << "] failed!" << std::endl;
            continue;
        }

        int ratio = 1;
        if (ref_ratio.size() > 0) 
            ratio = ref_ratio[level][0];
        
        if (level == finest_level) {
            ratio = 1;
        }
        CreateWriteHDF5AttrInt(grp, "ref_ratio", 1, &ratio);

        for (int k = 0; k < AMREX_SPACEDIM; ++k) {
            cellsizes[k] = (double)geom[level].CellSize()[k];
        }
        // Visit has issues with vec_dx, and is ok with a single "dx" value
        CreateWriteHDF5AttrDouble(grp, "Vec_dx", AMREX_SPACEDIM, cellsizes);
        // For VisIt Chombo plot
        CreateWriteHDF5AttrDouble(grp, "dx", 1, &cellsizes[0]);

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            lo[i] = (double)geom[level].ProbLo(i);
            hi[i] = (double)geom[level].ProbHi(i);
        }
        CreateWriteHDF5AttrDouble(grp, "prob_lo", AMREX_SPACEDIM, lo);
        CreateWriteHDF5AttrDouble(grp, "prob_hi", AMREX_SPACEDIM, hi);

        int domain[AMREX_SPACEDIM*2];
        Box tmp(geom[level].Domain());
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            domain[i] = tmp.smallEnd(i);
            domain[i+AMREX_SPACEDIM] = tmp.bigEnd(i);
        }

        hid_t aid = H5Screate(H5S_SCALAR);
        hid_t domain_attr = H5Acreate(grp, "prob_domain", comp_dtype, aid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(domain_attr, comp_dtype, domain);
        H5Aclose(domain_attr);
        H5Sclose(aid);

        int type[AMREX_SPACEDIM];
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            type[i] = (int)geom[level].Domain().ixType().test(i) ? 1 : 0;
        }
        CreateWriteHDF5AttrInt(grp, "domain_type", AMREX_SPACEDIM, type);

        CreateWriteHDF5AttrInt(grp, "steps", 1, &level_steps[level]);

        int ngrid = bArray[level].size();
        CreateWriteHDF5AttrInt(grp, "ngrid", 1, &ngrid);
        double cur_time = (double)time;
        CreateWriteHDF5AttrDouble(grp, "time", 1, &cur_time);

        int ngrow = mf[level]->nGrow();
        CreateWriteHDF5AttrInt(grp, "ngrow", 1, &ngrow);

        /* hsize_t npts = ngrid*AMREX_SPACEDIM*2; */
        /* double *realboxes = new double [npts]; */
        /* for (int i = 0; i < bArray[level].size(); ++i) */
        /* { */
        /*     const Box &b(bArray[level][i]); */
        /*     RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo()); */
        /*     for (int n = 0; n < AMREX_SPACEDIM; ++n) { */
        /*         /1* HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n'; *1/ */
        /*         realboxes[i*AMREX_SPACEDIM*2 + n] = loc.lo(n); */
        /*         realboxes[i*AMREX_SPACEDIM*2 + AMREX_SPACEDIM + n] = loc.hi(n); */
        /*     } */
        /* } */
        /* CreateWriteDsetDouble(grp, "Boxes", npts, realboxes); */
        /* delete [] realboxes; */

        H5Gclose(grp);
    }

    H5Tclose(comp_dtype);
}

void WriteMultiLevelPlotfileHDF5 (const std::string& plotfilename, 
				  int nlevels,
                         	  const Vector<const MultiFab*>& mf,
                         	  const Vector<std::string>& varnames,
                         	  const Vector<Geometry>& geom, 
				  Real time, 
				  const Vector<int>& level_steps,
                         	  const Vector<IntVect>& ref_ratio,
                         	  const std::string &versionName,
                         	  const std::string &levelPrefix,
                         	  const std::string &mfPrefix,
                         	  const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfileHDF5");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int myProc(ParallelDescriptor::MyProc());
    int nProcs(ParallelDescriptor::NProcs());
    

    herr_t  ret;
    int ndim = AMREX_SPACEDIM;
    int finest_level = nlevels-1;
    int ncomp = mf[0]->nComp();
    /* double total_write_start_time(ParallelDescriptor::second()); */
    std::string filename(plotfilename + ".h5");

    // Write out root level metadata
    hid_t fapl, dxpl, fid, grp;

    if(ParallelDescriptor::IOProcessor()) {
        // Have only one rank to create and write metadata (header)
        fapl = H5Pcreate (H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(fapl, MPI_COMM_SELF, MPI_INFO_NULL);
        H5Pset_coll_metadata_write(fapl, true);
        H5Pset_all_coll_metadata_ops(fapl, true);

        // Create the HDF5 file
        fid = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        if (fid < 0) 
            FileOpenFailed(filename.c_str());

        H5Pclose(fapl);

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        WriteGenericPlotfileHeaderHDF5(fid, nlevels, mf, boxArrays, varnames, geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
        H5Fclose(fid);
    }

    ParallelDescriptor::Barrier();

    hid_t babox_id;
    babox_id = H5Tcreate (H5T_COMPOUND, 2 * AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
	H5Tinsert (babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (babox_id, "hi_i", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
	H5Tinsert (babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (babox_id, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (babox_id, "hi_i", 2 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (babox_id, "hi_j", 3 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
	H5Tinsert (babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (babox_id, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (babox_id, "lo_k", 2 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (babox_id, "hi_i", 3 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (babox_id, "hi_j", 4 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (babox_id, "hi_k", 5 * sizeof(int), H5T_NATIVE_INT);
    }
 
    hid_t center_id = H5Tcreate (H5T_COMPOUND, AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
	H5Tinsert (center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
	H5Tinsert (center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (center_id, "j", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
	H5Tinsert (center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (center_id, "j", 1 * sizeof(int), H5T_NATIVE_INT);
	H5Tinsert (center_id, "k", 2 * sizeof(int), H5T_NATIVE_INT);
    }
 
    fapl = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl,  ParallelDescriptor::Communicator(), MPI_INFO_NULL);
    int alignment = 16 * 1024 * 1024;
    H5Pset_alignment(fapl, alignment, alignment);
    H5Pset_coll_metadata_write(fapl, true);
    H5Pset_all_coll_metadata_ops(fapl, true);

    dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);

    // Only use async for writing actual data
    #ifdef AMREX_USE_HDF5_ASYNC
    H5Pset_vol_async(fapl);
    H5Pset_dxpl_async(dxpl, true);
    #endif

    // All process open the file
    fid = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
    if (fid < 0) 
        FileOpenFailed(filename.c_str());

    // Write data for each level
    char level_name[32];
    for (int level = 0; level <= finest_level; ++level) {
        sprintf(level_name, "level_%d", level);
        grp = H5Gopen(fid, level_name, H5P_DEFAULT);
        if (grp < 0) {
            std::cout << "H5Gopen [" << level_name << "] failed!" << std::endl;
            continue;
        }

        /* const MultiFab* data; */
        /* std::unique_ptr<MultiFab> mf_tmp; */
        /* if (mf[level]->nGrow() > 0) { */
        /*     mf_tmp.reset(new MultiFab(mf[level]->boxArray(), */
        /*                               mf[level]->DistributionMap(), */
        /*                               mf[level]->nComp(), 0, MFInfo(), */
        /*                               mf[level]->Factory())); */
        /*     MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0); */
        /*     data = mf_tmp.get(); */
        /* } else { */
        /*     data = mf[level]; */
        /* } */
	/* VisMF::Write(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix)); */

        // Get the boxes assigned to all ranks and calculate their offsets and sizes
        Vector<int> procMap = mf[level]->DistributionMap().ProcessorMap();
        const BoxArray& grids = mf[level]->boxArray();
        hid_t boxdataset, boxdataspace;
        hid_t offsetdataset, offsetdataspace;
        hid_t centerdataset, centerdataspace;
        std::string bdsname("boxes");
        std::string odsname("data:offsets=0");
        std::string centername("boxcenter");
        std::string dataname("data:datatype=0");
        hsize_t  flatdims[1], count[1];
        flatdims[0] = grids.size();
        
        flatdims[0] = grids.size();
        boxdataspace = H5Screate_simple(1, flatdims, NULL);
       
   
        boxdataset = H5Dcreate(grp, bdsname.c_str(), babox_id, boxdataspace, H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
        
        // Create a boxarray sorted by rank
        std::map<int, Vector<Box> > gridMap;
        for(int i(0); i < grids.size(); ++i) {
            int gridProc(procMap[i]);
            Vector<Box> &boxesAtProc = gridMap[gridProc];
            boxesAtProc.push_back(grids[i]);
        }
        BoxArray sortedGrids(grids.size());
        Vector<int> sortedProcs(grids.size());
        int bIndex(0);
        for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
            int proc = it->first;
            Vector<Box> &boxesAtProc = it->second;
            for(int ii(0); ii < boxesAtProc.size(); ++ii) {
                sortedGrids.set(bIndex, boxesAtProc[ii]);
                sortedProcs[bIndex] = proc;
                ++bIndex;
            }
        }
        
        hsize_t  oflatdims[1];
        oflatdims[0] = sortedGrids.size() + 1;
        offsetdataspace = H5Screate_simple(1, oflatdims, NULL);
        offsetdataset   = H5Dcreate(grp, odsname.c_str(), H5T_NATIVE_LLONG, offsetdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        hsize_t centerdims[1];
        centerdims[0]   = sortedGrids.size() ;
        centerdataspace = H5Screate_simple(1, centerdims, NULL);
        centerdataset   = H5Dcreate(grp, centername.c_str(), center_id, centerdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        Vector<unsigned long long> offsets(sortedGrids.size() + 1);
        unsigned long long currentOffset(0L);
        for(int b(0); b < sortedGrids.size(); ++b) {
            offsets[b] = currentOffset;
            currentOffset += sortedGrids[b].numPts() * ncomp;
        }
        offsets[sortedGrids.size()] = currentOffset;
        
        Vector<unsigned long long> procOffsets(nProcs);
        int posCount(0);
        Vector<unsigned long long> procBufferSize(nProcs);
        unsigned long long totalOffset(0);
        for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
            int proc = it->first;
            Vector<Box> &boxesAtProc = it->second;
            BL_ASSERT(posCount == proc);
            procOffsets[posCount] = totalOffset;
            ++posCount;
            procBufferSize[proc] = 0L;
            for(int b(0); b < boxesAtProc.size(); ++b) {
                procBufferSize[proc] += boxesAtProc[b].numPts() * ncomp;
            }
            totalOffset += procBufferSize[proc];
        }
        
        if(ParallelDescriptor::IOProcessor()) {
            int vbCount(0);
            Vector<int> vbox(sortedGrids.size() * 2 * AMREX_SPACEDIM);
            Vector<int> centering(sortedGrids.size() * AMREX_SPACEDIM);
            count[0] = sortedGrids.size();
            for(int b(0); b < sortedGrids.size(); ++b) {
                for(int i(0); i < AMREX_SPACEDIM; ++i) {
                    vbox[(vbCount * 2 * AMREX_SPACEDIM) + i] = sortedGrids[b].smallEnd(i);
                    vbox[(vbCount * 2 * AMREX_SPACEDIM) + i + AMREX_SPACEDIM] = sortedGrids[b].bigEnd(i);
                    centering[vbCount * AMREX_SPACEDIM + i] = sortedGrids[b].ixType().test(i) ? 1 : 0;
                }
                ++vbCount;
            }
           
            // Only proc zero needs to write out this information
            ret = H5Dwrite(offsetdataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, dxpl, &(offsets[0]));
            if(ret < 0) { std::cout << "Write offset dataset failed! ret = " << ret << std::endl; }

            ret = H5Dwrite(centerdataset, center_id, H5S_ALL, H5S_ALL, dxpl, &(centering[0]));
            if(ret < 0) { std::cout << "Write center dataset failed! ret = " << ret << std::endl; }

            ret = H5Dwrite(boxdataset, babox_id, H5S_ALL, H5S_ALL, dxpl, &(vbox[0]));
            if(ret < 0) { std::cout << "Write box dataset failed! ret = " << ret << std::endl; }
        }
       
        
        BL_PROFILE_VAR("H5Dwritedata", h5dwd);
        hsize_t hs_procsize[1], hs_allprocsize[1], ch_offset[1];
        
        ch_offset[0]       = procOffsets[myProc];          // ---- offset on this proc
        hs_procsize[0]     = procBufferSize[myProc];       // ---- size of buffer on this proc
        hs_allprocsize[0]  = offsets[sortedGrids.size()];  // ---- size of buffer on all procs
        
        hid_t dataspace    = H5Screate_simple(1, hs_allprocsize, NULL);
        hid_t memdataspace = H5Screate_simple(1, hs_procsize, NULL);
        hid_t dataset      = H5Dcreate(grp, dataname.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, ch_offset, NULL, hs_procsize, NULL);
        
        Vector<Real> a_buffer(procBufferSize[myProc], -1.0);
        long dataCount(0);
        for(MFIter mfi(*mf[level]); mfi.isValid(); ++mfi) {
            const Box &vbox    = mfi.validbox();
            const Real *dataPtr = (*mf[level])[mfi].dataPtr();
            for(int i(0); i < vbox.numPts() * ncomp; ++i) {
                a_buffer[dataCount++] = dataPtr[i];
            }
        }
        
        /* ret = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE); */
           
        BL_PROFILE_VAR("H5DwriteGrids", h5dwg);

        ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace, dxpl, a_buffer.dataPtr());
        if(ret < 0) 
            std::cout << ParallelDescriptor::MyProc() << "Write data failed!  ret = " << ret << std::endl;

        BL_PROFILE_VAR_STOP(h5dwg);

        
        H5Sclose(memdataspace);
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Sclose(offsetdataspace);
        H5Dclose(offsetdataset);
        H5Sclose(centerdataspace);
        H5Dclose(centerdataset);
        H5Sclose(boxdataspace);
        H5Dclose(boxdataset);
        
        BL_PROFILE_VAR_STOP(h5dwd);

        H5Gclose(grp);
    } // For group

    H5Tclose(center_id);
    H5Tclose(babox_id);
    H5Pclose(fapl);
    H5Pclose(dxpl);
    H5Fclose(fid);

    /* double total_write_end_time(ParallelDescriptor::second()); */
    /* double total_write_time(total_write_end_time - total_write_start_time); */
    /* ParallelDescriptor::ReduceRealMax(total_write_time); */

    /* if(ParallelDescriptor::IOProcessor()) { */
    /*     std::cout << "WriteMultiLevelPlotfileHDF5 Time = " << total_write_time << "  seconds." << std::endl; */
    /* } */    

}

void
WriteSingleLevelPlotfileHDF5 (const std::string& plotfilename,
                          const MultiFab& mf, const Vector<std::string>& varnames,
                          const Geometry& geom, Real time, int level_step,
                          const std::string &versionName,
                          const std::string &levelPrefix,
                          const std::string &mfPrefix,
                          const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    WriteMultiLevelPlotfileHDF5(plotfilename, 1, mfarr, varnames, geomarr, time,
                            level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}

#endif
}
