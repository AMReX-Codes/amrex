
#include <fstream>
#include <iomanip>

#include <VisMF.H>
#include <PlotFileUtil.H>

std::string BoxLib::LevelPath (int level, const std::string &levelPrefix)
{
    return BoxLib::Concatenate(levelPrefix, level, 1);  // e.g., Level_5
}

std::string BoxLib::MultiFabHeaderPath (int level,
					const std::string &levelPrefix,
                                        const std::string &mfPrefix)
{
    return BoxLib::LevelPath(level, levelPrefix) + '/' + mfPrefix;  // e.g., Level_4/Cell
}

std::string BoxLib::LevelFullPath (int level,
                                   const std::string &plotfilename,
                                   const std::string &levelPrefix)
{
    std::string r(plotfilename);
    if ( ! r.empty() && r.back() != '/') {
	r += '/';
    }
    r += BoxLib::LevelPath(level, levelPrefix);  // e.g., plt00005/Level_5
    return r;
}

std::string BoxLib::MultiFabFileFullPrefix (int level,
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
BoxLib::PreBuildDirectorHierarchy (const std::string &dirName,
                                   const std::string &subDirPrefix,
                                   int nSubDirs, bool callBarrier)
{
  BoxLib::UtilCreateCleanDirectory(dirName, false);  // ---- dont call barrier
  for(int i(0); i < nSubDirs; ++i) {
    const std::string &fullpath = LevelFullPath(i, dirName);
    BoxLib::UtilCreateCleanDirectory(fullpath, false);  // ---- dont call barrier
  }

  if(callBarrier) {
    ParallelDescriptor::Barrier();
  }
}


void
BoxLib::WriteGenericPlotfileHeader (std::ostream &HeaderFile,
                                    int nlevels,
				    const Array<BoxArray> &bArray,
				    const Array<std::string> &varnames,
				    const Array<Geometry> &geom,
				    Real time,
				    const Array<int> &level_steps,
				    const Array<IntVect> &ref_ratio,
				    const std::string &versionName,
				    const std::string &levelPrefix,
				    const std::string &mfPrefix)
{
        BL_PROFILE("WriteGenericPlotfileHeader()");

        BL_ASSERT(nlevels <= bArray.size());
        BL_ASSERT(nlevels <= geom.size());
        BL_ASSERT(nlevels <= ref_ratio.size()+1);
        BL_ASSERT(nlevels <= level_steps.size());

        int finest_level(nlevels - 1);

	HeaderFile.precision(17);

	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
	HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

	// ---- this is the generic plot file type name
        HeaderFile << versionName << '\n';

        HeaderFile << varnames.size() << '\n';

        for (int ivar = 0; ivar < varnames.size(); ++ivar) {
	    HeaderFile << varnames[ivar] << "\n";
        }
        HeaderFile << BL_SPACEDIM << '\n';
        HeaderFile << time << '\n';
        HeaderFile << finest_level << '\n';
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            HeaderFile << Geometry::ProbLo(i) << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            HeaderFile << Geometry::ProbHi(i) << ' ';
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
            for (int k = 0; k < BL_SPACEDIM; ++k) {
                HeaderFile << geom[i].CellSize()[k] << ' ';
	    }
            HeaderFile << '\n';
        }
        HeaderFile << (int) Geometry::Coord() << '\n';
        HeaderFile << "0\n";

	for (int level = 0; level <= finest_level; ++level) {
	    HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
	    HeaderFile << level_steps[level] << '\n';
	    
	    for (int i = 0; i < bArray[level].size(); ++i)
	    {
		const Box &b(bArray[level][i]);
		RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
		for (int n = 0; n < BL_SPACEDIM; ++n) {
		    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
		}
	    }

	    HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
	}
}


void
BoxLib::WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
				 const Array<const MultiFab*>& mf,
				 const Array<std::string>& varnames,
				 const Array<Geometry>& geom, Real time, const Array<int>& level_steps,
				 const Array<IntVect>& ref_ratio)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int finest_level = nlevels-1;

    //
    // Only let 64 CPUs be writing at any one time.
    //
    int saveNFiles(VisMF::GetNOutFiles());
    VisMF::SetNOutFiles(64);

    const std::string versionName("HyperCLaw-V1.1");
    const std::string levelPrefix("Level_");
    const std::string mfPrefix("Cell");

    bool callBarrier(true);
    BoxLib::PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);

    if (ParallelDescriptor::IOProcessor()) {
      std::string HeaderFileName(plotfilename + "/Header");
      std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
	                                               std::ofstream::trunc |
						       std::ofstream::binary);
      if( ! HeaderFile.good()) {
        BoxLib::FileOpenFailed(HeaderFileName);
      }

      Array<BoxArray> boxArrays(nlevels);
      for(int level(0); level < boxArrays.size(); ++level) {
	boxArrays[level] = mf[level]->boxArray();
      }

      BoxLib::WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                         geom, time, level_steps, ref_ratio);
    }


    for (int level = 0; level <= finest_level; ++level)
    {
	VisMF::Write(*mf[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

    VisMF::SetNOutFiles(saveNFiles);
}

void
BoxLib::WriteSingleLevelPlotfile (const std::string& plotfilename,
				  const MultiFab& mf, const Array<std::string>& varnames,
				  const Geometry& geom, Real time, int level_step)
{
    Array<const MultiFab*> mfarr(1,&mf);
    Array<Geometry> geomarr(1,geom);
    Array<int> level_steps(1,level_step);
    Array<IntVect> ref_ratio;

    BoxLib::WriteMultiLevelPlotfile(plotfilename, 1, mfarr, varnames, geomarr, time,
				    level_steps, ref_ratio);
}


