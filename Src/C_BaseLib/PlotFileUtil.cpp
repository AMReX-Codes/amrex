
#include <fstream>
#include <iomanip>

#include <VisMF.H>
#include <PlotFileUtil.H>

static std::string LevelPath (int level)
{
    return BoxLib::Concatenate("Level_",level,1); // e.g., Level_5
}

static std::string MultiFabFileName (int level)
{
    return LevelPath(level) + "/Cell"; // e.g., Level_4/Cell
}

static std::string LevelFullPath (int level, const std::string& plotfilename)
{
    std::string r = plotfilename;
    if (!r.empty() && r[r.length()-1] != '/') {// could simply use back() in C++11
	r += '/';
    }
    r += LevelPath(level); // e.g., plt00005/Level_5
    return r;
}

static std::string MultiFabFileFullName (int level, const std::string& plotfilename)
{
    std::string r = plotfilename;
    if (!r.empty() && r[r.length()-1] != '/') {// could simply use back() in C++11
	r += '/';
    }
    r += MultiFabFileName(level);
    return r;
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
    VisMF::SetNOutFiles(64);

    if (ParallelDescriptor::IOProcessor())
    {
	// Only the I/O processor makes the directory if it doesn't already exist.
        if (!BoxLib::UtilCreateDirectory(plotfilename, 0755))
            BoxLib::CreateDirectoryFailed(plotfilename);
    
	std::string HeaderFileName = plotfilename + "/Header";

	std::ofstream HeaderFile;

        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out|std::ofstream::trunc|std::ofstream::binary);
        if (!HeaderFile.good()) {
            BoxLib::FileOpenFailed(HeaderFileName);
	}

	HeaderFile.precision(17);

	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
	HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        HeaderFile << "HyperCLaw-V1.1\n";

        HeaderFile << mf[0]->nComp() << '\n';

	// variable names
        for (int ivar = 0; ivar < mf[0]->nComp(); ++ivar) {
	    HeaderFile << varnames[ivar] << "\n";
        }
	// dimensionality
        HeaderFile << BL_SPACEDIM << '\n';
	// time
        HeaderFile << time << '\n';
	// maximum level number
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
        for (int i = 0; i <= finest_level; ++i)
        {
            for (int k = 0; k < BL_SPACEDIM; ++k) {
                HeaderFile << geom[i].CellSize()[k] << ' ';
	    }
            HeaderFile << '\n';
        }
        HeaderFile << (int) Geometry::Coord() << '\n';
        HeaderFile << "0\n";

	for (int level = 0; level <= finest_level; ++level)
	{
	    const std::string& fullpath = LevelFullPath(level, plotfilename);
	    //
	    // Only the I/O processor makes the directory if it doesn't already exist.
	    //
	    if (!BoxLib::UtilCreateDirectory(fullpath, 0755)) {
		BoxLib::CreateDirectoryFailed(fullpath);
	    }

	    const BoxArray& ba = mf[level]->boxArray();

	    HeaderFile << level << ' ' << ba.size() << ' ' << time << '\n';
	    HeaderFile << level_steps[level] << '\n';
	    
	    for (int i = 0; i < ba.size(); ++i)
	    {
		RealBox loc = RealBox(ba[i],geom[level].CellSize(),geom[level].ProbLo());
		for (int n = 0; n < BL_SPACEDIM; n++) {
		    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
		}
	    }

	    HeaderFile << MultiFabFileName(level) << '\n';
	}
    }

    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    for (int level = 0; level <= finest_level; ++level)
    {
	VisMF::Write(*mf[level], MultiFabFileFullName(level,plotfilename));
    }
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
