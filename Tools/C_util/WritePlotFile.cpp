// Lifted code from the HCAll version of HyperCLaw for writing plotfiles.
// Have to fake a level 0 data set, and then use incoming multifab as
// level 1 data. All this could be cleaner, but this is throw-away code
// anyway, so why bother?

#include <iostream.h>
#include <fstream.h>


#ifdef BL_ARCH_CRAY
#ifdef BL_USE_DOUBLE
#error "DOUBLE PRECISION NOT ALLOWED ON CRAY"
#endif
#endif

#ifndef        WIN32
#include <unistd.h>
#endif

#include <Utility.H>
#include <Tracer.H>
#include <ParallelDescriptor.H>
#include <VisMF.H>
#include <RealBox.H>
#include <Geometry.H>

#ifdef BL_USE_NEW_HFILES
#include <new>
using std::setprecision;
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <new.h>
#endif

#include <WritePlotFile.H>

aString
thePlotFileType ()
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const aString the_plot_file_type("HyperCLaw-V1.1");

    return the_plot_file_type;
}

void
writePlotFile (const aString&  dir,
	       ostream&        os,
	       int             level,
	       const MultiFab& mf,
	       const Geometry& geom,
	       const IntVect&  refRatio,
	       Real            bgVal)
{
    //
    // Faked data
    //
    const int NUM_STATE = mf.nComp();
    const Real cumTime = 0.0;
    const Real curTime = cumTime;
    const int finestLevel = 1;
    Array< Box > domain(finestLevel+1);
    const IndexType& ixType = mf.boxArray()[0].ixType();
    if (ixType != IndexType::TheCellType())
	BoxLib::Error("writePlotfile unable to handle non cell-centered data for now");
    Box tmpb = Box(geom.Domain()).convert(ixType);
    Array<int> corr(BL_SPACEDIM);
    for (int d = 0; d < BL_SPACEDIM; d++)
    {
	corr[d] = (ixType.ixType(d) == IndexType::CELL ? 1 : 0);
    }
    for (int M = 0; M < finestLevel; M++)
    {
	tmpb.coarsen(refRatio);
    }
    domain[0] = tmpb;
    const int levelSteps = 0;
    Array< Array< Real > > dx(finestLevel+1);
    Array< Array<RealBox> > grid_loc(finestLevel+1);
    const BoxArray& grids = mf.boxArray();
    for (int j = 0; j<= finestLevel; j++)
    {
	if (j!=0)
	    domain[j] = Box( domain[j-1] ).refine(refRatio);
	dx[j].resize(BL_SPACEDIM);
	for (int k = 0; k < BL_SPACEDIM; k++)
	{
	    dx[j][k] = (Geometry::ProbHi(k) - Geometry::ProbLo(k))/domain[j].length(k);
	}
	if (j==0)
	{
	    grid_loc[j].resize(1);
	    grid_loc[j][0] = RealBox(Geometry::ProbLo(),Geometry::ProbHi());
	} else {
	    grid_loc[j].resize(grids.length());
	    for (int L=0; L < grids.length(); L++)
	    {
		const Box bx = grids[L];
		grid_loc[j][L] = RealBox(D_DECL( dx[j][0]*bx.smallEnd(0),
						 dx[j][1]*bx.smallEnd(1),
						 dx[j][2]*bx.smallEnd(2) ),
					 D_DECL( dx[j][0]*(bx.bigEnd(0)+corr[0]),
						 dx[j][1]*(bx.bigEnd(1)+corr[1]),
						 dx[j][2]*(bx.bigEnd(2)+corr[2]) ));
	    }
	}	    
    }
    const int Coord = 0;
    BoxArray tba = BoxArray(&tmpb,1);
    MultiFab level0_dat(tba,mf.nComp(),mf.nGrow(),Fab_allocate);
    level0_dat.setVal(bgVal);
    
    int i, n;
    //
    // There is only one MultiFab written out at each level in HyperCLaw.
    //
    static const aString MultiFabBaseName("/MultiFab");
    //
    // Build the directory to hold the MultiFabs at this level.
    // The name is relative to the directory containing the Header file.
    //
    char buf[64];
    sprintf(buf, "Level_%d", level);
    aString Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    aString FullPath = dir;
    if (!FullPath.isNull() && FullPath[FullPath.length()-1] != '/')
	FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
	if (!Utility::UtilCreateDirectory(FullPath, 0755))
	    Utility::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
	if (level == 0)
	{
	    //
	    // The first thing we write out is the plot file type.
	    //
	    os << thePlotFileType() << '\n';
	    // Number of components, and names
	    int n_var = NUM_STATE;
	    os << n_var << '\n';
	    for (n = 0; n < NUM_STATE; n++)
	    {
		char buff[64];
		sprintf(buff, "state_%d", n);
		os << buff << '\n';
	    }
	    // dimensionality
	    os << BL_SPACEDIM << '\n';
	    // time
	    os << cumTime << '\n';
	    // finest amr level
	    os << finestLevel << '\n';
	    // prob domain
	    for (i = 0; i < BL_SPACEDIM; i++) os << Geometry::ProbLo(i) << ' ';
	    os << '\n';
	    for (i = 0; i < BL_SPACEDIM; i++) os << Geometry::ProbHi(i) << ' ';
	    os << '\n';
	    // refinement ratio
	    for (i = 0; i < finestLevel; i++) os << refRatio << ' ';
	    os << '\n';
	    // int domain
	    for (i = 0; i <= finestLevel; i++) os << domain[i] << ' ';
	    os << '\n';
	    // level steps
	    for (i = 0; i <= finestLevel; i++) os << levelSteps << ' ';
	    os << '\n';
	    for (i = 0; i <= finestLevel; i++)
	    {
		// cell size
		const Real* dx_lev = dx[i].dataPtr();
		for (int k = 0; k < BL_SPACEDIM; k++) os << dx_lev[k] << ' ';
		os << '\n';
	    }
	    // coordinate system
	    os << Coord << '\n';
	    os << "0\n"; // The bndry data.
	}
	//
	// Now write state data.
	//
	int ngrds          = (level==0 ? 1 : grids.length());
	Real cur_time      = curTime;
	const MultiFab& cell_dat = (level==0 ? level0_dat : mf);
	    
	os << level << ' ' << ngrds << ' ' << cur_time << '\n';
	// level steps
	os << levelSteps << '\n';
	
	for (i = 0; i < cell_dat.boxArray().length(); ++i)
	{
	    for (n = 0; n < BL_SPACEDIM; n++)
		// lo/hi position of this grid
		os <<grid_loc[level][i].lo(n) << ' ' << grid_loc[level][i].hi(n) << '\n';
	}
	
	//
	// Finally, the full relative pathname of the MultiFab.
	// The name is relative to the Header file containing this name.
	// It's the name that gets written into the Header.
	//
	aString PathNameInHeader = Level;
	PathNameInHeader += MultiFabBaseName;

	os << PathNameInHeader << '\n';
    }
    //
    // Now amend FullPath to contain full pathname of MF.
    //
    FullPath += MultiFabBaseName;

    //RunStats::addBytes(VisMF::Write(state[State_Type].newData(),FullPath,how));
    const MultiFab& cell_dat = (level==0 ? level0_dat : mf);
    VisMF::Write(cell_dat,FullPath,VisMF::OneFilePerCPU);
};

void
writePlotFile (const char*     name,
	       const MultiFab& mf,
	       const Geometry& geom,
	       const IntVect&  refRatio,
	       Real            bgVal)
{

    double dPlotFileTime0(ParallelDescriptor::second());
    //aString pltfile = Concatenate(root, num);
    aString pltfile(name);

    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!Utility::UtilCreateDirectory(pltfile, 0755))
            Utility::CreateDirectoryFailed(pltfile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    aString HeaderFileName = pltfile + "/Header";

#ifdef BL_USE_SETBUF
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
#endif
    ofstream HeaderFile;

#ifdef BL_USE_SETBUF
    HeaderFile.rdbuf()->setbuf(io_buffer.dataPtr(), io_buffer.length());
#endif
    int old_prec;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), ios::out|ios::trunc);

        if (!HeaderFile.good())
            Utility::FileOpenFailed(HeaderFileName);

        old_prec = HeaderFile.precision(15);
    }

    int finest_level = 1;
    for (int k = 0; k <= finest_level; k++)
    {
        RunStats write_pltfile_stats("write_pltfile", k);
        write_pltfile_stats.start();
        writePlotFile(pltfile, HeaderFile, k, mf, geom, refRatio, bgVal);
        write_pltfile_stats.end();
    }

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Accumulate # of bytes written to header file.
        //
        RunStats::addBytes(VisMF::FileOffset(HeaderFile));

        HeaderFile.precision(old_prec);

        if (!HeaderFile.good())
            BoxLib::Error("Amr::writePlotFile() failed");
    }

    double dPlotFileTime1(ParallelDescriptor::second());
    double dPlotFileTime(dPlotFileTime1 - dPlotFileTime0);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime);
    if(ParallelDescriptor::IOProcessor()) {
      cout << "Write plotfile time = " << dPlotFileTime << "  seconds." << endl;
    }
    
}

