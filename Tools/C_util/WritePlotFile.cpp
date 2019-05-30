// Lifted code from the HCAll version of HyperCLaw for writing plotfiles.
// Have to fake a level 0 data set, and then use incoming multifab as
// level 1 data. All this could be cleaner, but this is throw-away code
// anyway, so why bother?

#include <iostream>
#include <fstream>

#include <unistd.h>

#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <WritePlotFile.H>

using namespace amrex;
std::string
thePlotFileType ()
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");

    return the_plot_file_type;
}

void
writePlotFile (const std::string&        dir,
	       std::ostream&             os,
	       int                       level,
	       const MultiFab&           mf,
	       const Geometry&           geom,
	       const IntVect&            refRatio,
	       Real                      bgVal,
               const Vector<std::string>& names)
{
    //
    // Faked data
    //
    const int NUM_STATE = mf.nComp();
    const Real cumTime = 0.0;
    const Real curTime = cumTime;
    const int finestLevel = 1;
    Vector< Box > domain(finestLevel+1);
    const IndexType& ixType = mf.boxArray().ixType();
    if (ixType != IndexType::TheCellType())
	amrex::Error("writePlotfile unable to handle non cell-centered data for now");
    Box tmpb = Box(geom.Domain()).convert(ixType);
    Vector<int> corr(BL_SPACEDIM);
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
    Vector< Vector< Real > > dx(finestLevel+1);
    Vector< Vector<RealBox> > grid_loc(finestLevel+1);
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
	    grid_loc[j].resize(grids.size());
	    for (int L=0; L < grids.size(); L++)
	    {
		const Box& bx = grids[L];
		grid_loc[j][L] = RealBox(AMREX_D_DECL( dx[j][0]*bx.smallEnd(0),
						 dx[j][1]*bx.smallEnd(1),
						 dx[j][2]*bx.smallEnd(2) ),
					 AMREX_D_DECL( dx[j][0]*(bx.bigEnd(0)+corr[0]),
						 dx[j][1]*(bx.bigEnd(1)+corr[1]),
						 dx[j][2]*(bx.bigEnd(2)+corr[2]) ));
	    }
	}	    
    }
    const int Coord = 0;
    BoxArray tba = BoxArray(&tmpb,1);
    DistributionMapping dm {tba};
    MultiFab level0_dat(tba,dm,mf.nComp(),mf.nGrow());
    for (int j=0; j<mf.nComp(); ++j)
        level0_dat.setVal(0.5*(mf.min(j)+mf.max(j)),j,1);
    //level0_dat.setVal(bgVal);
    
    int i, n;
    //
    // There is only one MultiFab written out at each level in HyperCLaw.
    //
    static const std::string MultiFabBaseName("/MultiFab");
    //
    // Build the directory to hold the MultiFabs at this level.
    // The name is relative to the directory containing the Header file.
    //
    std::string Level = amrex::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
	FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
	if (!amrex::UtilCreateDirectory(FullPath, 0755))
	    amrex::CreateDirectoryFailed(FullPath);
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
                if (names.size()==0)
                {
                    std::string name = amrex::Concatenate("state_", n, 1);

                    os << name << '\n';
                }
                else
                {
                    os << names[n] << '\n';
                }
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
	int ngrds          = (level==0 ? 1 : grids.size());
	Real cur_time      = curTime;
	const MultiFab& cell_dat = (level==0 ? level0_dat : mf);
	    
	os << level << ' ' << ngrds << ' ' << cur_time << '\n';
	// level steps
	os << levelSteps << '\n';
	
	for (i = 0; i < cell_dat.boxArray().size(); ++i)
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
	std::string PathNameInHeader = Level;
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
}

void
writePlotFile (const char*               name,
	       const MultiFab&           mf,
	       const Geometry&           geom,
	       const IntVect&            refRatio,
	       Real                      bgVal,
               const Vector<std::string>& names)
{

    Real dPlotFileTime0(ParallelDescriptor::second());
    //std::string pltfile = Concatenate(root, num);
    std::string pltfile(name);

    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(pltfile, 0755))
            amrex::CreateDirectoryFailed(pltfile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    std::string HeaderFileName = pltfile + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc);

        if (!HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);

        old_prec = HeaderFile.precision(30);
    }

    int finest_level = 1;
    for (int k = 0; k <= finest_level; k++)
    {
        writePlotFile(pltfile, HeaderFile, k, mf, geom, refRatio, bgVal, names);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.precision(old_prec);

        if (!HeaderFile.good())
            amrex::Error("Amr::writePlotFile() failed");
    }

    Real dPlotFileTime1(ParallelDescriptor::second());
    Real dPlotFileTime(dPlotFileTime1 - dPlotFileTime0);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "Write plotfile time = " << dPlotFileTime << "  seconds." << std::endl;
    }
    
}

void WritePlotFile(const Vector<MultiFab*>&   mfa,
                   const Vector<Box>&         probDomain,
		   AmrData&                   amrdToMimic,
		   const std::string&         oFile,
		   bool                       verbose,
                   const Vector<std::string>& varNames)
{
    // If varnames not provided, use names in original plotfile
    const Vector<std::string>& derives = (varNames.size()==0 ? amrdToMimic.PlotVarNames() : varNames);
    AMREX_ASSERT(derives.size()==(*mfa[0]).nComp());
    int ntype = derives.size();
    int finestLevel = mfa.size() - 1;
    AMREX_ALWAYS_ASSERT(finestLevel >= 0);
    
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(oFile,0755))
            amrex::CreateDirectoryFailed(oFile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
  
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream os;
  
    os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
  
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Opening file = " << oFileHeader << '\n';

    os.open(oFileHeader.c_str(), std::ios::out|std::ios::binary);

    if (os.fail())
        amrex::FileOpenFailed(oFileHeader);
    //
    // Start writing plotfile.
    //
    os << amrdToMimic.PlotFileVersion() << '\n';
    int n_var = ntype;
    os << n_var << '\n';
    for (int n = 0; n < ntype; n++) os << derives[n] << '\n';
    os << BL_SPACEDIM << '\n';
    os << amrdToMimic.Time() << '\n';
    os << finestLevel << '\n';
    int i;
    for (i = 0; i < BL_SPACEDIM; i++) os << amrdToMimic.ProbLo()[i] << ' ';
    os << '\n';
    for (i = 0; i < BL_SPACEDIM; i++) os << amrdToMimic.ProbHi()[i] << ' ';
    os << '\n';
    for (i = 0; i < finestLevel; i++) os << amrdToMimic.RefRatio()[i] << ' ';
    os << '\n';
    for (i = 0; i <= finestLevel; i++) os << probDomain[i] << ' ';
    os << '\n';
    for (i = 0; i <= finestLevel; i++) os << 0 << ' ';
    os << '\n';
    for (i = 0; i <= finestLevel; i++)
    {
        for (int k = 0; k < BL_SPACEDIM; k++)
            os << amrdToMimic.DxLevel()[i][k] << ' ';
        os << '\n';
    }
    os << amrdToMimic.CoordSys() << '\n';
    os << "0\n"; // The bndry data width.
    //
    // Write out level by level.
    //
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        //
        // Write state data.
        //
        int nGrids = amrdToMimic.boxArray(iLevel).size();

        std::string LevelStr = amrex::Concatenate("Level_", iLevel, 1);
    
        if (ParallelDescriptor::IOProcessor())
        {
            os << iLevel << ' ' << nGrids << ' ' << amrdToMimic.Time() << '\n';
            os << 0 << '\n';
    
            for (i = 0; i < nGrids; ++i)
            {
                for (int n = 0; n < BL_SPACEDIM; n++)
                {
                    os << amrdToMimic.GridLocLo()[iLevel][i][n]
                       << ' '
                       << amrdToMimic.GridLocHi()[iLevel][i][n]
                       << '\n';
                }
            }
            //
            // Build the directory to hold the MultiFabs at this level.
            //
            std::string Level(oFile);
            Level += '/';
            Level += LevelStr;
    
            if (!amrex::UtilCreateDirectory(Level, 0755))
                amrex::CreateDirectoryFailed(Level);
        }
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        //
        // Now build the full relative pathname of the MultiFab.
        //
        static const std::string MultiFabBaseName("/MultiFab");
    
        std::string PathName(oFile);
        PathName += '/';
        PathName += LevelStr;
        PathName += MultiFabBaseName;
    
        if (ParallelDescriptor::IOProcessor())
        {
            //
            // The full name relative to the Header file.
            //
            std::string RelativePathName(LevelStr);
            RelativePathName += '/';
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';
        }
        VisMF::Write(*mfa[iLevel], PathName, VisMF::OneFilePerCPU);
    }

    os.close();
}

void WritePlotFile(const Vector<MultiFab*>&   mfa,
		   AmrData&                   amrdToMimic,
		   const std::string&         oFile,
		   bool                       verbose,
                   const Vector<std::string>& varNames)
{
    WritePlotFile(mfa,amrdToMimic.ProbDomain(),amrdToMimic,oFile,verbose,varNames);
}
