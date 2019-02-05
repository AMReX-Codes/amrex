
#include <fstream>
#include <iomanip>

#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

#include <writePlotFile.H>

using namespace amrex;

void writePlotFile (const std::string& dir, 
		    const Vector<MultiFab*>& soln,
		    const Vector<MultiFab*>& exac, 
		    const Vector<MultiFab*>& alph,
		    const Vector<MultiFab*>& beta, 
		    const Vector<MultiFab*>& rhs, 
		    const std::vector<Geometry>& geom, 
		    const std::vector<BoxArray>& grids,
		    int nsoln, int iCpp, int iHyp)
{
  int n_data_items = 4+nsoln;
  int nlevels = exac.size();
  int f_lev = nlevels-1;

    //
    // Only let 64 CPUs be writing at any one time.
    //
    VisMF::SetNOutFiles(64);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(dir, 0755))
            amrex::CreateDirectoryFailed(dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    std::string HeaderFileName = dir + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);
        HeaderFile << "HyperCLaw-V1.1\n";

        HeaderFile << n_data_items << '\n';

	if (iCpp >= 0) {
	  HeaderFile << "solnCpp\n";
	}
	if (iHyp >= 0) {
	  HeaderFile << "solnHyp\n";	  
	}
	HeaderFile << "exac\n";
	HeaderFile << "alph\n";
	HeaderFile << "beta\n";
	HeaderFile << "rhs\n";

        HeaderFile << BL_SPACEDIM << '\n';
        HeaderFile << 0 << '\n';
        HeaderFile << f_lev << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom[0].ProbLo(i) << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom[0].ProbHi(i) << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < f_lev; i++) {
	  HeaderFile << 2 <<' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i <= f_lev; i++) {
	  HeaderFile << geom[i].Domain() << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i <= f_lev; i++) {
	  HeaderFile << 0 << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i <= f_lev; i++) {
	  for (int k = 0; k < BL_SPACEDIM; k++) {
            HeaderFile << geom[i].CellSize()[k] << ' ';
	  }
	  HeaderFile << '\n';
	}
        HeaderFile << geom[0].Coord() << '\n';
        HeaderFile << "0\n";
    }

    
    for (int ilev=0; ilev < nlevels; ilev++) {
      // Build the directory to hold the MultiFab at this level.
      // The name is relative to the directory containing the Header file.
      //
      static const std::string BaseName = "/Cell";

      std::string Level = amrex::Concatenate("Level_", ilev, 1);
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

      if (ParallelDescriptor::IOProcessor()) {
	HeaderFile << ilev << ' ' << grids[ilev].size() << ' ' << 0 << '\n';
	HeaderFile << 0 << '\n';
	
	for (int i = 0; i < grids[ilev].size(); ++i) {
	  RealBox loc = RealBox(grids[ilev][i], geom[ilev].CellSize(), geom[ilev].ProbLo());
	  for (int n = 0; n < BL_SPACEDIM; n++) {
	    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
	  }
	}
	//
	// The full relative pathname of the MultiFabs at this level.
	// The name is relative to the Header file containing this name.
	// It's the name that gets written into the Header.
	//
	if (n_data_items > 0) {
	  std::string PathNameInHeader = Level;
	  PathNameInHeader += BaseName;
	  HeaderFile << PathNameInHeader << '\n';
	}
      }
      //
      // We combine all of the multifabs 
      //
      const DistributionMapping& dm = soln[ilev]->DistributionMap();
      MultiFab  plotMF(grids[ilev], dm, n_data_items, 0);
      //
      // Cull data -- use no ghost cells.
      //
      int cnt=0;
      for (int isoln=0; isoln < nsoln; isoln++) {
	MultiFab::Copy(plotMF, *soln[ilev], isoln, cnt, 1, 0);
	cnt++;
      }
      MultiFab::Copy(plotMF, *exac[ilev], 0, cnt, 1, 0);
      cnt++;
      MultiFab::Copy(plotMF, *alph[ilev], 0, cnt, 1, 0);
      cnt++;
      MultiFab::Copy(plotMF, *beta[ilev], 0, cnt, 1, 0);
      cnt++;
      MultiFab::Copy(plotMF,  *rhs[ilev], 0, cnt, 1, 0);

      //
      // Use the Full pathname when naming the MultiFab.
      //
      std::string TheFullPath = FullPath;
      TheFullPath += BaseName;
      VisMF::Write(plotMF,TheFullPath,VisMF::OneFilePerCPU,true);
    }
}

