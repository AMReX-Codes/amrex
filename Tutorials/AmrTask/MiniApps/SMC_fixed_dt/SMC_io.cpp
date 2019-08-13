#include <fstream>

#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <SMC.H>
#include <SMC_F.H>

using namespace amrex;

void
SMC::writePlotFile (int istep)
{
    const std::string& dir = amrex::Concatenate("plt",istep,5);

    MultiFab mf(U.boxArray(), U.DistributionMap(), nplot, 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox();
	make_plot_3d(bx.loVect(), bx.hiVect(),
		     BL_TO_FORTRAN_3D(mf[mfi]),
		     BL_TO_FORTRAN_3D( U[mfi]));
    }

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
        HeaderFile << "NavierStokes-V1.1\n";

        HeaderFile << mf.nComp() << '\n';

	// For now, set this manually
	HeaderFile << "rho" << "\n";
	HeaderFile << "vx" << "\n";
	HeaderFile << "vy" << "\n";
	HeaderFile << "vz" << "\n";
	HeaderFile << "T" << "\n";
	HeaderFile << "Y(H2)" << "\n";
	HeaderFile << "Y(O2)" << "\n";
	HeaderFile << "Y(H2O)" << "\n";
	HeaderFile << "Y(H)" << "\n";
	HeaderFile << "Y(O)" << "\n";
	HeaderFile << "Y(OH)" << "\n";
	HeaderFile << "Y(HO2)" << "\n";
	HeaderFile << "Y(H2O2)" << "\n";
	HeaderFile << "Y(N2)" << "\n";

        HeaderFile << BL_SPACEDIM << '\n';
        HeaderFile << t << '\n';
        HeaderFile << 0 << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom.ProbLo(i) << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom.ProbHi(i) << ' ';
        HeaderFile << '\n';
        HeaderFile << '\n';
        HeaderFile << geom.Domain() << ' ';
        HeaderFile << '\n';
        HeaderFile << 0 << ' ';
        HeaderFile << '\n';
        for (int k = 0; k < BL_SPACEDIM; k++)
            HeaderFile << geom.CellSize()[k] << ' ';
        HeaderFile << '\n';
        HeaderFile << geom.Coord() << '\n';
        HeaderFile << "0\n";
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";

    std::string Level = amrex::Concatenate("Level_", 0, 1);
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
        HeaderFile << 0 << ' ' << mf.boxArray().size() << ' ' << 0 << '\n';
        HeaderFile << 0 << '\n';

        for (int i = 0; i < mf.boxArray().size(); ++i)
        {
            RealBox loc(mf.boxArray()[i],geom.CellSize(),geom.ProbLo());
            for (int n = 0; n < BL_SPACEDIM; n++)
                HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
        }

        std::string PathNameInHeader = Level;
        PathNameInHeader += BaseName;
        HeaderFile << PathNameInHeader << '\n';
    }
    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;

    VisMF::Write(mf,TheFullPath);
}

