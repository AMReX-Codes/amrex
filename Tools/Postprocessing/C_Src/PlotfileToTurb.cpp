/*
  Read a plotfile containing velocity fields and generate a file suitable for time-dependent input into
  IAMR or PeleLM through the boundary.  The format of such files is a folder containing an ascii file, HDR,
  and a binary file, DAT. The data is uniformly spaced, grow-extended and suitable for componentwise parabolic
  interpolation.  In 3D, each ij plane will be parallel to the inflow face of the intended application,
  and k will serve as a time coordinate.  The HDR file includes the size and grid-spacing of the data,
  followed by a list of offsets into the DATA file for each of the catenated planar fabs of velocity data.
  In 2D the j axis represents time, and the data is written as a single plane.
 */
#include <string>
#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>

#include <AMReX_BLFort.H>

extern "C" {
  void fliprowsy(amrex_real* u, ARLIM_P(ulo), ARLIM_P(uhi));
}

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " ifile=<pltfile> ofile=<turbname>\n";
  exit(1);
}

static
void
Extend (FArrayBox& xfab,
        FArrayBox& vfab,
        const Box& domain)
{
  Box tbx = vfab.box();

  tbx.setBig(0, domain.bigEnd(0) + 3);

  const int ygrow = BL_SPACEDIM==3 ? 3 : 1;

  tbx.setBig(1, domain.bigEnd(1) + ygrow);

  xfab.resize(tbx,1);

  xfab.copy(vfab);
  vfab.shift(0, domain.length(0));
  xfab.copy(vfab);
  vfab.shift(1, domain.length(1));
  xfab.copy(vfab);
  vfab.shift(0, -domain.length(0));
  xfab.copy(vfab);
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string ifile;
    pp.get("ifile",ifile);

    std::string ofile;
    pp.get("ofile",ofile);

    std::cout << "Reading " << ifile << " ... " << std::flush;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(ifile, fileType);
    std::cout << "done" << std::endl;

    if (!dataServices.AmrDataOk())
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    AmrData& amrData = dataServices.AmrDataRef();

    std::string TurbDir = ofile;

    if (ParallelDescriptor::IOProcessor())
      if (!UtilCreateDirectory(TurbDir, 0755))
        CreateDirectoryFailed(TurbDir);

    std::string Hdr = TurbDir; Hdr += "/HDR";
    std::string Dat = TurbDir; Dat += "/DAT";

    std::ofstream ifsd, ifsh;

    ifsh.open(Hdr.c_str(), std::ios::out|std::ios::trunc);
    if (!ifsh.good())
      FileOpenFailed(Hdr);

    ifsd.open(Dat.c_str(), std::ios::out|std::ios::trunc);
    if (!ifsd.good())
      FileOpenFailed(Dat);

    const int dir                   = BL_SPACEDIM - 1;
    const int           finestLevel = amrData.FinestLevel();
    const Box           domain      = amrData.ProbDomain()[finestLevel];
    Box                 xdm         = domain;
    const Vector<Real>& dx          = amrData.DxLevel()[finestLevel];
    IntVect             sm          = domain.smallEnd();
    IntVect             bg          = domain.bigEnd();
    std::string         names[3]    = { "x_velocity", "y_velocity", "z_velocity" };
    FArrayBox           xfab;
    //
    // Write the first part of the header.
    // Note that this is solely for periodic style inflow files.
    //
#if BL_SPACEDIM==2
    xdm.setBig(0, domain.bigEnd(0) + 3);
    xdm.setBig(1, domain.bigEnd(1) + 1);

    ifsh << xdm.length(0) << ' '
         << xdm.length(1) << ' '
         << 1             << '\n';

    ifsh << amrData.ProbSize()[0] + 2*dx[0] << ' '
         << amrData.ProbSize()[1]           << ' '
         << 0                               << '\n';

    ifsh << 1 << ' ' << 1 << ' ' << 0 << '\n';
#elif BL_SPACEDIM==3
    xdm.setBig(0, domain.bigEnd(0) + 3);
    xdm.setBig(1, domain.bigEnd(1) + 3);
    xdm.setBig(2, domain.bigEnd(2) + 1);

    ifsh << xdm.length(0) << ' '
         << xdm.length(1) << ' '
         << xdm.length(2) << '\n';

    ifsh << amrData.ProbSize()[0] + 2*dx[0] << ' '
         << amrData.ProbSize()[1] + 2*dx[1] << ' '
         << amrData.ProbSize()[2]           << '\n';

    ifsh << 1 << ' ' << 1 << ' ' << 1 << '\n';
#else
    Abort("Only 2-D & 3-D supported");
#endif

    for (int d = 0; d < BL_SPACEDIM; ++d)
    {
      std::cout << "Loading component " << d << " ... " << std::flush;

      BoxArray ba(1);

#if BL_SPACEDIM==3
      //
      // In 3-D we work on one cell wide Z-planes.
      // We first do the lo BL_SPACEDIM plane.
      // And then all the other planes in xhi -> xlo order.
      //
      // In 2-D we only write a single x-y plane of data per component.
      //
      bg[dir] = sm[dir];
#endif
      Box bx(sm,bg);
      ba.set(0,bx);
      DistributionMapping dmap(ba);
      MultiFab TMP(ba,dmap,1,0);

      amrData.FillVar(TMP, amrData.FinestLevel(), names[d], 0);
      amrData.FlushGrids(amrData.StateNumber(names[d]));

      Extend(xfab, TMP[0], domain);
      //
      // Write current position of data file to header file.
      //
      ifsh << ifsd.tellp() << std::endl;
#if BL_SPACEDIM==2
      //
      // Write the FAB to the data file after flipping rows in Y direction.
      //
      fliprowsy(BL_TO_FORTRAN(xfab));
      xfab.writeOn(ifsd);
#elif BL_SPACEDIM==3

      xfab.writeOn(ifsd);

      std::cout << "xfab min,max: " << xfab.min() << ' ' << xfab.max() << ' ' << bx << '\n';
      //
      // Now do all the planes in zhi -> zlo order.
      //
      for (int i = domain.bigEnd(dir); i >= domain.smallEnd(dir); i--)
      {
        sm[dir] = i;
        bg[dir] = i;
        Box bx(sm,bg);
        ba.set(0,bx);
        DistributionMapping dmap(ba);
        MultiFab TMP(ba,dmap,1,0);

        amrData.FillVar(TMP, amrData.FinestLevel(), names[d], 0);
        amrData.FlushGrids(amrData.StateNumber(names[d]));

        Extend(xfab, TMP[0], domain);
        //
        // Write current position of data file to header file.
        //
        ifsh << ifsd.tellp() << std::endl;
        //
        // Write the FAB to the data file.
        //
        xfab.writeOn(ifsd);

        std::cout << "xfab i,min,max: " << i<< ' ' << xfab.min() << ' ' << xfab.max() << ' ' << bx << '\n';
      }
#endif
      std::cout << "done" << std::endl;
    }
  }
  Finalize();
}
