#include <string>
#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>

#include <AMReX_BLFort.H>

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
        const Box& domain, int idir, int jdir, int kdir, int comp)
{
  Box tbx;
  tbx.setSmall(IntVect::TheZeroVector());
  tbx.setBig(0, domain.bigEnd(idir) + 3);
  tbx.setBig(1, domain.bigEnd(jdir) + 3);
  tbx.setBig(2, 0);
  const Box& vbx = vfab.box();

  xfab.resize(tbx,1);
  for (int i=0; i<vbx.length(idir); ++i) {
    for (int j=0; j<vbx.length(jdir); ++j) {
      IntVect iv;
      iv[kdir] = vbx.smallEnd(kdir);
      iv[idir] = i + vbx.smallEnd(idir);
      iv[jdir] = j + vbx.smallEnd(jdir);
      xfab(IntVect(i,j,0),0) = vfab(iv,comp);
    }
  }

  for (int i=0; i<3; ++i) {
    for (int j=0; j<vbx.length(jdir); ++j) {
      IntVect iv;
      iv[kdir] = vbx.smallEnd(kdir);
      iv[idir] = i + vbx.smallEnd(idir);
      iv[jdir] = j + vbx.smallEnd(jdir);
      xfab(IntVect(i+vbx.length(idir),j,0),0) = vfab(iv,comp);
    }
  }

  for (int i=0; i<vbx.length(idir); ++i) {
    for (int j=0; j<3; ++j) {
      IntVect iv;
      iv[kdir] = vbx.smallEnd(kdir);
      iv[idir] = i + vbx.smallEnd(idir);
      iv[jdir] = j + vbx.smallEnd(jdir);
      xfab(IntVect(i,j+vbx.length(jdir),0),0) = vfab(iv,comp);
    }
  }

  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      IntVect iv;
      iv[kdir] = vbx.smallEnd(kdir);
      iv[idir] = i + vbx.smallEnd(idir);
      iv[jdir] = j + vbx.smallEnd(jdir);
      xfab(IntVect(i+vbx.length(idir),j+vbx.length(jdir),0),0) = vfab(iv,comp);
    }
  }
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

    const std::string farg = amrex::get_command_argument(1);
    if (farg == "-h" || farg == "--help")
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    int nf = pp.countval("ifiles");
    AMREX_ALWAYS_ASSERT(nf >= 3);
    Vector<std::string> ifiles(nf);
    pp.getarr("ifiles",ifiles,0,nf);

    std::string ofile;
    pp.get("ofile",ofile);

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

    Box domain;
    {
      std::ifstream ifs(ifiles[0]);
      FArrayBox fab;
      fab.readFrom(ifs);
      ifs.close();
      domain = fab.box();
    }
    int idir = 1;
    int jdir = 2;
    int kdir = 0;
    IntVect bg;
    bg[0] = domain.length(idir);
    bg[1] = domain.length(jdir);
    bg[2] = nf;
    //
    // Write the first part of the header.
    //
    ifsh << bg[0] + 3 << ' '
         << bg[1] + 3 << ' '
         << bg[2] << '\n';

    Vector<Real> planeSize(2), dx(2);
    pp.getarr("planeSize",planeSize,0,2);
    dx[0] = planeSize[0] / bg[0];
    dx[1] = planeSize[1] / bg[1];
    ifsh << planeSize[0] + 2*dx[0] << ' '
         << planeSize[1] + 2*dx[1] << ' '
         << nf << '\n';

    ifsh << 1 << ' ' << 1 << ' ' << 0 << '\n';

    FArrayBox TMP;
    Vector<Real> fileTimes(nf);
    for (int d = 0; d < BL_SPACEDIM; ++d)
    {
      std::cout << "Loading component " << d << " ... " << std::flush;

      for (int k = 0; k < nf; ++k)
      {
        FArrayBox fileFab;
        {
          std::ifstream ifs(ifiles[k]);
          fileFab.readFrom(ifs);
          if (d==0) {
            ifs >> fileTimes[k];
          }
          ifs.close();
        }
        Extend(TMP, fileFab, domain, idir, jdir, kdir, d);
        //
        // Write current position of data file to header file.
        //
        ifsh << ifsd.tellp() << std::endl;
        //
        // Write the FAB to the data file.
        //
        TMP.writeOn(ifsd);
      }
      std::cout << "done" << std::endl;
    }
    // Write plane times
    for (int k = 0; k < nf; ++k)
    {
      ifsh << fileTimes[k] << std::endl;
    }
  }
  Finalize();
}
