
//
// Write's `.mat' files for MatLab -- one per component.
//

#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include "AMReX_Utility.H"
#include <AMReX_VisMF.H>

using namespace amrex;

//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

std::list<std::string> plot_vars;

static
bool
isPlotVar (const std::string& name)
{
  for (std::list<std::string>::const_iterator li = plot_vars.begin();
       li != plot_vars.end();
       ++li)
  {
    if (*li == name)
      return true;
  }

  return false;
}

static
void
PrintUsage (char* progName)
{
  std::cout << "\nUsage:\n"
            << progName
            << "\n\tinfile = inputFileName"
            << "\n\tplot_vars = list of plot variables (none specified --> ALL)"
            << "\n\t[-help]"
            << "\n\n";
  exit(1);
}

//
// Special MatLab value specifying float or double.
//

#ifdef BL_USE_DOUBLE
const int RealType = 0;
#else
const int RealType = 10;
#endif

//
// Special MatLab value specifying byte order.
//
static
int
ByteOrder ()
{
  const int BigEndian   = 1000;
  const int SmallEndian = 0;

  union
  {
    long Long;
    char Char[sizeof(long)];
  }
  SwapTest;

  SwapTest.Long = 1;

  return SwapTest.Char[0] == 1 ? SmallEndian : BigEndian;
}

static
void
WriteFab (std::ostream&         os,
          const FArrayBox& fab,
          const char*      name)

{
  int nx = fab.box().length(0);
  int ny = fab.box().length(1);
  int dim = BL_SPACEDIM;


#if (BL_SPACEDIM == 2)
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
    {
      int index = j*nx + i;
      const Real * ptr = fab.dataPtr();
      os.write((char*)(ptr+index),sizeof(Real));
    }

#elif (BL_SPACEDIM == 3)
  int nz = fab.box().length(2);
  os.write((char*)&nz,sizeof(int));

  for (int k = 0; k < nz; k++)
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
      {
        int index = k*(nx*ny) + j*nx + i;
        const Real * ptr = fab.dataPtr();
        os.write((char*)(ptr+index),sizeof(Real));
      }
#endif
}

static
void
WriteLoc (std::ostream&           os,
          Real* lo, Real* hi)
{

  Real buf[2*BL_SPACEDIM];

  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    buf[i]             = lo[i];
    buf[i+BL_SPACEDIM] = hi[i];
  }
  for (int i = 0; i < 2*BL_SPACEDIM; i++)
    os.write((char*)&(buf[i]),sizeof(Real));
}

static
void
Write (AmrData&       amrData,
       const std::string& iFile_name,
       const std::list<std::string> plot_vars)
{
  int finest_level = amrData.FinestLevel();

  VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

  char buf[128];

  for (int icomp = 0; icomp < amrData.NComp(); icomp++)
  {
    //
    // Write one component per file.
    //
    const std::string& CompName = amrData.PlotVarNames()[icomp];
    if (isPlotVar(CompName))
    {

      std::string file = iFile_name;
      file += '_';
      file += CompName;
      file += ".mat";

      std::ofstream os;

      os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

      os.open(file.c_str(), std::ios::out|std::ios::binary);

      if (os.fail()) {
        amrex::FileOpenFailed(file);
      }

      int dim = BL_SPACEDIM;
      os.write((char*)&dim,sizeof(int));

      int num_levels = finest_level+1;
      os.write((char*)&num_levels,sizeof(int));

//      Write the number of grids at each level
      for (int iLevel = 0; iLevel <= finest_level; ++iLevel)
      {
        MultiFab& mf = amrData.GetGrids(iLevel,icomp);
        const BoxArray& ba = mf.boxArray();
        int num_grids = ba.size();
        os.write((char*)&num_grids,sizeof(int));
      }

//      Write the (Real) physical locations of each grid at each level
      for (int iLevel = 0; iLevel <= finest_level; ++iLevel)
      {
        MultiFab& mf = amrData.GetGrids(iLevel,icomp);
        const BoxArray& ba = mf.boxArray();
        int num_grids = ba.size();
        for (int i = 0; i < num_grids; ++i)
        {
          for (int idim = 0; idim < BL_SPACEDIM; ++idim)
          {
            Real xlo = amrData.GridLocLo()[iLevel][i][idim];
            Real xhi = amrData.GridLocHi()[iLevel][i][idim];
            os.write((char*)&xlo,sizeof(Real));
            os.write((char*)&xhi,sizeof(Real));
          }
        }
      }

//      Write the (integer) dimensions of each grid at each level
      for (int iLevel = 0; iLevel <= finest_level; ++iLevel)
      {
        MultiFab& mf = amrData.GetGrids(iLevel,icomp);
        const BoxArray& ba = mf.boxArray();
        int num_grids = ba.size();
        for (int i = 0; i < num_grids; ++i)
        {
          for (int idim = 0; idim < BL_SPACEDIM; ++idim)
          {
            int n = ba[i].length(idim);
            os.write((char*)&n,sizeof(int));
          }
        }
      }

//      Write the (Real) actual data of each grid at each level
      for (int iLevel = 0; iLevel <= finest_level; ++iLevel)
      {
        MultiFab& mf = amrData.GetGrids(iLevel,icomp);
        const BoxArray& ba = mf.boxArray();
        int num_grids = ba.size();
        for (int ig = 0; ig < num_grids; ++ig) 
        {
          sprintf(buf,
                  "%s_%d_%d",
                  CompName.c_str(),
                  iLevel,
                  ig);
          WriteFab(os, mf[ig], buf);
        }
      }

      os.close();

      if (os.fail()) {
        amrex::FileOpenFailed(file);
      }
    }
  }
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    if (argc == 1)
      PrintUsage(argv[0]);

    if (ParallelDescriptor::NProcs() > 1) {
      amrex::Error("This is an inherently serial program");
    }

    ParmParse pp;

    if (pp.contains("help"))
      PrintUsage(argv[0]);
    //
    // MatLab expects native floating-point format.
    //
    FArrayBox::setFormat(FABio::FAB_NATIVE);
    //
    // Scan the arguments.
    //
    std::string iFile;
    pp.query("infile", iFile);
    if (iFile.empty()) {
      amrex::Abort("You must specify `infile'");
    }

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    
    DataServices dataServices(iFile, fileType);

    if (!dataServices.AmrDataOk())
    {
      //
      // This calls ParallelDescriptor::EndParallel() and exit()
      //
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    std::string plot_var;
    pp.query("plot_vars", plot_var);

    int npv = pp.countval("plot_vars");

    if (npv == 0) {
      for (int i = 0; i < amrData.PlotVarNames().size(); i++)
      {
        plot_vars.push_back(amrData.PlotVarNames()[i]);
      }

    } else {

      for (int i = 0; i < npv; i++)
      {
        pp.get("plot_vars", plot_var, i);
        plot_vars.push_back(plot_var);
      }

    }

    Write(amrData, iFile, plot_vars);
  }
  Finalize();
}

