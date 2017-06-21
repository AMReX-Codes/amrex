#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

using namespace amrex;

void
writePlotFile (const std::string& dir,
               const MultiFab&    mf,
               const Geometry&    geom,
	       const Real&        time)
{

  // Only let 64 CPUs be writing at any one time.
  VisMF::SetNOutFiles(64);
  
  // Only the I/O processor makes the directory if it doesn't already exist.
  if (ParallelDescriptor::IOProcessor())
    if (!amrex::UtilCreateDirectory(dir, 0755))
      amrex::CreateDirectoryFailed(dir);

  // Force other processors to wait till directory is built.
  ParallelDescriptor::Barrier();

  std::string HeaderFileName = dir + "/Header";
  
  VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
  
  std::ofstream HeaderFile;
  
  HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
  
  if (ParallelDescriptor::IOProcessor()) {

    // Only the IOProcessor() writes to the header file.
    HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
    if (!HeaderFile.good())
      amrex::FileOpenFailed(HeaderFileName);
    HeaderFile << "NavierStokes-V1.1\n";
    
    HeaderFile << mf.nComp() << '\n';
    
    // variable names
    HeaderFile << "density\n";
	
    // dimensionality
    HeaderFile << BL_SPACEDIM << '\n';
    // time
    HeaderFile << time << '\n';
    // maximum level number (0=single level)
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
  static const std::string BaseName = "/Cell";
  
  std::string Level = amrex::Concatenate("Level_", 0, 1);

  // Now for the full pathname of that directory.
  std::string FullPath = dir;
  if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
    FullPath += '/';
  FullPath += Level;

  // Only the I/O processor makes the directory if it doesn't already exist.
  if (ParallelDescriptor::IOProcessor())
    if (!amrex::UtilCreateDirectory(FullPath, 0755))
      amrex::CreateDirectoryFailed(FullPath);

  // Force other processors to wait till directory is built.
  ParallelDescriptor::Barrier();
  
  if (ParallelDescriptor::IOProcessor()) {
    HeaderFile << 0 << ' ' << mf.boxArray().size() << ' ' << 0 << '\n';
    HeaderFile << 0 << '\n';
    
    for (int i = 0; i < mf.boxArray().size(); ++i) {
      RealBox loc = RealBox(mf.boxArray()[i],geom.CellSize(),geom.ProbLo());
      for (int n = 0; n < BL_SPACEDIM; n++)
	HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
    }
    
    std::string PathNameInHeader = Level;
    PathNameInHeader += BaseName;
    HeaderFile << PathNameInHeader << '\n';
  }
  
  // Use the Full pathname when naming the MultiFab.
  std::string TheFullPath = FullPath;
  TheFullPath += BaseName;
  
  VisMF::Write(mf, TheFullPath);
}

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  bool verbose;
};

void test_assign_density(TestParams& parms)
{

  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  IntVect domain_lo(0 , 0, 0); 
  IntVect domain_hi(parms.nx - 1, parms.ny - 1, parms.nz-1); 
  const Box domain(domain_lo, domain_hi);

  // This says we are using Cartesian coordinates
  int coord = 0;

  // This sets the boundary conditions to be doubly or triply periodic
  int is_per[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++) 
    is_per[i] = 1; 
  Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);

  BoxArray ba(domain);
  ba.maxSize(parms.max_grid_size);
  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << "Number of boxes              : " << ba[0].size() << '\n' << '\n';
  }

  DistributionMapping dmap(ba);

  MultiFab acceleration(ba, dmap, 3, 1);
  acceleration.setVal(5.0, 1);

  MultiFab density(ba, dmap, 1, 0);
  density.setVal(0.0);

  MultiFab partMF(ba, dmap, 1, 1);
  partMF.setVal(0.0);

  typedef ParticleContainer<1> MyParticleContainer;
  MyParticleContainer myPC(geom, dmap, ba);
  myPC.SetVerbose(false);

  int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Total number of particles    : " << num_particles << '\n' << '\n';

  bool serialize = true;
  int iseed = 451;
  Real mass = 10.0;

  MyParticleContainer::ParticleInitData pdata = {mass};
  myPC.InitRandom(num_particles, iseed, pdata, serialize);
  myPC.AssignCellDensitySingleLevelFort(0, partMF, 0, 1, 0);
  
  myPC.InterpolateSingleLevelFort(acceleration, 0);

  MultiFab::Copy(density, partMF, 0, 0, 1, 0);

  writePlotFile("plt00000", density, geom, 0.0);
  myPC.Checkpoint("plt00000", "particle0", true);
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  
  ParmParse pp;
  
  TestParams parms;
  
  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);
  pp.get("nppc", parms.nppc);
  if (parms.nppc < 1 && ParallelDescriptor::IOProcessor())
    amrex::Abort("Must specify at least one particle per cell");
  
  parms.verbose = false;
  pp.query("verbose", parms.verbose);
  
  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles per cell : ";
    std::cout << parms.nppc  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }
  
  test_assign_density(parms);
  
  amrex::Finalize();
}
