// ---------------------------------------------------------------
// sliceutils.cpp
// ---------------------------------------------------------------
#include <iostream>
#include <stdint.h>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;


// ---------------------------------------------------------------
int main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  MultiFab sliceMF;

  // ---- this example reads a multifab named State_x from the
  // ---- directory slice_00340 and makes a plotfile named
  // ---- plt_slice_00340

  std::string plotfileName("plt_slice_00340");
  std::string sliceDirName("slice_00340");
  std::string sliceMFName(sliceDirName + "/State_x");
  int level_step(340);

  Vector<std::string> varnames(6);
  varnames[0] = "var0";
  varnames[1] = "var1";
  varnames[2] = "var2";
  varnames[3] = "var3";
  varnames[4] = "var4";
  varnames[5] = "var5";

  amrex::VisMF::Read(sliceMF, sliceMFName);

  const BoxArray &ba = sliceMF.boxArray();
  Box smfDomain(ba.minimalBox());
  RealBox realDomain(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
  Geometry geom(smfDomain, &realDomain);
  Real sliceTime(123.4);

  amrex::WriteSingleLevelPlotfile(plotfileName, sliceMF, varnames,
                                  geom, sliceTime, level_step);
    
  amrex::Finalize();
}

// ---------------------------------------------------------------
// ---------------------------------------------------------------
