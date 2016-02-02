// --------------------------------------------------------------------------
// MKDir.cpp
// --------------------------------------------------------------------------
//   this file tests making directories.
// --------------------------------------------------------------------------
#include <winstd.H>

#include <new>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>

#include <ParallelDescriptor.H>
#include <Utility.H>

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    BoxLib::Initialize(argc,argv);    
    BL_PROFILE_VAR("main()", pmain);
    const Real tStart(ParallelDescriptor::second());
    int nprocs(ParallelDescriptor::NProcs());

    int ndirs(256), nlevels(4);

    for(int i(0); i < ndirs; ++i) {
      std::stringstream dirname;
      dirname << "dir" << i;
      if(ParallelDescriptor::IOProcessor()) {
        BL_PROFILE_VAR("mkdirs", mkdirs);
        if( ! BoxLib::UtilCreateDirectory(dirname.str(), 0755)) {
          BoxLib::CreateDirectoryFailed(dirname.str());
        }
        BL_PROFILE_VAR_STOP(mkdirs);
        for(int level(0); level < nlevels; ++level) {
          std::stringstream dirname;
          dirname << "dir" << i << "/Level_" << level;
          BL_PROFILE_VAR("mkdirs", mkdirs);
          if( ! BoxLib::UtilCreateDirectory(dirname.str(), 0755)) {
            BoxLib::CreateDirectoryFailed(dirname.str());
          }
          BL_PROFILE_VAR_STOP(mkdirs);
        }
      }
      ParallelDescriptor::Barrier("waitfordir");
    }

    for(int i(0); i < ndirs; ++i) {
      BL_PROFILE_VAR("renamedirs", renamedirs);
      if(ParallelDescriptor::IOProcessor()) {
        std::stringstream dirname;
        dirname << "dir" << i;
        std::string newdirname;
        newdirname = dirname.str() + ".old";
	std::rename(dirname.str().c_str(), newdirname.c_str());
      }
      BL_PROFILE_VAR_STOP(renamedirs);
    }


    Real runTime(ParallelDescriptor::second() - tStart);
    ParallelDescriptor::ReduceRealMax(runTime, ParallelDescriptor::IOProcessorNumber());

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << std::endl;
      std::cout << "Finished." << std::endl;
      std::cout << "Run time = " << runTime << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
    }

    BL_PROFILE_VAR_STOP(pmain);
    BoxLib::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
