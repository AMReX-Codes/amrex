// --------------------------------------------------------------------------
// MKDir.cpp
// --------------------------------------------------------------------------
//   this file tests making directories.
// --------------------------------------------------------------------------

#include <new>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

using namespace amrex;

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    amrex::Initialize(argc,argv);    
    BL_PROFILE_VAR("main()", pmain);
    const Real tStart(ParallelDescriptor::second());
    int nprocs(ParallelDescriptor::NProcs());

    int ndirs(256), nlevels(4);

    if(ParallelDescriptor::IOProcessor()) {
      errno = 0;
      mkdir("testdir", 0755);
      std::cout << "_here 0:  errno = " << strerror(errno) << std::endl;
      errno = 0;
      rmdir("testdir");
      std::cout << "_here 1:  errno = " << strerror(errno) << std::endl;
      errno = 0;
      mkdir("testnest/n0/n1", 0755);
      std::cout << "_here 2:  errno = " << strerror(errno) << std::endl;
      errno = 0;
    }

    BL_PROFILE_VAR("mkdirs", mkdirs);
    for(int i(0); i < ndirs; ++i) {
      std::stringstream dirname;
      dirname << "dir" << i;
      if(ParallelDescriptor::IOProcessor()) {
        if( ! amrex::UtilCreateDirectory(dirname.str(), 0755)) {
          amrex::CreateDirectoryFailed(dirname.str());
        }
        for(int level(0); level < nlevels; ++level) {
          std::stringstream dirname;
          dirname << "dir" << i << "/Level_" << level;
          if( ! amrex::UtilCreateDirectory(dirname.str(), 0755)) {
            amrex::CreateDirectoryFailed(dirname.str());
          }
        }
      }
    }
    ParallelDescriptor::Barrier("waitfordir");
    BL_PROFILE_VAR_STOP(mkdirs);

    BL_PROFILE_VAR("renamedirs", renamedirs);
    for(int i(0); i < ndirs; ++i) {
      if(ParallelDescriptor::IOProcessor()) {
        std::stringstream dirname;
        dirname << "dir" << i;
        std::string newdirname;
        newdirname = dirname.str() + ".old";
	std::rename(dirname.str().c_str(), newdirname.c_str());
      }
    }
    ParallelDescriptor::Barrier("renamedirs");
    BL_PROFILE_VAR_STOP(renamedirs);


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
    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
