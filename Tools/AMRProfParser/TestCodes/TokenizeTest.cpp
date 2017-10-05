// --------------------------------------------------------------
// TokenizeTest.cpp
// --------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include <math.h>
using std::cout;
using std::endl;

#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>

using namespace amrex;


// --------------------------------------------------------------
int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  BL_PROFILE_REGION_START("main()");
  BL_PROFILE_VAR("main()", pmain);

  int myProc(ParallelDescriptor::MyProc());

  std::string ctFileName("CallTrace.txt");
  std::string fileLine;
  std::ifstream ifs(ctFileName);
  if(ifs.fail()) {
    cout << "**** Error:  could not open:  " << ctFileName << endl;
  } else {
    std::getline(ifs, fileLine);
    cout << endl << "fileLine = " << fileLine << endl;
    const std::vector<std::string> &lineTokens = amrex::Tokenize(fileLine, "::::");
    cout << "lineTokens.size() = " << lineTokens.size() << endl;

    while( ! ifs.eof()) {
      std::getline(ifs, fileLine);
      if( ! ifs.eof()) {
        cout << endl << "fileLine = " << fileLine << endl;
	char delims[1];
	//delims[0] = '\011';
	delims[0] = ' ';
	//std::string delimString(delims);
	std::string delimString(" ");
	const std::vector<std::string> &lineTokens = amrex::Tokenize(fileLine, delimString);
	cout << "lineTokens.size() = " << lineTokens.size() << endl;
	for(int i(0); i < lineTokens.size(); ++i) {
	  cout << i << "::token = " << lineTokens[i] << endl;
	}
      }
    }
  }
  ifs.close();



  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_REGION_STOP("main()");

  amrex::Finalize();
}

// --------------------------------------------------------------
// --------------------------------------------------------------

