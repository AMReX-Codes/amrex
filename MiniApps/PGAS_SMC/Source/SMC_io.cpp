
// #include <unistd.h>

// #include <iomanip>
// #include <iostream>
// #include <string>
// #include <ctime>

#include <Utility.H>
#include "SMC.H"

using std::string;

namespace 
{
    int version = -1;
}

// I/O routines for SMC

void
SMC::restart (Amr&     papa,
	      istream& is,
	      bool     bReadSpecial)
{

}

void
SMC::checkPoint(const std::string& dir,
		std::ostream&  os,
		VisMF::How     how,
		bool dump_old_default)
{
  AmrLevel::checkPoint(dir, os, how, false);
}

std::string
SMC::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");

    return the_plot_file_type;
}

void
SMC::setPlotVariables ()
{
  AmrLevel::setPlotVariables();
}

void
SMC::writePlotFile (const std::string& dir,
		    ostream&       os,
		    VisMF::How     how)
{

}
