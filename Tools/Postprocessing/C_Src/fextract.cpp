//
// C++ version of fextract.f90.
// This does not rely on FBaselib and is based on AmrDeriveLineValues.cpp
//
//  Author: Michele Rosso
//  Date  : May, 2018
//
#include <iostream>
#include <string.h>
#include <AMReX.H>
#include <AMReX_REAL.H>
#include "AMReX_Box.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"
#include "AMReX_Utility.H"
#include "AMReX_VisMF.H"
#include <stdlib.h>


using namespace amrex;
using namespace std;


//
// Prototypes
// 
void get_input_arguments ( const int argc, char** argv,
			   string& pltfile, string& slicefile,
			   Vector<string>& var, int& format, bool& lleft,
			   int& dir, Real& tol );
void help ();

//
// 
// 
int main ( int argc, char* argv[] )
{
    // Use "false" as third arguments to prevent AMReX
    // to check for an input file as first command line argument
    amrex::Initialize ( argc, argv, false );
    {

    if ( ParallelDescriptor::NProcs() > 1 )
	Abort("Only run this with one processor");
    
    // Input arguments
    string  dirname;
    string  pltfile;
    string  slicefile;
    Vector<string>  var;
    int     format = 10; 
    bool    lleft = false;
    int     dir = 1;
    Real    tol = 0.0;          
	
    get_input_arguments ( argc, argv, pltfile, slicefile, var, format,
			  lleft, dir, tol );

    // Set the name of the chosen direction 
    if ( dir==1 )
    {
	dirname = "x";
    }
    else if ( dir==2 )
    {
	dirname = "y";
    }
    else if ( dir==3 )
    {
	dirname = "z";
    }

    // For consistency with the fortran version, set format=format-1
    // This is because C++ scientific format is x.xxxxExxx
    // rather than Fortran 0.xxxxExxx.
    --format;
    

    // Start dataservices (no clue why we need to do this )
    DataServices::SetBatchMode();

    // Define the type of file
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices (pltfile, fileType);

    if( ! dataServices.AmrDataOk())
    {
	DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    // set amrData to be the handler to the data in the plotfile
    AmrData &amrData = dataServices.AmrDataRef();
    
    int finestLevel(amrData.FinestLevel());
    Box domain(amrData.ProbDomain()[finestLevel]);
    Box line(domain);		// Define a box covering the whole domain

    // Shrink "Box line" to coincide with the line of interest 
    for ( int idir = 0; idir < AMREX_SPACEDIM; ++idir)
    {
	if ( idir != (dir-1)  )
	{
	    int index;

	    index = domain.smallEnd ( idir );

	    if ( !lleft )
		index = domain.smallEnd(idir) + ( domain.length(idir) / 2 );
	    
	    line.setSmall ( idir, index );
	    line.setBig ( idir, index );
	}
    }
    
    // Now extract the variables of interest
    Vector< unique_ptr<FArrayBox> > lineFabs( var.size() );

    for (int nv = 0; nv < lineFabs.size(); ++nv )
    {
  	lineFabs[nv].reset( new FArrayBox(line,1) );
        amrData.FillVar(lineFabs[nv].get(), lineFabs[nv]->box(), finestLevel, var[nv],
  			ParallelDescriptor::IOProcessorNumber());
    }

    // 
    // Now print to file
    // 
    ofstream  outfile;

    outfile.open( slicefile.c_str(), ofstream::out | ofstream::trunc );

    // Print header
    outfile <<  "# 1-d slice in " + dirname + "-direction, file: " + pltfile << endl;

    outfile << scientific << uppercase << setprecision(format);
    
    outfile << "# time = " << setw(20) << amrData.Time() << endl;
    outfile << "#" << setw(24) << dirname;

    for (int nv = 0; nv < var.size(); ++nv )
	outfile << " " << setw(24) << var[nv].c_str();

    outfile << endl;

    // Print data
    IntVect iv = lineFabs[0] -> smallEnd();
    Vector<Real>  x(3);

    for (int m = 0; m < lineFabs[0]->size (); ++m )
    {
	amrData.CellLoc ( finestLevel, iv.shift(dir-1, std::min(1,m) ), x );
	
	outfile << setw(25) << x[dir-1];

	for (int nv = 0; nv < lineFabs.size(); ++nv )
	{
	    Real  val = *( lineFabs[nv]->dataPtr() + m );

	    // Set to zero all the values lower than a certain tolerance
	    if ( std::abs(val) < tol )
		val = 0.;
	    
	    outfile << setw(25) << val; 
	}
	
	outfile	<< endl;
    }
  
    outfile.close();
    
    }
    amrex::Finalize();
    return 0;
}

//
// Parse command line arguments
//  
void get_input_arguments ( const int argc, char** argv,
			   string& pltfile, string& slicefile,
			   Vector<string>& var, int& format, bool& lleft,
			   int& dir, Real& tol )
{
    
    int i = 1; // skip program name

    while ( i < argc )
    {

	if ( !strcmp(argv[i],"-p") || !strcmp(argv[i],"--pltfile") )
	{
	    pltfile = argv[++i];
	}
	else if ( !strcmp(argv[i],"-s") || !strcmp(argv[i],"--slicefile") )
	{
	    slicefile = argv[++i];
	}
	else if ( !strcmp(argv[i],"-d") || !strcmp(argv[i],"--direction") )
	{
	    dir = atoi(argv[++i]);

	    if ( ( dir < 1 ) || ( dir > 3) )
	    {
		cout << "\nDirection parameter has value " << dir << endl;
		cout << "Valid values are 1, 2, 3. " << endl;
		help ();
		exit (EXIT_FAILURE);
	    }
	    
	}
	else if ( !strcmp(argv[i],"-v") || !strcmp(argv[i],"--variable") )
	{
	    string tmp = argv[++i];

	    while ( !tmp.empty() )
	    {
		int idx = tmp.find (",");

		if ( idx == string::npos )
		    idx = tmp.size() ;

		var.push_back ( tmp.substr(0,idx) );
		tmp.erase (0, idx+1);
	    }
	}
	else if ( !strcmp(argv[i],"-l") || !strcmp(argv[i],"--lower_left") )
	{
	    lleft = true; 
	}
	else if ( !strcmp(argv[i],"-f") || !strcmp(argv[i],"--format") )
	{
	    format = atoi(argv[++i]);
	}
	else if ( !strcmp(argv[i],"-t") || !strcmp(argv[i],"--tolerance") )
	{
	    tol = atof(argv[++i]);
	}
	else
	{
	    std::cout << "\n\nOption " << argv[i] << " not recognized" << std::endl;
	    help ();
	    exit ( EXIT_FAILURE );
	}

	// Go to the next parameter name
	++i;
    }

    if ( pltfile.empty () )
    {
	std::cout << "\n\nRequired option [-p|--pltfile] is missing" << std::endl;
    	help ();
    	exit (EXIT_FAILURE);
    }

    // Give default name to slicefile
    if ( slicefile.empty() )
    {
	string tmp = pltfile; 
	
	// Remove final / is present
	//  ( plotfile name are directories!)
	int idx = tmp.rfind("/");

	if ( idx == (tmp.length ()-1) )
	    tmp.erase(idx);

	idx = tmp.rfind("/");
	
	if ( idx == string::npos )
	{
	    slicefile = tmp + ".slice";   
	}
	else
	{
	    slicefile = tmp.substr(idx+1) + ".slice";
	}

	cout << "Slicefile set to " + slicefile << endl;
    }

    if ( var.size() == 0)
    {
	std::cout << "\n\nNeed to specify at least one variable" << std::endl;
    	help ();
    	exit (EXIT_FAILURE);
    }


    // Print summary
    cout << "\nPltfile is " + pltfile
	 << "\nSlicefile is set to " + slicefile
	 << "\nDirection to extract is " << dir
	 << "\nVariable(s) to extract is(are): ";

    for (int nv = 0; nv < var.size(); ++nv)
	cout << var[nv] + " ";
    
    cout<< "\nPrint lower left = " << lleft 
	<< "\nFormat is " << format
	<< "\nTolerance is " << tol << endl;
    
}



//
// Print usage info
// 
void help ()
{
    std::cout << "\n\nExtract at 1D slice through a plotfile in any coordinate direction."
	      << "\nWorks with 1-, 2-, or 3-d datasets."
	      << "\n "
	      << "\nUsage:"
	      << "\n   fextract [-p plotfile] [-s outfile] [-d dir] [-v variable] plotfile"
	      << "\n "
	      << "\nargs [-p|--pltfile]   plotfile   : plot file directory (depreciated, optional)"
	      << "\n     [-s|--slicefile] slice file : slice file          (optional)"
	      << "\n     [-d|--direction] idir       : slice direction {1 (default), 2, or 3}"
	      << "\n     [-v|--variable]  varname(s) : only output the values of variable varname"
	      << "\n                                  (comma separated string for multiple variables)"
	      << "\n     [-l|--lower_left]           : slice through lower left corner instead of center"
	      << "\n     [-f|--format]               : Output format, Number of sig-figs to write {10 (default)} "
	      << "\n     [-t|--tolerance]            : All values with smaller abs than tolerance are printed as 0"
	      << "\n "
	      << "\nNote the plotfile at the end of the commandline is only required if you do"
	      << "\nnot use the depreciated '-p' option"
	      << "\n "
	      << "\nIf a job_info file is present in the plotfile, that information is made"
	      << "\navailable at the end of the slice file (commented out), for reference.\n\n";
       
}
    
