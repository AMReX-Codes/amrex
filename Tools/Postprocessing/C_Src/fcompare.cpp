//
// C++ version of fcompare.f90
//
// Take 2 plotfiles as input and compare each level zone by zone for differences.
// The grids at each level may not be identical, but must have the same
// problem domain.
//
// This does not rely on FBaselib.
//
//  Author: Michele Rosso
//  Date  : May, 2018
//
#include "AMReX_DataServices.H"

using namespace amrex;
using namespace std;

//
// Prototypes
//
void GetInputArgs ( const int argc, char** argv,
		    bool& ghost, int& norm,
		    string& diffvar, string& zonevar,
		    string& file1, string& file2 );

void PrintHelp ( const char* progName );


//
//
//
int main ( int argc, char* argv[] )
{

    // Use "false" as third arguments to prevent AMReX
    // to check for an input file as first command line argument
    amrex::Initialize ( argc, argv, false );
    {

    // This runs in serial for now
    if ( ParallelDescriptor::NProcs() > 1 )
    	Abort("Only run this with one processor");

    // Input arguments
    bool    do_ghost=false;
    int     norm=0;
    string  diffVar;
    string  zoneVar;
    string  file1;
    string  file2;


    GetInputArgs ( argc, argv, do_ghost, norm, diffVar, zoneVar,
		   file1, file2 );

    // Start dataservices (no clue why we need to do this )
    DataServices::SetBatchMode();

    // Define the type of file
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices1 (file1, fileType);
    DataServices dataServices2 (file2, fileType);

    // NOTE:
    // fcompare.f90 checks that both files have the same dimensionality
    // and the dimensionality s the same as AMREX_SPACEDIM used to build fcompare itself.
    // fcompare.cpp does not need to do this since dataservices take care
    // of checking that when loading the files.
    if ( !dataServices1.AmrDataOk() || !dataServices2.AmrDataOk() )
	DataServices::Dispatch(DataServices::ExitRequest, NULL);

    // get data from plot files
    AmrData& data1 = dataServices1.AmrDataRef();
    AmrData& data2 = dataServices2.AmrDataRef();

    // Check if both plotfiles have the same numbers of levels
    if ( data1.FinestLevel () != data2.FinestLevel () )
	Abort ("number of levels do not match");

    int finestLevel = data1.FinestLevel ();

    // Check that the problem domain (bounding box) are the same size
    // check that boxarray and distribution maps are the same
    // for both plotfiles
    const Vector<Box>& bb1 = data1.ProbDomain ();
    const Vector<Box>& bb2 = data2.ProbDomain ();

    for (int nl = 0; nl <= finestLevel; ++nl )
    {
	if ( bb1[nl] != bb2[nl] )
    	    Abort ("grids do no match");

	const BoxArray& ba1 = data1.boxArray(nl);
	const BoxArray& ba2 = data2.boxArray(nl);
	const DistributionMapping& dm1 = data1.DistributionMap (nl);
	const DistributionMapping& dm2 = data2.DistributionMap (nl);

	AMREX_ALWAYS_ASSERT ( ba1 == ba2 );
	AMREX_ALWAYS_ASSERT ( dm1 == dm2 );
    }

    // Check if plotfiles have same numbers of variables
    if ( data1.NComp () != data2.NComp () )
	Warning ("WARNING: number of variables do not match");

    // Check grid spacing is the same for both plotfiles
    const Vector<Vector<Real> >& dxLevel1 = data1.DxLevel();
    const Vector<Vector<Real> >& dxLevel2 = data2.DxLevel();

    for (int nl = 0; nl <= finestLevel; ++nl )
	for (int nd = 0; nd < AMREX_SPACEDIM; ++nd )
	    if ( dxLevel1[nl][nd] != dxLevel2[nl][nd] )
		Abort ("grid dx does not match");

    // Check for the same number of ghost cells
    int ng = 0;

    if ( data1.NGrow() != data2.NGrow() )
    {
	Warning ("grids have different numbers of ghost cells");
    }
    else
    {
	if (do_ghost)
	    ng = data1.NGrow();
    }

    // Check if variables in file 1 are also in file 2
    // and print a warning if this is not the case
    const Vector<string>& varNames1 = data1.PlotVarNames ();
    const Vector<string>& varNames2 = data2.PlotVarNames ();

    bool         allVarsFound = true;
    Vector<int>  varMap(data1.NComp(),-1);

    for (int nv2 = 0; nv2 < data2.NComp (); ++nv2 )
    {
	bool varFound = false;

	for (int nv1 = 0; nv1 < data1.NComp (); ++nv1 )
	{
	    if ( varNames2[nv2] == varNames1[nv1] )
	    {
		varMap[nv1] = nv2;
		varFound = true;
		break;
	    }
	}


	if ( !varFound )
	{
	    string warnmsg = "WARNING: variable " + varNames2[nv2]
		+ " not found in plotfile " + file1;

	    Warning (warnmsg);

	    allVarsFound = false;
	}

    }

    //
    // Go level-by-level and patch-by-patch and compare the data
    //
    cout << endl;

    cout << setw(25) << " variable name  "
    	 << setw(25) << " absolute error "
    	 << setw(25) << " relative error "
    	 << endl;

    cout << setw(25) << " "
    	 << setw(22) << "(||A - B||)"
    	 << setw(29) << "(||A - B||/||A||)"
    	 << endl;

    cout << string(77, '-') << std::endl;

    vector<Real>  aerror(data1.NComp());
    vector<Real>  rerror(data1.NComp());
    vector<Real>  rerror_den(data1.NComp());
    vector<bool>  has_nan(data1.NComp());
    Real          globalError = 0.0;
    bool          anyNaNs = false;

    // Loop over levels
    for (int nl = 0; nl <= finestLevel; ++nl )
    {

	// Zero arrays
	fill ( aerror.begin(),     aerror.end(), 0.0 );
	fill ( rerror.begin(),     rerror.end(), 0.0 );
	fill ( rerror_den.begin(), rerror_den.end(), 0.0 );
	fill ( has_nan.begin(),   has_nan.end(), false );

	// We already checked that boxarray is the same for
	// both plotfiles
	const BoxArray& ba = data1.boxArray(nl);

	// Loop over boxes on same level
    	for (int nb = 0; nb < ba.size (); ++nb )
    	{
	    // Create a Multifab based on a boxarray with only
	    // one box.
	    // This way we only have in memory a box at a time, but we can
	    // still use the capabilities of MultiFab  class
	    BoxArray             sba(ba[nb]); // Single Box Array
	    DistributionMapping  dm(sba);
	    MultiFab             mf1( sba, dm, 1, ng );
	    MultiFab             mf2( sba, dm, 1, ng );
	    MultiFab             mfdiff( sba, dm, 1, ng );

            //
	    // In doing the calculation below, assume that the
	    // data are always cell-centered
	    //
	    for (int nv = 0; nv < data1.NComp(); ++nv)
	    {
		if ( varMap[nv] == -1 )
		    continue;

		data1.FillVar ( mf1, nl, varNames1[nv] );
		data2.FillVar ( mf2, nl, varNames2[varMap[nv]] );

		// check for NaNs -- comparisons don't work when they are present
		// note: regardless of do_ghost, this will check the ghostcells
		// too if they are present.
		if ( mf1.contains_nan() || mf2.contains_nan() )
		{
		    has_nan[nv] = true;
		    anyNaNs = true;
		    continue;
		}

		MultiFab::LinComb ( mfdiff, 1.0, mf1, 0, -1.0, mf2, 0, 0, 1, ng );

		if ( norm == 0 )
		{
		    aerror[nv]       = amrex::max ( mfdiff.norm0(0,ng), aerror[nv] );
		    rerror_den[nv]   = amrex::max ( mf1.norm0(0,ng),    rerror_den[nv] );
		}
		else if ( norm == 1 )
		{
		    aerror[nv]       = aerror[nv] + mfdiff.norm1(0,ng);
		    rerror_den[nv]   = rerror_den[nv] + mf1.norm1();
		}
		else if ( norm == 2 )
		{
		    aerror[nv]       = aerror[nv] + MultiFab::Dot ( mfdiff, 0, mfdiff, 0, 1, ng );
		    rerror_den[nv]   = rerror_den[nv] + MultiFab::Dot ( mf1, 0, mf1, 0, 1, ng );
		}

		rerror[nv] = aerror[nv];

    	    }

	    data1.FlushGrids();
	    data2.FlushGrids();
    	}

	// Compute normalized norms for this level and print out value
	// What  if nghost are included?
	// In this case the normalization by grid volume is not correct
	if ( norm > 0 )
	{
	    // Metric coefficient
	    Real dvol = 1.0;

	    for (int nd = 0; nd < AMREX_SPACEDIM; ++nd )
		dvol = dvol * dxLevel1[nl][nd];


	    for (int nv = 0; nv < data1.NComp(); ++nv )
	    {
		aerror[nv]  = pow ( aerror[nv]*dvol, 1.0/norm );
		rerror[nv]  = aerror[nv] / pow ( rerror_den[nv]*dvol, 1.0/norm );
	    }
	}
	else
	{
	    for (int nv = 0; nv < data1.NComp(); ++nv )
		rerror[nv]  = aerror[nv] / rerror_den[nv];
	}

	//
	// print out the comparison report for this level
	//
	Print() << " level = " << setw(2) << nl << endl;

	for (int nv = 0; nv < data1.NComp(); ++nv )
	{
	   if ( varMap[nv] == -1 )
	   {
	       Print() << " " << setw(22) << varNames1[nv] << "  "
		       << setw(50) << "< variable not present in both files > "
		       << endl;
	   }
	   else if ( has_nan[nv] == true )
	   {
	       Print() << " " << setw(24) << varNames1[nv] << "  "
	   	       << setw(50) << "< NaN present > "
	   	       << endl;
	   }
	   else
	   {
	       Real aerr = 0.0;
	       Real rerr = 0.0;

	       if ( aerror[nv] > 0.0 )
		   aerr = amrex::min ( amrex::max ( aerror[nv], 1e-99 ), 1e+98 );

	       if ( rerror[nv] > 0.0 )
		   rerr = amrex::min ( amrex::max ( rerror[nv], 1e-99 ), 1e+98 );

	       Print() << scientific;
	       Print() << " " << setw(22) << varNames1[nv] << "  "
		       << setw(24) << setprecision(10)  << uppercase << aerr
		       << "  " << setw(24) << setprecision(10)  << rerr
		       << endl;
	   }


	}

	// compute global error
	for ( unsigned int nv = 0; nv < aerror.size(); ++nv )
	    globalError = amrex::max (globalError, aerror[nv] );

    }

    if ( globalError == 0.0 && !anyNaNs)
    {
	if ( !allVarsFound )
	    Warning ("WARNING: not all variables present in both files");

	Print() << "PLOTFILES AGREE" << endl;

	exit (EXIT_SUCCESS);
    }
    else
    {
	exit (EXIT_FAILURE);
    }

    }
    Finalize ();
}



//
// Parse command line arguments
//
void GetInputArgs ( const int argc, char** argv,
		    bool& ghost, int& norm,
		    string& diffvar, string& zonevar,
		    string& file1, string& file2 )
{

    int i = 1; // skip program name

    while ( i < argc - 2 ) // Skip last two: those are the files names
    {

	if ( !strcmp(argv[i],"-g") || !strcmp(argv[i],"--ghost") )
	{
	    ghost = true;
	}
	else if ( !strcmp(argv[i],"-n") || !strcmp(argv[i],"--norm") )
	{
	    norm = atoi(argv[++i]);

	    if ( ( norm < 0 ) || ( norm > 2) )
	    {
		cout << "\nNorm type is " << norm << endl;
		cout << "Valid values are 0, 1, 2. " << endl;
		PrintHelp ( argv[0] );
		exit (EXIT_FAILURE);
	    }

	}
	else if ( !strcmp(argv[i],"-d") || !strcmp(argv[i],"--diffvar") )
	{
	    Warning ("Option -d/--diffvar is not supported");
	    //diffvar = argv[++i];
	}
	else if ( !strcmp(argv[i],"-z") || !strcmp(argv[i],"--zone_info") )
	{
	    Warning ("Option -z/--zone_info is not supported");
	    //zonevar = argv[++i];
	}
	else
	{
	    std::cout << "\n\nOption " << argv[i] << " not recognized" << std::endl;
	    PrintHelp ( argv[0] );
	    exit ( EXIT_FAILURE );
	}

	// Go to the next parameter name
	++i;
    }

    // get name of input files
    if (argc > 2)
    {
	file1 = argv[argc-2];
	file2 = argv[argc-1];
    }

    if ( file1.empty() || file2.empty() )
    {
	PrintHelp (argv[0]);
	Abort("Missing one or both input files");
    }

    // Print summary
    cout << "\nFile 1 is " + file1
	 << "\nFile 2 is " + file2;

    if (ghost)
	cout << "\nGhost cells will be compared (if stored)";

    if (!diffvar.empty())
	cout << "\nPlotfile with differences for variable " + diffvar
	     << " will be output";

    if (!zonevar.empty())
	cout << "\nInfo on max error zone for var " + zonevar
	     << " will be given";

    cout << endl << endl;
}



//
// Print usage info
//
void PrintHelp (const char* progName)
{
    cout  << "\nCompare two plotfiles, zone by zone, to machine precision"
	  << "\nand report the maximum absolute and relative errors for each"
	  << "\nvariable."
	  << "\n "
	  << "\nUsage:"
	  << "\n " << progName << " [-g|--ghost] [-n|--norm num] file1 file2"
	// << "\n " << progName << " [-g|--ghost] [-n|--norm num] [--diffvar var] [-z|--zone_info var] file1 file2"
	  << "\n "
	  << "\nOptional arguments:"
	  << "\n   -g|--ghost         : compare the ghost cells too (if stored)"
	  << "\n   -n|--norm num      : what norm to use (default is 0 for inf norm)"
	  // << "\n   -d|--diffvar var   : output a plotfile showing the differences for "
	  // << "\n                        variable var"
	  // << "\n   -z|--zone_info var : output the information for a zone corresponding"
	  // << "\n                        to the maximum error for the given variable"
	  << "\n\n" << endl;

}
