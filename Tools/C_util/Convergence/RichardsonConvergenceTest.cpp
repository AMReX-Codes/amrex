
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
using std::ios;

#include <unistd.h>

#include <WritePlotFile.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_AVGDOWN_F.H>
#include "AMReX_ArrayLim.H"
#include "DebugDump.H"
#include <iomanip>

#ifdef AMREX_DEBUG
#include <TV_TempWrite.H>
#endif

#define GARBAGE 666.e+40

using namespace amrex;
using std::setw;
using std::setprecision;

static
void
PrintUsage (const char* progName)
{
  std::cout << "This utility performs a convergence test among three" << std::endl
            << "plotfiles which have the same geometrical domain"        << std::endl
            << "but a factor of refinement between"                      << std::endl
            << "the cells from each plotfile at the same level."         << std::endl
            << "For instance, it works to diff two plotfiles having"     << std::endl
            << "  Plotfile 1: 25x25  base grid, Ref_Ratio = 2"            << std::endl
            << "  Plotfile 2: 50x50  base grid, Ref_Ratio = 2"            << std::endl
            << "  Plotfile 3:100x100 base grid, Ref_Ratio = 2"            << std::endl
            << "In both cases, the geometrical region which is refined"  << std::endl
            << "should be the same.  So, this is generally good for"       << std::endl
            << "comparing cases ran using a fixed grid file."            << std::endl;
  std::cout << '\n';
  std::cout << "Usage:" << '\n';
  std::cout << progName << '\n';
  std::cout << "    coarFile= coarse file name  (input)" << '\n';
  std::cout << "    mediFile= medium plot file  (input)" << '\n';
  std::cout << "    fineFile= finest plot file  (input)" << '\n';
  std::cout << "    mediError = medium error file (optional output)" << '\n';
  std::cout << "    coarError = coarse error file (optional output)" << '\n';
  std::cout << "   [-help]" << '\n';
  std::cout << "   [-verbose]" << '\n';
  std::cout << '\n';
  exit(1);
}
/**********/
bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2)
{
  const Vector<std::string>& derives1 = amrd1.PlotVarNames();
  const Vector<std::string>& derives2 = amrd2.PlotVarNames();
  int length = derives1.size();
  if (length != derives2.size())
    return false;
  for (int i=0; i<length; ++i)
    if (derives1[i] != derives2[i])
      return false;
  return true;
}
/**********/
Vector<string>
getVarNames(const AmrData& amrd1,
            const AmrData& amrd2)
{
  Vector<string> names;
  const Vector<std::string>& derives1 = amrd1.PlotVarNames();
  const Vector<std::string>& derives2 = amrd2.PlotVarNames();
  int length = derives1.size();
  if (length == derives2.size())
  {
    names = derives1;
  }
  else
  {
    int minlen = std::min(derives1.size(), derives2.size());
    names.resize(minlen);
    for(int ilen = 0; ilen < minlen; ilen++)
    {
      names[ilen] = derives1[ilen];
    }
  }

  return names;
}
/**********/
IntVect
getRefRatio(const Box& crse,
            const Box& fine)
{
  // Compute refinement ratio between crse and fine boxes, return invalid
  // IntVect if there is none suitable
  ParmParse pp("");
  Vector<int> rr_in(BL_SPACEDIM,-1);

  int Nrr = pp.countval("ref_ratio");

  BL_ASSERT(Nrr==0 || Nrr==BL_SPACEDIM || Nrr==1);
  if (Nrr>0)
  {
    pp.queryarr("ref_ratio",rr_in,0,Nrr);
    return IntVect(rr_in);
  }

  IntVect ref_ratio;
  for (int i=0; i<BL_SPACEDIM; ++i)
    ref_ratio[i] = fine.length(i) / crse.length(i);

  // Check results
  Box test1 = ::Box(fine).coarsen(ref_ratio);
  Box test2 = ::Box(test1).refine(ref_ratio);
  if (test1 != crse  ||  test2 != fine)
    ref_ratio = IntVect();
  return ref_ratio;
}

/**********/
void
getErrorNorms(Vector<Real>& a_norms, //one for each comp
              Vector<string>& a_names,
              const string& a_fineFile,
              const string& a_coarFile,
              const string& a_errFile,
              const int& a_norm,
              bool verbose)
{
  //
  // Scan the arguments.
  //
  std::string iFile1, iFile2, difFile;
  int norm;
  iFile1  = a_coarFile;
  iFile2  = a_fineFile;
  difFile = a_errFile;
  norm   = a_norm;

  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  DataServices dataServices1(iFile1, fileType);
  DataServices dataServices2(iFile2, fileType);

  if (!dataServices1.AmrDataOk() || !dataServices2.AmrDataOk())
    amrex::Abort("ERROR: Dataservices not OK");

  //
  // Generate AmrData Objects
  //
  AmrData& amrData1 = dataServices1.AmrDataRef();
  AmrData& amrData2 = dataServices2.AmrDataRef();

  //
  // Initial Tests
  //
  if (amrData1.FinestLevel() != amrData2.FinestLevel())
    amrex::Abort("ERROR: Finest level is not the same in the two plotfiles");

  int nComp1      = amrData1.NComp();
  int nComp2      = amrData2.NComp();

  int nComp = std::min(nComp1,nComp2);

  if (!amrDatasHaveSameDerives(amrData1,amrData2))
  {
    std::cout << "Warning: Plotfiles do not have the same state variables" << std::endl;
    std::cout << "         Setting nComp to be the minimum number of variables: " << nComp << std::endl;
  }

  a_names = getVarNames(amrData1, amrData2);

  int finestLevel = amrData1.FinestLevel();

  Vector<int> destComps(nComp);
  for (int i = 0; i < nComp; i++)
    destComps[i] = i;


  //
  // Compute the error
  //
  Vector<MultiFab*> error(finestLevel+1);

  if (ParallelDescriptor::IOProcessor())
    std::cout << "Level  L"<< norm << " norm of Error in Each Component" << std::endl
              << "-----------------------------------------------" << std::endl;

  Vector<Real> sum_norms(nComp);
  for (int iComp = 0; iComp < nComp; iComp++)
    sum_norms[iComp] = 0.0;
  for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
  {
    //
    // Construct level box array and check that the lengths are the same
    //
    const BoxArray& ba1 = amrData1.boxArray(iLevel);
    const BoxArray& ba2 = amrData2.boxArray(iLevel);

    if (ba1.size() != ba2.size())
      std::cout << "Warning: BoxArray lengths are not the same at level " << iLevel << std::endl;

    //
    // Construct refinement ratio, build the coarsened boxarray
    // (amrData2 is the refined plotfile)
    //
    const Box& domain1     = amrData1.ProbDomain()[iLevel];
    const Box& domain2     = amrData2.ProbDomain()[iLevel];
    IntVect refine_ratio   = getRefRatio(domain1, domain2);
    if (refine_ratio == IntVect())
      amrex::Error("Cannot find refinement ratio from data to exact");

    if (verbose)
      std::cerr << "level = " << iLevel << "  Ref_Ratio = " << refine_ratio
                << std::endl;

    BoxArray ba2Coarse(ba2);
    ba2Coarse.coarsen(refine_ratio);

    // Define new_data1 in case the boxarrays are not the same
    DistributionMapping dm2Coarse(ba2Coarse);
    MultiFab new_data1(ba2Coarse,dm2Coarse, 1,0);

    //
    // Construct MultiFab for errors
    //
    error[iLevel] = new MultiFab(ba2Coarse, dm2Coarse, nComp, 0);
    error[iLevel]->setVal(GARBAGE);

    //
    // For each component, average the fine fields down and calculate
    // the errors
    //
    Vector<Real> norms(nComp);
    for (int iComp = 0; iComp < nComp; iComp++)
      norms[iComp] = 0.0;

    for (int iComp = 0; iComp < nComp; ++iComp)
    {
      MultiFab& data1     = amrData1.GetGrids(iLevel, iComp);
      MultiFab& data2Fine = amrData2.GetGrids(iLevel, iComp);

      BoxArray baf;
      if (iLevel<finestLevel)
      {
        baf = BoxArray(amrData1.boxArray(iLevel+1)).coarsen(refine_ratio);
      }


      // Copy the data at the coarse level one component at a time
      new_data1.copy(data1,0,0,1);

      //
      // Calculate the errors  for each FAB in the MultiFab
      //
      for (MFIter mfi(new_data1); mfi.isValid(); ++mfi)
      {
        //
        // Create the Coarsened version of data2
        //
        int index = mfi.index();

        const Box& bx = ba2Coarse[index];
        FArrayBox data2Coarse(bx, 1);
        int ncCoarse = 1;



        FORT_CV_AVGDOWN(data2Coarse.dataPtr(),
                        ARLIM(bx.loVect()), ARLIM(bx.hiVect()),
                        &ncCoarse,
                        data2Fine[mfi].dataPtr(),
                        ARLIM(data2Fine[mfi].loVect()),
                        ARLIM(data2Fine[mfi].hiVect()),
                        bx.loVect(), bx.hiVect(),
                        refine_ratio.getVect());


        //
        // Calculate the errors on this FAB for this component
        //

        (*error[iLevel])[mfi].copy(new_data1[mfi], 0, iComp, 1);
        (*error[iLevel])[mfi].minus(data2Coarse  , 0, iComp, 1);

        if (iLevel<finestLevel)
        {
          std::vector< std::pair<int,Box> > isects = baf.intersections(bx);

          for (int ii = 0; ii < isects.size(); ii++)
            (*error[iLevel])[mfi].setVal(0,isects[ii].second,iComp,1);
        }

        Real grdL2 = (*error[iLevel])[mfi].norm(norm, iComp, 1);

        if (norm != 0)
        {
          norms[iComp] = norms[iComp] + pow(grdL2, norm);
        }
        else
        {
          norms[iComp] = std::max(norms[iComp], grdL2);
        }
      }
    }


    //
    // Output Statistics
    //
    if (ParallelDescriptor::IOProcessor())
      std::cout << "  " << iLevel << "    ";


#ifdef BL_USE_MPI
    MPI_Datatype datatype = ParallelDescriptor::Mpi_typemap<Real>::type();
    if (ParallelDescriptor::IOProcessor())
    {
      Vector<Real> tmp(nComp);
      for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
        if (proc != ParallelDescriptor::IOProcessorNumber())
        {
          MPI_Status stat;
          int rc = MPI_Recv(tmp.dataPtr(), nComp, datatype,
                            MPI_ANY_SOURCE, proc, ParallelDescriptor::Communicator(),
                            &stat);

          if (rc != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);

          for (int iComp = 0; iComp < nComp; iComp++)
            if (norm != 0)
            {
              norms[iComp] = norms[iComp] + tmp[iComp];
            }
            else
            {
              norms[iComp] = std::max(norms[iComp], tmp[iComp]);
            }
        }
    }
    else
    {
      int rc = MPI_Send(norms.dataPtr(), nComp, datatype,
                        ParallelDescriptor::IOProcessorNumber(),
                        ParallelDescriptor::MyProc(),
                        ParallelDescriptor::Communicator());

      if (rc != MPI_SUCCESS)
        ParallelDescriptor::Abort(rc);
    }
#endif


    Real vol = 1.0;
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
      vol *= amrData1.DxLevel()[iLevel][dir];

    if (ParallelDescriptor::IOProcessor())
    {
      for (int iComp = 0; iComp < nComp; iComp++)
      {
        if (norm != 0)
        {
          norms[iComp] = norms[iComp] * vol;
          norms[iComp] = pow(norms[iComp], (1.0/norm));
        }

        if (verbose)
          std::cout << norms[iComp] << " ";

        sum_norms[iComp] = sum_norms[iComp] + norms[iComp];
      }

      if (verbose)
        std::cout << std::endl;
    }
  }

  if (verbose) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "  Sum    ";
      for (int iComp = 0; iComp < nComp; iComp++)
        {
          std::cout << sum_norms[iComp] << " ";
        }
    }
  }

  a_norms= sum_norms;
  //only want to write this once
  if (!difFile.empty() && (a_norm == 0))
  {
    WritePlotFile(error, amrData1, difFile, verbose);
  }

  for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
  {
    delete error[iLevel];
  }
}


int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);

  if (argc == 1)
    PrintUsage(argv[0]);

  ParmParse pp;

  if (pp.contains("help"))
    PrintUsage(argv[0]);

  bool verbose = false;
  if (pp.contains("verbose"))
    verbose = true;

  FArrayBox::setFormat(FABio::FAB_IEEE_32);
  //
  // Scan the arguments.
  //
  std::string coarFile, mediFile, fineFile;
  std::string coarError, mediError;


  pp.query("fineFile", fineFile);
  if (fineFile.empty())
    amrex::Abort("You must specify fineFile");

  pp.query("mediFile", mediFile);
  if (mediFile.empty())
    amrex::Abort("You must specify mediFile");

  pp.query("coarFile", coarFile);
  if (coarFile.empty())
    amrex::Abort("You must specify coarFile");


  pp.query("mediError", mediError);
  pp.query("coarError", coarError);

  //l2 is not supported
  for(int inorm = 0; inorm <=1; inorm++)
  {
    Vector<Real> normsMedi, normsCoar;
    Vector<string> namesMedi, namesCoar;
    getErrorNorms(normsMedi, namesMedi, fineFile, mediFile, mediError, inorm, verbose);
    getErrorNorms(normsCoar, namesCoar, mediFile, coarFile, coarError, inorm, verbose);
    int ncompMedi = normsMedi.size();
    int ncompCoar = normsCoar.size();
    int ncomp = std::min(ncompMedi, ncompCoar);
    amrex::Print() << "\\begin{table}[p]" << std::endl;
    amrex::Print() << "\\begin{center}" << std::endl;
    amrex::Print() << "\\begin{tabular}{|cccc|} \\hline" << std::endl;
    amrex::Print() << "Variable & $e_{4h \\rightarrow 2h}$ & Order & $e_{2h \\rightarrow h}$\\\\" << std::endl;;
    amrex::Print() << "\\hline " << std::endl;

    for (int icomp = 0; icomp < ncomp; icomp++)
    {
      bool printOrder = false;
      Real order;
      Real finenorm = normsMedi[icomp];
      Real coarnorm = normsCoar[icomp];
      if(std::abs(finenorm) > 1.0e-40)
      {
        order = log(std::abs(coarnorm/finenorm))/log(2.0);
        printOrder = true;
      }

      amrex::Print() << namesMedi[icomp] << "&    \t ";
      amrex::Print() << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << normsCoar[icomp]  << " & ";
      if(printOrder)
      {
        amrex::Print()      << setw(12)
                            << setprecision(3)
                            << setiosflags(ios::showpoint)
                            << setiosflags(ios::scientific);
        amrex::Print()  << order << " & ";
      }
      else
      {
        amrex::Print() << "------------ &";
      }

      amrex::Print() << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << normsMedi[icomp];
      amrex::Print() << " \\\\ " << std::endl;

    }
    amrex::Print() << "\\hline " << std::endl;
    amrex::Print() << "\\end{tabular}" << std::endl;
    amrex::Print() << "\\end{center}" << std::endl;
    amrex::Print() << "\\caption{Solution error convergence rates using $L_";
    if(inorm == 0)
    {
      amrex::Print() <<"\\infty$ norm." << std::endl;;
    }
    else
    {
      amrex::Print() << inorm << "$ norm.";
    }
    amrex::Print() << "}" << std::endl;
    amrex::Print() << "\\end{table}" << std::endl;
    amrex::Print() << std::endl << std::endl;
  }


  amrex::Finalize();
  return 0;
}
