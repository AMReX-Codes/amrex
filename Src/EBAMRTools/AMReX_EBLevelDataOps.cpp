#include <cmath>
#include <iomanip>
#include "AMReX_SPMD.H"

#include "AMReX_EBLevelDataOps.H"
#include "AMReX_FaceIterator.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_DistributionMapping.H"
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include "AMReX_PlotFileUtil.H"
#include <cstdlib>
namespace amrex
{
  void 
  EBLevelDataOps::
  viewEBLevel(const FabArray<EBCellFAB> * a_data, const EBLevelGrid* a_eblg)
  {
    Vector<string> names(a_data->nComp());
    for(int icomp = 0; icomp < a_data->nComp(); icomp++)
    {
      names[icomp] = string("eb_var_") + EBArith::convertInt(icomp);
    }
    string filename("debug_file.plt");
    writeSingleLevelEBPlotFile(filename, *a_data, *a_eblg, names);

    string command = "visit -o " + filename + string("/Header");
    int ret = std::system(command.c_str());
    amrex::Print() << "data output to " << filename << ".  Visit was called and got return value " << ret << endl;
  }
  ///writes plotfile that visit can eat.   Just single-valued stuff
  void 
  EBLevelDataOps::
  writeSingleLevelEBPlotFile(const std::string         & a_filename,
                             const FabArray<EBCellFAB> & a_data,
                             const EBLevelGrid         & a_eblg,
                             const Vector<string>      & a_varNames)
  {
    MultiFab mfdata;
    makePlotMultiFab(mfdata, a_data, a_eblg);
    Vector<string> varNameArr(a_varNames.size());
    Vector<string>& varNameArrCast = static_cast<Vector<string>& >(varNameArr);
    varNameArrCast = a_varNames;
    varNameArr.push_back(string("vfrac"));

    Geometry geom(a_eblg.getDomain());
    Real time = 0; int level_step = 0;  

    WriteSingleLevelPlotfile(a_filename, mfdata, varNameArr, geom, time, level_step);
  }

  ///writes plotfile that visit can eat.   Just single-valued stuff
  void 
  EBLevelDataOps::
  writeEBAMRPlotFile(const std::string                   & a_filename,
                     const Vector<FabArray<EBCellFAB>* > & a_data,
                     const Vector<EBLevelGrid>           & a_eblg,
                     const Vector<int>                   & a_refRat,
                     const Vector<string>                & a_varNames)
  {
    int nlevels = a_data.size();
    Vector<IntVect>   refRatArr(nlevels, 2*IntVect::Unit);
    for(int ilev = 0; ilev < a_refRat.size(); ilev++)
    {
      refRatArr[ilev] = a_refRat[ilev]*IntVect::Unit;
    }


    Vector<string>   varNameArr(a_varNames.size());
    Vector<string>& varNameArrCast = static_cast<Vector<string>& >(varNameArr);
    varNameArrCast = a_varNames;
    varNameArr.push_back(string("vfrac"));

    Vector<const MultiFab*> mfdata(nlevels);
    Vector<Geometry> geom(nlevels);
    for(int ilev = 0; ilev < nlevels; ilev++)
    {
      geom[ilev].define(a_eblg[ilev].getDomain());
      mfdata[ilev] = new MultiFab();
      MultiFab* castfab = const_cast<MultiFab*>(mfdata[ilev]);
      makePlotMultiFab(*castfab, *a_data[ilev], a_eblg[ilev]);
    }

    Real time = 0; Vector<int> level_step(nlevels, 0);  
    WriteMultiLevelPlotfile(a_filename, nlevels, mfdata, varNameArr, geom, time, level_step, refRatArr);

    for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete mfdata[ilev];
    }
  }

  /// tacks on volume fraction as the last variable
  void 
  EBLevelDataOps::
  makePlotMultiFab(MultiFab                              & a_mfdata,
                   const FabArray<EBCellFAB>             & a_ebdata,
                   const EBLevelGrid                     & a_eblg)

  {
    DistributionMapping dm = a_ebdata.DistributionMap();
    BoxArray            ba = a_ebdata.boxArray();
    int ngrow = 0; //simplifies setting vfrac
    int ncomp = a_ebdata.nComp() + 1 ; // + 1 for vol fraction
    a_mfdata.define(ba, dm, ncomp, ngrow, MFInfo(), FArrayBoxFactory());
    
    for(MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
      BaseFab<Real>       & mfbf = static_cast<BaseFab<Real>&>(a_mfdata[mfi]);
      const BaseFab<Real> & ebbf = a_ebdata[mfi].getSingleValuedFAB();
      int isrc = 0; int idst = 0; int inco = ncomp-1; //last component gets vfrac.
      mfbf.copy(ebbf, isrc, idst, inco);
      Box valid = ba[mfi];
      EBISBox ebis = a_eblg.getEBISL()[mfi];
      for(BoxIterator bit(valid); bit.ok(); ++bit)
      {
        Real vfrac;
        if(ebis.isRegular(bit()))
        {
          vfrac = 1;
        }
        else if(ebis.isCovered(bit()))
        {
          vfrac = 0;
        }
        else
        {
          VolIndex vof(bit(), 0);
          vfrac = ebis.volFrac(vof);
        }
        mfbf(bit(), ncomp-1) = vfrac;
      }
    }
  }

  //-----------------------------------------------------------------------
  Real 
  EBLevelDataOps::
  kappaDotProduct(Real&                       a_volume,
                  const FabArray<EBCellFAB> & a_data1,
                  const FabArray<EBCellFAB> & a_data2,
                  const EBLevelGrid         & a_eblg)
  {
    BL_PROFILE("EBLevelDataOps::kappaDotProduct");
    Real accum = 0.0;

    a_volume = 0.0;

    for (MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      const Box   & box = a_eblg.getDBL()[mfi];
      EBISBox ebisbox = a_eblg.getEBISL()[mfi];
      IntVectSet ivs(box);

      Real cur      = 0;
      Real curVolume= 0;
      
      for(VoFIterator vofit(ivs, ebisbox.getEBGraph()); vofit.ok(); ++vofit)
      {
        Real kappa = ebisbox.volFrac(vofit());
        curVolume += kappa;
        for(int icomp = 0; icomp < a_data1.nComp(); icomp++)
        {
          Real value1 = a_data1[mfi](vofit(), icomp);
          Real value2 = a_data2[mfi](vofit(), icomp);
          cur += value1*value2*kappa*kappa; //makes no sense to me either but that is how it is coded.
        }

        a_volume += curVolume;
        accum += cur;
      }

    }
    gatherBroadCast(accum, a_volume, 1);

    if (a_volume > 0.0)
    {
      accum = accum / a_volume;
    }

    return accum;
  }
  //-----------------------------------------------------------------------

  void 
  EBLevelDataOps::
  aliasIntoMF(shared_ptr<MultiFab>      & a_regData,
              const FabArray<EBCellFAB> & a_ebcData,
              const EBLevelGrid         & a_eblg)
  {
    int ngrow               = a_ebcData.nGrow();
    int nvar                = a_ebcData.nComp();
    BoxArray ba             = a_eblg.getDBL();
    DistributionMapping  dm = a_eblg.getDM();
    Vector<Real*> ptrs(a_ebcData.local_size(), NULL);
    for(MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
      const BaseFab<Real>& bf = a_ebcData[mfi].getSingleValuedFAB();
      int ibox = mfi.LocalIndex();
      ptrs[ibox] = (Real*)(bf.dataPtr(0));
    }
    a_regData = shared_ptr<MultiFab>(new MultiFab(ba, dm, nvar, ngrow, ptrs));
  }

  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  aliasIntoMF(shared_ptr<MultiFab>      & a_regData,
              const FabArray<EBFluxFAB> & a_ebfData,
              const int                 & a_faceDir,
              const EBLevelGrid         & a_eblg)
  {
    int ngrow               = a_ebfData.nGrow();
    int nvar                = a_ebfData.nComp();
    BoxArray ba             = a_eblg.getDBL();
    IntVect indexType = IntVect::Zero;
    indexType[a_faceDir] = 1;
    BoxArray baFace = convert(ba, indexType);
    DistributionMapping  dm = a_eblg.getDM();
    Vector<Real*> ptrs(a_ebfData.local_size(), NULL);
    for(MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
      const BaseFab<Real>& bf = a_ebfData[mfi][a_faceDir].getSingleValuedFAB();
      int ibox = mfi.LocalIndex();
      ptrs[ibox] = (Real*)(bf.dataPtr(0));
    }
    a_regData = shared_ptr<MultiFab>(new MultiFab(baFace, dm, nvar, ngrow, ptrs));
  }
  //-----------------------------------------------------------------------
  void
  EBLevelDataOps::
  checkData(const FabArray<EBCellFAB>&a_data, const string& label)
  {
    barrier();
    amrex::Print() << "==== checking " << label << " for nans and infs =====" << "\n";
    checkForBogusNumbers(a_data);
    barrier();
    amrex::Print() << "Getting maxs and mins for " << label << "\n";
    for(int icomp = 0; icomp < a_data.nComp(); icomp++)
    {
      Real maxVal, minVal;
      getMaxMin(maxVal,minVal,a_data,icomp);
      amrex::Print() << "max value = " << maxVal << " for comp = " << icomp << "\n";
      amrex::Print() << "min value = " << minVal << " for comp = " << icomp << "\n";
    }
    amrex::Print() << "==================================================================== " << "\n";
    barrier();
  }
  //-----------------------------------------------------------------------
  Real 
  EBLevelDataOps::
  parallelSum(const Real& a_value)
  {
    // Find the sum of all a_value's
    Real sum = a_value;
#ifdef BL_USE_MPI
    Real sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &sum, 1, ParallelDescriptor::Mpi_typemap<Real>::type(),
                               MPI_SUM, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelSum");
    }
#endif
    return sum;
  }
  //-----------------------------------------------------------------------
  int  
  EBLevelDataOps::
  parallelSum(const int& a_value)
  {
    // Find the sum of all a_value's
    int sum = a_value;
#ifdef BL_USE_MPI
    int sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &sum, 1, MPI_INT,MPI_SUM, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelSum");
    }
#endif
    return sum;
  }
  //-----------------------------------------------------------------------
  long long 
  EBLevelDataOps::
  parallelSum(const long long& a_value)
  {
    // Find the sum of all a_value's
    long sum = a_value;
#ifdef BL_USE_MPI
    long long sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &sum, 1, MPI_LONG_LONG,MPI_SUM, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelSum");
    }
#endif
    return sum;
  }
  //-----------------------------------------------------------------------
  int 
  EBLevelDataOps::
  parallelMin(const int& a_value)
  {
    // Find the minimum of a_value's
    int val = a_value;
#ifdef BL_USE_MPI
    int sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &val, 1, MPI_INT,MPI_MIN, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMin");
    }
#endif
    return val;
  }
  //-----------------------------------------------------------------------
  int  
  EBLevelDataOps::
  parallelMax(const int& a_value)
  {
    // Find the maximum of a_value's
    int val = a_value;
#ifdef BL_USE_MPI
    int sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &val, 1, MPI_INT,MPI_MAX, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMax");
    }
#endif
    return val;
  }
  //-----------------------------------------------------------------------
  Real 
  EBLevelDataOps::
  parallelMin(const Real& a_value)
  {
    // Find the minimum of a_value's
    Real val = a_value;
#ifdef BL_USE_MPI
    Real sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &val, 1, ParallelDescriptor::Mpi_typemap<Real>::type(),
                               MPI_MIN, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMin");
    }
#endif
    return val;
  }
  //-----------------------------------------------------------------------
  Real 
  EBLevelDataOps::
  parallelMax(const Real& a_value)
  {
    // Find the maximum of a_value's
    Real val = a_value;
#ifdef BL_USE_MPI
    Real sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &val, 1, ParallelDescriptor::Mpi_typemap<Real>::type(),
                               MPI_MAX, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMax");
    }
#endif
    return val;
  }


  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  setCoveredVal(FabArray<EBCellFAB>&a_data,
                const int&           a_comp,
                const Real&          a_value)
  {
    BL_PROFILE("EBLevelDataOps::setCoveredVal(cell,comp)");
    for(MFIter mfi(a_data); mfi.isValid(); ++mfi)
    {
      a_data[mfi].setCoveredCellVal(a_value,a_comp);
    }
  }
  //-----------------------------------------------------------------------
  bool 
  EBLevelDataOps::
  checkForBogusNumbers(const FabArray<EBCellFAB>&a_data)
  {
    //this function checks for nans and infs
    bool dataIsNANINF = false;
    int ncomp = a_data.nComp();
    VolIndex badvof;
    Real badval;
    int badcomp;
    BoxArray dbl = a_data.boxArray();
    for(MFIter mfi(a_data); mfi.isValid(); ++mfi)
    {
      const EBCellFAB& dataEBFAB = a_data[mfi];
      const Box& region = dbl[mfi];
      IntVectSet ivsBox(region);
      const EBISBox& ebisBox = dataEBFAB.getEBISBox();
      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        for (int icomp = 0; icomp < ncomp;++icomp)
        {
          Real val = dataEBFAB(vof,icomp);
          if (std::isnan(val) || std::isinf(val) || std::abs(val)>1.e16)
          {
            badvof = vof;
            badval = val;
            badcomp = icomp;
            //                  amrex::Print() << "      icomp = " << icomp << " vof = " << vof << " val = " << val << "\n";
            dataIsNANINF = true;
            //                  amrex::Error();
          }
        }
      }
    }
    if (dataIsNANINF)
    {
      amrex::Print() << "first bad vof cell = "  << badvof.gridIndex() << ", cellindex = " << badvof.cellIndex() << "\n";
      amrex::Print() << "bad val = "  << badval << " at comp " << badcomp << "\n";
      amrex::Warning("Found a NaN or Infinity.");
    }
    return dataIsNANINF;
  }
  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  getMaxMin(Real&                       a_maxVal,
            Real&                       a_minVal,
            const FabArray<EBCellFAB>& a_data,
            const int&                  a_comp,
            const bool&                 a_doAbs)
  {
    BL_PROFILE("EBLevelDataOps::getMaxMin");
    //this function gets the max and min (valid) values
    a_maxVal = -1.e99;
    a_minVal =  1.e99;

    BoxArray dbl = a_data.boxArray();
    for(MFIter mfi(a_data); mfi.isValid(); ++mfi)
    {
      const EBCellFAB& dataEBFAB = a_data[mfi];
      const Box&        dblBox   = dbl[mfi];
      const IntVectSet ivsBox(dblBox);
      const EBISBox& ebisBox = dataEBFAB.getEBISBox();
      VoFIterator vofit(ivsBox, ebisBox.getEBGraph());
      for (vofit.reset();vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        const Real& val = dataEBFAB(vof,a_comp);
        if (a_doAbs)
        {
          a_maxVal = std::max(a_maxVal,std::abs(val));
          a_minVal = std::min(a_minVal,std::abs(val));
        }
        else
        {
          a_maxVal = std::max(a_maxVal,val);
          a_minVal = std::min(a_minVal,val);
        }
      }
    }
    a_minVal = EBLevelDataOps::parallelMin(a_minVal);
    a_maxVal = EBLevelDataOps::parallelMax(a_maxVal);
  }


  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  setVal(FabArray<EBCellFAB> & a_result,
         const Real          & a_value,
         int                   a_comp)
  {
    for(MFIter mfi(a_result); mfi.isValid(); ++mfi)
    {
      if(a_comp >= 0)
      {
        a_result[mfi].setVal(a_comp,a_value);
      }
      else
      {
        a_result[mfi].setVal(a_value);
      }
    }
  }


  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  setVal(FabArray<EBFluxFAB> & a_result,
         const Real          & a_value)
  {
    for(MFIter mfi(a_result); mfi.isValid(); ++mfi)
    {
      a_result[mfi].setVal(a_value);
    }
  }


  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  axby( FabArray<EBCellFAB>&       a_lhs,
        const FabArray<EBCellFAB>& a_x,
        const FabArray<EBCellFAB>& a_y,
        const Real& a_a,
        const Real& a_b)
  {
    for(MFIter mfi(a_lhs); mfi.isValid(); ++mfi)
    {
      a_lhs[mfi].axby(a_x[mfi], a_y[mfi], a_a, a_b);
    }
  }

  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  assign(FabArray<EBCellFAB>& a_lhs,
         const FabArray<EBCellFAB>& a_rhs,
         const Real& a_scale)
  {
    BL_PROFILE("EBLevelDataOps::assign(to,from)");
    setVal(a_lhs, 0.);
    incr(a_lhs, a_rhs, a_scale);
  }

  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  incr( FabArray<EBCellFAB>& a_lhs,
        const FabArray<EBCellFAB>& a_rhs,
        const Real& a_scale)
  {
    BL_PROFILE("EBLevelDataOps::incr");
    BL_ASSERT(a_lhs.nComp()  == a_rhs.nComp());

    for(MFIter mfi(a_lhs); mfi.isValid(); ++mfi)
    {
      a_lhs[mfi].plus(a_rhs[mfi], 0, 0, a_lhs.nComp(), a_scale);
    }
  }
  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  scale(FabArray<EBCellFAB>& a_result,
        const Real&           a_value)
  {
    BL_PROFILE("EBLevelDataOps::scale");
    for(MFIter mfi(a_result); mfi.isValid(); ++mfi)
    {
      a_result[mfi] *= a_value;
    }
  }

  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  scale(FabArray<EBCellFAB>       & a_result,
        const FabArray<EBCellFAB> & a_value)
  {
    BL_PROFILE("EBLevelDataOps::scale_2");
    for(MFIter mfi(a_result); mfi.isValid(); ++mfi)
    {
      a_result[mfi] *= a_value[mfi];
    }
  }
  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  kappaWeight(FabArray<EBCellFAB>& a_data)
  {
    BL_PROFILE("EBLevelDataOps::kappaWeight");
    for(MFIter mfi(a_data); mfi.isValid(); ++mfi)
    {
      a_data[mfi].kappaWeight();
    }
  }

  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  gatherBroadCast(Real& a_accum, Real& a_volume, const int& a_p)
  {
    BL_PROFILE("EBLevelDataOps::gatherBroadcast1");
    //   Vector<Real> accum(1,a_accum);
//   gatherBroadCast(accum, a_volume, a_p);
//   a_accum = accum[0];
#ifdef BL_USE_MPI
    Real tmp=a_volume;
    MPI_Allreduce(&tmp, &a_volume, 1, ParallelDescriptor::Mpi_typemap<Real>::type(),
                  MPI_SUM, MPI_COMM_WORLD);
    tmp = a_accum;
    if (a_p==0)
    {
      MPI_Allreduce(&tmp, &a_accum, 1, ParallelDescriptor::Mpi_typemap<Real>::type(),
                    MPI_MAX, MPI_COMM_WORLD);
    }
    else
    {
      MPI_Allreduce(&tmp, &a_accum, 1, ParallelDescriptor::Mpi_typemap<Real>::type(),
                    MPI_SUM, MPI_COMM_WORLD);
    }
#endif

  }
  //-----------------------------------------------------------------------
  Real 
  EBLevelDataOps::
  norm(Real&                       a_volume,
       const FabArray<EBCellFAB>&  a_data,
       const EBLevelGrid        &  a_eblg,
       int                         a_p,
       int                         a_comp,
       bool                        a_useGhost)
  {
    Real integral  = 0;
    Real volume    = 0;
    for(MFIter mfi(a_eblg.getDBL(), a_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      const EBCellFAB& data =    a_data[mfi];
      Box     grid  = a_eblg.getDBL()  [mfi];
      EBISBox ebis  = a_eblg.getEBISL()[mfi];
      Box region = grid;
      if(a_useGhost)
      {
        region = data.box();
        region &= a_eblg.getDomain();
      }
      for(VoFIterator vofit(IntVectSet(region), ebis.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex & vof = vofit();
        Real value   = std::abs(data(vof, a_comp));
        Real volFrac = ebis.volFrac(vof);
        volume += volFrac;
        if(a_p == 0)
        {
          //p = 0 --- max norm
          if(value > integral)
          {
            integral = value;
          }
        }
        else if(a_p == 1)
        {
          // p = 1 L1 norm (integral |phi| dV)
          integral += volFrac*value;
        }
        else if(a_p == 2)
        {
          // p = 1 L2 norm (integral phi^2 dV)
          integral += volFrac*value*value;
        }
        else
        {
          amrex::Error("bogus norm value sent into EBLevelDataOps::norm");
        }
      }
    }
    //the values above are local to this proc.  now gather-broadcast to make them global
    gatherBroadCast(integral, volume, a_p);
    
    Real norm;
    if(a_p == 0)
    {
      //p = 0 --- max norm
      norm = integral;
    }
    else if(a_p == 1)
    {
      // p = 1 L1 norm (integral |phi| dV)
      //this check should be OK since there is no grid spacing here 
      // any single cell should have a volfrac > 1.0e-10
      if(volume > 1.0e-10)
      {
        norm = integral/volume;
      }
      else
      {
        norm = 0;
      }
    }
    else if(a_p == 2)
    {
      // p = 1 L2 norm (integral phi^2 dV)
      //this check should be OK since there is no grid spacing here 
      // any single cell should have a volfrac > 1.0e-10
      if(volume > 1.0e-10)
      {
        norm = sqrt(integral/volume);
      }
      else
      {
        norm = 0;
      }
    }
    else
    {
      amrex::Error("bogus norm value sent into EBLevelDataOps::norm");
    }

    return norm;
  }
  //-----------------------------------------------------------------------
  void 
  EBLevelDataOps::
  compareError(const FabArray<EBCellFAB>       &   a_errorFine,
               const FabArray<EBCellFAB>       &   a_errorCoar,
               const EBLevelGrid               &   a_eblgFine,
               const EBLevelGrid               &   a_eblgCoar,
               Vector<string>                 a_names,
               bool                                a_useGhost)
  {
    BL_ASSERT(a_errorFine.nComp() == a_errorCoar.nComp());

    Box domCoar = a_eblgCoar.getDomain();
    Box domCoFi = a_eblgFine.getDomain();
    domCoFi.coarsen(2);
    //this assumes fine is a on a domain refined by two of the coarse domain
    if(domCoar != domCoFi)
    {
      amrex::Error("domains do not match for ebleveldataops::compareerror");
    }
    const int ncomp = a_errorFine.nComp();
    bool useDefaultNames = (a_names.size() < ncomp);


    Vector<Real> coarNorm[3], fineNorm[3], orders[3];  //one for each type of norm;
    for (int p = 0; p <=2;  p++)
    {
      coarNorm[p].resize(ncomp, 0.);
      fineNorm[p].resize(ncomp, 0.);
      orders  [p].resize(ncomp, 0.);
    }
  
    for (int comp = 0; comp < ncomp; comp++)
    {
      for (int p = 0; p <=2;  p++)
      {
        Real coarVolu, fineVolu;
        coarNorm[p][comp] = norm(coarVolu, a_errorCoar, a_eblgCoar, p, comp, a_useGhost);
        fineNorm[p][comp] = norm(fineVolu, a_errorFine, a_eblgFine, p, comp, a_useGhost);

        if ((std::abs(fineNorm[p][comp]) > 1.0e-10) && (std::abs(coarNorm[p][comp]) > 1.0e-10))
        {
          orders[p][comp] = log(std::abs(coarNorm[p][comp]/fineNorm[p][comp]))/log(2.0);
        }
      }
    }

    int ncoar = a_eblgCoar.getDomain().size()[0];
    int nmedi = 2*ncoar;

    amrex::Print()    << setw(12)
                      << setprecision(6)
                      << setiosflags(ios::showpoint)
                      << setiosflags(ios::scientific);

    amrex::Print() << "\\begin{table}[p]" << "\n";
    amrex::Print() << "\\begin{center}" << "\n";
    amrex::Print() << "\\begin{tabular}{|cccc|} \\hline" << "\n";
    amrex::Print() << "Variable & \t norm \t \t & $e_{c}$ \t& Order \t& $e_{f}$\\\\" << "\n";;
    amrex::Print() << "\\hline " << "\n";

    for (int inorm = 0; inorm <= 2; inorm++)
    {

      for (int icomp = 0; icomp < ncomp; icomp++)
      {
        Real normCoar = coarNorm[inorm][icomp];
        Real normFine = fineNorm[inorm][icomp];
        Real order    =   orders[inorm][icomp];

        if (!useDefaultNames)
        {
          amrex::Print() << setw(14) << a_names[icomp] << " &    \t ";
        }
        else
        {
          amrex::Print() << "var" << icomp << " &    \t ";
        }

        if(inorm == 0)
        {
          amrex::Print() << "L_\\infty  & \t ";
        }
        else
        {
          amrex::Print() << "L_" << inorm << "       & \t ";
        }

        amrex::Print() << setw(12)
                       << setprecision(6)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << normCoar  << " & "
                       << setw(12)
                       << setprecision(3)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << order << " & "
                       << setw(12)
                       << setprecision(6)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << normFine;
        amrex::Print() << " \\\\ " << "\n";
      }

    }
    amrex::Print() << "\\hline " << "\n";
    amrex::Print() << "\\end{tabular}" << "\n";
    amrex::Print() << "\\end{center}" << "\n";
    amrex::Print() << "\\caption{Convergence rates for $nx_f = " << nmedi << ", nx_c = nx_f/2$.} " << "\n";
    amrex::Print() << "\\end{table}" << "\n";
    amrex::Print() << "\n" << "\n";

  }
}
