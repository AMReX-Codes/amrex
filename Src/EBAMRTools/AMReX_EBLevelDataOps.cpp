/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#include <cmath>

#include "SPMD.H"

#include "AMReX_EBLevelDataOps.H"
#include "AMReX_FaceIterator.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_BoxIterator.H"

namespace amrex
{

  void
  EBLevelDataOps::checkData(const FabArray<EBCellFAB>&a_data, const string& label)
  {
    barrier();
    amrex::Print() << "==== checking " << label << " for nans and infs =====" << std::"\n";
    checkForBogusNumbers(a_data);
    barrier();
    amrex::Print() << "Getting maxs and mins for " << label << std::"\n";
    for(int icomp = 0; icomp < a_data.nComp(); icomp++)
    {
      Real maxVal, minVal;
      getMaxMin(maxVal,minVal,a_data,icomp);
      amrex::Print() << "max value = " << maxVal << " for comp = " << icomp << std::"\n";
      amrex::Print() << "min value = " << minVal << " for comp = " << icomp << std::"\n";
    }
    amrex::Print() << "==================================================================== " << std::"\n";
    barrier();
  }

/*****/
  Real EBLevelDataOps::parallelSum(const Real& a_value)
  {
    // Find the sum of all a_value's
    Real sum = a_value;
#ifdef BL_USE_MPI
    Real sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &sum, 1, MPI_CH_REAL,MPI_SUM, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelSum");
    }
#endif
    return sum;
  }

/*****/
  int  EBLevelDataOps::parallelSum(const int& a_value)
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

/*****/
  long long EBLevelDataOps::parallelSum(const long long& a_value)
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

/*****/
  int EBLevelDataOps::parallelMin(const int& a_value)
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

/*****/
  int  EBLevelDataOps::parallelMax(const int& a_value)
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
/*****/
  Real EBLevelDataOps::parallelMin(const Real& a_value)
  {
    // Find the minimum of a_value's
    Real val = a_value;
#ifdef BL_USE_MPI
    Real sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &val, 1, MPI_CH_REAL,MPI_MIN, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMin");
    }
#endif
    return val;
  }

/*****/
  Real EBLevelDataOps::parallelMax(const Real& a_value)
  {
    // Find the maximum of a_value's
    Real val = a_value;
#ifdef BL_USE_MPI
    Real sendBuf = a_value;
    int result = MPI_Allreduce(&sendBuf, &val, 1, MPI_CH_REAL,MPI_MAX, MPI_COMM_WORLD);

    if (result != MPI_SUCCESS)
    {
      amrex::Error("Sorry, but I had a communication error in EBLevelDataOps::parallelMax");
    }
#endif
    return val;
  }


  void EBLevelDataOps::setCoveredVal(FabArray<EBCellFAB>&a_data,
                                     const int&           a_comp,
                                     const Real&          a_value)
  {
    CH_TIME("EBLevelDataOps::setCoveredVal(cell,comp)");
    for(MFIter mfi(a_data); mfi.isValid(); ++mfi)
    {
      a_data[mfi].setCoveredCellVal(a_value,a_comp);
    }
  }



  bool 
  EBLevelDataOps::checkForBogusNumbers(const FabArray<EBCellFAB>&a_data)
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
          if (std::isnan(val) || std::isinf(val) || Abs(val)>1.e16)
          {
            badvof = vof;
            badval = val;
            badcomp = icomp;
            //                  amrex::Print() << "      icomp = " << icomp << " vof = " << vof << " val = " << val << std::"\n";
            dataIsNANINF = true;
            //                  amrex::Error();
          }
        }
      }
    }
    if (dataIsNANINF)
    {
      amrex::Print() << "first bad vof = "  << badvof << "\n";
      amrex::Print() << "bad val = "  << badval << " at comp " << badcomp << "\n";
      amrex::Warning("Found a NaN or Infinity.");
    }
    return dataIsNANINF;
  }
  void EBLevelDataOps::getMaxMin(Real&                       a_maxVal,
                                 Real&                       a_minVal,
                                 const FabArray<EBCellFAB>& a_data,
                                 const int&                  a_comp,
                                 const bool&                 a_doAbs)
  {
    CH_TIME("EBLevelDataOps::getMaxMin");
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
          a_maxVal = Max(a_maxVal,Abs(val));
          a_minVal = Min(a_minVal,Abs(val));
        }
        else
        {
          a_maxVal = Max(a_maxVal,val);
          a_minVal = Min(a_minVal,val);
        }
      }
    }
    a_minVal = EBLevelDataOps::parallelMin(a_minVal);
    a_maxVal = EBLevelDataOps::parallelMax(a_maxVal);
  }


  void EBLevelDataOps::setVal(FabArray<EBCellFAB>& a_result,
                              const Real&           a_value,
                              const int&            a_comp)
  {
    for(MFIter mfi(a_result); mfi.isValid(); ++mfi)
    {
      a_result[mfi].result.setVal(a_comp,a_value);
    }
  }


  void EBLevelDataOps::axby( FabArray<EBCellFAB>&       a_lhs,
                             const FabArray<EBCellFAB>& a_x,
                             const FabArray<EBCellFAB>& a_y,
                             const Real& a,
                             const Real& b)
  {
    //  CH_assert(a_lhs.boxArray() == a_x.boxArray());

    DataIterator dit = a_lhs.dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
    {
      DataIndex d = dit[mybox];

      EBCellFAB& data = a_lhs[d];
      data.axby(a_x[d], a_y[d], a, b);
    }
  }

  void EBLevelDataOps::assign(FabArray<EBCellFAB>& a_lhs,
                              const FabArray<EBCellFAB>& a_rhs,
                              const Real& a_scale)
  {
    CH_TIME("EBLevelDataOps::assign(to,from)");
    setVal(a_lhs, 0.);
    incr(a_lhs, a_rhs, a_scale);
  }

  void EBLevelDataOps::incr( FabArray<EBCellFAB>& a_lhs,
                             const FabArray<EBCellFAB>& a_rhs,
                             const Real& a_scale)
  {
    //  CH_assert(a_lhs.boxArray() == a_rhs.boxArray());
    DataIterator dit = a_lhs.dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
    {
      DataIndex d = dit[mybox];

      a_lhs[d].plus(a_rhs[d], a_scale);
    }
  }


  void EBLevelDataOps::scale(FabArray<EBCellFAB>& a_result,
                             const Real&           a_value,
                             const int&            a_comp)
  {
    DataIterator dit = a_result.dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
    {
      DataIndex d = dit[mybox];

      EBCellFAB& result = a_result[d];
  
      result.mult(a_value,a_comp, 1);
    }
  }

//-----------------------------------------------------------------------
  void EBLevelDataOps::kappaWeight(FabArray<EBCellFAB>& a_data)
  {
    DataIterator dit = a_data.dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
    {
      DataIndex d = dit[mybox];

      EBCellFAB& data = a_data[d];

      kappaWeight(data);
    }
  }
//-----------------------------------------------------------------------


  void EBLevelDataOps::gatherBroadCast(Real& a_accum, Real& a_volume, const int& a_p)
  {
    //   Vector<Real> accum(1,a_accum);
//   gatherBroadCast(accum, a_volume, a_p);
//   a_accum = accum[0];
#ifdef BL_USE_MPI
    Real tmp=a_volume;
    MPI_Allreduce(&tmp, &a_volume, 1, MPI_CH_REAL, MPI_SUM, MPI_COMM_WORLD);
    tmp = a_accum;
    if (a_p==0)
    {
      MPI_Allreduce(&tmp, &a_accum, 1, MPI_CH_REAL,
                    MPI_MAX, MPI_COMM_WORLD);
    }
    else
    {
      MPI_Allreduce(&tmp, &a_accum, 1, MPI_CH_REAL,
                    MPI_SUM, MPI_COMM_WORLD);
    }
#endif

  }

  void EBLevelDataOps::gatherBroadCast(Vector<Real>& a_accum, Real& a_volume, const int& a_p)
  {
#ifdef BL_USE_MPI
    {
      CH_TIME("MPI_Barrier EBLevelDataOps::gatherBroadcast");
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    //gather what each processor thinks is the accum and volume
    int ncomp = a_accum.size();
    Real volumeLoc = a_volume;
    Vector<Real> accumLoc  = a_accum;

    Vector<Real> volumeVec;
    Vector<Vector<Real> > accumVec;
    int baseproc = 0;
    gather(volumeVec, volumeLoc, baseproc);
    gather( accumVec,  accumLoc, baseproc);

    a_volume = 0.;
    for (int i=0; i<ncomp; i++) a_accum[i]=0.0;
    if (procID() == baseproc)
    {
      for (int ivec = 0; ivec < numProc(); ivec++)
      {
        a_volume += volumeVec[ivec];
        Vector<Real> cur =   accumVec[ivec];
        for (int i=0; i<ncomp; i++)
        {
          if (a_p == 0)
          {
            if (cur[i] > a_accum[i])
            {
              a_accum[i] = cur[i];
            }
          }
          else
          {
            a_accum[i]  += cur[i];
          }
        }
      }
    }

    //broadcast the sum to all processors.
    broadcast( a_accum, baseproc);
    broadcast(a_volume, baseproc);
  }
  /*********/
  Real 
  EBLevelDataOps::
  norm(Real&                       a_volume,
       const FabArray<EBCellFAB>& a_data,
       const ProblemDomain&        a_domain,
       int                         a_p,
       int                         a_comp)
  {
  }
  /*********/
  void 
  EBLevelDataOps::
  compareError(std::vector<Real>               &   a_orders,
               const FabArray<EBCellFAB>       &   a_errorFine,
               const FabArray<EBCellFAB>       &   a_errorCoar,
               const Box                       &   a_domain,
               std::vector<string>                 a_names)
  {
    amrex::Error("not implemented");
  }
}
