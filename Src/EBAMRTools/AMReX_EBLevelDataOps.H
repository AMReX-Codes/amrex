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

#ifndef _EBLEVELDATAOPS_H___
#define _EBLEVELDATAOPS_H___

#include "AMReX_IntVect.H"
#include "AMReX_BoxArray.H"
#include "AMReX_FabArray.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBISLayout.H"


namespace amrex
{
///
/**
   Repository of useful EB functions over a level.  All functions are static so just
   use EBFabArrayOps::blah(...);
*/
  class EBFabArrayOps
  {
  public:
    EBFabArrayOps()
      {
      }

    ~EBFabArrayOps()
      {
      }

    static bool checkForBogusNumbers(const FabArray<EBCellFAB> & a_data);


    static void  checkData(const FabArray<EBCellFAB>&a_data, const string& label);

    ///
    /**
     */
    static int  parallelSum(const int&  a_value);

    ///
    /**
     */
    static long long parallelSum(const long long&  a_value);

    ///
    /**
     */
    static Real parallelSum(const Real& a_value);

    ///
    /**
     */
    static int  parallelMin(const int&  a_value);

    ///
    /**
     */
    static int  parallelMax(const int&  a_value);

    ///
    /**
     */
    static Real parallelMin(const Real& a_value);

    ///
    /**
     */
    static Real parallelMax(const Real& a_value);


    ///
    /**
     */
    static void getMaxMin(Real&                       a_maxVal,
                          Real&                       a_minVal,
                          const FabArray<EBCellFAB>& a_data,
                          const int&                  a_comp,
                          const bool&                 a_doAbs=false);

    ///
    /**
     */
    static void setCoveredVal(FabArray<EBCellFAB>&   a_lData);


    ///
    /**
     */
    static void scale(FabArray<EBCellFAB>& a_lhs,
                      const Real& a_scale);


    ///
    /**
     */
    static void setVal(FabArray<EBCellFAB>& a_result,
                       const Real&           a_value);



    ///
    /**
     */
    static void incr( FabArray<EBCellFAB>& a_lhs,
                      const FabArray<EBCellFAB>&   a_rhs,
                      const Real& a_scale);

    ///
    /**
     */
    static void assign( FabArray<EBCellFAB>& a_lhs,
                        const FabArray<EBCellFAB>&   a_rhs,
                        const Real& a_scale);


    //! Multiply each datum in \a a_data by the volume fraction of the
    //! corresponding cell.
    //! \param a_data The data to be multiplied by the volume fraction.
    static  void kappaWeight(FabArray<EBCellFAB>& a_data);


    ///
    /**
     */
    static void gatherBroadCast(Real& a_accum, Real& a_volume, const int& a_p);

    ///
    /**
     */
    static void gatherBroadCast(std::vector<Real>& a_accum, Real& a_volume, const int& a_p);


    ///
    static void axby( FabArray<EBCellFAB>&       a_lhs,
                      const FabArray<EBCellFAB>& a_x,
                      const FabArray<EBCellFAB>& a_y,
                      const Real& a,
                      const Real& b);
    ///
    /**
     */
    static  Real norm(Real&                       a_volume,
                      const FabArray<EBCellFAB>& a_data,
                      const ProblemDomain&        a_domain,
                      int                         a_p,
                      int                         a_comp);


    ///
    /**
     */
    static void compareError(std::vector<Real>               &   a_orders,
                             const FabArray<EBCellFAB>       &   a_errorFine,
                             const FabArray<EBCellFAB>       &   a_errorCoar,
                             const Box                       &   a_domain,
                             std::vector<string> a_names = std::vector<string>());
  private:

  };
}

#endif
