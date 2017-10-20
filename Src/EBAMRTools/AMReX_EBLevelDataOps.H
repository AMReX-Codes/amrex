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
#include "AMReX_MultiFab.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBFluxFAB.H"
#include "AMReX_EBFluxFactory.H"
#include "AMReX_MultiFab.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_EBArith.H"


//for creating simple grids for tests
namespace amrex
{

  ///
  /**
     Repository of useful EB functions over a level.  All functions are static so just
     use EBLevelDataOps::blah(...);
  */
  class EBLevelDataOps
  {
  public:
    EBLevelDataOps()
      {
      }

    ~EBLevelDataOps()
      {
      }

    ///quick and dirty debugging vis.
    static void  viewEBLevel(const FabArray<EBCellFAB> * a_data, const EBLevelGrid* a_eblg);

    ///writes plotfile that visit can eat.   Just single-valued stuff
    static void writeSingleLevelEBPlotFile(const std::string         & a_filename,
                                           const FabArray<EBCellFAB> & a_data,
                                           const EBLevelGrid         & a_eblg,
                                           const Vector<string>      & a_varNames);


    ///writes plotfile that visit can eat.   Just single-valued stuff
    static void writeEBAMRPlotFile(const std::string                   & a_filename,
                                   const Vector<FabArray<EBCellFAB>* > & a_data,
                                   const Vector<EBLevelGrid>           & a_eblg,
                                   const Vector<int>                   & a_refRat,
                                   const Vector<string>                & a_varNames);

    /// tacks on volume fraction as the last variable
    static void makePlotMultiFab(MultiFab                              & a_mfSingleValued,
                                 const FabArray<EBCellFAB>             & a_ebdata,
                                 const EBLevelGrid                     & a_eblg);

    ///
    /**
       Take the regular pointers of the eb cell data and put it into a MF
     */
    static void aliasIntoMF(shared_ptr<MultiFab>      & a_regData,
                            const FabArray<EBCellFAB> & a_ebcData,
                            const EBLevelGrid         & a_eblg);

    ///
    /**
       Take the regular pointers of the eb face data and put it into a MF
     */
    static void aliasIntoMF(shared_ptr<MultiFab>      & a_regData,
                            const FabArray<EBFluxFAB> & a_ebfData,
                            const int                 & a_faceDir,
                            const EBLevelGrid         & a_eblg);

    ///
    /**
     */
    static bool checkForBogusNumbers(const FabArray<EBCellFAB> & a_data);


    ///
    /**
     */
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
    static void setCoveredVal(FabArray<EBCellFAB> &  a_lData,
                              const int           &  a_comp,
                              const Real          &  a_value);


    ///
    /**
     */
    static void scale(FabArray<EBCellFAB>& a_lhs,
                      const Real& a_scale);


    ///
    /**
    */
    static void scale(FabArray<EBCellFAB>       & a_lhs,
                      const FabArray<EBCellFAB> & a_scale);


    ///
    /**
       If comp < 0, set all components
     */
    static void setVal(FabArray<EBCellFAB> & a_result,
                       const Real          & a_value,
                       int  a_comp = -1);

    ///
    /**
       If comp < 0, set all components
     */
    static void setVal(FabArray<EBFluxFAB> & a_result,
                       const Real          & a_value);



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
    static void axby( FabArray<EBCellFAB>&       a_lhs,
                      const FabArray<EBCellFAB>& a_x,
                      const FabArray<EBCellFAB>& a_y,
                      const Real& a,
                      const Real& b);
    ///
    /**
       set useGhost = false unless you are testing cf interp or some such
     */
    static  Real norm(Real&                       a_volume,
                      const FabArray<EBCellFAB>& a_data,
                      const EBLevelGrid        &  a_eblg,
                      int                         a_p,
                      int                         a_comp,
                      bool a_useGhost = false);


    /// prints out latex summary of error comparison.  useful for convergence tests
    /**
       set useGhost = false unless you are testing cf interp or some such
     */
    static void compareError(const FabArray<EBCellFAB>       &   a_errorFine,
                             const FabArray<EBCellFAB>       &   a_errorCoar,
                             const EBLevelGrid               &   a_eblgFine,
                             const EBLevelGrid               &   a_eblgCoar,
                             Vector<string> a_names = Vector<string>(),
                             bool a_useGhost = false);
    ///
    /**
     */
    static  Real kappaDotProduct(Real&                       a_volume,
                                 const FabArray<EBCellFAB> & a_data1,
                                 const FabArray<EBCellFAB> & a_data2, 
                                 const EBLevelGrid         & a_eblg);

  private:

  };
}

#endif
