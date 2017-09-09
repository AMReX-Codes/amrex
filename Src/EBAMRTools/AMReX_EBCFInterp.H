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

#ifndef _NWOEBQUADCFINTERP_H___
#define _NWOEBQUADCFINTERP_H___

#include <iostream>
#include <math.h>
#include "AMReX_IntVect.H"
#include "AMReX_Box.H"
#include "AMReX_BaseFab.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_FabArray.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_LayoutData.H"
#include "AMReX_AggStencil.H"
#include "AMReX_MFIter.H"


namespace  amrex
{

  ///
  /**
     Class to fill the ghost cells over coarse-fine interfaces.  The 
     order of the polynomial used in the interpolation is adjustable.
     Second derivatives are as far as I go, for now.
   */
  class EBCFInterp
  {
  public:

    ///  
    EBCFInterp()
      {
        m_isDefined = false;
      }

    /// 
    /** makes fully defined object.   
       If orderOfPolynomial = 0,  no first or second derivative stencils included
       If orderOfPolynomial = 1,  no second derivative stencils included.
       Second derivatives are as far as I go, for now.
       Slow mode is slow.  Don't use it unless you are stepping through a debugger.
    */
    EBCFInterp(const EBLevelGrid      & a_eblgFine,
               const EBLevelGrid      & a_eblgCoar,
               const int              & a_refRat,
               const int              & a_ghostCellsInData,
               const int              & a_ghostCellsToFill,
               int  a_orderOfPolynomial = 2,
               bool a_slowMode = false)
      {
        define(a_eblgFine, a_eblgCoar, a_refRat, a_ghostCellsInData, a_ghostCellsToFill, a_slowMode);
      }


    ///used for defining a null-defined object.  all ghost cells at the coarse-fine interface filled
    /** makes fully defined object.   
       If orderOfPolynomial = 0,  no first or second derivative stencils included
       If orderOfPolynomial = 1,  no second derivative stencils included.
       Second derivatives are as far as I go, for now.
       Slow mode is slow.  Don't use it unless you are stepping through a debugger.
    */
    void  define(const EBLevelGrid      & a_eblgFine,
                 const EBLevelGrid      & a_eblgCoar,
                 const int              & a_refRat,
                 const int              & a_ghostCellsInData,
                 const int              & a_ghostCellsToFill,
                 int  a_orderOfPolynomial = 2,
                 bool a_slowMode = false);

    ///
    ~EBCFInterp()
      {  }

    ///
    /**
       Coarse / Fine interpolation operator.
    */
    void coarseFineInterp(FabArray<EBCellFAB>&       a_phif,
                          const FabArray<EBCellFAB>& a_phic,
                          int isrc, int idst, int inco);


    ///
    /**
       Coarse-fine interpolation operator with homogeneous BCs.
       Does same interpolation as coarseFineInterp function,
       but with all coarse-level data set to 0.
    */
    void coarseFineInterpH(FabArray<EBCellFAB>& a_phif,
                           int isrc, int idst, int inco);


    /// true if the define function or the full constructor has been called
    bool isDefined() const
      {
        return m_isDefined;
      }
  protected:
    void defineInternals();
  
    IntVectSet getCFIVS(const MFIter& a_mfi);

    void
    getStencil(VoFStencil           & a_stencil,
               const VolIndex       & a_vofFine,
               const EBISBox        & a_ebisFine,
               const EBISBox        & a_ebisCoFi);


    EBLevelGrid m_eblgCoar;
    EBLevelGrid m_eblgFine;
    EBLevelGrid m_eblgCoFi;

    int               m_refRat;
    int               m_dataGhost;
    int               m_fillGhost;
    bool              m_isDefined;
    LayoutData<std::shared_ptr< AggStencil<EBCellFAB, EBCellFAB> > > m_stencil;
    int m_orderOfPolynomial;
    int m_buffGhost;
    //for slow mode, we are just debugging so do stuff pointwise
    bool m_slowMode; //usually false
    LayoutData<Array<VoFStencil> > m_slowStencils;//usually left undefined
    LayoutData<Array<VolIndex  > >  m_slowVoFs;   //usually left undefined

  private:


    //disallowed for all the usual reasons
    EBCFInterp(EBCFInterp& a_input);

    //disallowed for all the usual reasons
    void operator=(EBCFInterp& a_input);


  };
}

#endif
