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

/// Quadratic coarse-fine interpolation utility 
  class NWOEBQuadCFInterp
  {
  public:

    ///  
    NWOEBQuadCFInterp()
      {
        m_isDefined = false;
      }

    /// makes fully defined object.   all ghost cells at the coarse-fine interface filled
    NWOEBQuadCFInterp(const EBLevelGrid      & a_eblgFine,
                      const EBLevelGrid      & a_eblgCoar,
                      const int              & a_refRat,
                      const int              & a_ghostCellsInData)
      {
        define(a_eblgFine, a_eblgCoar, a_refRat, a_ghostCellsInData);
      }


    ///used for defining a null-defined object.  all ghost cells at the coarse-fine interface filled
    void  define(const EBLevelGrid      & a_eblgFine,
                 const EBLevelGrid      & a_eblgCoar,
                 const int              & a_refRat,
                 const int              & a_ghostCellsInData);

    ///
    ~NWOEBQuadCFInterp()
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
  
    void
    getStencils(std::vector<std::shared_ptr< BaseStencil> >  & a_stencils, 
                std::vector<std::shared_ptr< BaseIndex  > >  & a_baseDstVoFs,
                const MFIter                                 & a_mfi);


    void
    getStencil(VoFStencil           & a_stencil,
               const VolIndex       & a_vofFine,
               const EBISBox        & a_ebisFine,
               const EBISBox        & a_ebisCoFi);


    EBLevelGrid m_eblgCoar;
    EBLevelGrid m_eblgFine;
    EBLevelGrid m_eblgCoFi;

    int               m_refRat;
    int               m_ghost;
    bool              m_isDefined;
    LayoutData<std::shared_ptr< AggStencil<EBCellFAB, EBCellFAB> > > m_stencil;
    FabArray<EBCellFAB> m_bufferCoFi;
  

  private:


    //disallowed for all the usual reasons
    NWOEBQuadCFInterp(NWOEBQuadCFInterp& a_input);

    //disallowed for all the usual reasons
    void operator=(NWOEBQuadCFInterp& a_input);


  };
}

#endif
