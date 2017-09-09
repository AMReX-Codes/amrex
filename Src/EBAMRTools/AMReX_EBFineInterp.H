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

#ifndef _EBPWLFINEINTERP_H___
#define _EBPWLFINEINTERP_H___

#include "AMReX_FabArray.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBArith.H"
#include "AMReX_LayoutData.H"
#include "AMReX_AggStencil.H"


namespace amrex
{

///
/**
   Class to interpolate data from a coarser grid to a finer grid.   The order
   of the polynomial used in the interpolation is adjustable.
*/
  class EBFineInterp
  {
  public:
    ///
    /**
       Default constructor.  User must subsequently call define().
    */
    EBFineInterp()
      {
        m_isDefined = false;
      }

    ///
    ~EBFineInterp()
      {
      }

    ///
    /**
       Defining constructor.  Constructs a valid object.
       Equivalent to default construction followed by define().
       If orderOfPolynomial = 0,  no first or second derivative stencils included
       If orderOfPolynomial = 1,  no second derivative stencils included
       This does not do slope limiting.   Second derivatives are as far as I go, for now.
       Slow mode is slow.  Don't use it unless you are stepping through a debugger.
    */
    EBFineInterp(const EBLevelGrid   & a_eblgFine,
                 const EBLevelGrid   & a_eblgCoar,
                 const int           & a_nref,
                 const int           & a_ghostCellsInData,
                 int a_orderOfPolynomial = 2,
                 bool a_slowMode = false)
      {

        define(a_eblgFine, a_eblgCoar, a_nref, a_ghostCellsInData, a_orderOfPolynomial);
      }

    ///
    /**
       Defines this object.  Existing information is overriden.
       If orderOfPolynomial = 0,  no first or second derivative stencils included
       If orderOfPolynomial = 1,  no second derivative stencils included
       This does not do slope limiting.
       Slow mode is slow.  Don't use it unless you are stepping through a debugger.
    */
    void define(const EBLevelGrid   & a_eblgFine,
                const EBLevelGrid   & a_eblgCoar,
                const int           & a_nref,
                const int           & a_ghostCellsInData,
                int a_orderOfPolynomial = 2,
                bool a_slowMode = false);

    ///
    /**
       Returns true if this object was created with the defining
       constructor or if define() has been called.
    */
    bool isDefined() const
      {
        return m_isDefined;
      }

    ///
    /**
       Interpolate the fine data from the coarse data
       over the intersection to the order of accuracy requested.
       NOTE: This function sends these arrays
       to Fortran so you must set the covered
       cell values to something.  Unitialized memory
       might cause IEEE problems with some compilers.
       {\bf Arguments:}\\
       coarData: The data over the coarse layout.\\
       fineData: The data over the fine layout.  \\
       The variables to interpolate.  Those not
       in this range will be left alone.   This range of variables
       must be in both the coarse and fine data.
    */
    void
    interpolate(FabArray<EBCellFAB>& a_fineData,
                const FabArray<EBCellFAB>& a_coarseData,
                int isrc, int idst, int inco);

  protected:

    void defineInternals();

    void
    getStencil(VoFStencil           & a_stencil,
               const VolIndex       & a_vofFine,
               const EBISBox        & a_ebisFine,
               const EBISBox        & a_ebisCoFi,
               const MFIter         & a_mfi);
    bool m_isDefined;

    EBLevelGrid m_eblgCoar;
    EBLevelGrid m_eblgFine;
    EBLevelGrid m_eblgCoFi;
    int m_refRat;
    int m_buffGhost;
    int m_dataGhost;
    LayoutData<shared_ptr<AggStencil<EBCellFAB, EBCellFAB> > > m_stencil;
    int m_orderOfPolynomial;

    //for slow mode, we are just debugging so do stuff pointwise
    bool m_slowMode; //usually false
    LayoutData<Array<VoFStencil> > m_slowStencils;//usually left undefined
    LayoutData<Array<VolIndex  > >  m_slowVoFs;   //usually left undefined
  private:
    //disallowed for all the usual reasons
    EBFineInterp(const EBFineInterp& ebcin);
    void operator=(const EBFineInterp& fabin);

  };
}

#endif
