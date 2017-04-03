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


#ifndef _AGGEBPWLFILLPATCH_H_
#define _AGGEBPWLFILLPATCH_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Stencils.H"
#include "AggStencil.H"
#include "EBIndexSpace.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"

///
/**
   Fills one layer of fine level ghost cells
   by piecewise linear interpolation
   from the coarse level.
*/
class AggEBPWLFillPatch
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  AggEBPWLFillPatch();

  ///
  virtual  ~AggEBPWLFillPatch();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.\\
     radius: the total number of
     ghost cells filled.
     forceNoEBCF = true promises that the embedded boundary 
     will not cross the coarse-fine interface.   This will save you  
     a  lot of compute time if you are certain you are tagging all irregulkar cells.
  */
  AggEBPWLFillPatch(const DisjointBoxLayout& a_dblFine,
                    const DisjointBoxLayout& a_dblCoar,
                    const EBISLayout& a_ebislFine,
                    const EBISLayout& a_ebislCoar,
                    const ProblemDomain& a_domainCoar,
                    const int& a_nref,
                    const int& a_nvar,
                    const int& a_radius,
                    const IntVect& a_ghost,
                    const bool& a_forceNoEBCF = false,
                    const EBIndexSpace* const a_eb = Chombo_EBIS::instance());


  /// saner interface
  AggEBPWLFillPatch(const EBLevelGrid & a_eblgFine,
                    const EBLevelGrid & a_eblgCoar,
                    const int& a_nref,
                    const int& a_nvar,
                    const int& a_radius,
                    const IntVect& a_ghost,
                    const bool& a_forceNoEBCF = false,
                    const EBIndexSpace* const a_eb = Chombo_EBIS::instance())
  {
    setDefaultValues();
    define(a_eblgFine.getDBL(),
           a_eblgCoar.getDBL(),
           a_eblgFine.getEBISL(),
           a_eblgCoar.getEBISL(),
           a_eblgCoar.getDomain(),
           a_nref, a_nvar, a_radius, a_ghost, a_forceNoEBCF, a_eb);
  }
                    
                    
  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
     radius: the total number of
     ghost cells filled.
     forceNoEBCF = true promises that the embedded boundary 
     will not cross the coarse-fine interface.   This will save you  
     a  lot of compute time if you are certain you are tagging all irregulkar cells.
  */
  void define(const DisjointBoxLayout& a_dblFine,
              const DisjointBoxLayout& a_dblCoar,
              const EBISLayout& a_ebislFine,
              const EBISLayout& a_ebislCoar,
              const ProblemDomain& a_domainCoar,
              const int& a_nref,
              const int& a_nvar,
              const int& a_radius,
              const IntVect& a_ghost,
              const bool& a_forceNoEBCF = false,
              const EBIndexSpace* const a_eb = Chombo_EBIS::instance());

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
     over the intersection
     of the fine layout with the refined coarse layout. \\
     {\bf Arguments:}\\
     coarDataOld: The data over the coarse layout at coarse time old.\\
     coarDatanew: The data over the coarse layout at coarse time new.\\
     fineData: The data over the fine layout.  \\
     Error occurs unless  coarTimeOld <= fineTime <= coarTimeNew \\
     Fine and coarse data must
     have the same number of variables.\\
     variables:  The variables to interpolate.  Those not
     in this range will be left alone.   This range of variables
     must be in both the coarse and fine data.
  */
  virtual void
  interpolate(LevelData<EBCellFAB>& a_fineData,
              const LevelData<EBCellFAB>& a_coarDataOld,
              const LevelData<EBCellFAB>& a_coarDataNew,
              const Real& a_coarTimeOld,
              const Real& a_coarTimeNew,
              const Real& a_fineTime,
              const Interval& a_variables) const;

protected:

  //internal use only functions

  virtual void
  definePieceWiseLinearFillPatch(const DisjointBoxLayout& a_dblfine,
                                 const DisjointBoxLayout& a_dblCoar);

  void getOffsets(const LayoutData< Vector<VolIndex> >&  a_srcVoFsCoar,
                  const LayoutData<IntVectSet>&          a_irregRegionsFine,
                  const LayoutData<Vector<VoFStencil> >  a_loStencils[SpaceDim],
                  const LayoutData<Vector<VoFStencil> >  a_hiStencils[SpaceDim],
                  const LayoutData<IntVectSet>           a_coarLoInterp[SpaceDim],
                  const LayoutData<IntVectSet>           a_coarHiInterp[SpaceDim],
                  const LayoutData<IntVectSet>           a_coarCeInterp[SpaceDim]);


  void
  interpolateFAB(EBCellFAB& a_fine,
                 const EBCellFAB& a_coarOld,
                 const EBCellFAB& a_coarNew,
                 const Real& a_coarTimeOld,
                 const Real& a_coarTimeNew,
                 const Real& a_fineTime,
                 const DataIndex& a_datInd,
                 const Interval& a_variables) const;

  void
  makeStencils();

  void defineSlopeHolders(const LayoutData<IntVectSet>& a_irregRegionsCoFi);

  void
  getIVS(LayoutData<IntVectSet>& a_irregRegionsFine,
         LayoutData<IntVectSet>& a_irregRegionsCoFi,
         LayoutData< Vector<VolIndex> >&   a_srcVoFs);

  void
  defineAggStencils(LayoutData<Vector<VoFStencil> >  a_loStencils[SpaceDim],
                    LayoutData<Vector<VoFStencil> >  a_hiStencils[SpaceDim],
                    const LayoutData< Vector<VolIndex> >&   a_srcVoFs);
  void
  getLoHiCenIVS(LayoutData<IntVectSet>  a_coarLoInterp[SpaceDim],
                LayoutData<IntVectSet>  a_coarHiInterp[SpaceDim],
                LayoutData<IntVectSet>  a_coarCeInterp[SpaceDim]);

  void
  getSten(LayoutData<Vector<VoFStencil> >  a_loStencils[SpaceDim],
          LayoutData<Vector<VoFStencil> >  a_hiStencils[SpaceDim],
          LayoutData<IntVectSet>           a_coarLoInterp[SpaceDim],
          LayoutData<IntVectSet>           a_coarHiInterp[SpaceDim],
          LayoutData<IntVectSet>           a_coarCeInterp[SpaceDim],
          const LayoutData< Vector<VolIndex> >&   a_srcVoFs);

  void getSlopes(const EBCellFAB& a_coarDataOld,
                 const EBCellFAB& a_coarDataNew,
                 const DataIndex& a_dit,
                 const Interval& a_variables) const;

  void deltaMinMod(Real& a_slope, Real& a_slopeLo, Real& a_slopeHi) const;

  void setDefaultValues();

  struct
  {
    size_t offset;
    int    dataID;
  } typedef access_t;

  struct
  {
    access_t dest_access; //offsets for fine data
    int      slop_index;  //index into vector of slope offsetsf
    IntVect  coariv;
    IntVect  fineiv;
  } typedef fine_logic_t;

  struct
  {
    access_t slop_access;      //offsets for slope data
    access_t data_access;      //offsets for coarsened fine data
    bool     has_lo[SpaceDim]; //whether there is a lo slope stencil
    bool     has_hi[SpaceDim]; //whether there is a hi slope stencil
    bool     ivs_lo[SpaceDim]; //whether inside ivs centered
    bool     ivs_hi[SpaceDim]; //whether inside ivs centered
    bool     ivs_ce[SpaceDim]; //whether inside ivs centered
  } typedef coar_logic_t;

  //data to provide fast access to fine points and their associated slopes
  LayoutData< Vector<fine_logic_t> > m_fineOffsets;
  //data so we can do the wacky delta minmod thing fast.
  LayoutData< Vector<coar_logic_t> > m_coarOffsets;

  bool m_forceNoEBCF;
  bool m_isDefined;

  //non-EB fill patch
  PiecewiseLinearFillPatch* m_patcher;

  //ghost cells in solution
  IntVect m_ghost;
  //domain at the coarse levle
  ProblemDomain m_coarDomain;
  //layouts
  DisjointBoxLayout m_coarsenedFineGrids;
  DisjointBoxLayout m_fineGrids;
  DisjointBoxLayout m_coarGrids;
  EBISLayout m_coarsenedFineEBISL;
  EBISLayout m_fineEBISL;


  mutable LevelData<EBCellFAB> m_coarOnFDataOld;
  mutable LevelData<EBCellFAB> m_coarOnFDataNew;

  int m_refRat;
  int m_nComp;
  int m_radius;
  int m_coarGhostRad;

  //building this pig for speed only.   this is scratch data for slopes
  mutable LayoutData< BaseIVFAB<Real> > m_slopeLoOld[SpaceDim];
  mutable LayoutData< BaseIVFAB<Real> > m_slopeHiOld[SpaceDim];
  mutable LayoutData< BaseIVFAB<Real> > m_slopeCeOld[SpaceDim];
  mutable LayoutData< BaseIVFAB<Real> > m_slopeLoNew[SpaceDim];
  mutable LayoutData< BaseIVFAB<Real> > m_slopeHiNew[SpaceDim];
  mutable LayoutData< BaseIVFAB<Real> > m_slopeCeNew[SpaceDim];

  LayoutData< RefCountedPtr<AggStencil <EBCellFAB, BaseIVFAB<Real> > > > m_stenLo[SpaceDim];
  LayoutData< RefCountedPtr<AggStencil <EBCellFAB, BaseIVFAB<Real> > > > m_stenHi[SpaceDim];

private:
  //disallowed for all the usual reasons.  No code specified because Jeff says that is better.
  AggEBPWLFillPatch(const AggEBPWLFillPatch& ebcin);

  void operator=(const AggEBPWLFillPatch& fabin);

};

#include "NamespaceFooter.H"
#endif
