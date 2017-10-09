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


#ifndef _FLATPLATEGEOM_H_
#define _FLATPLATEGEOM_H_

#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"
#include "AMReX_Box.H"
#include "AMReX_IntVect.H"
#include "AMReX_BaseIF.H"
#include "AMReX_BaseFab.H"
#include "AMReX_Moments.H"
#include "AMReX_IrregNode.H"
#include "AMReX_GeometryService.H"


namespace amrex
{
  ///
  /**
     This defines a geometry that is an infinitely thin, finite size plate
     aligned with a coordinate plane.
  */
  class FlatPlateGeom: public GeometryService
  {
  public:

    ///
    /**
    */
    FlatPlateGeom(const int&      a_normalDir,
                  const Real&     a_plateLocation, //where plate is in the normal direciton.
                  const RealVect& a_plateLo,       //only tangential locations used
                  const RealVect& a_plateHi)      //only tangential locations used
    {
      m_normalDir      = a_normalDir;
      m_plateLocation  = a_plateLocation;
      m_plateLo        = a_plateLo;
      m_plateHi        = a_plateHi;
      m_plateLo[a_normalDir] = a_plateLocation;
      m_plateHi[a_normalDir] = a_plateLocation;
    }

    ///
    ~FlatPlateGeom()
    {
    }

    bool
    isCellCut(const IntVect            & a_iv,
              const Box                & a_domain,
              const RealVect           & a_origin,
              const Real               & a_dx) const;

    /**
       Return true if every cell in region is regular at the
       refinement described by dx.
    */
    virtual bool isRegular(const Box&           a_region,
                           const Box& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;
    ///
    /**
       Return true if every cell in region is covered at the
       refinement described by dx.
    */
    virtual bool isCovered(const Box&           a_region,
                           const Box& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;

    ///return true if the cell face is cut by the plate
    bool isFaceCut(Real                     & a_areaFracLo,
                   Real                     & a_areaFracHi,
                   const IntVect            & a_iv, 
                   const int                & a_faceDir,
                   const Side::LoHiSide     & a_sd, 
                   const Box                & a_domain,
                   const RealVect           & a_origin,
                   const Real               & a_dx) const;


    int 
    getNumVolumes(const IntVect  & a_iv, 
                  const Box      & a_domain, 
                  const RealVect & a_origin, 
                  const Real     & a_dx) const;

    void
    addIrregularNodes(Vector<IrregNode>   & a_nodes,
                      const BaseFab<int>       & a_numVolumes,
                      const IntVect            & a_iv,
                      const Box                & a_domain,
                      const RealVect           & a_origin,
                      const Real               & a_dx) const;
    /**
    */
    virtual void fillGraph(BaseFab<int>&              a_regIrregCovered,
                           Vector<IrregNode>&    a_nodes,
                           const Box&                 a_validRegion,
                           const Box&                 a_ghostRegion,
                           const Box&                 a_domain,
                           const RealVect&            a_origin,
                           const Real&                a_dx) const;

  private:

    int      m_normalDir;    
    Real     m_plateLocation;
    RealVect m_plateLo;      
    RealVect m_plateHi;          
  };
}
#endif
