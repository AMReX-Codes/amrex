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

#ifndef _GEOMETRYSERVICE_H_
#define _GEOMETRYSERVICE_H_


#include <cmath>
#include <cstdlib>

#include "AMReX_LoHiSide.H"
#include "AMReX_RealVect.H"
#include "AMReX_Box.H"
#include "AMReX_IrregNode.H"


namespace amrex
{
  ///
  /**
     This is a pure base class that EBIndexSpace uses
     to interface with geometry-generation classes.
  */
  class GeometryService
  {
  public:
      
    enum InOut
    {
      Regular,
      Covered,
      Irregular
    };
      
    GeometryService();
      
    virtual ~GeometryService();
      
    ///
    /**
       Return true if every cell in region is regular at the
       refinement described by dx.
    */
    virtual bool isRegular(const Box&           a_region,
                           const Box&           a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const = 0;
      
    ///
    /**
       Return true if every cell in region is covered at the
       refinement described by dx.
    */
    virtual bool isCovered(const Box&           a_region,
                           const Box&           a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const = 0;
      
    virtual bool isIrregular(const Box&           a_region,
                             const Box&           a_domain,
                             const RealVect&      a_origin,
                             const Real&          a_dx) const ;
    ///
    /**
       irregGraph contains a complete list of irregular vofs. \\
       validRegion is the the region of the graph \\
       domain is the domain of compuation \\
       regIrregCovered = 1 for regular, 0 for irregular, -1 for covered.
       regIrregCovered must contain the valid region grown by 1. \\
       In chombospeak, \\
       ghostregion = (grow(validRegion, 1) & domain);
    */
    virtual void fillGraph(BaseFab<int>&        a_regIrregCovered,
                           Array<IrregNode>&    a_nodes,
                           const Box&           a_validRegion,
                           const Box&           a_ghostRegion,
                           const Box&           a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const = 0;
      
      
      
    virtual bool canGenerateMultiCells() const;
      
    virtual InOut InsideOutside(const Box&           a_region,
                                const Box&           a_domain,
                                const RealVect&      a_origin,
                                const Real&          a_dx) const ;
      
      
    ///handy functions to do rectangle intersections in real space
    static bool intersection(const RealVect& a_lo1, const RealVect& a_hi1,
                             const RealVect& a_lo2, const RealVect& a_hi2);
      
    static bool intersection(const Box&           a_region,
                             const RealVect&      a_origin,
                             const Real&          a_dx,
                             const RealVect&      a_lower,
                             const RealVect&      a_upper);
      
  };
}

#endif
