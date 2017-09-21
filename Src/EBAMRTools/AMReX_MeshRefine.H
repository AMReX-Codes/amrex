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


#ifndef _MeshRefine_Wrapper_
#define _MeshRefine_Wrapper_

#include <cmath>
#include <cstdlib>
#include <AMReX_Box.H>
#include <AMReX_IntVectSet.H>
#include <AMReX_BoxArray.H>
#include <AMReX_RealVect.H>
#include <AMReX_EBLevelGrid.H>


namespace amrex
{
  class GridParameters
  {
  public:
    int           whichGeom    ;
    IntVect       nCells;
    int           maxGridSize;
    int           properNest;
    int           blockFactor;
    int           bufferSize;
    Real          fillRatio;
    int           maxLevel;
    int           numLevels;

    // number of ghost cells (in each direction) for the EBISBoxes
    int           ghostEBISBox;

    Array<int>   refRatio;
    Box           coarsestDomain;
    Real          coarsestDx;
    Real          domainLength;
    RealVect      probLo;
    RealVect      probHi;

    virtual void coarsen(int a_factor);
    virtual void refine(int a_factor);
    virtual void pout() const;
  
    GridParameters()
      {
      }
    virtual ~GridParameters()
      {
      }
  };

  /// function to make grids
  /**
     This makes a set of grids from a set of tags.   
     This is not as general as Chombo's meshrefine (it has
     none of the stuff necessary for lbase != 0).
     This is useful for  writing tests, however.  
  */
  extern void
  MeshRefine(Array<BoxArray>           &   a_grids,
             const Array<IntVectSet>   &   a_tags,
             const Array<int>          &   a_refRat,
             const int                       &   a_maxLev,
             const int                       &   a_blockingFactor,
             const int                       &   a_properNestingRadius,
             const int                       &   a_maxGridSize,
             const Real                      &   a_gridEfficiency,
             const Box                       &   a_coarsestDomain);


  ///for tests---get grids where the irregular cells are all refined
  extern void getAllIrregEBLG(Array<EBLevelGrid>   & a_eblg,
                              const GridParameters  & a_params);

  extern void
  getGridParameters(GridParameters&  a_params,
                    int forceMaxLevel = -1,
                    bool a_verbose = false);

}
#endif
