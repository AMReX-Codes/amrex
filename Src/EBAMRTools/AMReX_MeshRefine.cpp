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


#include <AMReX_MeshRefine.H>
#include <AMReX.H>
#include <AMReX_RealBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_BoxIterator.H>




/// Factory class to produce EBCellFABs.
/**
   This makes a set of grids from a set of tags.   
   This is not as general as Chombo's meshrefine--it is useful for 
   writing tests, however.
*/
namespace amrex
{
  class MRAmr: public AmrMesh
  {
  public:
    using AmrMesh::AmrMesh;  // inherit AmrMesh constructors

    virtual ~MRAmr() {}

    //! Manually tag.  Note that tags is built on level lev grids coarsened by bf_lev.
    virtual void ManualTagsPlacement (int lev, amrex::TagBoxArray& tags,
                                      const amrex::Array<amrex::IntVect>& bf_lev) override
      {
        IntVectSet coarsenedTagsLev = (*m_vec_ivs_ptr)[lev];
        coarsenedTagsLev.coarsen(bf_lev[lev][0]);

        for (MFIter mfi(tags); mfi.isValid(); ++mfi)
        {
          TagBox& tag = tags[mfi];
          const Box& bx = tag.box();
          for(BoxIterator boxit(bx); boxit.ok(); ++boxit)
          {
            if(coarsenedTagsLev.contains(boxit()))
            {
              tag(boxit()) = TagBox::SET;
            }
          }
        }
      }

    //not safe but this is a private class hidden away.  shhh.
    void setTags(const std::vector<IntVectSet>* a_vec_ivs_ptr)
      {
        m_vec_ivs_ptr = a_vec_ivs_ptr; 
      }
      
  private:
    
    const std::vector<IntVectSet>* m_vec_ivs_ptr;

  };

  void MeshRefine(std::vector<BoxArray>           &   a_grids,
                  const std::vector<IntVectSet>   &   a_tags,
                  const std::vector<int>          &   a_refRat,
                  const int                       &   a_maxLev,
                  const int                       &   a_blockingFactor,
                  const int                       &   a_properNestingRadius,
                  const int                       &   a_maxGridSize,
                  const Real                      &   a_gridEfficiency,
                  const Box                       &   a_coarsestDomain)
  {
    RealBox prob_domain(D_DECL(0.,0.,0.), D_DECL(1.,1.,1.));

    Array<int> n_cell(SpaceDim);
    IntVect domsize =  a_coarsestDomain.size();
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      n_cell[idir] = domsize[idir];
    }

    MRAmr amr(&prob_domain, a_maxLev, n_cell, -1, a_refRat);

    amr.setTags(&a_tags);
    amr.SetMaxGridSize(a_maxGridSize);
    amr.SetBlockingFactor(a_blockingFactor);
    amr.SetGridEff(a_gridEfficiency);
    amr.SetNProper(a_properNestingRadius);
    amr.MakeNewGrids();

    a_grids.resize(amr.finestLevel() + 1);
    for (int lev = 0; lev <= amr.finestLevel(); ++lev)
    {
      a_grids[lev] = amr.boxArray(lev);
    }
    
  }
}

