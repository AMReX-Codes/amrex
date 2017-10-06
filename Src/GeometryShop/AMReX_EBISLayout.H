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

#ifndef _EBISLAYOUT_H_
#define _EBISLAYOUT_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_EBISBox.H"
#include "AMReX_EBGraph.H"
#include "AMReX_EBData.H"
#include <AMReX_DistributionMapping.H>
#include <AMReX_MFIter.H>

namespace amrex
{

  class EBIndexSpace;
  class EBISLayout;
  template <class FAB> class FabArray;  // avoid circular depencies
///
/**
   EBISLayoutImplem represents the geometric information
   on a union of rectangles.  This is not really part
   of the public interface.  Users should use EBISLayout.
*/
  class EBISLayoutImplem
  {
  public:

    ///
    EBISBox operator[](const MFIter& a_datInd) const;

    ///
    EBISBox operator[](const int& a_boxindex) const;

    ///
    EBISLayoutImplem();

    ///
    ~EBISLayoutImplem();

    ///
    /**
       Returns the index of the VoF corresponding to coarsening
       the input VoF by the input ratio.  It is an error if
       the ratio is greater than the maximum coarsening ratio
       or if the vof does not exist at the input data index.
    */
    VolIndex coarsen(const VolIndex& a_vof,
                     const int     & a_ratio,
                     const MFIter  & a_datInd) const;

    ///
    /**
       Returns the index of the face corresponding to coarsening
       the input face by the input ratio.  It is an error if
       the ratio is greater than the maximum coarsening ratio
       or if the face does not exist at the input data index.
    */
    FaceIndex coarsen(const FaceIndex& a_face,
                      const int      & a_ratio,
                      const MFIter   & a_datInd) const;

    ///
    /**
       Returns the indices of the VoFs corresponding to refining
       the input VoF by the input ratio.  It is an error if
       the ratio is greater than the maximum refinement ratio
       or if the vof does not exist at the input data index.
    */
    Vector<VolIndex> refine(const VolIndex &  a_vof,
                                 const int      &  a_ratio,
                                 const MFIter   &  a_datInd) const;

    ///
    /**
       Returns the indices of the Faces corresponding to refining
       the input Face by the input ratio.  It is an error if
       the ratio is greater than the maximum refinement ratio
       or if the face does not exist at the input data index.
    */
    Vector<FaceIndex> refine(const FaceIndex &  a_face,
                                  const int       &  a_ratio,
                                  const MFIter    &  a_datInd) const;

    ///
    /**
       Sets the maximum level of refinement that this {\tt EBISLayout}
       will have to perform.  Creates and holds new {\tt EBISLayout}s
       at intermediate levels of refinement.
       Default is one (no refinement done).
    */
    void setMaxRefinementRatio(const int& a_maxRefine);


    ///
    /**
       Sets the maximum level of coarsening that this {\tt EBISLayout}
       will have to perform.  Creates and holds new {\tt EBISLayout}s
       at intermediate levels of refinement.
       Default is one (no coarsening done).
    */
    void setMaxCoarseningRatio(const int& a_maxCoarsen);

    int getMaxCoarseningRatio() const
      {
        return m_maxCoarseningRatio;
      }

    
    int getMaxRefinementRatio() const
      {
        return m_maxRefinementRatio;
      }

    ///
    void define(const Box               & a_domain,
                const BoxArray          & a_grids,
                const DistributionMapping & a_dm,
                const int               & a_nghost,
                const FabArray<EBGraph> & a_graph,
                const FabArray<EBData>  & a_data);


    bool isDefined() const {return m_defined;}

    shared_ptr<FabArray<EBGraph> > getAllGraphs() const
      {
        return m_ebGraph;
      } 

    shared_ptr<FabArray<EBData> > getAllData() const
      {
        return m_ebData;
      }
 
    EBISLayoutImplem(const Box               & a_domain,
                     const BoxArray          & a_grids,
                     const DistributionMapping & a_dm,
                     const int               & a_nghost,
                     const FabArray<EBGraph> & a_graph,
                     const FabArray<EBData>  & a_data)
      {
        define(a_domain, a_grids, a_dm, a_nghost, a_graph, a_data);
      }

    int getGhost() const
      {
        return m_nghost;
      }

    DistributionMapping getDM() const
      {
        return m_dm;
      }

    BoxArray getDBL() const
      {
        return m_dblInputDom;
      }
  private:
    ///
    Box m_domain;

    ///
    BoxArray  m_dblInputDom;

    DistributionMapping m_dm;
    ///
    Vector<EBISLayout> m_fineLevels;

    ///
    Vector<EBISLayout> m_coarLevels;

    ///
    int m_nghost;

    ///
    int m_maxCoarseningRatio;

    ///
    int m_maxRefinementRatio;


    ///
    shared_ptr<FabArray<EBGraph> > m_ebGraph;
    shared_ptr<FabArray<EBData > > m_ebData;

    bool m_defined;

  private:

    void operator=(const EBISLayoutImplem& ebiin)
      {;}
    EBISLayoutImplem(const EBISLayoutImplem& ebiin)
      {;}

  };

///
/**
   EBISLayout represents the geometric information
   on a union of rectangles.    This is
   another ref-counted pointer object so copying
   it is cheap but has the ref-counted pointer semantic.
*/
  class EBISLayout
  {
  public:
    ///
    /**
       Access the {\tt EBISBox} associated with the input DataIndex.
       Only constant access is permitted.
    */
    EBISBox operator[](const MFIter& a_index) const
      {
        return (*m_implem)[a_index];
      }
    

    ///
    EBISBox operator[](const int& a_boxindex) const
      {
        return (*m_implem)[a_boxindex];
      }

    ///
    /**
       Sets the maximum level of refinement that this {\tt EBISLayout}
       will have to perform.  Creates and holds new {\tt EBISLayout}s
       at intermediate levels of refinement.
       Default is one (no refinement done).
    */
    void setMaxRefinementRatio(const int& a_maxRefine)
      {
        m_implem->setMaxRefinementRatio(a_maxRefine);
      }



    ///
    /**
       Sets the maximum level of coarsening that this {\tt EBISLayout}
       will have to perform.  Creates and holds new {\tt EBISLayout}s
       at intermediate levels of refinement.
       Default is one (no coarsening done).
    */
    void setMaxCoarseningRatio(const int& a_maxCoarsen)
      {
        m_implem->setMaxCoarseningRatio(a_maxCoarsen);
      }


    int getMaxCoarseningRatio() const
      {
        return m_implem->getMaxCoarseningRatio();
      }

    int getMaxRefinementRatio() const
      {
        return m_implem->getMaxRefinementRatio();
      }

    ///
    EBISLayout()
      {
        m_implem = std::shared_ptr<EBISLayoutImplem>(new EBISLayoutImplem());
      }

      EBISLayout (const EBISLayout& rhs) = default;
      EBISLayout (EBISLayout&& rhs) noexcept = default;
      EBISLayout& operator= (const EBISLayout& rhs) = default;

      EBISLayout& operator= (EBISLayout&& rhs) = delete;

    ///
    ~EBISLayout()
      {
      }

    ///
    /**
       Returns the index of the VoF corresponding to coarsening
       the input VoF by the input ratio.  It is an error if
       the ratio is greater than the maximum coarsening ratio
       or if the vof does not exist at the input data index.
    */
    VolIndex coarsen(const VolIndex& a_vof,
                     const int     & a_ratio,
                     const MFIter  & a_datInd) const
      {
        return m_implem->coarsen(a_vof, a_ratio, a_datInd);
      }

    ///
    /**
       Returns the indices of the VoFs corresponding to refining
       the input VoF by the input ratio.  It is an error if
       the ratio is greater than the maximum refinement ratio
       or if the vof does not exist at the input data index.
    */
    Vector<VolIndex> refine(const VolIndex& a_vof,
                                 const int     & a_ratio,
                                 const MFIter  & a_datInd) const
      {
        return m_implem->refine(a_vof, a_ratio, a_datInd);
      }

    ///
    /**
       Returns the indices of the Faces corresponding to refining
       the input Face by the input ratio.  It is an error if
       the ratio is greater than the maximum refinement ratio
       or if the face does not exist at the input data index.
    */
    Vector<FaceIndex> refine(const FaceIndex & a_face,
                                  const int       & a_ratio,
                                  const MFIter    & a_datInd) const
      {
        return m_implem->refine(a_face, a_ratio, a_datInd);
      }

    ///
    void define(const Box               & a_domain,
                const BoxArray          & a_grids,
                const DistributionMapping & a_dm,
                const int               & a_nghost,
                const FabArray<EBGraph> & a_graph,
                const FabArray<EBData>  & a_data)
      {
        m_implem->define(a_domain, a_grids, a_dm, a_nghost, a_graph, a_data);
      }


    ///
    EBISLayout(const Box               & a_domain,
               const BoxArray          & a_grids,
               const DistributionMapping & a_dm,
               const int               & a_nghost,
               const FabArray<EBGraph> & a_graph,
               const FabArray<EBData>  & a_data)
      {
        m_implem = std::shared_ptr<EBISLayoutImplem>(new EBISLayoutImplem(a_domain, a_grids, a_dm, a_nghost, a_graph, a_data));
      }

    bool isDefined() const
      {
        return m_implem->isDefined();
      }

     shared_ptr<FabArray<EBGraph> > getAllGraphs() const
      {
        return m_implem->getAllGraphs();
      }

     shared_ptr<FabArray<EBData> > getAllData() const
      {
        return m_implem->getAllData();
      }

    int getGhost() const
      {
        return m_implem->getGhost();
      }
    DistributionMapping getDM() const
      {
        return m_implem->getDM();
      }

    BoxArray getDBL() const
      {
        return m_implem->getDBL();
      }
  private:
    ///
    std::shared_ptr<EBISLayoutImplem> m_implem;
  };


}
#endif
