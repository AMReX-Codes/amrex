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


#ifndef _EBGRAPH_H_
#define _EBGRAPH_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_RealVect.H"
#include "AMReX_BaseFab.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_IrregNode.H"
#include "AMReX_GraphNode.H"
#include "AMReX_EBCellFlag.H"

namespace amrex
{
  class EBGraph;
  /// Geometric description within a box
  /**
     EBGraphImplem represents the geometric information
     of the domain at a given refinement
     within the boundaries of a particular box.
  */
  class EBGraphImplem
  {
  public:
        
    ///sure hope this does not matter
    typedef int value_type;

    
    ///
    /**
       The EBGraphImplem can be all regular or all
       covered or contain irregular cells.
    */
    enum TAG
    {
      AllCovered = 0,
      AllRegular, 
      HasIrregular
    };

    
    ///
    /**
       Makes an EBGraphImplem whose contents are undefined.
    */
    EBGraphImplem();
        
    ///
    ~EBGraphImplem();
        
    ///
    /**
       changing this definition to make things run faster/better
    */
    const BaseFab<char>& getMask(int& a_regIrregCovered) const;
        
    ///
    /**
       Return the VoFs in the cell.
    */
    Array<VolIndex> getVoFs(const IntVect& a_iv) const;
        
    ///
    /**
       Return true if every cell in the EBGraphImplem is a regular cell.
    */
    bool isAllRegular() const;
        
    ///
    /**
       Return true neither all regular nor all covered
    */
    bool hasIrregular() const;
        
    ///
    /**
       Return true if  every cell in the EBGraphImplem is a covered cell.
    */
    bool isAllCovered() const;
        
    ///
    /**
       Return true if a_iv is an irregular cell.
    */
    bool isIrregular(const IntVect& a_iv) const;

        
    ///
    /**
       Return true if every cell in a_box is a
       covered cell.
    */
    bool isCovered(const Box& a_box) const;
        
    ///
    /**
       Return true if every cell in a_box is a
       regular cell.
    */
    bool isRegular(const Box& a_box) const;
        
    ///
    bool isRegular(const IntVect& a_iv) const;
        
    ///
    bool isCovered(const IntVect& a_iv) const;
        
    ///
    bool isMultiValued(const IntVect& a_iv) const;
        
    ///
    bool isConnected(const VolIndex& a_vof1,
                     const VolIndex& a_vof2) const;
    ///
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    Array<FaceIndex> getFaces(const VolIndex& a_vof,
                                    const int& a_idir,
                                    const Side::LoHiSide& a_sd) const;
        
    ///
    Array<Array<VolIndex> >  getVoFSets(const Box& a_region) const;
        
    ///
    /**
       Returns the corresponding set of VoFs from the next finer
       EBGraphImplem (factor of two refinement).  The result is only
       defined if this {\tt EBGraphImplem} was defined by coarsening.
    */
    Array<VolIndex> refine(const VolIndex& a_coarVoF) const;
        
    ///
    /**
       Returns the corresponding set of faces from the next finer
       EBGraphImplem (factor of two refinement).  The result is only
       defined if this {\tt EBGraphImplem} was defined by coarsening.
    */
    Array<FaceIndex> refine(const FaceIndex& a_coarFace, const EBGraphImplem& a_fineGraph) const;
        
    ///
    /**
       Returns the corresponding face from the next coarser
       EBGraphImplem (factor of two refinement).  This only makes
       sense if there is a coarser level.
    */
    FaceIndex coarsen(const FaceIndex& a_coarFace) const;
        
    ///
    /**
       Returns the corresponding  VoF from the next coarser
       EBGraphImplem (same solution location, different index space, factor
       of two refinement ratio).
    */
    VolIndex coarsen(const VolIndex& a_fineVoF) const;
        
    ///
    /**
       Completely construct the graph from IrregNodes.
       This makes a DENSE representation of the graph.  To make
       sparse representations of all regular or all covered graphs,
       use setToAllRegular or setToAllCovered.\\
       irregGraph contains a complete list of irregular vofs. \\
       validRegion is the the region of the graph \\
       domain is the domain of compuation \\
       regIrregCovered = 1 for regular, 0 for irregular, -1 for covered.
       regIrregCovered must contain the valid region grown by 1. \\
       In chombospeak, \\
       regIrregCovered.box().contains(grow(validRegion, 1) & domain); \\
       otherwise an error will result
    */
    void
    buildGraph(const BaseFab<int>&      a_regIrregCovered,
               const Array<IrregNode>& a_irregGraph,
               const Box&               a_validRegion,
               const Box&     a_domain);
        
        
    ///
    /**
       Set the graph to all regular cells (sparse);
    */
    void setToAllRegular();
        
    ///
    /**
       Set the graph to all covered cells (sparse);
    */
    void setToAllCovered();
        
    ///
    /**
       Define an all-regular EBGraphImplem with the input box as its
       region.
    */
    void define(const Box& box);
        
    EBGraphImplem(const Box& box);
        
    ///
    /**
       Copy the information from a_source to the over the intersection
       of the box a_region, the box of the current EBGraphImplem and
       the box of a_source.  The Interval arguments are ignored.
       This function is required by LevelData.
    */
    ///
    const Box& getRegion() const;
        
    ///
    const Box& getDomain() const;
        
    EBGraphImplem& copy(const EBGraphImplem& src,
                        const Box&        srcbox,
                        int               srccomp,
                        const Box&        destbox,
                        int               destcomp,
                        int               numcomp);
        
        
    //null op because this does not make any sense in this context
    EBGraphImplem& plus(const EBGraphImplem& src,
                        const Box&        srcbox,
                        const Box&        destbox,
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1)
      {
        amrex::Error("this function should not be called");
        return *this;
      }
        
        
    //this includes all meta data
    std::size_t nBytesFull() const;


    //this includes all meta data
    std::size_t copyToMemFull(void*      dst) const;

    //this includes all meta data
    std::size_t copyFromMemFull(const void* src);

    //this is for copy and fillboundary where all the meta data has been filled
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const;
        
    //this is for copy and fillboundary where all the meta data has been filled
    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;
        
        
    //this is for copy and fillboundary where all the meta data has been filled
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);
    ///
    bool isDefined() const;
        
    ///
    bool isDomainSet() const;
        
    //stuff below is not part of the public API
        
    //define vofs to be the coarsened vofs of the inputs
    //this must go first in coarsening.  defines happen here
    void coarsenVoFs(const EBGraphImplem&   a_fineGraph,
                     const Box& a_coarRegion);
        
    //define the faces to be the coarsened faces of the input
    void coarsenFaces(const EBGraphImplem& a_fineEBIS,
                      const Box          & a_coarRegion);
        
        
    //fix fineVoF->coarseVoF information of the input
    void fixFineToCoarse(EBGraphImplem& a_fineEBIS,
                         const Box    & a_coarRegion) const;
        
    //don't ask
    Array<int>
    coarsenFaces(const VolIndex&       a_coarVoF,
                 const EBGraphImplem&        a_fineGraph,
                 const int&            a_idir,
                 const Side::LoHiSide& a_sd);
        
    ///
    void setDomain(const Box& a_domain);
        
    ///
    IntVectSet getIrregCells(const Box& a_subbox) const;
        
    ///
    IntVectSet getMultiCells(const Box& a_subbox) const;
        
    ///
    Array<FaceIndex>
    getAllFaces(const IntVect& a_iv,
                const int& a_idir,
                const Side::LoHiSide& a_sd) const;
        
    ///
    /**
       0       outside problem domain, multi-valued, covered
       1       regular
    */
    void fillMask(BaseFab<char>& a_mask) const;
        
    ///
    /**
       0       outside problem domain
       1       covered
       2       regular
       3       irregular
       handy for debugging, can write out to hdf5 files.
    */
    void fillCellTypeMask(BaseFab<char>& a_mask) const;
        
        
        
    ///
    void addFullIrregularVoFs(const IntVectSet& a_vofsToChange);
        
        
    ///
    void
    getRegNextToMultiValued(IntVectSet&      a_vofsToChange,
                            const Box &      a_region) const;
        
    ///
    Array<FaceIndex> getMultiValuedFaces(const int& a_idir, const Box& a_box, const EBGraph& a_ebgraph) const;
        
    //bool operator<(const EBGraphImplem& ebiin) const;
        
        
    int numVoFs(const IntVect& a_iv) const
      {
        Array<VolIndex> vofs = getVoFs(a_iv);
        return vofs.size();
      }
    const IntVectSet& getAllIrregCells() const
      {
        return m_irregIVS;
      }
    const IntVectSet& getAllMultiCells() const
      {
        return m_multiIVS;
      }
    
    Box getFullRegion() const
      {
        return m_fullRegion;
      }

    const EBCellFlagFab& getEBCellFlagFab() const
      {
        return m_cellFlags;
      }

    //for debugging
    static bool s_verbose;

  ///    
  Array<VolIndex>
  getVoFs (const VolIndex& a_vof,
           const int& a_dir,
           const Side::LoHiSide& a_sd,
           const int& a_steps) const;
  private:

    void setCellFlags();

    void fixCellFlagType();

    ///
    EBCellFlagFab m_cellFlags;
    
    ///for serialization
    GraphNode getGraphNode(const IntVect& a_iv) const;

    ///box over which this graph is defined
    Box m_region;

    ///box NOT intersected  with domain
    Box m_fullRegion;

    ///computational domain at this level of refinement
    Box m_domain;
        
    ///
    /**
       If this is allregular or allcovered,
       the BaseFab below is undefined.
    */
    TAG m_tag;

    ///
    /**
       the graph description of the box.
       If the pointer GraphNode.m_thisNode == 0, the cell is covered.
       If the pointer GraphNode.m_thisNode == 1, the cell is regular.
       Otherwise, the pointer goes to the node of the graph.
    */
    BaseFab<GraphNode> m_graph;
        
    ///
    bool m_isDefined;
        
    ///
    bool m_isDomainSet;
        
    ///
    IntVectSet m_irregIVS;
        
    ///
    IntVectSet m_multiIVS;
        
        
    //mutable because it gets built when needed
    mutable BaseFab<char> m_mask;
    mutable bool m_isMaskBuilt;
        
  private:
        
    void operator=(const EBGraphImplem& ebiin) = delete;
    EBGraphImplem(const EBGraphImplem& ebiin) = delete;
        
    friend class EBIndexSpace;
    friend class EBISLevel;
  };
        
  /// Geometric description within a box
  /**
     EBGraph represents the geometric information
     of the domain at a given refinement
     within the boundaries of a particular box.
     This is a refcounted object with all
     that implies.
  */
  class EBGraph
  {
  public:
    ///sure hope this does not matter
    typedef int value_type;
        
    ///
    /**
       The EBGraph can be all regular or all
       covered or contain irregular cells.
    */
    enum TAG
    {
      AllRegular = 0,
      AllCovered,
      HasIrregular
    };
        
    ///
    /**
       Makes an EBGraph whose contents are undefined.
    */
    EBGraph();
        
    ///
    ~EBGraph();

    EBGraph (EBGraph&& rhs) noexcept = default;
        
    ///
    /**
       If the graph has irregular cells
       get a mask over the size of the graph.
       -1 for covered cells or multivalued cells, 0 for singlevalued irregular, and 1 for regular.
       If the graph is all regular, regIrregCovered returns 1 and
       the mask is undefined.
       If the graph is all regular, regIrregCovered returns -1 and
       the mask is undefined.
    */
    const  BaseFab<char>& getMask(int& a_regIrregCovered) const;
        
    ///
    /**
       Return the VoFs in the cell.
    */
    Array<VolIndex> getVoFs(const IntVect& a_iv) const;
        
    int numVoFs(const IntVect& a_iv) const
      {
        return m_implem->numVoFs(a_iv);
      }

    ///
    Box getFullRegion() const
      {
        return m_implem->getFullRegion();
      }

    ///
    const EBCellFlagFab& getEBCellFlagFab() const
      {
        return m_implem->getEBCellFlagFab();
      }

    ///
    /**
       Return true if every cell in the EBGraph is a regular cell.
    */
    bool isAllRegular() const;
        
    ///
    /**
       Return true neither all regular nor all covered
    */
    bool hasIrregular() const;
        
    ///
    /**
       Return true if  every cell in the EBGraph is a covered cell.
    */
    bool isAllCovered() const;
        
    ///
    /**
       Return true if a_iv is an irregular cell.
    */
    bool isIrregular(const IntVect& a_iv) const;
        
    ///
    /**
       Return true if every cell in a_box is a
       covered cell.
    */
    bool isCovered(const Box& a_box) const;
        
    ///
    /**
       Return true if every cell in a_box is a
       regular cell.
    */
    bool isRegular(const Box& a_box) const;
        
    ///
    bool isRegular(const IntVect& a_iv) const;
        
    ///
    bool isCovered(const IntVect& a_iv) const;
        
    ///
    bool isMultiValued(const IntVect& a_iv) const;
        
    ///
    bool isConnected(const VolIndex& a_vof1,
                     const VolIndex& a_vof2) const;
    ///
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    Array<FaceIndex> getFaces(const VolIndex& a_vof,
                                    const int& a_idir,
                                    const Side::LoHiSide& a_sd) const;
        
    ///
    Array<Array<VolIndex> >  getVoFSets(const Box& a_region) const;
    ///
    /**
       Returns the corresponding set of VoFs from the next finer
       EBGraph (factor of two refinement).  The result is only
       defined if this {\tt EBGraph} was defined by coarsening.
    */
    Array<VolIndex> refine(const VolIndex& a_coarVoF) const;
        
    ///
    /**
       Returns the corresponding set of faces from the next finer
       EBGraph (factor of two refinement).  The result is only
       defined if this {\tt EBGraph} was defined by coarsening.
    */
    Array<FaceIndex> refine(const FaceIndex& a_coarFace, const EBGraph& a_fineGraph) const;
        
    ///
    /**
       Returns the corresponding face from the next coarser
       EBGraph (factor of two refinement).  This only makes
       sense if there is a coarser level.
    */
    FaceIndex coarsen(const FaceIndex& a_coarFace) const;
        
    ///
    /**
       Returns the corresponding  VoF from the next coarser
       EBGraph (same solution location, different index space, factor
       of two refinement ratio).
    */
    VolIndex coarsen(const VolIndex& a_fineVoF) const;
        
    ///
    /**
       Define an all-regular EBGraph with the input box as its
       region.  The component argument is ignored.
       This function is required by LevelData.
    */
    EBGraph(const Box& a_box, int a_comps, 
            bool       alloc = true,
            bool       shared = false);
        
        
    ///
    /**
       Completely construct the graph from IrregNodes.
       This makes a DENSE representation of the graph.  To make
       sparse representations of all regular or all covered graphs,
       use setToAllRegular or setToAllCovered.\\
       irregGraph contains a complete list of irregular vofs. \\
       validRegion is the the region of the graph \\
       domain is the domain of compuation \\
       regIrregCovered = 1 for regular, 0 for irregular, -1 for covered.
       regIrregCovered must contain the valid region grown by 1. \\
       In chombospeak, \\
       regIrregCovered.box().contains(grow(validRegion, 1) & domain); \\
       otherwise an error will result
    */
    void
    buildGraph(const BaseFab<int>&      a_regIrregCovered,
               const Array<IrregNode>& a_irregGraph,
               const Box&               a_validRegion,
               const Box&     a_domain);
        
    ///
    /**
       Set the graph to all regular cells (sparse);
    */
    void setToAllRegular();
        
    ///
    /**
       Set the graph to all covered cells (sparse);
    */
    void setToAllCovered();
        
    ///
    /**
       Define an all-regular EBGraph with the input box as its
       region.
    */
    void define(const Box& box);
        
    ///need for fabarray
    void resize(const Box& box, const int & _ncomp)
      {
        define(box);
      }
    ///
    /**
       Copy the information from a_source to the over the intersection
       of the box a_region, the box of the current EBGraph and
       the box of a_source.  The Interval arguments are ignored.
       This function is required by LevelData.
    */
        
    EBGraph& copy(const EBGraph   & src,
                  const Box&        srcbox,
                  int               srccomp,
                  const Box&        destbox,
                  int               destcomp,
                  int               numcomp)
      {
        m_implem->copy(*src.m_implem, srcbox, srccomp, destbox, destcomp, numcomp);
        return *this;
      }
        
        
    //null op because this does not make any sense in this context
    EBGraph& plus(const EBGraph   & src,
                  const Box&        srcbox,
                  const Box&        destbox,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp=1)
      {
        amrex::Error("this function should not be called");
        return *this;
      }
        

    /// this is for copy and fillboundary where the meta data has already been filled
    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const
      {
        return m_implem->copyToMem(srcbox, srccomp, numcomp, dst);
      }
        
        
    /// this is for copy and fillboundary where the meta data has already been filled
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src)
      {
        return m_implem->copyFromMem(dstbox, dstcomp, numcomp, src);
      }
        

    ///    
    Array<VolIndex>
    getVoFs (const VolIndex& a_vof,
             const int& a_dir,
             const Side::LoHiSide& a_sd,
             const int& a_steps) const
      {
        return m_implem->getVoFs(a_vof, a_dir, a_sd, a_steps);
      }

    /// this is for copy and fillboundary where the meta data has already been filled
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const
      {
        return m_implem->nBytes(bx, start_comp, ncomps);
      }


    ///this includes all meta data
    std::size_t copyToMemFull(void*      dst) const
      {
        return m_implem->copyToMemFull(dst);
      }
        
        

    ///this includes all meta data
    std::size_t copyFromMemFull(const void* src)
      {
        return m_implem->copyFromMemFull(src);
      }
        

    const IntVectSet& getAllIrregCells() const
      {
        return m_implem->getAllIrregCells();
      }
    const IntVectSet& getAllMultiCells() const
      {
        return m_implem->getAllMultiCells();
      }

    ///this includes all meta data
    std::size_t nBytesFull() const
      {
        return m_implem->nBytesFull();
      }

    ///
    bool hasIrregular(const Box& a_box) const
      {
        bool retval = true;
        if(isRegular(a_box) || isCovered(a_box))
        {
          retval = false;
        }
        return retval;
      }

    ///
    const Box& getRegion() const;
        
    ///
    const Box& getDomain() const;
        
    ///
    bool isDefined() const;
        
    ///
    bool isDomainSet() const;
        
    //stuff below is not part of the public API
        
    //define vofs to be the coarsened vofs of the inputs
    //this must go first in coarsening.  defines happen here
    void coarsenVoFs(const EBGraph&   a_fineGraph,
                     const Box& a_coarRegion);
        
    //define the faces to be the coarsened faces of the input
    void coarsenFaces(const EBGraph& a_coarGhostGraph,
                      const Box    & a_region);
        
    //fix fineVoF->coarseVoF information of the input
    void fixFineToCoarse(EBGraph      & a_fineEBIS,
                         const Box    & a_coarRegion) const;
        
    //don't ask
    Array<int>
    coarsenFaces(const VolIndex&       a_coarVoF,
                 const EBGraph&        a_fineGraph,
                 const int&            a_idir,
                 const Side::LoHiSide& a_sd);
        
    ///
    void setDomain(const Box& a_domain);
        
    ///
    IntVectSet getIrregCells(const Box& a_subbox) const;
        
    ///
    IntVectSet getMultiCells(const Box& a_subbox) const;
        
    ///
    EBGraph& operator=(const EBGraph& ebiin);
        
    ///
    EBGraph(const EBGraph& a_ebiin);
        
    ///
    Array<FaceIndex>
    getAllFaces(const IntVect& a_iv,
                const int& a_idir,
                const Side::LoHiSide& a_sd) const;
        
    void fillMask(BaseFab<char>& a_mask) const;
        
    ///
    /**
       0       outside priblem domain
       1       covered
       2       regular
       3       irregular
    */
    void fillCellTypeMask(BaseFab<char>& a_mask) const;
        
    ///
    /**
       This sets the  BaseFab on the input region to
       1 = regular
       0 = irregular
       -1 = covered
       2 outside the domain
       This does NOT resize the basefab (it was not thread safe to do so).
    */
    void fillIntMask(BaseFab<int>& a_mask,
                     const Box   & a_subbox) const;
        
        
    ///
    EBGraph(const Box& box);
    ///pointer comparison
    bool operator==(const EBGraph& a_ebiin);
        
    /// pointer comparison
    bool operator<(const EBGraph& a_ebiin) const;
        
    ///
    /**
       This stuff required by LevelData in parallel:
    */
        
    ///
    void addFullIrregularVoFs(const IntVectSet& a_vofsToChange);
        
    ///
    void
    getRegNextToMultiValued(IntVectSet&    a_vofsToChange,
                            const Box&    a_region) const;
        
    ///
    Array<FaceIndex> getMultiValuedFaces(const int& a_idir, const Box& a_box) const;
        
    ///
    static bool preAllocatable () { return false; }
    static bool isCopyOMPSafe () { return false; }
        
  private:
        
    ///
    std::shared_ptr<EBGraphImplem> m_implem;
        
    friend class EBIndexSpace;
    friend class EBISLevel;
  };
        
        
        
  /*******************************/
  /*  inline EBGraph functions */
  /*******************************/
        
  inline
  const BaseFab<char>&
  EBGraph::getMask(int& a_regIrregCovered) const
  {
    return m_implem->getMask(a_regIrregCovered);
  }
        
  inline void
  EBGraph::buildGraph(const BaseFab<int>&      a_regIrregCovered,
                      const Array<IrregNode>& a_irregGraph,
                      const Box&               a_validRegion,
                      const Box&     a_domain)
  {
    m_implem->buildGraph(a_regIrregCovered, a_irregGraph, a_validRegion, a_domain);
  }
  /*******************************/
  inline bool EBGraph::isDefined() const
  {
    return m_implem->isDefined();
  }
  /*******************************/
  inline bool EBGraph::isConnected(const VolIndex& a_vof1,
                                   const VolIndex& a_vof2) const
  {
    return m_implem->isConnected(a_vof1, a_vof2);
  }
  /*******************************/
  inline const Box& EBGraph::getRegion() const
  {
    return m_implem->getRegion();
  }
  /*******************************/
  inline EBGraph& EBGraph::operator=(const EBGraph& a_ebiin)
  {
    if (&a_ebiin != this)
    {
      m_implem = a_ebiin.m_implem;
    }
    return *this;
  }
  /*******************************/
  inline IntVectSet EBGraph::getMultiCells(const Box& a_subbox) const
  {
    return m_implem->getMultiCells(a_subbox);
  }
  /*******************************/
  inline IntVectSet EBGraph::getIrregCells(const Box& a_subbox) const
  {
    return m_implem->getIrregCells(a_subbox);
  }
  /*******************************/
  inline Array<VolIndex> EBGraph::getVoFs(const IntVect& a_iv) const
  {
    return m_implem->getVoFs(a_iv);
  }
  /*******************************/
  inline bool EBGraph::isRegular(const IntVect& a_iv) const
  {
    return m_implem->isRegular(a_iv);
  }
  /*******************************/
  inline bool EBGraph::isAllRegular() const
  {
    return m_implem->isAllRegular();
  }
  /*******************************/
  inline bool EBGraph::isAllCovered() const
  {
    return m_implem->isAllCovered();
  }
  /*******************************/
  inline bool EBGraph::isIrregular(const IntVect& a_iv) const
  {
    return m_implem->isIrregular(a_iv);
  }
  /*******************************/
  inline bool EBGraph::isCovered(const IntVect& a_iv) const
  {
    return m_implem->isCovered(a_iv);
  }
  /*******************************/
  inline bool EBGraph::isMultiValued(const IntVect& a_iv) const
  {
    return m_implem->isMultiValued(a_iv);
  }
  /*******************************/
  inline bool EBGraph::isCovered(const Box& a_box) const
  {
    return m_implem->isCovered(a_box);
  }
  /*******************************/
  inline bool EBGraph::isRegular(const Box& a_box) const
  {
    return m_implem->isRegular(a_box);
  }
  /*******************************/
  inline Array<FaceIndex>
  EBGraph::getFaces(const VolIndex& a_vof,
                    const int& a_idir,
                    const Side::LoHiSide& a_sd) const
  {
    return m_implem->getFaces(a_vof, a_idir, a_sd);
  }
  /*******************************/
  inline Array<VolIndex>
  EBGraph::refine(const VolIndex& a_coarVoF) const
  {
    return m_implem->refine(a_coarVoF);
  }
  /*******************************/
  inline Array<FaceIndex>
  EBGraph::refine(const FaceIndex& a_coarFace, const EBGraph& a_fineGraph) const
  {
    return m_implem->refine(a_coarFace, *a_fineGraph.m_implem);
  }
  /*******************************/
  inline FaceIndex EBGraph::coarsen(const FaceIndex& a_fineFace) const
  {
    return m_implem->coarsen(a_fineFace);
  }
  /*******************************/
  inline VolIndex EBGraph::coarsen(const VolIndex& a_fineVoF) const
  {
    return m_implem->coarsen(a_fineVoF);
  }
  /*******************************/
  inline bool EBGraph::operator==(const EBGraph& a_ebiin)
  {
    return (m_implem == a_ebiin.m_implem);
  }
  /*******************************/
  inline bool EBGraph::operator<(const EBGraph& a_ebiin) const
  {
    return (m_implem < a_ebiin.m_implem);
  }
  /*******************************/
  inline const Box& EBGraph::getDomain() const
  {
    return (m_implem->getDomain());
  }
  /*******************************/
  inline void EBGraph::setDomain(const Box& a_domain)
  {
    m_implem->setDomain(a_domain);
  }
  /*******************************/
  inline void EBGraph::setToAllRegular()
  {
    m_implem->setToAllRegular();
  }
  /*******************************/
  inline void EBGraph::setToAllCovered()
  {
    m_implem->setToAllCovered();
  }
  /*******************************/
  inline void EBGraph::
  coarsenVoFs(const EBGraph& a_fineEBIS,
              const Box& a_coarRegion)
  {
    m_implem->coarsenVoFs(*a_fineEBIS.m_implem, a_coarRegion);
  }
  /*******************************/
  inline void EBGraph::
  coarsenFaces(const EBGraph& a_fineGraph, const Box& a_coarRegion)
  {
    m_implem->coarsenFaces(*a_fineGraph.m_implem, a_coarRegion);
  }
  /*******************************/
  inline void EBGraph::
  fixFineToCoarse(EBGraph  & a_fineEBIS,
                  const Box& a_coarRegion) const
  {
    m_implem->fixFineToCoarse(*a_fineEBIS.m_implem, a_coarRegion);
  }
  /*******************************/
  inline Array<FaceIndex>
  EBGraph::getAllFaces(const IntVect& a_iv,
                       const int& a_idir,
                       const Side::LoHiSide& a_sd) const
  {
    return m_implem->getAllFaces(a_iv, a_idir, a_sd);
  }
        
  /*******************************/
        
}
#endif
