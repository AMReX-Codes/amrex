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

#ifndef _GRAPHNODE_H_
#define _GRAPHNODE_H_

#include "AMReX_REAL.H"
#include "AMReX_Box.H"
#include "AMReX_IntVect.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include <cassert>


namespace amrex
{
  ///
  /**
     This is the graph data held for each irregular node.
     Warning: each GraphNode holds a <EM>Vrray</EM> of GraphNodeImplems.
  */
  class GraphNodeImplem
  {
  public:
    ///
    inline
    GraphNodeImplem();
                  
    ///
    inline
    ~GraphNodeImplem();
                  
    ///
    GraphNodeImplem& operator=(const GraphNodeImplem& a_impin);
                  
    ///
    GraphNodeImplem(const GraphNodeImplem& a_impin);
                  
    /// If true, this represents a regular cell and only m_coarserNode will be valid
    bool m_isRegular;
                  
    /// Used to mark a GraphNode as invalid - used in getting connected components of the graph
    bool m_isValid;
                  
    ///
    /**
       Cell indicies of neighboring cells.
       If the arc is a boundary face, the int = -1
    */
    std::array<Vector<int>, 2*SpaceDim> m_arc;
                  
    /// Cell index of this vof is the index into the vector
                  
    /// Cell index of next coarser vof
    int m_coarserNode;
                  
    /// Index into node vector (for construction).
    int m_nodeInd;

    mutable bool m_verbose;

    /// List of finer vofs
    Vector<VolIndex> m_finerNodes;
                  
    /// Return the index into the arc vector
    inline int index(int a_idir, Side::LoHiSide a_side) const;
                  
    ///
    int linearSize() const;
                  
    ///
    void linearOut(void* buffer ) const;
                  
    ///
    void linearIn(void* buffer );
                  
    friend class EBISLevel;
  };
                  
  ///
  /**
     This is a list showing the connectivity of a given cell.  This
     has the special property of also being able to set itself to regular
     or covered by setting the current nodeimplem pointer to 1 or 0.
  */
  class GraphNode
  {
  public:
    ///
    /**
       Constructor sets node to regular as default.
    */
    inline GraphNode();
                  
    ///
    /**
     */
    inline ~GraphNode();
                  
    ///
    int size() const;
                  
    ///
    /**
       Deletes memory if it has been allocated and sets
       members to 0
    */
    inline void clear();
                  
    ///
    /**
       Return true if the node is covered (m_cellList==0).
    */
    inline bool isCovered() const;
                  
    ///
    /**
       Return true if the node is regular.
    */
    inline bool isRegular() const;
                  
    ///
    /**
       Return true if the node is regular and has a single-valued parent
       (m_cellList==1).
    */
    inline bool isRegularWithSingleValuedParent() const;
                  
    ///
    /**
       Return true if the node is regular and has a multi-valued parent
       (m_cellList is a valid pointer)
    */
    inline bool isRegularWithMultiValuedParent() const;
                  
    ///
    /**
       Return true if the node is neither regular or covered.
    */
    inline bool isIrregular() const;
                  
    ///
    /**
       Return true if the node has a valid m_cellList pointer
    */
    inline bool hasValidCellList() const;
                  
    ///
    /**
       Set Node to regular.  If previously set to irregular,
       deletes their memory.
    */
    inline void defineAsRegular();
                  
    ///
    /**
       Set Node to covered.  If previously set to irregular,
       deletes their memory.
    */
    inline void defineAsCovered();
                  
    ///
    /**
       Get the faces in the direction and side
       for the vof in the list.  if the vof's cell
       index is not found in the list, abort.
       Use the input vof's grid index for a the grid index
       of the list.
    */
    Vector<FaceIndex>
    getFaces(const VolIndex&       a_vof,
             const int&            a_idir,
             const Side::LoHiSide& a_sd,
             const Box&  a_domain) const;
                  
    ///
    Vector<FaceIndex>
    getFaces(const IntVect&        a_iv,
             const int&            a_idir,
             const Side::LoHiSide& a_sd,
             const Box&  a_domain) const;
                  
    ///
    /**
       Return all the vofs in the list, using the input
       intvect for the gridIndex
    */
    Vector<VolIndex>
    getVoFs(const IntVect& a_iv) const;
                  
    ///
    const GraphNode& operator=(const GraphNode& ebiin);
                  
    ///
    GraphNode(const GraphNode& ebiin);
                  
    ///
    /**
       Returns the corresponding set of VoFs from the next finer
       EBGraph (factor of two refinement).  The result is only
       defined if this {\tt EBGraph} was defined by coarsening.
    */
    Vector<VolIndex> refine(const VolIndex& a_coarVoF) const;
                  
    ///
    /**
       Returns the corresponding  VoF from the next coarser
       EBGraph (same solution location, different index space, factor
       of two refinement ratio).
    */
    VolIndex coarsen(const VolIndex& a_fineVoF) const;
                  
    ///
    /**
     */
    void addIrregularNode(const GraphNodeImplem& a_nodein, int cellIndex);
    ///
    /**
     */
    void addIrregularNode(const GraphNodeImplem& a_nodein);
                  
    ///
    /**
       Get all sets of connected vofs within the box (used in coarsening)
    */
    Vector<Vector<VolIndex> >  getVoFSets(const Box& a_box) const;
                  
    ///
    int linearSize() const;
                  
    ///
    void linearOut(void* buffer ) const;
                  
    ///
    void linearIn(void* buffer );
                  
    /// internal use only
    inline void setDefaults();
                  
    ///
    /**
       The connectivity data at this point.
       If m_cellList == 0, node is covered.
       If m_cellList == 1, node is regular and parent is single-valued.
       otherwise, it is a real list (node may still be regular).
    */
    Vector<GraphNodeImplem>* m_cellList;
    mutable bool m_verbose;
    friend class EBISLevel;
  };

/* Inline functions */
/*******************************/
  inline GraphNode::~GraphNode()
  {
    clear();
  }

/*******************************/
  inline bool GraphNode::isCovered() const
  {
    return (m_cellList == ((Vector<GraphNodeImplem>*) 0));
  }

/*******************************/
  inline bool GraphNode::isRegular() const
  {
    return (isRegularWithSingleValuedParent() ||
            isRegularWithMultiValuedParent());
  }

/*******************************/
  inline bool GraphNode::isRegularWithSingleValuedParent() const
  {
    return (m_cellList == ((Vector<GraphNodeImplem>*) 1));
  }

/*******************************/
  inline bool GraphNode::isRegularWithMultiValuedParent() const
  {
    return (((m_cellList != ((Vector<GraphNodeImplem>*) 0)) &&
             ((*m_cellList).size() == 1) &&
             ((*m_cellList)[0]).m_isRegular));
  }

/*******************************/
  inline bool GraphNode::hasValidCellList() const
  {
    return (!isCovered() && !isRegularWithSingleValuedParent());
  }

/*******************************/
  inline bool GraphNode::isIrregular() const
  {
    return ((!isRegular()) && (!isCovered()));
  }

/*******************************/
  inline void GraphNode::setDefaults()
  {
    m_verbose = false;
    m_cellList = (Vector<GraphNodeImplem>*)1;
  }

/*******************************/
  inline GraphNode::GraphNode()
  {
    setDefaults();
  }

/*******************************/
  inline void GraphNode::clear()
  {
    if (hasValidCellList())
    {
      delete m_cellList;
    }
    setDefaults();
  }

/*******************************/
  inline void GraphNode::defineAsRegular()
  {
    clear();
    m_cellList = (Vector<GraphNodeImplem>*) 1;
  }

/*******************************/
  inline void GraphNode::defineAsCovered()
  {
    clear();
    m_cellList = (Vector<GraphNodeImplem>*) 0;
  }

/*******************************/
  inline GraphNodeImplem::GraphNodeImplem()
    : m_isRegular(false),
      m_isValid(true),
      m_coarserNode(-1),
      m_nodeInd(-1)
  {
    m_verbose = false;
  }

/*******************************/
  inline GraphNodeImplem::~GraphNodeImplem()
  {
  }

/*******************************/
  inline int GraphNodeImplem::index(int            a_idir,
                                    Side::LoHiSide a_sd) const
  {
    assert(a_idir >= 0 && a_idir < SpaceDim);
    int retval;
    if (a_sd == Side::Lo)
    {
      retval = a_idir;
    }
    else
    {
      retval = a_idir + SpaceDim;
    }
    return retval;
  }


}

#endif
