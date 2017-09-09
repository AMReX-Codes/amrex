
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


#ifndef _EBISBOX_H___
#define _EBISBOX_H___

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_RealVect.H"
#include "AMReX_BaseFab.H"
#include "AMReX_SPMD.H"
#include "AMReX_EBGraph.H"
#include "AMReX_EBData.H"


namespace amrex
{

/// Geometric description within a box
///
/**
   EBISBox represents the geometric information
   of the domain at a given refinement
   within the boundaries of a particular box.
   It is implemented as a reference-counted
   pointer to an EBISBox object.  This means
   that assignment and copying are computationally
   inexpensive but it also means that these functions
   have the reference-counted pointer semantic.  If
   one copies an EBISBox and changes the copy,
   the original one changes too.
*/
  class EBISBox
  {
  public:

    ///
    /**
       Makes an EBISBox whose contents are undefined.
    */
    EBISBox();

    ///
    ~EBISBox();

    ///
    /**
       Returns all the multi-valued cells
       within the input Box a_subbox.
    */
    IntVectSet getMultiCells(const Box& a_subbox) const;

    ///
    /**
       Returns the irregular cells of the EBISBox that are
       within the input subbox.
    */
    IntVectSet getIrregIVS(const Box& a_subbox) const;

    ///
    /**
       Returns the irregular cells that have non-zero boundary area
    */
    IntVectSet boundaryIVS(const Box& a_subbox) const;

    ///
    /**
       Returns the number of VoFs in the cell
    */
    int numVoFs(const IntVect& a_iv) const;

    ///
    /**
       Return the VoFs in the cell.
    */
    Array<VolIndex> getVoFs(const IntVect& a_iv) const;

    ///
    /**
       Return the VoFs in the direction and distance given.
       Chases the connectivity the number of steps given.
       a_steps = 0 just returns the vof input.
       a_steps < 0 is an error.  Stops at the domain
       boundary without error.  Will be an error if the
       number of steps steps out of the EBISBox region
       without stepping out of the domain.  If this happens,
       you probably need to define your EBISLayout with
       more ghost cells.
    */
    Array<VolIndex> getVoFs(const VolIndex& a_vof,
                             const int& a_dir,
                             const Side::LoHiSide& a_sd,
                             const int& a_steps) const;

    ///
    /**
       Return true if every cell in the EBISBox is a regular cell.
    */
    bool isAllRegular() const;

    ///
    /**
       Return true if  every cell in the EBISBox is a covered cell.
    */
    bool isAllCovered() const;

    ///
    /**
       Return true if a_iv is a regular cell.
    */
    bool isRegular(const IntVect& a_iv) const;

    ///
    /**
       Return true if a_iv is an irregular cell.
    */
    bool isIrregular(const IntVect& a_iv) const;

    ///
    /**
       Return true if a_iv is a covered cell.
    */
    bool isCovered(const IntVect& a_iv) const;

    ///
    /**
       Return true if this cell is multi-valued.
    */
    bool isMultiValued(const IntVect& a_iv) const;

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
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    Array<FaceIndex> getFaces(const VolIndex& a_vof,
                               const int& a_idir,
                               const Side::LoHiSide& a_sd) const;

    ///
    Array<FaceIndex> getAllFaces(const IntVect& a_iv,
                                  const int& a_idir,
                                  const Side::LoHiSide& a_sd) const;
    ///
    /**
       Returns a RealVect whose component in the uninteresting
       direction normal to the face is undefined.
       Returns the centroid of input face in the (one or two)
       interesting directions.  Return the
       zero vector if the face is covered or regular.
       The answer is given as a normalized (by grid spacing)
       offset from the center of the cell face
       (all numbers range from -0.5 to 0.5).
    */
    RealVect centroid(const FaceIndex& facein) const;

    ///
    /**
       Return the number of  faces on the side and direction
       of the input VoF.
    */
    int numFaces(const VolIndex& a_vof,
                 const int& a_idir,
                 const Side::LoHiSide& a_sd) const;

    ///
    /**
       Return the volume fraction of teh input VoF.
    */
    Real volFrac(const VolIndex& a_vof) const;

    ///
    /**
       Return the area fraction scaling of the input VoF.
       Currently, this is 1/Max_faces(areaFrac).
    */
    Real areaFracScaling(const VolIndex& a_vof) const;

    ///
    /**
       Return true if the two VoFs are connected.
    */
    bool isConnected(const VolIndex& a_vof1,
                     const VolIndex& a_vof2) const;

    ///
    /**
       Return the area fraction of the face.
       Returns zero if the two vofs in the face
       are not actually connected.
    */
    Real areaFrac(const FaceIndex& a_face1) const;

    ///
    /**
       Return the sum of area fractions of the high
       side of the vof in the given direction.
    */
    Real sumArea(const VolIndex& a_vof,
                 const int& a_idir,
                 const Side::LoHiSide& a_sd) const;

    ///
    /**
       Returns the centroid of the input VoF. Return the zero
       vector if the VoF is regular or covered.  The answer is given
       as a normalized (by grid spacing) offset from the center
       of the cell (all numbers range from -0.5 to 0.5).
    */
    RealVect centroid(const VolIndex& a_vof) const;

    ///
    /**
       Returns the corresponding set of VoFs from the next finer
       EBISBox (factor of two refinement).  The result is only
       defined if this {\tt EBISBox} was defined by coarsening.
    */
    Array<VolIndex> refine(const VolIndex& a_coarVoF) const;

    ///
    /**
       Returns the corresponding  VoF from the next coarser
       EBISBox (same solution location, different index space, factor
       of two refinement ratio).
    */
    VolIndex coarsen(const VolIndex& a_fineVoF) const;

    ///
    /**
     */
    RealVect bndryCentroid(const VolIndex& a_vof) const;

    ///
    /**
     */
    Real bndryArea(const VolIndex& a_vof) const;

    ///
    /**
       Return the normal to the boundary at the
       input VoF.  If said normal is undefined,
       returns the zero vector.
    */
    RealVect normal(const VolIndex& a_vof) const;

    ///
    const Box& getRegion() const;

    ///
    const Box& getDomain() const;

    ///
    void setDomain(const Box& a_domain);

    ///
    /**
       Returns the corresponding set of faces from the next finer
       EBISBox (factor of two refinement).  The result is only
       defined if this  EBISBox was defined by coarsening.
    */
    Array<FaceIndex> refine(const FaceIndex& a_coarFace,const EBISBox& a_fineEBISBox) const;

    ///
    /**
       Returns the corresponding  face from the next coarser
       {\tt EBISLevel} (same solution location,
       different index space, factor
       of two refinement ratio).
    */
    FaceIndex coarsen(const FaceIndex& a_fineFace) const;

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

    /// define from scratch using irregular graph
    void
    define(const BaseFab<int>&      a_regIrregCovered,
           const Array<IrregNode>& a_irregGraph,
           const Box&               a_validRegion,
           const Box&               a_domain);

    ///
    /**
     */
    const EBGraph& getEBGraph() const;

    ///
    /**
     */
    const EBData& getEBData() const;

    /// define from scratch using a graph and an ebdata
    void
    define(const EBGraph&  a_graph,
           const EBData&   a_data);

    EBISBox(const EBGraph&  a_graph,
            const EBData&   a_data)
      {
        define(a_graph, a_data);
      }

    ///pointer copy
    EBISBox& operator=(const EBISBox& ebiin);

    ///pointer copy
    EBISBox(const EBISBox& ebiin);

    ///pointer comparison
    bool operator==(const EBISBox& ebiin);

    EBISBox (EBISBox&& rhs) noexcept = default;

  private:

    ///this is already a ref-counted object
    EBGraph  m_graph;

    ///this is already a ref-counted objecto
    EBData   m_data;

  };



/*******************************/
  inline bool
  EBISBox::isRegular(const IntVect& a_iv) const
  {
    return m_graph.isRegular(a_iv);
  }
/*******************************/
  inline bool
  EBISBox::isIrregular(const IntVect& a_iv) const
  {
    return m_graph.isIrregular(a_iv);
  }
/*******************************/
  inline bool
  EBISBox::isAllCovered() const
  {
    return m_graph.isAllCovered();
  }
/*******************************/
  inline bool
  EBISBox::isAllRegular() const
  {
    return m_graph.isAllRegular();
  }
/*******************************/
  inline bool
  EBISBox::isCovered(const IntVect& a_iv) const
  {
    return m_graph.isCovered(a_iv);
  }

  inline bool
  EBISBox::isMultiValued(const IntVect& a_iv) const
  {
    return m_graph.isMultiValued(a_iv);
  }

/*******************************/
  inline bool
  EBISBox::isCovered(const Box& a_box) const
  {
    return m_graph.isCovered(a_box);
  }
/*******************************/
  inline bool
  EBISBox::isRegular(const Box& a_box) const
  {
    return m_graph.isRegular(a_box);
  }

  inline Array<VolIndex>
  EBISBox::getVoFs(const IntVect& a_iv) const
  {
    return m_graph.getVoFs(a_iv);
  }

/*******************************/
  inline Array<FaceIndex>
  EBISBox::getAllFaces(const IntVect& a_iv,
                       const int& a_idir,
                       const Side::LoHiSide& a_sd) const
  {
    return m_graph.getAllFaces(a_iv, a_idir, a_sd);
  }
/*******************************/
  inline Array<FaceIndex>
  EBISBox::getFaces(const VolIndex& a_vof,
                    const int& a_idir,
                    const Side::LoHiSide& a_sd) const
  {
    return m_graph.getFaces(a_vof, a_idir, a_sd);
  }

  inline int
  EBISBox::numVoFs(const IntVect& a_iv) const
  {
    return m_graph.numVoFs(a_iv);
  }
}
#endif
