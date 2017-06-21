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

#include "AMReX_EBGraph.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_FaceIterator.H"
#include "AMReX_parstream.H"

namespace amrex
{
  Box ebg_debbox1(IntVect(D_DECL(30,14,0)), IntVect(D_DECL(49, 33,0)));
  Box ebg_debbox2(IntVect(D_DECL(14,14,0)), IntVect(D_DECL(33, 33,0)));
  IntVect ebg_debiv(D_DECL(33, 13, 0));
  /*******************************/
  std::vector<FaceIndex> EBGraph::getMultiValuedFaces(const int&  a_idir,
                                                      const Box&  a_box) const
  {
    return
      m_implem->getMultiValuedFaces(a_idir, a_box, *this);
  }
        
        
  /*******************************/
  std::vector<FaceIndex> 
  EBGraphImplem::
  getMultiValuedFaces(const int&     a_idir,
                      const Box&     a_box,
                      const EBGraph& a_ebgraph) const
  {
    std::vector<FaceIndex> multiValuedFaces;
    Box ghostRegion = a_box;
    ghostRegion.grow(a_idir, 1);
    const IntVectSet ivsMulti = a_ebgraph.getMultiCells(ghostRegion);
    //Use faceiterator to only stop at faces once
    FaceIterator faceit(ivsMulti, a_ebgraph, a_idir, FaceStop::SurroundingWithBoundary);
    for (faceit.reset(); faceit.ok(); ++faceit)
    {
      const IntVect& ivHi = faceit().gridIndex(Side::Hi);
      const IntVect& ivLo = faceit().gridIndex(Side::Lo);
      if (a_box.contains(ivLo) || a_box.contains(ivHi))
      {
        bool isMulti = false;
        if (a_ebgraph.getDomain().contains(ivHi))
        {
          const std::vector<FaceIndex> faces = getAllFaces(ivHi, a_idir, Side::Lo);
          isMulti = (faces.size() > 1);
        }
        else if (a_ebgraph.getDomain().contains(ivLo))
        {
          const std::vector<FaceIndex> faces = getAllFaces(ivLo, a_idir, Side::Hi);
          isMulti = (faces.size() > 1);
        }
        else
        {
          amrex::Error("EBGraph::getMultiValuedFaces --  domain does not contain either ivHi or ivLo");
        }
        if (isMulti)
        {
          multiValuedFaces.push_back(faceit());
        }
      }
    }
    return multiValuedFaces;
  }
        
  /*******************************/
  void 
  EBGraph::
  getRegNextToMultiValued(IntVectSet&    a_vofsToChange,
                          const Box &    a_region) const
  {
    m_implem->getRegNextToMultiValued(a_vofsToChange, a_region);
  }
        
  /*******************************/
  void 
  EBGraphImplem::
  getRegNextToMultiValued(IntVectSet&    a_vofsToChange,
                          const Box &    a_region) const
  {
    Box ghostRegion = m_region;
    Box region = getRegion();
    region &= a_region;
    BL_ASSERT(ghostRegion.contains(region));
        
    //loop through multiply-valued vofs in the grown region
    //if any of the vofs next to the multiply-valued vofs
    //are regular, collect them so they can be changed to
    //irregular vofs with unit vol fracs and area fracs.
    a_vofsToChange = IntVectSet();
    IntVectSet multiIVS = getMultiCells(ghostRegion);
    for (VoFIterator vofit(multiIVS, *this); vofit.ok(); ++vofit)
    {
      const VolIndex& multiVoF = vofit();
      const IntVect& multiIV = multiVoF.gridIndex();
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        for (SideIterator sit; sit.ok(); ++sit)
        {
          IntVect otherIV = multiIV + sign(sit())*BASISV(idir);
          if ((region.contains(otherIV) && isRegular(otherIV)))
          {
            a_vofsToChange |= otherIV;
          }
        }
      }
    }
  }
        
  /*******************************/
  void EBGraph::addFullIrregularVoFs(const IntVectSet& a_vofsToChange)
  {
    m_implem->addFullIrregularVoFs(a_vofsToChange);
  }
        
  /*******************************/
  void EBGraphImplem::addFullIrregularVoFs(const IntVectSet& a_vofsToChange)
  {
    if (!a_vofsToChange.isEmpty())
    {
      BL_ASSERT(isDefined());
      BL_ASSERT(isDomainSet());
      //this is all supposed to be called for regular vofs
      //to be changed to full irregular vofs
      BL_ASSERT(!isAllCovered());
      if (isAllRegular())
      {
        m_tag = HasIrregular;
        m_multiIVS = IntVectSet();
        m_irregIVS = IntVectSet();
        m_graph.resize(m_region, 1);
        for (BoxIterator bit(m_region); bit.ok(); ++bit)
        {
          m_graph(bit(), 0).defineAsRegular();
        }
      }
        
      //  //now for changing vofs
      for (IVSIterator ivsit(a_vofsToChange); ivsit.ok(); ++ivsit)
      {
        const IntVect&  iv = ivsit();
        //needs to be a regular cell to start with
        BL_ASSERT(isRegular(iv));
        VolIndex vof(iv, 0);
        
        GraphNode regularNode;
        regularNode.defineAsRegular();
        
        //create node and its arcs
        //the coarse-fine info is created later.
        //this operation must be done before all that
        //the finer ones can be set trivially but not the
        //coarse ones.
        GraphNodeImplem node;
        node.m_finerNodes = regularNode.refine(vof);
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
          {
        
            int gNodeIndex = IrregNode::index(idir, sit());
            std::vector<int>& nodeArcsDir = node.m_arc[gNodeIndex];
            nodeArcsDir.resize(0);
            //find which vof node is connected to in each direction.
            //cannot use isConnected here because it will always return
            //true since one vof is still regular
            IntVect  otherIV = iv + sign(sit())*BASISV(idir);
            if (m_domain.contains(otherIV))
            {
              std::vector<VolIndex> otherVoFs = getVoFs(otherIV);
              bool found = false;
              for (int iother = 0; iother < otherVoFs.size(); iother++)
              {
                const VolIndex& otherVoF = otherVoFs[iother];
                std::vector<FaceIndex> otherFaces = getFaces(otherVoF, idir, flip(sit()));
                //there are rare cases where the number of other faces is greater than 1
                for (int iface = 0; iface < otherFaces.size(); iface++)
                {
                  nodeArcsDir.push_back(otherVoF.cellIndex());
                  found = true;
                }
              }
              if (!found)
              {
                amrex::Error("former regular vof not connected to anything");
              }
            }
            else
            {
              //boundary arc
              nodeArcsDir.resize(1, -1);
            }
        
          }
        }
        //finally add node into graph.
        //again coarse and fine info have to be added on later.
        m_graph(iv, 0).addIrregularNode(node);
        (m_irregIVS) |= iv;
        if (m_graph(iv, 0).size() > 1)
        {
          amrex::Error("that vof was already irregular");
        }
      }
    }
  }
        
  /*******************************/
  IntVectSet EBGraphImplem::getIrregCells(const Box& a_subbox) const
  {
    IntVectSet retval;
    retval = m_irregIVS;
    retval &= a_subbox;
    return retval;
  }
        
  /*******************************/
  IntVectSet EBGraphImplem::getMultiCells(const Box& a_subbox) const
  {
    IntVectSet retval;

    retval = m_multiIVS;
    retval &= a_subbox;
    return retval;
  }
        
  bool EBGraphImplem::isMultiValued(const IntVect& a_iv) const
  {
    bool retval = (m_multiIVS.contains(a_iv));
    return retval;
  }
        
  /*******************************/
  void EBGraphImplem::setToAllRegular()
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    m_tag = AllRegular;
    m_irregIVS = IntVectSet();
    m_multiIVS = IntVectSet();
  }
        
  /*******************************/
  void EBGraphImplem::setToAllCovered()
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    m_tag = AllCovered;
    m_irregIVS = IntVectSet();
    m_multiIVS = IntVectSet();

  }
        
        
  /*******************************/
  void 
  EBGraphImplem::
  buildGraph(const BaseFab<int>          & a_regIrregCovered,
             const std::vector<IrregNode>& a_irregGraph,
             const Box                   & a_validRegion,
             const Box                   & a_domain)
  {
    define(a_validRegion);
    setDomain(a_domain);
        
    //remember that this forces a dense representation.
    m_tag = HasIrregular;
    m_multiIVS = IntVectSet();
    m_irregIVS = IntVectSet();
    m_graph.resize(m_region, 1);
        
    //set regular and covered cells
    for (BoxIterator bit(m_region); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      if (a_regIrregCovered(iv, 0) == 1) //regular cell
      {
        m_graph(iv, 0).defineAsRegular();
      }
      else if (a_regIrregCovered(iv, 0) == -1) //covered cell
      {
        m_graph(iv, 0).defineAsCovered();
      }
      else if (a_regIrregCovered(iv, 0) != 0)
      {
        amrex::Error("invalid flag");
      }
    }
        
    //now for irregular cells
    //add the vofs
    for (int ivecIrreg = 0; ivecIrreg < a_irregGraph.size(); ivecIrreg++)
    {
      const IrregNode& inputNode = a_irregGraph[ivecIrreg];
      const IntVect& iv =inputNode.m_cell;
        
      GraphNodeImplem newImplem;
      newImplem.m_nodeInd = ivecIrreg;
        
      m_graph(iv, 0).addIrregularNode(newImplem, inputNode.m_cellIndex);
      (m_irregIVS) |= iv;
      if (m_graph(iv, 0).size() > 1)
      {
        (m_multiIVS) |= iv;
      }
    }
        
    //add the faces
    for (int ivecIrreg = 0; ivecIrreg < a_irregGraph.size(); ivecIrreg++)
    {
      const IrregNode& inputNode = a_irregGraph[ivecIrreg];
      const IntVect& iv =inputNode.m_cell;
      std::vector<GraphNodeImplem>& vecNodes =
        *(m_graph(iv, 0).m_cellList);
      //pick out which node we are talking about
      //by maching its nodeInd with ivecIrreg
      bool foundNode = false;
      GraphNodeImplem* nodePtr = NULL;
      for (int ivecGraph = 0; ivecGraph < vecNodes.size(); ivecGraph++)
      {
        if (vecNodes[ivecGraph].m_nodeInd == ivecIrreg)
        {
          foundNode = true;
          nodePtr = &(vecNodes[ivecGraph]);
        }
      }
      if (!foundNode)
      {
        amrex::Error("EBGraph: internal error in construction");
      }
      //now add the arcs in the input to the node
      GraphNodeImplem& node = *nodePtr;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        for (SideIterator sit; sit.ok(); ++sit)
        {
          int irregIndex = IrregNode::index(idir, sit());
          int gNodeIndex = IrregNode::index(idir, sit());
          const std::vector<int>& irregArcs = inputNode.m_arc[irregIndex];
          std::vector<int>& nodeArcs = node.m_arc[gNodeIndex];
          nodeArcs.resize(irregArcs.size());
          for (int iarc = 0; iarc < irregArcs.size(); iarc++)
          {
            int otherNodeInd = irregArcs[iarc];
            if (otherNodeInd == -1)
            {
              //if otherNodeInd == -1, boundary arc.
              //just make the arc in our node = -1
              //to signify the same
              nodeArcs[iarc] = -1;
            }
            else if (otherNodeInd == -2)
            {
              //this means that the vof is connected
              //to a regular vof,
              //which always have a cell index of 0
              nodeArcs[iarc] =  0;
            }
            else
            {
              nodeArcs[iarc] =  otherNodeInd;
            }
          }
        }
      }
        
    }
  }
        
  /*******************************/
  const Box& EBGraphImplem::getRegion() const
  {
    return m_region;
  }
        
  /*******************************/
  const Box& EBGraphImplem::getDomain() const
  {
    return m_domain;
  }
        
  /*******************************/
  EBGraphImplem::EBGraphImplem(const Box& a_box)
    :m_isMaskBuilt(false)
  {
    define(a_box);
  }
        
  /*******************************/
  void EBGraphImplem::define(const Box& a_region)
  {
    BL_ASSERT(!a_region.isEmpty());
        
    m_tag = AllRegular;
    m_region = a_region;
    m_irregIVS = IntVectSet();
    m_multiIVS = IntVectSet();
    m_mask.clear();
    m_isMaskBuilt = false;
    m_isDefined= true;
  }
        
  /*******************************/
  void EBGraphImplem::setDomain(const Box& a_domain)
  {
    m_isDomainSet = true;
    m_domain = a_domain;
  }
        
  /*******************************/
  EBGraphImplem::EBGraphImplem()
  {
    m_isDefined = false;
    m_isDomainSet = false;
  }
        
  /*******************************/
  EBGraphImplem::~EBGraphImplem()
  {
  }
        
  /*******************************/
  bool EBGraphImplem::isDefined() const
  {
    return m_isDefined;
  }
        
  /*******************************/
  bool EBGraphImplem::isDomainSet() const
  {
    return m_isDomainSet;
  }
        
  /*******************************/
  const BaseFab<char>& 
  EBGraphImplem::
  getMask(int& a_regIrregCovered) const
  {
    if (this->isAllRegular())
    {
      a_regIrregCovered = 1;
    }
    else if (this->isAllCovered())
    {
      a_regIrregCovered = -1;
    }
    else
    {
      a_regIrregCovered = 0;
      if (!m_isMaskBuilt)
      {
        Box maskBox = m_region & m_domain;
        m_mask.resize(maskBox, 1);
        fillMask(m_mask);
        m_isMaskBuilt = true;
      }
    }
    return m_mask;
  }
        
  /*******************************/
  std::vector<VolIndex> 
  EBGraphImplem::
  getVoFs(const IntVect& a_iv) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    std::vector<VolIndex> retvec;
    if (m_tag == AllRegular)
    {
      retvec.push_back(VolIndex(a_iv, 0));
    }
    else if (m_tag == AllCovered)
    {
      //return an empty vector
    }
    else if (m_tag == HasIrregular)
    {
      BL_ASSERT(m_region.contains(a_iv));
      BL_ASSERT(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retvec = node.getVoFs(a_iv);
    }
    return retvec;
  }
        
  /*******************************/
  bool 
  EBGraphImplem::
  isRegular(const IntVect& a_iv) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    bool retval;
    if (m_tag == AllRegular)
    {
      retval = true;
    }
    else if (m_tag == AllCovered)
    {
      retval = false;
    }
    else if (m_tag == HasIrregular)
    {
      //BL_ASSERT(m_region.contains(a_iv)); //picked up my m_graph already
      //BL_ASSERT(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retval = node.isRegular();
    }
    else
    {
      retval = false;
      amrex::Error("EBGraphImplem::isRegular:Bogus Tag");
    }
    return retval;
  }
        
  /*******************************/
  std::vector<FaceIndex> 
  EBGraphImplem::
  getAllFaces(const IntVect&        a_iv,
              const int&            a_idir,
              const Side::LoHiSide& a_sd) const
  {
    std::vector<FaceIndex> retval(0);
    std::vector<VolIndex> vofs = getVoFs(a_iv);
    for (int ivof= 0; ivof < vofs.size(); ivof++)
    {
      std::vector<FaceIndex> faces = getFaces(vofs[ivof], a_idir, a_sd);
      for(int iface = 0; iface < faces.size(); iface++)
      {
        retval.push_back(faces[iface]);
      }
    }
    return retval;
  }
        
  /*******************************/
  bool 
  EBGraphImplem::
  isIrregular(const IntVect& a_iv) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    bool retval;
    if (m_tag == AllRegular)
    {
      retval = false;
    }
    else if (m_tag == AllCovered)
    {
      retval = false;
    }
    else if (m_tag == HasIrregular)
    {
      BL_ASSERT(m_region.contains(a_iv));
      BL_ASSERT(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retval = node.isIrregular();
    }
    else
    {
      retval = false;
      amrex::Error("EBGraphImplem::isIrregular:Bogus Tag");
    }
    return retval;
  }
        
  /*******************************/
  bool EBGraphImplem::isAllCovered() const
  {
    return (m_tag == AllCovered);
  }
        
  /*******************************/
  bool EBGraphImplem::isAllRegular() const
  {
    return m_tag == AllRegular;
  }
        
  /*******************************/
  bool EBGraph::hasIrregular() const
  {
    return m_implem->hasIrregular();
  }
        
  /*******************************/
  bool EBGraphImplem::hasIrregular() const
  {
    return m_tag == HasIrregular;
  }
        
  /*******************************/
  bool 
  EBGraphImplem::
  isCovered(const IntVect& a_iv) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    bool retval;
    if (m_tag == AllRegular)
    {
      retval = false;
    }
    else if (m_tag == AllCovered)
    {
      retval = true;
    }
    else if (m_tag == HasIrregular)
    {
      //BL_ASSERT(m_region.contains(a_iv)); this check picked up by m_graph
      //BL_ASSERT(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retval = node.isCovered();
    }
    else
    {
      retval = false;
      amrex::Error("EBGraphImplem::isIrregular:Bogus Tag");
    }
    return retval;
  }
        
  /*******************************/
  bool 
  EBGraphImplem::
  isCovered(const Box& a_box) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    bool retval;
    if (m_tag == AllRegular)
    {
      retval = false;
    }
    else if (m_tag == AllCovered)
    {
      retval = true;
    }
    else if (m_tag == HasIrregular)
    {
      BL_ASSERT(m_region.contains(a_box));
      BL_ASSERT(m_domain.contains(a_box));
      retval = true;
      BoxIterator bit(a_box);
      for (bit.reset(); bit.ok() && retval; ++bit)
      {
        if (!isCovered(bit()))
        {
          retval = false;
        }
      }
    }
    else
    {
      retval = false;
      amrex::Error("EBGraphImplem::isCovered:Bogus Tag");
    }
    return retval;
  }
        
  /*******************************/
  bool 
  EBGraphImplem::
  isRegular(const Box& a_box) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    bool retval;
    if (m_tag == AllRegular)
    {
      retval = true;
    }
    else if (m_tag == AllCovered)
    {
      retval = false;
    }
    else if (m_tag == HasIrregular)
    {
      Box region = m_region & a_box;
      region &= m_domain;
      
      retval = true;
      for (BoxIterator bit(region); bit.ok() ; ++bit)
      {
        if (!isRegular(bit()))
        {
          retval = false;
          break;
        }
      }
    }
    else
    {
      retval = false;
      amrex::Error("EBGraphImplem::isRegular:Bogus Tag");
    }
    return retval;
  }
        
  /*******************************/
  std::vector<FaceIndex> 
  EBGraphImplem::
  getFaces(const VolIndex&       a_vof,
           const int&            a_idir,
           const Side::LoHiSide& a_sd) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    BL_ASSERT(m_region.contains(a_vof.gridIndex()));
    BL_ASSERT(m_domain.contains(a_vof.gridIndex()));
    BL_ASSERT((a_idir >= 0) && (a_idir < SpaceDim));
    BL_ASSERT((a_sd == Side::Lo) || (a_sd == Side::Hi));
        
    std::vector<FaceIndex> retvec;
    if (m_tag == AllRegular)
    {
      IntVect otherIV = a_vof.gridIndex()
        + sign(a_sd)*BASISV(a_idir);
      int otherCellInd = 0;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
      }
      VolIndex otherVoF(otherIV, otherCellInd);
      retvec.push_back(FaceIndex(a_vof, otherVoF, a_idir));
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const IntVect& iv = a_vof.gridIndex();
      const GraphNode& node = m_graph(iv, 0);
      retvec = node.getFaces(a_vof, a_idir, a_sd, m_domain);
    }
        
    return retvec;
  }
  /*******************************/
  void EBGraph::fillMask(BaseFab<char>& a_mask) const
  {
    m_implem->fillMask(a_mask);
  }
  /*******************************/
  void EBGraph::fillIntMask(BaseFab<int>& a_mask) const
  {
    Box interiorRegion = getRegion();
    interiorRegion &= getDomain();
    a_mask.resize(interiorRegion, 1);
    for(BoxIterator boxit(interiorRegion); boxit.ok(); ++boxit)
    {
      const IntVect& iv = boxit();
      if(isRegular(iv))
      {
        a_mask(iv, 0) = 1;
      }
      else if(isCovered(iv))
      {
        a_mask(iv, 0) = -1;
      }
      else
      {
        a_mask(iv, 0) = 0;
      }
    }
  }
        
  /*******************************/
  void 
  EBGraph::
  fillCellTypeMask(BaseFab<char>& a_mask) const
  {
    m_implem->fillCellTypeMask(a_mask);
  }
        
        
  /*******************************/
  void 
  EBGraphImplem::
  fillMask(BaseFab<char>& a_mask) const
  {
    if (m_tag == AllRegular )
    {
      a_mask.setVal(1);
    }
    else if (m_tag == AllCovered )
    {
      a_mask.setVal(0);
    }
    else
    {
      Box b = a_mask.box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        if (isCovered(bit())|| isMultiValued(bit()))
        {
          a_mask(bit(), 0) = 0;
        }
        else 
        {
          a_mask(bit(), 0) = 1;
        }
      }
                  
    }
  }
  /*******************************/
  void 
  EBGraphImplem::
  fillCellTypeMask(BaseFab<char>& a_mask) const
  {
    // This is added for VisIt support.
    //  0       outside priblem domain
    //  1       covered
    //  2       regular
    //  3       irregular
    if (m_tag == AllCovered )
    {
      a_mask.setVal(1);
    }
    else
      if (m_tag == AllRegular )
      {
        a_mask.setVal(2);
      }
      else
      {
        Box b = a_mask.box() & m_graph.box();
        if (b.isEmpty()) return;
        for (BoxIterator bit(b); bit.ok(); ++bit)
        {
          if (isCovered(bit()))
          {
            a_mask(bit(), 0) = 1;
          }
          else
            if (isRegular(bit()))
            {
              a_mask(bit(), 0) = 2;
            }
            else
              if (isIrregular(bit()))
              {
                a_mask(bit(), 0) = 3;
              }
              else
              {
                a_mask(bit(), 0) = 0;
              }
        }
      }
  }
        
  std::vector<FaceIndex> EBGraph:: getIrregFaces(const Box& a_box,
                                                 int        a_dir) const
  {
    return m_implem->getIrregFaces(a_box, a_dir);
  }
  /***********/      
  std::vector<FaceIndex>  
  EBGraphImplem::
  getIrregFaces(const Box& a_box,
                int        a_dir) const
  {
    std::vector<FaceIndex> faces;
    if (m_tag == AllRegular || m_tag == AllCovered)
    {
      // do nothing
    }
    else
    {
      Box b = a_box & m_domain;
      IntVectSet ivs = this->getIrregCells(b);
      for (IVSIterator it(ivs); it.ok(); ++it)
      {
        const GraphNode& node = m_graph(it(), 0);
        IntVect hi = it() + BASISV(a_dir);
        BL_ASSERT(node.isIrregular());
        if (b.contains(hi) && m_graph(hi, 0).isIrregular())
        {
          std::vector<FaceIndex> faces = node.getFaces(it(), a_dir, Side::Hi, m_domain);
          for(int iface = 0; iface < faces.size(); iface++)
          {
            faces.push_back(faces[iface]);
          }
        }
      }
    }
    return faces;
  }
        
  /*******************************/
  FaceIndex 
  EBGraphImplem::
  coarsen(const FaceIndex& a_fineFace) const
  {
    VolIndex loVoFCoar = coarsen(a_fineFace.getVoF(Side::Lo));
    VolIndex hiVoFCoar = coarsen(a_fineFace.getVoF(Side::Hi));
    return FaceIndex(loVoFCoar, hiVoFCoar, a_fineFace.direction());
  }
        
  /*******************************/
  std::vector<FaceIndex> 
  EBGraphImplem::
  refine(const FaceIndex&     a_coarFace,
         const EBGraphImplem& a_fineGraph) const
  {
    std::vector<FaceIndex> retval;
    IntVect ivLoCoar = a_coarFace.gridIndex(Side::Lo);
    IntVect ivHiCoar = a_coarFace.gridIndex(Side::Hi);
    int direction = a_coarFace.direction();
    if (m_region.contains(ivLoCoar) && m_region.contains(ivHiCoar))
    {
      //interior face
      std::vector<VolIndex> loVoFFine = refine(a_coarFace.getVoF(Side::Lo));
      std::vector<VolIndex> hiVoFFine = refine(a_coarFace.getVoF(Side::Hi));
      int idir = a_coarFace.direction();
      for (int ilo = 0; ilo < loVoFFine.size(); ilo++)
      {
        for (int ihi = 0; ihi < hiVoFFine.size(); ihi++)
        {
          if (a_fineGraph.isConnected(loVoFFine[ilo], hiVoFFine[ihi]))
          {
            FaceIndex newFace(loVoFFine[ilo], hiVoFFine[ihi], idir);
            retval.push_back(newFace);
          }
        }
      }
    }
    else if (m_region.contains(ivLoCoar))
    {
      //boundary face on the hi side of the domain
      std::vector<VolIndex> loVoFsFine = refine(a_coarFace.getVoF(Side::Lo));
      Box fineRegion = m_region;
      fineRegion.refine(2);
      for (int ivof = 0; ivof < loVoFsFine.size(); ivof++)
      {
        VolIndex& loVoFFine = loVoFsFine[ivof];
        IntVect ivHiFine = loVoFFine.gridIndex() + BASISV(direction);
        if (!fineRegion.contains(ivHiFine))
        {
          std::vector<FaceIndex> fineFaces = a_fineGraph.getFaces(loVoFFine, direction, Side::Hi);
          retval.insert(retval.end(), fineFaces.begin(), fineFaces.end());
        }
      }
    }
    else if (m_region.contains(ivHiCoar))
    {
      //boundary face on the low side of the domain
      std::vector<VolIndex> hiVoFsFine = refine(a_coarFace.getVoF(Side::Hi));
      Box fineRegion = m_region;
      fineRegion.refine(2);
        
      for (int ivof = 0; ivof < hiVoFsFine.size(); ivof++)
      {
        VolIndex& hiVoFFine = hiVoFsFine[ivof];
        IntVect ivLoFine = hiVoFFine.gridIndex() - BASISV(direction);
        if (!fineRegion.contains(ivLoFine) )
        {
          std::vector<FaceIndex> fineFaces = a_fineGraph.getFaces(hiVoFFine, direction, Side::Lo);
          retval.insert(retval.end(), fineFaces.begin(), fineFaces.end());
        }
      }
    }
    else
    {
      amrex::Error("neither vof of the face is inside the domain");
    }
    return retval;
  }
        
  /*******************************/
  bool 
  EBGraphImplem::
  isConnected(const VolIndex& a_vof1,
              const VolIndex& a_vof2) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(m_isDomainSet);
    const IntVect& iv1 = a_vof1.gridIndex();
    const IntVect& iv2 = a_vof2.gridIndex();
        
    BL_ASSERT(m_region.contains(iv1));
    BL_ASSERT(m_region.contains(iv2));
    VolIndex vofLo, vofHi;
        
    if (iv1 > iv2)
    {
      vofLo = a_vof2;
      vofHi = a_vof1;
    }
    else
    {
      vofLo = a_vof1;
      vofHi = a_vof2;
    }
    int direction;
    bool dirfound;
    const IntVect& ivLo = vofLo.gridIndex();
    const IntVect& ivHi = vofHi.gridIndex();
        
    dirfound = false;
    for (int idir = 0; ((idir<SpaceDim) && !dirfound); idir++)
    {
      if ((ivHi - ivLo) == BASISV(idir))
      {
        direction = idir;
        dirfound = true;
      }
    }
    //if not neigboring intvects, no way it can be connected
    if (!dirfound) return false;
        
    std::vector<FaceIndex> faces = getFaces(vofLo, direction, Side::Hi);
    bool voffound = false;
    for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      if (face.getVoF(Side::Hi) == vofHi)
      {
        voffound = true;
      }
    }
    return voffound;
  }
        
        
  /*******************************/
  VolIndex 
  EBGraphImplem::
  coarsen(const VolIndex& a_fineVoF) const
  {
    VolIndex retval;
    IntVect ivfine = a_fineVoF.gridIndex();
    IntVect ivcoar = ivfine;
    ivcoar.coarsen(2);
    if (!m_domain.contains(ivfine))
    {
      //boundary vof
      retval = VolIndex(ivcoar, -1);
    }
    else if ((m_tag == AllRegular ) || (m_tag == AllCovered ))
    {
      retval = VolIndex(ivcoar, 0);
    }
    else
    {
      BL_ASSERT(m_tag == HasIrregular);
      const IntVect& iv = a_fineVoF.gridIndex();
      const GraphNode& node = m_graph(iv, 0);
      retval = node.coarsen(a_fineVoF);
    }
    return retval;
  }
        
  /*******************************/
  std::vector<VolIndex> 
  EBGraphImplem::
  refine(const VolIndex& a_coarVoF) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
        
    std::vector<VolIndex> retval(0);
    IntVect ivCoar = a_coarVoF.gridIndex();
    BL_ASSERT(m_domain.contains(ivCoar));
        
    if (m_tag == AllRegular)
    {
      const IntVect& iv = a_coarVoF.gridIndex();
      Box refbox = Box(iv,iv);
      refbox.refine(2);
      BoxIterator bit(refbox);
      for (bit.reset(); bit.ok(); ++bit)
      {
        retval.push_back(VolIndex(bit(), 0));
      }
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else
    {
      BL_ASSERT(m_tag == HasIrregular);
      const IntVect& iv = a_coarVoF.gridIndex();
      const GraphNode& node = m_graph(iv, 0);
      retval = node.refine(a_coarVoF);
    }
    return retval;
  }
        
        
  /*******************************/
  EBGraphImplem& 
  EBGraphImplem::
  copy(const EBGraphImplem& a_source,
       const Box&           a_srcbox,
       int                  a_srccomp,
       const Box&           a_destbox,
       int                  a_destcomp,
       int                  a_numcomp)
  {
    BL_ASSERT(isDefined());

    Box testbox = m_region & a_srcbox;
    if(!testbox.isEmpty())
    {
      setDomain(a_source.m_domain);
      m_region &= m_domain;

      Box regionTo  = testbox;
      Box regionFrom= testbox;
      if (isRegular(regionTo) && a_source.isRegular(regionFrom))
      {
        return *this;
      }
      else if (isCovered(regionTo) && a_source.isCovered(regionFrom))
      {
        return *this;
      }
      else if (a_source.isCovered(regionFrom) && regionTo.contains(m_region))
      {
        setToAllCovered();
        return *this;
      }
      else if (a_source.isRegular(regionFrom) && regionTo.contains(m_region))
      {
        setToAllRegular();
        return *this;
      }
      else if (isAllRegular() && a_source.isAllCovered())
      {
        //define the basefab as all regular and set the region to
        //covered in the intersection
        m_tag = HasIrregular;
        m_multiIVS = IntVectSet();
        m_irregIVS = IntVectSet();
        m_graph.resize(m_region, 1);
        GraphNode regularNode;
        regularNode.defineAsRegular();
        m_graph.setVal(regularNode);
        Box interBox = m_region & regionTo;
        for (BoxIterator bit(interBox); bit.ok(); ++bit)
        {
          m_graph(bit(), 0).defineAsCovered();
        }
      }
      else if (isAllCovered() && a_source.isAllRegular())
      {
        //define the basefab as all covered and set the region to
        //regular in the intersection
        m_tag = HasIrregular;
        m_multiIVS = IntVectSet();
        m_irregIVS = IntVectSet();
        m_graph.resize(m_region, 1);
        GraphNode  coveredNode;
        coveredNode.defineAsCovered();
        m_graph.setVal(coveredNode);
        Box interBox = m_region & regionTo;
        for (BoxIterator bit(interBox); bit.ok(); ++bit)
        {
          m_graph(bit(), 0).defineAsRegular();
        }
      }
      else
      {
        //one or both has irregular cells.
        //because i am sick of combinatorics,
        //use basefab copy to transfer the data.
        
        //see if we need to generate a source fab
        BaseFab<GraphNode>* srcFabPtr;
        bool needToDelete;
        if (a_source.hasIrregular())
        {
          srcFabPtr = (BaseFab<GraphNode>*)&a_source.m_graph;
          needToDelete = false;
        }
        else
        {
          needToDelete = true;
          srcFabPtr = new BaseFab<GraphNode>(regionFrom, 1);
          GraphNode srcVal;
          if (a_source.isAllRegular())
          {
            srcVal.defineAsRegular();
          }
          else
          {
            //this really has to be true
            BL_ASSERT(a_source.isAllCovered());
            srcVal.defineAsCovered();
          }
          srcFabPtr->setVal(srcVal);
        }
        //define our graph if i need to. leave alone otherwise
        if (isAllRegular() || isAllCovered())
        {
          m_multiIVS = IntVectSet();
          m_irregIVS = IntVectSet();
          m_graph.resize(m_region, 1);
        }
        
        //copy the data
        m_tag = HasIrregular;
///do the slow one for debugging
///      m_graph.slowCopy(*srcFabPtr, regionFrom, 0,regionTo, 0, 1);
        m_graph.copy(*srcFabPtr, regionFrom, 0,regionTo, 0, 1);
        
        //if we needed to new the basefab, clean it up
        if (needToDelete)
        {
          delete srcFabPtr;
        }
      }
      //  now fix up the IntVectSets to match the information
      if (a_source.hasIrregular())
      {
        IntVectSet ivsInterIrreg = (a_source.m_irregIVS);
        ivsInterIrreg &= regionTo;
        ivsInterIrreg &= m_region;
        
        if (!ivsInterIrreg.isEmpty())
        {
          for (IVSIterator it(ivsInterIrreg); it.ok(); ++it)
          {
            IntVect iv = it();
            (m_irregIVS) |= iv;
            if (numVoFs(iv) > 1) // this will be correct since we already
              // did a m_graph copy operation
            {
              (m_multiIVS) |= iv;
            }
          }
        }
      }
    }
    return *this;
  }
        
  /*******************************/
  void 
  EBGraphImplem::
  coarsenVoFs(const EBGraphImplem& a_fineGraph,
              const Box&           a_coarRegion)
  {
    //this also defines the boxes
//    m_region = a_coarRegion;
    m_domain = a_fineGraph.getDomain();
    m_domain.coarsen(2);
    m_region &= m_domain;
    m_isDomainSet = true;
    m_isDefined = true;
    Box refRegion = a_coarRegion;
    refRegion.refine(2);
    if (a_fineGraph.isCovered(refRegion))
    {
      m_tag = AllCovered;
    }
    else if (a_fineGraph.isRegular(refRegion))
    {
      m_tag = AllRegular;
    }
    else
    {
      m_tag = HasIrregular;
      m_multiIVS = IntVectSet();
      m_irregIVS = IntVectSet();
      m_graph.resize(m_region, 1);
      for (BoxIterator bit(a_coarRegion); bit.ok(); ++bit)
      {
        const IntVect& iv = bit();
        Box fineBox = Box(iv, iv);
        fineBox.refine(2);

        if (a_fineGraph.isRegular(fineBox))
        {
          m_graph(bit(), 0).defineAsRegular();
        }
        else if (a_fineGraph.isCovered(fineBox))
        {
          m_graph(bit(), 0).defineAsCovered();
        }
        else
        {
          //get sets of all connected vofs in the box
          //the number of sets will be the number of
          //vofs in the coarse cell
        
          std::vector<std::vector<VolIndex> > fineVoFSets
            = a_fineGraph.getVoFSets(fineBox);
        
          for (int iset = 0; iset < fineVoFSets.size(); iset++)
          {
            GraphNodeImplem newImplem;
            newImplem.m_finerNodes = fineVoFSets[iset];

            m_graph(bit(), 0).addIrregularNode(newImplem);
            (m_irregIVS) |= bit();
            if (m_graph(bit(), 0).size() > 1)
            {
              (m_multiIVS) |= bit();
            }
          }
        }
      }
    }
  }
        
  /*******************************/
  std::vector<std::vector<VolIndex> >  
  EBGraphImplem::
  getVoFSets(const Box& a_region) const
  {
    std::vector<std::vector<VolIndex> > retval;
    //gather all  vofs
    std::vector<VolIndex> allVoFs;
    for (BoxIterator bit(a_region); bit.ok(); ++bit)
    {
      std::vector<VolIndex> newVoFs = getVoFs(bit());
      allVoFs.insert(allVoFs.begin(), newVoFs.begin(), newVoFs.end());
    }
        
    std::vector<bool> beenAdded(allVoFs.size(), false);
    for (int ivof = 0; ivof < allVoFs.size(); ivof++)
    {
      if (!beenAdded[ivof])
      {
        const VolIndex& thisVoF = allVoFs[ivof];
        //we have a vof to start with.
        //now find all vofs connected to this inside
        //this cell----not necessarily directly---
        //hence the inner loop of kvof.
        std::vector<VolIndex> thisVoFSet(1, thisVoF);
        beenAdded[ivof] = true;
        bool doneAdding = false;
        while (!doneAdding)
        {
          int oldSetSize = thisVoFSet.size();
          //i can start at ivof+1 here because 0 to ivof
          //has always been added.
          for (int jvof = ivof+1; jvof < allVoFs.size(); jvof++)
          {
            const VolIndex& testVoF = allVoFs[jvof];
            //this might be a tad confusing to people because
            //the length of the vector can be changing inside
            //the loop.  fortran would puke horribly here.
            for (int kvof = 0; kvof < thisVoFSet.size(); kvof++)
            {
              //because of the nature of push_back,
              //the vector's pointer can change inside so beware.
              //need to hold the vof by value here
              VolIndex innerVoF = thisVoFSet[kvof];
              if (!beenAdded[jvof] && isConnected(innerVoF, testVoF))
              {
                thisVoFSet.push_back(testVoF);
                beenAdded[jvof] = true;
              }
            }
          }
          doneAdding = (thisVoFSet.size() == oldSetSize);
        } //end while !doneadding
        //add the finished bunch to the list of  vofs.
        retval.push_back(thisVoFSet);
      } //end this bunch of connected vofs
    }
        
    //can't hurt to check in case some wacky case was missed.
    for (int ivof = 0; ivof < allVoFs.size(); ivof++)
    {
      if (!beenAdded[ivof])
      {
        amrex::Error("We seem to have missed a vof in getVoFSets");
      }
    }
        
    return retval;
  }
        
  /*******************************/
  std::vector<int> 
  EBGraphImplem::
  coarsenFaces(const VolIndex&       a_coarVoF,
               const EBGraphImplem&  a_fineGraph,
               const int&            a_idir,
               const Side::LoHiSide& a_sd)
  {
    BL_ASSERT(m_isDomainSet);
    BL_ASSERT(a_fineGraph.isDomainSet());
    std::vector<int> retval;
        
    IntVect coarIV = a_coarVoF.gridIndex();
    IntVect otherIV= coarIV + sign(a_sd)*BASISV(a_idir);
    std::vector<VolIndex> theseFineVoFs = refine(a_coarVoF);
    EBGraphImplem&  coarGhostGraph  = *this;  //trying to reduce communcation in this
        
    if (m_domain.contains(otherIV))
    {
        
      //interior faces.
      //just get all possible vofs to connect to and
      //check connectivity
      std::vector<VolIndex> otherCoarVoFs = coarGhostGraph.getVoFs(otherIV);
      for (int iotherCoar = 0; iotherCoar < otherCoarVoFs.size(); iotherCoar++)
      {
        const VolIndex& otherCoarVoF = otherCoarVoFs[iotherCoar];
        std::vector<VolIndex> otherFineVoFs = coarGhostGraph.refine(otherCoarVoF);
        bool addThisFace = false;
        for (int iotherFine = 0; iotherFine < otherFineVoFs.size(); iotherFine++)
        {
          const VolIndex& otherFineVoF = otherFineVoFs[iotherFine];
          for (int ithisFine = 0; ithisFine < theseFineVoFs.size(); ithisFine++)
          {
            const VolIndex& thisFineVoF = theseFineVoFs[ithisFine];
            //                  int thesetwocon = -1;
            if (a_fineGraph.isConnected(thisFineVoF, otherFineVoF))
            {
              addThisFace = true;
              //                      thesetwocon = 1;
            }
          }
        }
        if (addThisFace)
        {
          retval.push_back(iotherCoar);
        }
      }
    }
    else
    {
        
      //boundary faces.
      //if there are any boundary faces on the fine level,
      //make one here too
      bool hasBoundaryFaceFine = false;
      for (int ithis = 0; ithis < theseFineVoFs.size() && !hasBoundaryFaceFine ; ithis++)
      {
        const VolIndex& thisVoF = theseFineVoFs[ithis];
        std::vector<FaceIndex> fineFaces =
          a_fineGraph.getFaces(thisVoF,a_idir, a_sd);
        for (int iface = 0; iface < fineFaces.size() && !hasBoundaryFaceFine; iface++)
        {
          if (fineFaces[iface].isBoundary())
          {
            hasBoundaryFaceFine = true;
          }
        }
      }
      if (hasBoundaryFaceFine)
      {
        //remember that -1 is the secret code for boundary face arcs
        retval.push_back(-1);
      }
    }
        
    return retval;
  }
        
  /*******************************/
  void 
  EBGraphImplem::
  coarsenFaces(const EBGraphImplem& a_fineGraph,
               const Box& a_coarRegion)
  {
    BL_ASSERT(m_isDomainSet);
    BL_ASSERT(a_fineGraph.isDomainSet());
    if (hasIrregular())
    {
      Box region = m_region;
      region &= a_coarRegion;
      for (BoxIterator bit(region); bit.ok(); ++bit)
      {
        if (isIrregular(bit()))
        {
          std::vector<VolIndex> vofsCoar = getVoFs(bit());
          std::vector<GraphNodeImplem>& nodes =
            *(m_graph(bit(), 0).m_cellList);
          for (int ivof = 0; ivof < vofsCoar.size(); ivof++)
          {
            const VolIndex& vofCoar= vofsCoar[ivof];
            GraphNodeImplem& node = nodes[vofCoar.cellIndex()];
            for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
              {
                std::vector<int> coarArcs =
                  coarsenFaces(vofsCoar[ivof],
                               a_fineGraph,
                               idir, sit());
        
                int nodeind = IrregNode::index(idir, sit());
                node.m_arc[nodeind] = coarArcs;
              }
            }
          }
        }
      }
    }
  }
        
  /*******************************/
  void EBGraphImplem::fixFineToCoarse(EBGraphImplem& a_fineGraph, const Box& a_coarRegion) const
  {
    if (hasIrregular())
    {
      for (BoxIterator bit(a_coarRegion); bit.ok(); ++bit)
      {
        if (isIrregular(bit()))
        {
          const IntVect& ivCoar = bit();
        
          const std::vector<GraphNodeImplem>&
            vofsCoar = *(m_graph(ivCoar, 0).m_cellList);
        
          int numVofsCoar = vofsCoar.size();
        
          for (int icoar = 0; icoar < numVofsCoar; icoar++)
          {
            VolIndex vofCoar(ivCoar, icoar);
            std::vector<VolIndex> vofsFine = refine(vofCoar);
        
            for (int ifine = 0; ifine < vofsFine.size(); ifine++)
            {
              const IntVect& ivFine = vofsFine[ifine].gridIndex();
        
              if (a_fineGraph.isIrregular(ivFine))
              {
                int cellIndexFine = vofsFine[ifine].cellIndex();
        
                std::vector<GraphNodeImplem>&
                  nodesFine = *(a_fineGraph.m_graph(ivFine, 0).m_cellList);
        
                nodesFine[cellIndexFine].m_coarserNode = icoar;
              }
              else if ((numVofsCoar > 1) && (a_fineGraph.isRegular(ivFine)))
              {
                GraphNode& nodeFine = a_fineGraph.m_graph(ivFine, 0);
        
                nodeFine.m_cellList = new std::vector<GraphNodeImplem>(1);
        
                (*(nodeFine.m_cellList))[0].m_isRegular   = true;
                (*(nodeFine.m_cellList))[0].m_coarserNode = icoar;
              }
            }
          }
        }
      }
    }
  }
        
  /*******************************/
  EBGraph::EBGraph(const Box& a_box, int a_comps, bool alloc, bool shared)
    : m_implem( std::shared_ptr<EBGraphImplem>( new EBGraphImplem(a_box) ) )
  {
  }
        
  /*******************************/
  EBGraph::EBGraph(const Box& a_box)
    :  m_implem( std::shared_ptr<EBGraphImplem>( new EBGraphImplem(a_box) ) )
  {
  }
        
  /*******************************/
  EBGraph::EBGraph(const EBGraph& a_ebiin)
    :  m_implem( a_ebiin.m_implem )
  {
  }
        
  /*******************************/
  EBGraph::EBGraph()
    : m_implem( std::shared_ptr<EBGraphImplem>( new EBGraphImplem() ) )
  {
  }
        
  /*******************************/
  EBGraph::~EBGraph()
  {
  }
        
        
  /*******************************/
  void EBGraph::define(const Box& a_box)
  {
    m_implem->define(a_box);
  }
  /*******************************/
  /// below lies serialization land.  enter at thy own risk
  ///management is not responsible for any gibbering madness resulting 
  ///from ignoring this warning.
  /*******************************/
  std::size_t 
  EBGraphImplem::
  nBytesFull() const
  {
    //first the region and domain
    BL_PROFILE("EBGraphImplem::nbytesFull");
    size_t retval = 0;
    //the tag
    retval += sizeof(TAG);

    //region, domain, graphnode box
    retval +=  3*Box::linearSize();
    if(m_tag == HasIrregular)
    {
      for (BoxIterator bit(m_graph.box()); bit.ok(); ++bit)
      {
        const GraphNode& node = m_graph(bit(), 0);
        int nodeSize = node.linearSize();
        retval += nodeSize;
      }
      retval += m_irregIVS.linearSize();
      retval += m_multiIVS.linearSize();
    }
    return retval;
  }
  /*******************************/
  std::size_t
  EBGraphImplem::
  copyToMemFull(void*      a_buf) const
  {
    //first the region and domain
    BL_PROFILE("EBGraphImplem::copyToMemFulll");
    size_t retval = 0;
    size_t incrval = 0;
    //the tag
    TAG* intbuf = (TAG *) a_buf;
    *intbuf = m_tag;
    retval += sizeof(TAG);
    intbuf++;

    //now into byte mode
    //region, domain, graphnode box
    unsigned char* buf = (unsigned char*) intbuf;
    m_region.linearOut(buf);
    incrval = m_region.linearSize();
    buf    += incrval;
    retval += incrval;

    m_domain.linearOut(buf);
    incrval = m_region.linearSize();
    buf    += incrval;
    retval += incrval;

    Box graphbox = m_graph.box();
    graphbox.linearOut(buf);
    incrval = m_graph.box().linearSize();
    buf    += incrval;
    retval += incrval;

    if(m_tag == HasIrregular)
    {
      for (BoxIterator bit(m_graph.box()); bit.ok(); ++bit)
      {
        //IntVect ivdeb(D_DECL(10, 7, 0));
        //int ideb = 0;
        //if(bit() == ivdeb)
        //{
        //  ideb = 1;
        //}
        const GraphNode& node = m_graph(bit(), 0);
        node.linearOut(buf);
        incrval = node.linearSize();
        buf    += incrval;
        retval += incrval;

      }
      m_irregIVS.linearOut(buf);
      incrval =  m_irregIVS.linearSize();
      buf    += incrval;
      retval += incrval;
      m_multiIVS.linearOut(buf);
      incrval =  m_multiIVS.linearSize();
      buf    += incrval;
      retval += incrval;
    }
    return retval;
  }
        
  /*******************************/
  std::size_t 
  EBGraphImplem::
  copyFromMemFull(const void* a_buf)
  {
    //first the region and domain
    BL_PROFILE("EBGraphImplem::copyFromMemFulll");
    size_t retval = 0;
    size_t incrval = 0;
    //the tag
    TAG* intbuf = (TAG*) a_buf;
    m_tag = *intbuf;
    retval += sizeof(TAG);
    intbuf++;

    //now into byte mode
    //region, domain, graphnode box
    unsigned char* buf = (unsigned char*) intbuf;
    m_region.linearIn(buf);
    incrval = m_region.linearSize();
    buf    += incrval;
    retval += incrval;

    m_domain.linearIn(buf);
    incrval = m_region.linearSize();
    buf    += incrval;
    retval += incrval;

    Box graphbox;
    graphbox.linearIn(buf);
    incrval = m_graph.box().linearSize();
    buf    += incrval;
    retval += incrval;

    if(m_tag == HasIrregular)
    {
      m_graph.resize(graphbox, 1);
      for (BoxIterator bit(m_graph.box()); bit.ok(); ++bit)
      {
        //IntVect ivdeb(D_DECL(10, 7, 0));
        //int ideb = 0;
        //if(bit() == ivdeb)
        //{
        //  ideb = 1;
        //}
        GraphNode& node = m_graph(bit(), 0);
        node.linearIn(buf);
        incrval = node.linearSize();
        buf    += incrval;
        retval += incrval;

      }
      m_irregIVS.linearIn(buf);
      incrval =  m_irregIVS.linearSize();
      buf    += incrval;
      retval += incrval;

      m_multiIVS.linearIn(buf);
      incrval =  m_multiIVS.linearSize();
      buf    += incrval;
      retval += incrval;
    }

    m_isDefined   = true;
    m_isDomainSet = true;
    m_isMaskBuilt = false;
    return retval;
  }
  /*******************************/

  GraphNode
  EBGraphImplem::
  getGraphNode(const IntVect& a_iv) const
  {
    GraphNode retval;
    if(isAllRegular())
    {
      retval.defineAsRegular();
    }
    else if(isAllCovered())
    {
      retval.defineAsCovered();
    }
    else
    {
      retval = m_graph(a_iv, 0);
    }
    return retval;
  }
  /*******************************/
  std::size_t 
  EBGraphImplem::
  nBytes(const Box& a_region, int start_comp, int ncomps) const
  {
    std::size_t linearSize = 0;
    //domain
    linearSize +=  Box::linearSize();
    for (BoxIterator bit(a_region); bit.ok(); ++bit)
    {
      GraphNode node = getGraphNode(bit());
      size_t nodeSize = node.linearSize();
      linearSize += nodeSize;
    }

    return linearSize;
  }
        
  /*******************************/
  std::size_t
  EBGraphImplem::
  copyToMem (const Box& a_region,
             int        srccomp,
             int        numcomp,
             void*      a_buf) const
  {
    BL_ASSERT(isDefined());

    std::size_t retval = 0;
    unsigned char* buffer = (unsigned char*) a_buf;

    m_domain.linearOut(buffer);
    size_t incrval = Box::linearSize();
    buffer    += incrval;
    retval += incrval;

    for (BoxIterator bit(a_region); bit.ok(); ++bit)
    {
      GraphNode node = getGraphNode(bit());
      size_t nodeSize = node.linearSize();
      node.linearOut(buffer);
      buffer += nodeSize;
      retval += nodeSize;
    }
    return retval;
  }
        
  /*******************************/
  std::size_t 
  EBGraphImplem::
  copyFromMem (const Box&  a_region,
               int         dstcomp,
               int         numcomp,
               const void* a_buf)
  {
    BL_ASSERT(isDefined());

    std::size_t retval = 0;
    unsigned char* buffer = (unsigned char*) a_buf;

    m_domain.linearIn(buffer);
    m_isDomainSet = true;

    size_t incrval = Box::linearSize();
    buffer    += incrval;
    retval += incrval;

    m_region &= m_domain;


    if (isAllRegular() || isAllCovered())
    {
      m_tag = HasIrregular;
      m_multiIVS = IntVectSet();
      m_irregIVS = IntVectSet();
      m_graph.resize(m_region, 1);
    }
    for (BoxIterator bit(a_region); bit.ok(); ++bit)
    {

      GraphNode& node = m_graph(bit(), 0);

      node.linearIn(buffer);
      if (node.isIrregular())
      {
        (m_irregIVS)|=bit();
      }
      if (node.size()>1)
      {
        (m_multiIVS)|=bit();
      }
      size_t nodeSize = node.linearSize();
      buffer += nodeSize;
      retval += nodeSize;

    }

    return retval;
  }
}
