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
#include "AMReX_Print.H"


namespace amrex
{
  static const IntVect ebg_debiv(D_DECL(0, 0, 0));
  bool EBGraphImplem::s_verbose = false;
  /*******************************/
  Vector<FaceIndex> EBGraph::getMultiValuedFaces(const int&  a_idir,
                                                      const Box&  a_box) const
  {
    return
      m_implem->getMultiValuedFaces(a_idir, a_box, *this);
  }
        
        
  /*******************************/
  Vector<FaceIndex> 
  EBGraphImplem::
  getMultiValuedFaces(const int&     a_idir,
                      const Box&     a_box,
                      const EBGraph& a_ebgraph) const
  {
    Vector<FaceIndex> multiValuedFaces;
    Box ghostRegion = a_box;
    ghostRegion.grow(a_idir, 1);
    const IntVectSet& ivsMulti = a_ebgraph.getMultiCells(ghostRegion);
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
          const Vector<FaceIndex> faces = getAllFaces(ivHi, a_idir, Side::Lo);
          isMulti = (faces.size() > 1);
        }
        else if (a_ebgraph.getDomain().contains(ivLo))
        {
          const Vector<FaceIndex> faces = getAllFaces(ivLo, a_idir, Side::Hi);
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
            Vector<int>& nodeArcsDir = node.m_arc[gNodeIndex];
            nodeArcsDir.resize(0);
            //find which vof node is connected to in each direction.
            //cannot use isConnected here because it will always return
            //true since one vof is still regular
            IntVect  otherIV = iv + sign(sit())*BASISV(idir);
            if (m_domain.contains(otherIV))
            {
              Vector<VolIndex> otherVoFs = getVoFs(otherIV);
              bool found = false;
              for (int iother = 0; iother < otherVoFs.size(); iother++)
              {
                const VolIndex& otherVoF = otherVoFs[iother];
                Vector<FaceIndex> otherFaces = getFaces(otherVoF, idir, flip(sit()));
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
    IntVectSet retval = m_irregIVS;
    retval &= a_subbox;
    return retval;
  }
        
  /*******************************/
  IntVectSet EBGraphImplem::getMultiCells(const Box& a_subbox) const
  {
    IntVectSet retval = m_multiIVS;
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
    EBCellFlag flag;
    flag.setRegular();
    m_cellFlags.setVal(flag);
    m_cellFlags.setType(FabType::regular);
  }
        
  /*******************************/
  void EBGraphImplem::setToAllCovered()
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    m_tag = AllCovered;
    m_irregIVS = IntVectSet();
    m_multiIVS = IntVectSet();

    EBCellFlag flag;
    flag.setCovered();
    m_cellFlags.setVal(flag);
    m_cellFlags.setType(FabType::covered);
  }
/*******************************/
  Vector<VolIndex>
  EBGraphImplem::
  getVoFs (const VolIndex& a_vof,
           const int& a_dir,
           const Side::LoHiSide& a_sd,
           const int& a_steps) const
  {
    assert((a_dir >= 0) && (a_dir < SpaceDim));
    assert(a_steps >= 0);

    Vector<VolIndex> retVoFs(1, a_vof);
    for (int irad = 1; irad <= a_steps; irad++)
    {
      Vector<VolIndex> tempVoFs(0);
      for (int ivof = 0; ivof < retVoFs.size(); ivof++)
      {
        const VolIndex& stepVoF = retVoFs[ivof];
        Vector<FaceIndex> faces = getFaces(stepVoF, a_dir, a_sd);
        for (int iface = 0; iface < faces.size(); iface++)
        {
          const FaceIndex& face = faces[iface];
          if (!face.isBoundary())
          {
            const VolIndex& flipVoF = face.getVoF(a_sd);
            tempVoFs.push_back(flipVoF);
          }
        }
      }
      retVoFs = tempVoFs;
    }
    return retVoFs;
  }
        
  /*******************************/
  void 
  EBGraphImplem::
  setCellFlags()
  {
    BL_PROFILE("ebgraph::setCellFlags()");
    if(isAllRegular())
    {
      EBCellFlag flag;
      flag.setRegular();
      m_cellFlags.setVal(flag);
      m_cellFlags.setType(FabType::regular);
    }
    else if(isAllCovered())
    {
      EBCellFlag flag;
      flag.setCovered();
      m_cellFlags.setVal(flag);
      m_cellFlags.setType(FabType::covered);
    }
    else
    {
      Box bx = m_region & m_domain;
      bx &= m_cellFlags.box();
      int ipt = 0;
      for (BoxIterator bi(bx); bi.ok(); ++bi)
      {
        const IntVect& iv = bi();
        auto& cellflag = m_cellFlags(iv);
        if (isRegular(iv)) 
        {
          cellflag.setRegular();
        } 
        else if (isCovered(iv)) 
        {
          cellflag.setCovered();
        } 
        else if (isMultiValued(iv)) 
        {
          cellflag.setMultiValued(numVoFs(iv));
        } 
        else 
        {
          cellflag.setSingleValued();
        }
        ipt++;
      }
      fixCellFlagType();
      Box domain = m_domain;

      const Box& ibx = amrex::grow(m_cellFlags.box(),-1) & domain;
    //begin debug
    //if(m_region.contains(ebg_debiv))
    //{
    //  amrex::Print() << "EBGraph::setCellFlags initial phase flag at " << ebg_debiv << " = " << m_cellFlags(ebg_debiv, 0).getValue();
    //  amrex::Print() << ", fullregion = " << m_fullRegion << ", ibx = " << ibx << endl;
    //}
    // end debug
      for (BoxIterator bi(ibx); bi.ok(); ++bi)
      {
        const IntVect& iv = bi();
        EBCellFlag& cellflag = m_cellFlags(iv, 0);
        cellflag.setDisconnected();
        cellflag.setConnected(IntVect::TheZeroVector());

        const auto& vofs = getVoFs(iv);
        for (const auto& vi : vofs)
        {
          std::array<int,AMREX_SPACEDIM> dirs = {AMREX_D_DECL(0,1,2)};
          do
          {
            IntVect offset_0 = IntVect::TheZeroVector();
            for (SideIterator sit_0; sit_0.ok(); ++sit_0)
            {
              offset_0[dirs[0]] = amrex::sign(sit_0());

              const Vector<VolIndex>& vofs_0 = getVoFs(vi, dirs[0], sit_0(), 1);
              for (const auto& vi_0 : vofs_0)
              {
                cellflag.setConnected(offset_0);

#if (AMREX_SPACEDIM >= 2)
                IntVect offset_1 = offset_0;
                for (SideIterator sit_1; sit_1.ok(); ++sit_1)
                {
                  offset_1[dirs[1]] = amrex::sign(sit_1());

                  const auto& vofs_1 = getVoFs(vi_0, dirs[1], sit_1(), 1);
                  for (const auto& vi_1 : vofs_1)
                  {
                    cellflag.setConnected(offset_1);

#if (AMREX_SPACEDIM == 3)
                    IntVect offset_2 = offset_1;
                    for (SideIterator sit_2; sit_2.ok(); ++sit_2)
                    {
                      offset_2[dirs[2]] = amrex::sign(sit_2());

                      const auto& vofs_2 = getVoFs(vi_1, dirs[2], sit_2(), 1);
                      for (const auto& vi_2 : vofs_2)
                      {
                        cellflag.setConnected(offset_2);
                      }
                    }
#endif
                  }
                }
#endif
              }
            }
          } while (std::next_permutation(dirs.begin(), dirs.end()));
        }
      }        
    }
    //begin debug
    //if(m_region.contains(ebg_debiv))
    //{
    //  amrex::Print() << "EBGraph::setCellFlags flag at " << ebg_debiv << " = " << m_cellFlags(ebg_debiv, 0).getValue() << endl;
    //}
    //end debug
  }        
  /*******************************/
  void 
  EBGraphImplem::
  buildGraph(const BaseFab<int>          & a_regIrregCovered,
             const Vector<IrregNode>& a_irregGraph,
             const Box                   & a_validRegion,
             const Box                   & a_domain)
  {
    define(a_validRegion);
    setDomain(a_domain);
    m_region  &= m_domain;

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
      Vector<GraphNodeImplem>& vecNodes = *(m_graph(iv, 0).m_cellList);
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
          const Vector<int>& irregArcs = inputNode.m_arc[irregIndex];
          Vector<int>& nodeArcs = node.m_arc[gNodeIndex];
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
    setCellFlags();
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
    m_fullRegion = a_region;
    m_irregIVS = IntVectSet();
    m_multiIVS = IntVectSet();
    m_mask.clear();
    m_isMaskBuilt = false;
    m_isDefined= true;
 
    m_cellFlags.resize(m_fullRegion, 1);
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
  Vector<VolIndex> 
  EBGraphImplem::
  getVoFs(const IntVect& a_iv) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
    Vector<VolIndex> retvec;
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
  Vector<FaceIndex> 
  EBGraphImplem::
  getAllFaces(const IntVect&        a_iv,
              const int&            a_idir,
              const Side::LoHiSide& a_sd) const
  {
    Vector<FaceIndex> retval(0);
    const Vector<VolIndex>& vofs = getVoFs(a_iv);
    for (int ivof= 0; ivof < vofs.size(); ivof++)
    {
      const Vector<FaceIndex>& faces = getFaces(vofs[ivof], a_idir, a_sd);
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
  Vector<FaceIndex> 
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
        
    Vector<FaceIndex> retvec;
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
  void EBGraph::fillIntMask(BaseFab<int>& a_mask,
                            const Box   & a_subbox) const
  {
    Box fullRegion = getFullRegion();

    fullRegion &= a_subbox;
    Box interiorRegion = fullRegion;
    interiorRegion &= a_mask.box();
    const Box& domain = getDomain();

    for(BoxIterator boxit(interiorRegion); boxit.ok(); ++boxit)
    {
      const IntVect& iv = boxit();
      if(!domain.contains(iv))
      {
        //set values outside domain to 2
        a_mask(iv, 0) = 2;
      }
      else if(isRegular(iv))
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
      const Box& b = a_mask.box();
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

  /*******************************/
  FaceIndex 
  EBGraphImplem::
  coarsen(const FaceIndex& a_fineFace) const
  {
    const VolIndex& loVoFCoar = coarsen(a_fineFace.getVoF(Side::Lo));
    const VolIndex& hiVoFCoar = coarsen(a_fineFace.getVoF(Side::Hi));
    return FaceIndex(loVoFCoar, hiVoFCoar, a_fineFace.direction());
  }
        
  /*******************************/
  Vector<FaceIndex> 
  EBGraphImplem::
  refine(const FaceIndex&     a_coarFace,
         const EBGraphImplem& a_fineGraph) const
  {
    Vector<FaceIndex> retval;
    const IntVect& ivLoCoar = a_coarFace.gridIndex(Side::Lo);
    const IntVect& ivHiCoar = a_coarFace.gridIndex(Side::Hi);
    int direction = a_coarFace.direction();
    if (m_region.contains(ivLoCoar) && m_region.contains(ivHiCoar))
    {
      //interior face
      const Vector<VolIndex>& loVoFFine = refine(a_coarFace.getVoF(Side::Lo));
      const Vector<VolIndex>& hiVoFFine = refine(a_coarFace.getVoF(Side::Hi));
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
      const Vector<VolIndex>& loVoFsFine = refine(a_coarFace.getVoF(Side::Lo));
      Box fineRegion = m_region;
      fineRegion.refine(2);
      for (int ivof = 0; ivof < loVoFsFine.size(); ivof++)
      {
        const VolIndex& loVoFFine = loVoFsFine[ivof];
        IntVect ivHiFine = loVoFFine.gridIndex() + BASISV(direction);
        if (!fineRegion.contains(ivHiFine))
        {
          const Vector<FaceIndex>& fineFaces = a_fineGraph.getFaces(loVoFFine, direction, Side::Hi);
          retval.insert(retval.end(), fineFaces.begin(), fineFaces.end());
        }
      }
    }
    else if (m_region.contains(ivHiCoar))
    {
      //boundary face on the low side of the domain
      const Vector<VolIndex> hiVoFsFine = refine(a_coarFace.getVoF(Side::Hi));
      Box fineRegion = m_region;
      fineRegion.refine(2);
        
      for (int ivof = 0; ivof < hiVoFsFine.size(); ivof++)
      {
        const VolIndex& hiVoFFine = hiVoFsFine[ivof];
        IntVect ivLoFine = hiVoFFine.gridIndex() - BASISV(direction);
        if (!fineRegion.contains(ivLoFine) )
        {
          const Vector<FaceIndex>& fineFaces = a_fineGraph.getFaces(hiVoFFine, direction, Side::Lo);
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
        
    const Vector<FaceIndex>& faces = getFaces(vofLo, direction, Side::Hi);
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
  Vector<VolIndex> 
  EBGraphImplem::
  refine(const VolIndex& a_coarVoF) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(isDomainSet());
        
    Vector<VolIndex> retval(0);
    IntVect ivCoar = a_coarVoF.gridIndex();
    BL_ASSERT(m_domain.contains(ivCoar));
        
    if (m_tag == AllRegular)
    {
      const IntVect& iv = a_coarVoF.gridIndex();
      Box refbox(iv,iv);
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
    BL_PROFILE("EBGraphImplem::copy");
    BL_ASSERT(isDefined());
    const Box& testbox = m_region & a_srcbox;

    
    if(!testbox.isEmpty())
    {
      setDomain(a_source.m_domain);
      m_region &= m_domain;

      Box regionTo  = testbox;
      Box regionFrom= testbox;
      if (isRegular(regionTo) && a_source.isRegular(regionFrom))
      {
      }
      else if (isCovered(regionTo) && a_source.isCovered(regionFrom))
      {
      }
      else if (a_source.isCovered(regionFrom) && regionTo.contains(m_region))
      {
        setToAllCovered();
      }
      else if (a_source.isRegular(regionFrom) && regionTo.contains(m_region))
      {
        setToAllRegular();

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
        std::unique_ptr<BaseFab<GraphNode> > raii;
        if (a_source.hasIrregular())
        {
          srcFabPtr = (BaseFab<GraphNode>*)&a_source.m_graph;
        }
        else
        {
          raii.reset(new BaseFab<GraphNode>(regionFrom, 1));
          srcFabPtr = raii.get();
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
        //these conditionals only depend on m_tag
        if (isAllRegular())
        {
          GraphNode regularNode;
          regularNode.defineAsRegular();
          m_graph.setVal(regularNode);
        }
        else if(isAllCovered())
        {
          GraphNode  coveredNode;
          coveredNode.defineAsCovered();
          m_graph.setVal(coveredNode);
        }
        //copy the data and set the tag properly
        m_tag = HasIrregular;

        m_graph.copy(*srcFabPtr, regionFrom, 0,regionTo, 0, 1);
        
        //  now fix up the IntVectSets to match the information
        if (a_source.hasIrregular())
        {
          IntVectSet ivsInterIrreg = (a_source.m_irregIVS);
          ivsInterIrreg &= regionTo;
          ivsInterIrreg &= m_region;
          (m_irregIVS) |= ivsInterIrreg;
        
          IntVectSet ivsInterMulti = (a_source.m_multiIVS);
          ivsInterMulti &= regionTo;
          ivsInterMulti &= m_region;
          m_multiIVS |= ivsInterMulti;
        }
      }
    }

    m_cellFlags.copy(a_source.m_cellFlags, a_srcbox, 0, a_destbox, 0, 1);
    fixCellFlagType();
    return *this;
  }
  /*******************************/
  void 
  EBGraphImplem::
  fixCellFlagType ()
  {
    Box bx = m_region;
    int nregular=0, nsingle=0, nmulti=0, ncovered=0;
    int ncells = bx.numPts();
    for (BoxIterator bi(bx); bi.ok(); ++bi)
    {
      const IntVect& iv = bi();
      if (isRegular(iv)) 
      {
        ++nregular;
      } 
      else if (isCovered(iv)) 
      {
        ++ncovered;
      } 
      else if (isMultiValued(iv)) 
      {
        ++nmulti;
//          amrex::Abort("EBLevel: multi-value cell not supported yet");
      } 
      else 
      {
        ++nsingle;
      }
    }
        
    if (nregular == ncells) 
    {
      m_cellFlags.setType(FabType::regular);
    } 
    else if (ncovered == ncells) 
    {
      m_cellFlags.setType(FabType::covered);
    } 
    else if (nmulti > 0) 
    {
      m_cellFlags.setType(FabType::multivalued);
    } 
    else 
    {
      m_cellFlags.setType(FabType::singlevalued);
    }
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
        
          const Vector<Vector<VolIndex> >& fineVoFSets
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
  Vector<Vector<VolIndex> >  
  EBGraphImplem::
  getVoFSets(const Box& a_region) const
  {
    Vector<Vector<VolIndex> > retval;
    //gather all  vofs
    Vector<VolIndex> allVoFs;
    for (BoxIterator bit(a_region); bit.ok(); ++bit)
    {
      const Vector<VolIndex>& newVoFs = getVoFs(bit());
      allVoFs.insert(allVoFs.end(), newVoFs.begin(), newVoFs.end());
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
        Vector<VolIndex> thisVoFSet(1, thisVoF);
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
  Vector<int> 
  EBGraphImplem::
  coarsenFaces(const VolIndex&       a_coarVoF,
               const EBGraphImplem&  a_fineGraph,
               const int&            a_idir,
               const Side::LoHiSide& a_sd)
  {
    BL_ASSERT(m_isDomainSet);
    BL_ASSERT(a_fineGraph.isDomainSet());
    Vector<int> retval;
        
    IntVect coarIV = a_coarVoF.gridIndex();
    IntVect otherIV= coarIV + sign(a_sd)*BASISV(a_idir);
    const Vector<VolIndex>& theseFineVoFs = refine(a_coarVoF);
    EBGraphImplem&  coarGhostGraph  = *this;  //trying to reduce communcation in this
        
    if (m_domain.contains(otherIV))
    {
        
      //interior faces.
      //just get all possible vofs to connect to and
      //check connectivity
      const Vector<VolIndex>& otherCoarVoFs = coarGhostGraph.getVoFs(otherIV);
      for (int iotherCoar = 0; iotherCoar < otherCoarVoFs.size(); iotherCoar++)
      {
        const VolIndex& otherCoarVoF = otherCoarVoFs[iotherCoar];
        const Vector<VolIndex>& otherFineVoFs = coarGhostGraph.refine(otherCoarVoF);
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
        const Vector<FaceIndex>& fineFaces =
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
          const Vector<VolIndex>& vofsCoar = getVoFs(bit());
          Vector<GraphNodeImplem>& nodes =
            *(m_graph(bit(), 0).m_cellList);
          for (int ivof = 0; ivof < vofsCoar.size(); ivof++)
          {
            const VolIndex& vofCoar= vofsCoar[ivof];
            GraphNodeImplem& node = nodes[vofCoar.cellIndex()];
            for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
              {
                const Vector<int>&& coarArcs =
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
    setCellFlags();
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
        
          const Vector<GraphNodeImplem>&
            vofsCoar = *(m_graph(ivCoar, 0).m_cellList);
        
          int numVofsCoar = vofsCoar.size();
        
          for (int icoar = 0; icoar < numVofsCoar; icoar++)
          {
            VolIndex vofCoar(ivCoar, icoar);
            const Vector<VolIndex>& vofsFine = refine(vofCoar);
        
            for (int ifine = 0; ifine < vofsFine.size(); ifine++)
            {
              const IntVect& ivFine = vofsFine[ifine].gridIndex();
        
              if (a_fineGraph.isIrregular(ivFine))
              {
                int cellIndexFine = vofsFine[ifine].cellIndex();
        
                Vector<GraphNodeImplem>&
                  nodesFine = *(a_fineGraph.m_graph(ivFine, 0).m_cellList);
        
                nodesFine[cellIndexFine].m_coarserNode = icoar;
              }
              else if ((numVofsCoar > 1) && (a_fineGraph.isRegular(ivFine)))
              {
                GraphNode& nodeFine = a_fineGraph.m_graph(ivFine, 0);
        
                nodeFine.m_cellList = new Vector<GraphNodeImplem>(1);
        
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

    //region, domain, fullregion, graphnode box
    retval +=  4*Box::linearSize();
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
    incrval = m_domain.linearSize();
    buf    += incrval;
    retval += incrval;

    m_fullRegion.linearOut(buf);
    incrval = m_fullRegion.linearSize();
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
    incrval = m_domain.linearSize();
    buf    += incrval;
    retval += incrval;

    m_fullRegion.linearIn(buf);
    incrval = m_fullRegion.linearSize();
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

    m_cellFlags.resize(m_fullRegion,1);
    setCellFlags();

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
    BL_PROFILE("EBGraphImplem::nBytes");
    std::size_t linearSize = 0;
    //domain
    linearSize +=  Box::linearSize();
    //for the cell flags
    linearSize += m_cellFlags.nBytes(a_region, 0, 1);
    //flag for reg or covered
    linearSize += sizeof(int);
    bool isRegBox = isRegular(a_region);
    bool isCovBox = isCovered(a_region);
    if((!isRegBox) && (!isCovBox))
    {
      for (BoxIterator bit(a_region); bit.ok(); ++bit)
      {
        const GraphNode& node = getGraphNode(bit());
        size_t nodeSize = node.linearSize();
        linearSize += nodeSize;
      }
    }

    //amrex::AllPrint() << "linearsize: a_region = " << a_region << ", retval = " << linearSize << endl;;
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
    BL_PROFILE("EBGraphImplem::copyToMem");
    BL_ASSERT(isDefined());

    std::size_t retval = 0;
    unsigned char* buffer = (unsigned char*) a_buf;

    m_domain.linearOut(buffer);
    size_t incrval = Box::linearSize();
    buffer += incrval;
    retval += incrval;

    incrval = m_cellFlags.copyToMem(a_region, 0, 1, buffer);
    buffer += incrval;
    retval += incrval;

    int regIrregCovCode = 0;
    if(isRegular(a_region))
    {
      regIrregCovCode = -1;
    }
    if(isCovered(a_region))
    {
      regIrregCovCode =  1;
    }
    //amrex::AllPrint() << "copytomem   region = " << a_region << ", regirregcov = " << regIrregCovCode << endl;

    int* intbuf = (int*) buffer;
    *intbuf = regIrregCovCode;
    incrval = sizeof(int);
    buffer += incrval;
    retval += incrval;

    if(regIrregCovCode == 0)
    {
      for (BoxIterator bit(a_region); bit.ok(); ++bit)
      {
        const GraphNode& node = getGraphNode(bit());
        size_t nodeSize = node.linearSize();
        node.linearOut(buffer);
        buffer += nodeSize;
        retval += nodeSize;
      }

    }

    //amrex::AllPrint() << "copytomem:   a_region = " << a_region << ", retval = " << retval << endl;;
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
    BL_PROFILE("EBGraphImplem::copyFromMem");
    BL_ASSERT(isDefined());

    std::size_t retval = 0;
    unsigned char* buffer = (unsigned char*) a_buf;

    m_domain.linearIn(buffer);
    m_isDomainSet = true;

    size_t incrval = Box::linearSize();
    buffer    += incrval;
    retval += incrval;

    incrval = m_cellFlags.copyFromMem(a_region, 0, 1, buffer);
    buffer += incrval;
    retval += incrval;

    m_region &= m_domain;

    int regIrregCovCode = 0;
    int* intbuf = (int*) buffer;
    regIrregCovCode = *intbuf;

    // amrex::AllPrint() << "copyfrommem region = " << a_region << ", regirregcov = " << regIrregCovCode << endl;

    incrval = sizeof(int);
    buffer += incrval;
    retval += incrval;

    bool allRegInput = (regIrregCovCode == -1);
    bool allCovInput = (regIrregCovCode ==  1);

    if(allRegInput && isRegular(a_region))
    {
      //amrex::AllPrint() << "copyfrommem: a_region = " << a_region << ", retval = " << retval << endl;;
      return retval;
    }
    else if(allCovInput && isCovered(a_region))
    {
      //amrex::AllPrint() << "copyfrommem: a_region = " << a_region << ", retval = " << retval << endl;;
      return retval;
    }

    //to get to this point, something is going to have to change in this object
    if (isAllRegular() || isAllCovered())
    {
      m_multiIVS = IntVectSet();
      m_irregIVS = IntVectSet();
      m_graph.resize(m_region, 1);
      if(isAllRegular())
      {
        GraphNode regularNode;
        regularNode.defineAsRegular();
        m_graph.setVal(regularNode);
      }
      if(isAllCovered())
      {
        GraphNode coveredNode;
        coveredNode.defineAsCovered();
        m_graph.setVal(coveredNode);
      }
      m_tag = HasIrregular;
    }
    for (BoxIterator bit(a_region); bit.ok(); ++bit)
    {

      GraphNode& node = m_graph(bit(), 0);
      if(allRegInput)
      {
        node.defineAsRegular();
      }
      else if(allCovInput)
      {
        node.defineAsCovered();
      }
      else
      {
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
    }


    fixCellFlagType();

    //amrex::AllPrint() << "copyfrommem: a_region = " << a_region << ", retval = " << retval << endl;;
    return retval;
  }
}
