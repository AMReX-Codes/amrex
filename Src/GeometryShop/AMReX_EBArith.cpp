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

#include "AMReX_EBArith.H"


namespace amrex
{
  //-----
  ///returns true if coarsenable (either by agglomeration or otherwise)
   bool
   EBArith::
   createCoarserEBLG(EBLevelGrid       &  a_eblgCoar,
                     const EBLevelGrid &  a_eblgFine,
                     const int         &  a_refRat,
                     const int         &  a_minBoxSize,
                     const int         &  a_maxBoxSize)
   {
     bool coarsenable;
     int testRef = a_minBoxSize*a_refRat;
     if(a_eblgFine.coarsenable(testRef))
     {
       coarsen(a_eblgCoar, a_eblgFine, a_refRat);
       coarsenable = true;
     }
     else
     {
       coarsenable = false;

       //turning this stuff off so I can make some progress
//       //see if we can agglomerate.
//       long numPtsRemaining = a_eblgFine.getDomain().numPts();
//       for(int ibox = 0; ibox < a_eblgFine.getDBL().size(); ibox++)
//       {
//         numPtsRemaining -= a_eblgFine.getDBL()[ibox].numPts();
//       }
//       coarsenable = (numPtsRemaining == 0);
//       if(coarsenable)
//       {
//         Box coarDom = coarsen(a_eblgFine.getDomain(), a_refRat);
//         BoxArray ba(coarDom);
//         ba.maxSize(a_maxBoxSize);
//         DistributionMapping dm(ba);

//         a_eblgCoar = EBLevelGrid(ba, dm, coarDom, a_eblgFine.getGhost());
//       }
     }
     return coarsenable;
   }
  //-----
  
  RealVect
  EBArith::
  getDomainNormal(int a_idir, Side::LoHiSide a_side)
  {
    RealVect normal = BASISREALV(a_idir);
    if (a_side == Side::Hi)
    {
      normal *= -1.0;
    }
    return normal;
  }
  //-----
  Box EBArith::
  adjCellBox(const Box            & a_valid, 
             const int            & a_idir, 
             const Side::LoHiSide & a_side, 
             const int            & a_len)
  {
    Box retval;
    if(a_side == Side::Lo)
    {
      retval = adjCellLo(a_valid, a_idir, a_len);
    }
    else
    {
      retval = adjCellHi(a_valid, a_idir, a_len);
    }
    return retval;
  }
  //-----
  void EBArith::
  ExtrapolateBC(BaseFab<Real>&    a_state,
                const Box&      a_valid,
                Real            a_dx,
                int             a_dir,
                Side::LoHiSide  a_side)
  {
    BL_PROFILE("EBArith::ExtrapolateBC");

    int isign = sign(a_side);

    Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
    toRegion &= a_state.box();

    for (BoxIterator bit(toRegion); bit.ok(); ++bit)
    {
      const IntVect& ivTo = bit();

      IntVect ivClose = ivTo -   isign*BASISV(a_dir);
      IntVect ivFar   = ivTo - 2*isign*BASISV(a_dir);

      for (int icomp = 0; icomp < a_state.nComp(); icomp++)
      {
        Real nearVal = a_state(ivClose, icomp);
        Real farVal  = a_state(ivFar,   icomp);

        Real ghostVal = 2.*nearVal - farVal;
        a_state(ivTo, icomp) = ghostVal;
      }
    }
  }

  //---------
  Real
  EBArith::
  getDiagWeight(  VoFStencil&     a_vofStencil,
                  const VolIndex& a_vof,
                  int             a_ivar)
  {
    //has to be zero because adding to this later
    Real retval = 0;
    bool found = false;
    for (int ivof = 0; ivof  < a_vofStencil.size(); ivof++)
    {
      if ((a_vofStencil.vof(ivof) == a_vof) && (a_vofStencil.variable(ivof) == a_ivar))
      {
        found = true;
        //additive in case there are more than one entry with vof == a_vof
        retval += a_vofStencil.weight(ivof);
      }
    }
    if (!found)
    {
      //      MayDay::Warning("no diagonal weight, probably an empty cell");
      retval = 1;
    }
    return retval;
  }
  //---------
  void
  EBArith::
  getMultiColors(Vector<IntVect>& a_colors)
  {

#if AMREX_SPACEDIM==2
    a_colors.resize(4);
    a_colors[0] = IntVect::Zero;//(0,0)
    a_colors[1] = IntVect::Unit;//(1,1)
    a_colors[2] = IntVect::Zero + BASISV(1);//(0,1)
    a_colors[3] = IntVect::Zero + BASISV(0);//(1,0)
#elif AMREX_SPACEDIM==3
    a_colors.resize(8);
    a_colors[0] = IntVect::Zero;//(0,0,0)
    a_colors[1] = IntVect::Zero + BASISV(0) + BASISV(1);//(1,1,0)
    a_colors[2] = IntVect::Zero + BASISV(1) + BASISV(2);//(0,1,1)
    a_colors[3] = IntVect::Zero + BASISV(0) + BASISV(2);//(1,0,1)
    a_colors[4] = IntVect::Zero + BASISV(1);//(0,1,0)
    a_colors[5] = IntVect::Zero + BASISV(0);//(1,0,0)
    a_colors[6] = IntVect::Zero + BASISV(2);//(0,0,1)
    a_colors[7] = IntVect::Unit;//(1,1,1)
#endif
  }
  void
  EBArith::
  computeCoveredFaces(Vector<VolIndex>&     a_coveredFace,
                      IntVectSet&           a_coveredSets,
                      IntVectSet&           a_irregIVS,
                      const int&            a_idir,
                      const Side::LoHiSide& a_sd,
                      const EBISBox&        a_ebisBox,
                      const Box&            a_region)
  {
    BL_PROFILE("EBArith::computeCoveredFaces");
    //first compute the sets where where the covered faces exist
    //start with all irregular cells and subtract off  cells
    //whose vofs all have faces in the given direction
    a_irregIVS =  a_ebisBox.getIrregIVS(a_region);
    a_coveredSets = a_irregIVS;
    a_coveredFace.resize(0);
    for (IVSIterator ivsit(a_irregIVS); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      Vector<VolIndex> vofs = a_ebisBox.getVoFs(iv);
      bool allVoFsHaveFaces = true;
      for (int ivof = 0; ivof < vofs.size(); ivof++)
      {
        const VolIndex& vof = vofs[ivof];
        Vector<FaceIndex> faces = a_ebisBox.getFaces(vof, a_idir, a_sd);
        if (faces.size() == 0)
        {
          allVoFsHaveFaces = false;
          a_coveredFace.push_back(vof);
        }
      }
      if (allVoFsHaveFaces)
        a_coveredSets -= iv;
    }
  }

  /*****************************/
  FaceStencil
  EBArith::
  getInterpStencil(const FaceIndex&     a_face,
                   const IntVectSet&    a_cfivs,
                   const EBISBox&       a_ebisBox,
                   const Box           & a_domainBox)
  {
    BL_PROFILE("EBArith::getInterpStencil");

    FaceStencil sten;
#if AMREX_SPACEDIM==2
    getInterpStencil2D(sten, a_face, a_cfivs, a_ebisBox, a_domainBox);

#elif AMREX_SPACEDIM==3
    getInterpStencil3D(sten, a_face, a_cfivs, a_ebisBox, a_domainBox);
#else
    bogus_bl_spacedim_macro();
#endif
    //EBArith::computeInterpStencil(sten, a_face, a_ebisBox, a_domainBox, a_face.direction());
    return sten;
  }
  void
  EBArith::
  getInterpStencil2D(FaceStencil&         a_sten,
                     const FaceIndex&     a_face,
                     const IntVectSet&    a_coarseFineIVS,
                     const EBISBox&       a_ebisBox,
                     const Box           & a_domainBox)
  {
    BL_PROFILE("EBArith::getInterpStencil2D");
    BL_ASSERT(SpaceDim == 2);
    a_sten.clear();
    RealVect faceCentroid = a_ebisBox.centroid(a_face);
    int faceDir = a_face.direction();
    int tanDir = 1 - faceDir;
    Side::LoHiSide tanFaceSide;
    if (faceCentroid[tanDir] > 0.)
    {
      tanFaceSide = Side::Hi;
    }
    else if (faceCentroid[tanDir] < 0.)
    {
      tanFaceSide = Side::Lo;
    }
    else
    {
      tanFaceSide = Side::Lo;
      // MayDay::Error("EBArith: getInterpStencil2D faceCentroid[tanDir] = 0");
    }

    FaceIndex otherFace;
    bool uniqueFace = EBArith::getAdjacentFace(otherFace, a_face,
                                               a_ebisBox, a_domainBox,
                                               tanDir, tanFaceSide);

    //drop order of interpolation to zero if the other face is not unique
    //or if this face or the otherface is on the coarse fine interface
    bool dropOrder = false;
    dropOrder = !uniqueFace;
    if (!dropOrder)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        if (a_coarseFineIVS.contains(   a_face.gridIndex(sit())) ||
            a_coarseFineIVS.contains(otherFace.gridIndex(sit())))
        {
          dropOrder = true;
        }
      }
    }

    if (dropOrder)
    {
      a_sten.add(a_face, 1.0);
    }
    else
    {
      BL_ASSERT(a_face.isDefined());
      BL_ASSERT(otherFace.isDefined());
      Real dist = std::abs(faceCentroid[tanDir]);
      Real thisWeight = 1.0 - dist;
      Real otherWeight =  dist;
      a_sten.add(a_face, thisWeight);
      a_sten.add(otherFace, otherWeight);
    }
  }
  /*****************************/
  bool
  EBArith::
  getAdjacentFace(FaceIndex           & a_adjacentFace,
                  const FaceIndex     & a_face,
                  const EBISBox       & a_ebisBox,
                  const Box           & a_domain,
                  const int           & a_idir,
                  const Side::LoHiSide& a_side)
  {
    BL_PROFILE("EBArith::getAdjacentFace");
    bool uniqueFace = true;
    int faceDir = a_face.direction();
    VolIndex loFaceVoF = a_face.getVoF(Side::Lo);
    VolIndex hiFaceVoF = a_face.getVoF(Side::Hi);
    //figure out if we have the connectivity to find a face
    //with which to interpolate.  If so,
    //find the other face
    if (!a_face.isBoundary())
    {
      Vector<FaceIndex> loTanFaces = a_ebisBox.getFaces(loFaceVoF, a_idir, a_side);
      Vector<FaceIndex> hiTanFaces = a_ebisBox.getFaces(hiFaceVoF, a_idir, a_side);
      if ((loTanFaces.size() != 1) || (hiTanFaces.size() != 1))
      {
        uniqueFace = false;
      }
      else if ((loTanFaces[0].isBoundary()) && (hiTanFaces[0].isBoundary()))
      {
        uniqueFace = false;
      }
      else
      {
        const FaceIndex& loTanFace = loTanFaces[0];
        const FaceIndex& hiTanFace = hiTanFaces[0];
        VolIndex loOtherVoF = loTanFace.getVoF(a_side);
        VolIndex hiOtherVoF = hiTanFace.getVoF(a_side);
        if (!(a_ebisBox.isConnected(loOtherVoF, hiOtherVoF)))
        {
          uniqueFace = false;
        }
        if (uniqueFace)
        {
          Vector<FaceIndex> otherFaces = a_ebisBox.getFaces(loOtherVoF, faceDir, Side::Hi);
          if (otherFaces.size() != 1)
          {
            uniqueFace = false;
          }
          else
          {
            a_adjacentFace = otherFaces[0];
          }
        }
      }
    } //end if !face.isBoundary
    else
    {
      //boundary face.
      IntVect loVoFiv = loFaceVoF.gridIndex();
      IntVect hiVoFiv = hiFaceVoF.gridIndex();
      IntVect loDomiv = a_domain.smallEnd();
      IntVect hiDomiv = a_domain.bigEnd();

      if (hiVoFiv[faceDir] == loDomiv[faceDir])
      {
        Vector<FaceIndex> hiTanFaces = a_ebisBox.getFaces(hiFaceVoF, a_idir, a_side);
        if (hiTanFaces.size() != 1)
        {
          uniqueFace = false;
        }
        else if (hiTanFaces[0].isBoundary())
        {
          //in 3d,  can be a boundary in two directions
          uniqueFace = false;
        }
        else
        {
          VolIndex hiOtherVoF = hiTanFaces[0].getVoF(a_side);
          const Vector<FaceIndex>& otherFaces = a_ebisBox.getFaces(hiOtherVoF, faceDir, Side::Lo);
          if (otherFaces.size() != 1)
          {
            uniqueFace = false;
          }
          else
          {
            a_adjacentFace = otherFaces[0];
          }
        }

      }
      else if (loVoFiv[faceDir] == hiDomiv[faceDir])
      {
        Vector<FaceIndex> loTanFaces = a_ebisBox.getFaces(loFaceVoF, a_idir, a_side);
        if (loTanFaces.size() != 1)
        {
          uniqueFace = false;
        }
        else if (loTanFaces[0].isBoundary())
        {
          //in 3d,  can be a boundary in two directions
          uniqueFace = false;
        }
        else
        {
          VolIndex loOtherVoF = loTanFaces[0].getVoF(a_side);
          const Vector<FaceIndex>& otherFaces = a_ebisBox.getFaces(loOtherVoF, faceDir, Side::Hi);
          if (otherFaces.size() != 1)
          {
            uniqueFace = false;
          }
          else
          {
            a_adjacentFace = otherFaces[0];
          }
        }

      }
      else
      {
        amrex::Error("boundary face confusion");
      }

    }
    return uniqueFace;
  }

  /*****************************/
  void
  EBArith::
  getInterpStencil3D(FaceStencil&         a_sten,
                     const FaceIndex&     a_face,
                     const IntVectSet&    a_coarseFineIVS,
                     const EBISBox&       a_ebisBox,
                     const Box           & a_domainBox)
  {
    BL_PROFILE("EBArith::getInterpStencil3D");
    BL_ASSERT(SpaceDim == 3);
    a_sten.clear();

    int faceDir = a_face.direction();
    RealVect faceCentroid = a_ebisBox.centroid(a_face);
    //first check to see if this is a genuinely 3d face
    //centroid > 0 in all directions
    bool really3D = true;
    int izerod = 7;
    Real tol = 1.0e-8;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != faceDir)
      {
        if (std::abs(faceCentroid[idir]) < tol)
        {
          really3D= false;
          izerod = idir;
        }
      }
    }
    if (!really3D)
    {
      //use 2D stencil because we have a centroid that is zero
      //in the given direction
      int tanDir = 3 - izerod - faceDir;
      Side::LoHiSide tanFaceSide;
      if (faceCentroid[tanDir] > 0)
        tanFaceSide = Side::Hi;
      else
        tanFaceSide = Side::Lo;

      bool dropOrder = false;
      FaceIndex otherFace;
      bool uniqueFace = EBArith::getAdjacentFace(otherFace, a_face,
                                                 a_ebisBox, a_domainBox,
                                                 tanDir, tanFaceSide);
      dropOrder = !uniqueFace;

      if ( !dropOrder )
      {
        for (SideIterator sit; sit.ok(); ++sit)
        {
          if (a_coarseFineIVS.contains(   a_face.gridIndex(sit())) ||
              a_coarseFineIVS.contains(otherFace.gridIndex(sit())))
          {
            dropOrder = true;
          }
        }
      }

      if (dropOrder)
      {
        a_sten.add(a_face, 1.0);
      }
      else
      {
        BL_ASSERT(a_face.isDefined());
        BL_ASSERT(otherFace.isDefined());
        Real dist = std::abs(faceCentroid[tanDir]);
        Real thisWeight = 1 - dist;
        Real otherWeight =  dist;
        a_sten.add(a_face, thisWeight);
        a_sten.add(otherFace, otherWeight);
      }
    }
    else
    {
      int tanDirs[2];
      Side::LoHiSide tanFaceSides[2];
      {
        int itan = 0;
        for (int idir = 0; idir < 3; idir++)
        {
          if (idir != faceDir)
          {
            tanDirs[itan] = idir;
            if (faceCentroid[tanDirs[itan]] > 0)
              tanFaceSides[itan] = Side::Hi;
            else
              tanFaceSides[itan] = Side::Lo;
            itan++;
          }
        }
      }

      bool dropOrder = false;
      FaceIndex face01, face10, face11;
      int ixDir = tanDirs[0];
      int iyDir = tanDirs[1];
      Side::LoHiSide xSide = tanFaceSides[0];
      Side::LoHiSide ySide = tanFaceSides[1];
      //figure out whether to drop order.
      //if not,  get the faces involved
      VolIndex vofLo = a_face.getVoF(Side::Lo);
      VolIndex vofHi = a_face.getVoF(Side::Hi);

      //first get face10 which is the in
      //the  ixdir direction from the input face
      if (!dropOrder)
      {
        bool uniqueFace10
          = EBArith::getAdjacentFace(face10, a_face, a_ebisBox, a_domainBox, ixDir, xSide);
        dropOrder = !uniqueFace10;
        //now get face11 which is the in
        //the  ixdir,iydir direction from the input face (diagonal)
        if (uniqueFace10)
        {
          bool uniqueFace11
            = EBArith::getAdjacentFace(face11, face10, a_ebisBox, a_domainBox, iyDir, ySide);
          dropOrder = !uniqueFace11;
        }
      }
      //the  ixdir,iydir direction from the input face (diagonal)
      //now get face01 which is the in
      //the  iydir direction from the input face
      if (!dropOrder)
      {
        //first get face01 which is the in
        //the  iydir direction from the input face
        bool uniqueFace01
          = EBArith::getAdjacentFace(face01, a_face, a_ebisBox, a_domainBox, iyDir, ySide);
        dropOrder = !uniqueFace01;
        //now get face11 which is the in
        //the  ixdir,iydir direction from the input face (diagonal)
        //compute temp face and see if it is the same as
        //the one computed from the other direction.
        //if not , drop order
        FaceIndex face11temp;
        if (uniqueFace01)
        {
          bool uniqueFace11
            = EBArith::getAdjacentFace(face11temp, face01, a_ebisBox, a_domainBox, ixDir, xSide);
          dropOrder = !uniqueFace11;
          if ((!dropOrder) && !(face11temp == face11) )
          {
            dropOrder = true;
          }
        }

        //finally if any of the stencil faces are in the coarse-fine interface, drop order
        if (!dropOrder)
        {
          for (SideIterator sit; sit.ok(); ++sit)
          {
            if ((a_coarseFineIVS.contains(a_face.gridIndex(sit()))) ||
                (a_coarseFineIVS.contains(face01.gridIndex(sit()))) ||
                (a_coarseFineIVS.contains(face10.gridIndex(sit()))) ||
                (a_coarseFineIVS.contains(face11.gridIndex(sit()))))
            {
              dropOrder = true;
            }
          }
        }
      }
      ///////////////////////
      //construct the stencils
      ///////////////////////
      if (dropOrder)
      {
        a_sten.add(a_face, 1.0);
      }
      else
      {
        FaceIndex face00 = a_face;
        Real xbar = std::abs(faceCentroid[ixDir]);
        Real ybar = std::abs(faceCentroid[iyDir]);
        Real f00coef = 1.0 - xbar - ybar + xbar*ybar;
        Real f10coef = xbar - xbar*ybar;
        Real f01coef = ybar - xbar*ybar;
        Real f11coef = xbar*ybar;

        for (SideIterator sit; sit.ok(); ++sit)
        {
          if (!face00.isBoundary())
            BL_ASSERT(face00.cellIndex(sit()) >= 0);
          if (!face01.isBoundary())
            BL_ASSERT(face01.cellIndex(sit()) >= 0);
          if (!face10.isBoundary())
            BL_ASSERT(face10.cellIndex(sit()) >= 0);
          if (!face11.isBoundary())
            BL_ASSERT(face11.cellIndex(sit()) >= 0);
        }
        BL_ASSERT(face00.isDefined());
        BL_ASSERT(face01.isDefined());
        BL_ASSERT(face10.isDefined());
        BL_ASSERT(face11.isDefined());

        a_sten.add(face00, f00coef);
        a_sten.add(face01, f01coef);
        a_sten.add(face10, f10coef);
        a_sten.add(face11, f11coef);
      }
    }
  }
  /*******/
  void
  EBArith::
  getAllVoFsInMonotonePath(Vector<VolIndex>& a_vofList,
                           const VolIndex&   a_vof,
                           const EBISBox&    a_ebisBox,
                           const int&        a_redistRad)
  {
    Vector<VolIndex> vofsStencil;
    IntVect timesMoved = IntVect::TheZeroVector();
    IntVect pathSign   = IntVect::TheZeroVector();
    EBArith::getAllVoFsInMonotonePath(vofsStencil, timesMoved,
                                      pathSign, a_vof, a_ebisBox,
                                      a_redistRad);
    a_vofList = vofsStencil;
  }
  /*******/
  void
  EBArith::
  getAllVoFsInMonotonePath(Vector<VolIndex>& a_vofList,
                           const IntVect&    a_timesMoved,
                           const IntVect&    a_pathSign,
                           const VolIndex&   a_vof,
                           const EBISBox&    a_ebisBox,
                           const int&        a_redistRad)
  {
    BL_PROFILE("EBArith::getAllVoFsInMonotonePath");
    const Box& domain = a_ebisBox.getDomain();
    if (domain.contains(a_vof.gridIndex()))
    {
      //check to see if we have already added it.
      //if not, add it
      bool found = false;
      for (int ivof = 0; ivof < a_vofList.size(); ivof++)
      {
        if (a_vofList[ivof] == a_vof)
        {
          found = true;
        }
      }
      if (!found)
      {
        a_vofList.push_back(a_vof);
      }
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        //only move redist radius times in a direction
        if (a_timesMoved[idir] < a_redistRad)
        {
          IntVect newTimesMoved = a_timesMoved;
          newTimesMoved[idir]++;
          //pathSign preserves monotonicity
          //we set pathsign to -1 once you start in the low
          //direction, 1 in the high direction.  It starts at zero

          //if we started in the low direction or have not started
          //add vofs on low side
          if ((a_pathSign[idir] == -1) || (a_pathSign[idir] == 0))
          {
            IntVect newSign = a_pathSign;
            newSign[idir] = -1;
            Vector<FaceIndex> facesLo =
              a_ebisBox.getFaces(a_vof, idir, Side::Lo);
            for (int iface = 0; iface < facesLo.size(); iface++)
            {
              VolIndex newVoF = facesLo[iface].getVoF(Side::Lo);
              getAllVoFsInMonotonePath(a_vofList, newTimesMoved, newSign,
                                       newVoF, a_ebisBox, a_redistRad);
            }
          }
          //if we started in the high direction or have not started
          //add vofs on high side
          if ((a_pathSign[idir] == 1) || (a_pathSign[idir] == 0))
          {
            IntVect newSign = a_pathSign;
            newSign[idir] = 1;
            Vector<FaceIndex> facesHi =
              a_ebisBox.getFaces(a_vof, idir, Side::Hi);
            for (int iface = 0; iface < facesHi.size(); iface++)
            {
              VolIndex newVoF = facesHi[iface].getVoF(Side::Hi);
              getAllVoFsInMonotonePath(a_vofList, newTimesMoved, newSign,
                                       newVoF, a_ebisBox, a_redistRad);
            }
          }
        } //end if (we are less than redist radius away)
      }//end loop over directions
    }
  }
  /*******/
  void
  EBArith::
  getVoFsDir(bool& a_hasClose, VolIndex& a_closeVoF,
             bool& a_hasFar,   VolIndex& a_farVoF,
             const EBISBox& a_ebisBox,
             const VolIndex& a_vof,
             int a_idir, Side::LoHiSide a_sd,
             IntVectSet*    a_cfivsPtr)
  {
    a_hasClose = false;
    a_hasFar   = false;
    bool checkCFIVS = false;
    if (a_cfivsPtr != NULL)
    {
      checkCFIVS = true;
    }
    //get faces on both sides to see in which direction we can do diffs
    Vector<FaceIndex> closeFaces = a_ebisBox.getFaces(a_vof, a_idir, a_sd);

    //boundary faces and multi-valued faces are to be one-sided away from
    a_hasClose = ((closeFaces.size() == 1) && (!closeFaces[0].isBoundary()));
    if (a_hasClose)
    {
      a_closeVoF = closeFaces[0].getVoF(a_sd);
      if (checkCFIVS && (*a_cfivsPtr).contains(a_closeVoF.gridIndex()))
      {
        a_hasClose = false;
      }
      Vector<FaceIndex> farFaces = a_ebisBox.getFaces(a_closeVoF, a_idir, a_sd);
      a_hasFar = ((farFaces.size() == 1) && (!farFaces[0].isBoundary()));
      if (a_hasFar)
      {
        a_farVoF =  farFaces[0].getVoF(a_sd);
        if (checkCFIVS && (*a_cfivsPtr).contains(a_farVoF.gridIndex()))
        {
          a_hasFar = false;
        }
      }
    }
  }
  /*******/
  void
  EBArith::
  getAllVoFsWithinRadius(Vector<VolIndex>& a_vofList,
                         const VolIndex&   a_vof,
                         const EBISBox&    a_ebisBox,
                         const int&        a_redistRad)
  {
    BL_PROFILE("EBArith::getAllVoFsWithinRadius(inner)");
    a_vofList.clear();
    Box grownBox(a_vof.gridIndex(), a_vof.gridIndex());
    grownBox.grow(a_redistRad);
    grownBox &= a_ebisBox.getDomain();
    IntVectSet vis(grownBox);
    VoFIterator vofit(vis, a_ebisBox.getEBGraph());
    a_vofList = vofit.getVector();
  }
  int
  EBArith::
  getFirstDerivStencil(VoFStencil&      a_sten,
                       const VolIndex&  a_vof,
                       const EBISBox&   a_ebisBox,
                       const int&       a_idir,
                       const Real&      a_dx,
                       IntVectSet*    a_cfivsPtr,
                       int ivar)

  {
    BL_ASSERT(a_dx > 0.0);
    BL_ASSERT(a_idir >= 0);
    BL_ASSERT(a_idir <  SpaceDim);
    int order;
    bool hasLo, hasLower, hasHi, hasHigher;
    VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
    EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF,   a_ebisBox, a_vof, a_idir, Side::Lo, a_cfivsPtr);
    EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisBox, a_vof, a_idir, Side::Hi, a_cfivsPtr);

    //clear any residual stencil info
    a_sten.clear();

    if (hasHi && hasLo)
    {
      //if we have vofs on both sides, we do our friend the centered difference
      order = 2;
      a_sten.add(hiVoF,  1.0, ivar);
      a_sten.add(loVoF, -1.0, ivar);
      a_sten *= (1.0/(2.*a_dx));
    }
    else if (hasHi)
    {
      //if we only have vofs on the high side, we see if we can take a higher
      //order one sided diff.   If not, we take simple one-sided diff and drop order
      if (hasHigher)
      {
        a_sten.add(hiVoF,       4.0, ivar);
        a_sten.add(higherVoF,  -1.0, ivar);
        a_sten.add(a_vof,      -3.0, ivar);
        a_sten *= (1.0/(2.*a_dx));
        order = 2;
      }
      else
      {
        a_sten.add(hiVoF,       1.0, ivar);
        a_sten.add(a_vof,      -1.0, ivar);
        a_sten *= (1.0/(a_dx));
        order = 1;
      }
    }
    else if (hasLo)
    {
      //if we only have vofs on the low side, we see if we can take a higher
      //order one sided diff.   If not, we take simple one-sided diff and drop order
      if (hasLower)
      {
        a_sten.add(loVoF,      -4.0, ivar);
        a_sten.add(lowerVoF,    1.0, ivar);
        a_sten.add(a_vof,       3.0, ivar);
        a_sten *= (1.0/(2.*a_dx));
        order = 2;
      }
      else
      {
        a_sten.add(a_vof,       1.0, ivar);
        a_sten.add(loVoF,      -1.0, ivar);
        a_sten *= (1.0/(a_dx));
        order = 1;
      }
    }
    else
    {
      //no vofs on either side.   return cleared stencil and order=0
      order = 0;
    }

    return order;
  }
///
/**
   Gets the stencil to take the first derivative of  cell centered data.
*/
  int
  EBArith::
  getSecondDerivStencil(VoFStencil&      a_sten,
                        const VolIndex&  a_vof,
                        const EBISBox&   a_ebisBox,
                        const int&       a_idir,
                        const Real&      a_dx,
                        IntVectSet*    a_cfivsPtr,
                        int ivar)
  {
    BL_PROFILE("EBArith::getSecondDerivStencil");
    BL_ASSERT(a_dx > 0.0);
    BL_ASSERT(a_idir >= 0);
    BL_ASSERT(a_idir <  SpaceDim);
    bool hasLo, hasLower, hasHi, hasHigher;
    VolIndex loVoF, lowerVoF, hiVoF, higherVoF;
    EBArith::getVoFsDir(hasLo, loVoF, hasLower, lowerVoF,   a_ebisBox, a_vof, a_idir, Side::Lo, a_cfivsPtr);
    EBArith::getVoFsDir(hasHi, hiVoF, hasHigher, higherVoF, a_ebisBox, a_vof, a_idir, Side::Hi, a_cfivsPtr);
    int order;

    //clear any residual stencil info
    a_sten.clear();


    if (hasHi && hasLo)
    {
      //if we have vofs on both sides, we do our friend the centered difference
      order = 2;
      a_sten.add(hiVoF,  1.0, ivar);
      a_sten.add(loVoF,  1.0, ivar);
      a_sten.add(a_vof, -2.0, ivar);
      a_sten *= (1.0/(a_dx*a_dx));
    }
    else if (hasHi)
    {
      //no vof on the low side
      //if we can, shift stencil one to the high side
      if (hasHigher)
      {
        a_sten.add(higherVoF,   1.0, ivar);
        a_sten.add(a_vof,       1.0, ivar);
        a_sten.add(hiVoF,      -2.0, ivar);
        a_sten *= (1.0/(a_dx*a_dx));
        order = 1;
      }
      else
      {
        //need 3 points
        order = 0;
      }
    }
    else if (hasLo)
    {
      //no vof on the high side
      //if we can, shift stencil one to the low side
      if (hasLower)
      {
        a_sten.add(lowerVoF,    1.0, ivar);
        a_sten.add(a_vof,       1.0, ivar);
        a_sten.add(loVoF,      -2.0, ivar);
        a_sten *= (1.0/(a_dx*a_dx));
        order = 1;
      }
      else
      {
        //need 3 points
        order = 0;
      }
    }
    else
    {
      //no vofs on either side.   return cleared stencil and order=0
      order = 0;
    }

    return order;
  }
  /***************/
  bool
  EBArith::isVoFHere(VolIndex& a_vof2,
                     const Vector<VolIndex>& a_vofsStencil,
                     const IntVect& a_cell2)
  {
    int whichVoF;
    bool found = isVoFHere(a_vof2, whichVoF, a_vofsStencil, a_cell2);

    return found;
  }
  /***********/
  bool
  EBArith::isVoFHere(VolIndex& a_vof2, int& a_whichVoF,
                     const Vector<VolIndex>& a_vofsStencil,
                     const IntVect& a_cell2)
  {
    BL_PROFILE("EBArith::isVoFHere");
    bool found = false;
    for (int isten = 0; isten < a_vofsStencil.size(); isten++)
    {
      if (a_vofsStencil[isten].gridIndex() == a_cell2)
      {
        if (found == true)
        {
          //vof not unique.  return false;
          return false;
        }
        else
        {
          found = true;
          a_vof2 = a_vofsStencil[isten];
          a_whichVoF = isten;
        }
      }
    }

    return found;
  }
  /***********/
  int
  EBArith::
  getMixedDerivStencil(VoFStencil&      a_sten,
                       const VolIndex&  a_vof,
                       const EBISBox&   a_ebisBox,
                       const int&       a_dir1,
                       const int&       a_dir2,
                       const Real&      a_dx1,
                       const Real&      a_dx2,
                       IntVectSet*    a_cfivsPtr,
                       int ivar)
  {
    BL_PROFILE("EBArith::getMixedDerivStencil");
    BL_ASSERT(a_dx1 > 0.0);
    BL_ASSERT(a_dx2 > 0.0);
    BL_ASSERT(a_dir1 >= 0);
    BL_ASSERT(a_dir2 >= 0);
    BL_ASSERT(a_dir1 <  SpaceDim);
    BL_ASSERT(a_dir2 <  SpaceDim);
    bool checkCFIVS = false;
    if (a_cfivsPtr != NULL)
    {
      checkCFIVS = true;
    }
    int radius = 1;
    Vector<VolIndex> vofList;

    EBArith::getAllVoFsInMonotonePath(vofList, a_vof, a_ebisBox, radius);


    //clear any residual stencil info
    a_sten.clear();

    //see how many corners have all vofs available for a mixed stencil derivative
    //and average the available stencils
    int numSten = 0;
    IntVect iv = a_vof.gridIndex();
    IntVect ivOne, ivTwo, ivCor;
    VolIndex vofOne, vofTwo, vofCor;
    for (SideIterator sit1; sit1.ok(); ++sit1)
    {
      int sign1 = sign(sit1());
      ivOne = iv + sign1*BASISV(a_dir1);
      bool isOneHere = EBArith::isVoFHere(vofOne, vofList, ivOne);
      if (isOneHere)
      {
        for (SideIterator sit2; sit2.ok(); ++sit2)
        {
          int sign2 = sign(sit2());
          ivTwo = iv + sign2*BASISV(a_dir2);
          bool isTwoHere = EBArith::isVoFHere(vofTwo, vofList, ivTwo);
          if (isTwoHere)
          {
            ivCor = iv + sign2*BASISV(a_dir2) +  sign1*BASISV(a_dir1);
            bool isCorHere = EBArith::isVoFHere(vofCor, vofList, ivCor);

            if (isCorHere && checkCFIVS)
            {
              if ( (*a_cfivsPtr).contains(vofCor.gridIndex()) ||
                   (*a_cfivsPtr).contains(vofOne.gridIndex()) ||
                   (*a_cfivsPtr).contains(vofTwo.gridIndex()))
              {
                isCorHere = false;
              }
            }
            if (isCorHere)
            {
              //if we have all vofs, construct the stencil.
              Real rsign = Real(sign1*sign2);

              VoFStencil mixedStencil;
              mixedStencil.add(a_vof,   1.0, ivar);
              mixedStencil.add(vofCor,  1.0, ivar);
              mixedStencil.add(vofOne, -1.0, ivar);
              mixedStencil.add(vofTwo, -1.0, ivar);

              mixedStencil *= rsign/(a_dx1*a_dx2);

              a_sten += mixedStencil;
              numSten++;
            }
          }
        }
      }
    }
    int order = 0;
    if (numSten > 0)
    {
      a_sten *= 1.0/numSten;

      if (numSten == 4)
      {
        order = 2;
      }
      else
      {
        order = 1;
      }
    }

    return order;
  }
  /******/
  RealVect
  EBArith::
  getFaceLocation(const FaceIndex& a_face,
                  const RealVect&  a_dx,
                  const RealVect&  a_probLo)
  {
    BL_PROFILE("EBArith::getFaceLocation");
    const IntVect& iv = a_face.gridIndex(Side::Hi);
    RealVect loc = a_probLo;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir == a_face.direction())
      {
        loc[idir] += a_dx[idir]*Real(iv[idir]);
      }
      else
      {
        loc[idir] += a_dx[idir]*(Real(iv[idir]) + 0.5);
      }
    }
    return loc;
  }
  /***/
  RealVect
  EBArith::
  getVofLocation(const VolIndex& a_vof,
                 const RealVect& a_dx,
                 const RealVect& a_probLo)
  {
    BL_PROFILE("EBArith::getVofLocation");
    const IntVect& iv = a_vof.gridIndex();
    RealVect loc = getIVLocation(iv,a_dx,a_probLo);
    return loc;
  }
  RealVect
  EBArith::
  getIVLocation(const IntVect&  a_iv,
                const RealVect& a_dx,
                const RealVect& a_probLo)
  {
    BL_PROFILE("EBArith::getIVLocation");
    RealVect loc = a_probLo;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      loc[idir] += a_dx[idir]*(Real(a_iv[idir]) + 0.5);
    }
    return loc;
  }
  /****/
  // This is the original way of doing the least squares stencil when eb normal points out of domain
  // Only use stencil of all available vofs when symmetry is needed
  void
  EBArith::
  getLeastSquaresGradSten(VoFStencil&     a_stencil,
                          Real&           a_weight,
                          const RealVect& a_normal  ,
                          const RealVect& a_centroid,
                          const VolIndex& a_vof,
                          const EBISBox&  a_ebisBox,
                          const RealVect& a_dx,
                          const Box& a_domain,
                          int a_ivar)
  {
    BL_PROFILE("EBArith::getLeastSquaresGradSten");
    bool needSymStencil = true;
    for (int idir=0; idir<SpaceDim; idir++)
    {
      if (needSymStencil)
      {
        if (abs(a_normal[idir]) != 1. && a_normal[idir] != 0.)
        {
          needSymStencil = false;
        }
      }
    }

    if (needSymStencil)
    {
      getLeastSquaresGradStenAllQuad(a_stencil, a_weight, a_normal,
                                     a_centroid, a_vof, a_ebisBox,
                                     a_dx,a_domain, a_ivar, true);
      // getLeastSquaresGradStenAllVoFs(a_stencil, a_weight, a_normal,
      //                                a_centroid, a_vof, a_ebisBox,
      //                                a_dx,a_domain, a_ivar);
    }
    else
    {
      IntVect quadrant;
      const RealVect& normal = a_ebisBox.normal(a_vof);
      //IntVect iv0 = a_vof.gridIndex();

      //Box domain = a_domain;
      for (int idir=0; idir<SpaceDim; idir++)
      {
        // if (iv0[idir]==domain.smallEnd(idir))
        //   {
        //     quadrant[idir]=1;
        //   }
        // else if (iv0[idir]==domain.bigEnd(idir))
        //   {
        //     quadrant[idir]=-1;
        //   }
        // else
        //   {
        if (normal[idir] < 0)
        {
          quadrant[idir]=-1;
        }
        else
        {
          quadrant[idir]=1;
        }
        // }
      }

      getLeastSquaresGradSten(a_stencil, a_weight, a_normal, a_centroid,
                              quadrant, a_vof, a_ebisBox, a_dx, a_domain,
                              a_ivar);

      if (a_stencil.size()==0)
      {
        getLeastSquaresGradStenAllQuad(a_stencil, a_weight, a_normal,
                                       a_centroid, a_vof, a_ebisBox,
                                       a_dx,a_domain, a_ivar);
      }
    }
  }
  /****/
  void
  EBArith::
  getLeastSquaresGradStenAllQuad(VoFStencil&          a_stencil,
                                 Real&                a_weight,
                                 const RealVect&      a_normal,
                                 const RealVect&      a_centroid,
                                 const VolIndex&      a_vof,
                                 const EBISBox&       a_ebisBox,
                                 const RealVect&      a_dx,
                                 const Box& a_domain,
                                 int                  a_ivar,
                                 bool                 a_doSymmetric)
  {
    BL_PROFILE("EBArith::getLeastSquaresGradStenAllQuad");
    // if we've gotten here, then we shouldn't have a stencil
    a_stencil.clear();
    a_weight = 0.;
    // number of valid stencils so far
    int numValidStencils = 0;

    // calculate the quadrant
    IntVect quadrant;
    //IntVect iv0 = a_vof.gridIndex();

    for (int idir=0; idir<SpaceDim; idir++)
    {
      if (a_normal[idir] < 0)
      {
        quadrant[idir]=-1;
      }
      else
      {
        quadrant[idir]=1;
      }
    }

    // next stencil candidate
    VoFStencil curStencil;
    IntVect curQuadrant;
    Real curWeight = 0.;
    bool validCurQuadrant;

    // in symmetric case and near domain boundaries this quadrant will not
    //  have been tried, so we try it here
    getLeastSquaresGradSten(curStencil, curWeight, a_normal, a_centroid,
                            quadrant, a_vof, a_ebisBox, a_dx, a_domain,
                            a_ivar);

    if (curStencil.size() != 0)
    {
      numValidStencils += 1;
      a_stencil += curStencil;
      a_weight += curWeight;
    }

    // cycle through directions, switching the sign of quadrant's
    //  components looking for a better one
    for (int idir=0; idir<SpaceDim; idir++)
    {
      curQuadrant = IntVect::TheUnitVector();
      // first try flipping this direction
      curQuadrant[idir] *= -1;
      validCurQuadrant = curQuadrant != quadrant;
      if (!a_doSymmetric)
      {
        validCurQuadrant = validCurQuadrant && (curQuadrant != -quadrant);
      }
      if (validCurQuadrant)
      {
        curStencil.clear();
        curWeight = 0.;
        getLeastSquaresGradSten(curStencil, curWeight, a_normal, a_centroid,
                                curQuadrant, a_vof, a_ebisBox, a_dx, a_domain,
                                a_ivar);

        if (curStencil.size() != 0)
        {
          numValidStencils += 1;
          a_stencil += curStencil;
          a_weight += curWeight;
        }
      }
#if AMREX_SPACEDIM == 3
      // now try flipping a second direction
      int jdir = (idir+1)%SpaceDim;
      curQuadrant[jdir] *= -1;
      validCurQuadrant = curQuadrant != quadrant;
      if (!a_doSymmetric)
      {
        validCurQuadrant = validCurQuadrant && (curQuadrant != -quadrant);
      }
      if (validCurQuadrant)
      {
        curStencil.clear();
        curWeight = 0.;
        getLeastSquaresGradSten(curStencil, curWeight, a_normal, a_centroid,
                                curQuadrant, a_vof, a_ebisBox, a_dx, a_domain,
                                a_ivar);

        if (curStencil.size() != 0)
        {
          numValidStencils += 1;
          a_stencil += curStencil;
          a_weight += curWeight;
        }
      }
#endif
    }
    // lastly, try flipping all directions
    curQuadrant = IntVect::TheUnitVector();
    curQuadrant *= -1;
    validCurQuadrant = curQuadrant != quadrant;
    if (!a_doSymmetric)
    {
      validCurQuadrant = validCurQuadrant && (curQuadrant != -quadrant);
    }
    if (validCurQuadrant)
    {
      curStencil.clear();
      curWeight = 0.;
      getLeastSquaresGradSten(curStencil, curWeight, a_normal, a_centroid,
                              curQuadrant, a_vof, a_ebisBox, a_dx, a_domain,
                              a_ivar);

      if (curStencil.size() != 0)
      {
        numValidStencils += 1;
        a_stencil += curStencil;
        a_weight += curWeight;
      }
    }

    if (numValidStencils > 1)
    {
      a_stencil *= 1./numValidStencils;
      a_weight *= 1./numValidStencils;
    }
  }
  /***/
  void
  EBArith::calculateWeightingMatrix(RealVect           x0,
                                    Vector<RealVect>&  xp,
                                    Vector<RealVect>&  weightMatrix,
                                    bool&              detZero)
  {
    BL_PROFILE("EBArith::calculateWeightingMatrix");
    int stenSize = xp.size();

    Vector<RealVect> deltaX = xp;
    for (int isten = 0; isten < stenSize; isten++)
    {
      deltaX[isten] -= x0;
    }

    Vector<RealVect> aTransA(SpaceDim,RealVect::Zero), invATransA(SpaceDim,RealVect::Zero);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (int jdir = 0; jdir < SpaceDim; jdir++)
      {
        for (int isten = 0; isten < stenSize; isten++)
        {
          aTransA[idir][jdir] = aTransA[idir][jdir]
            + deltaX[isten][idir]*deltaX[isten][jdir];
        }
      }
    }

    Real det;
#if   AMREX_SPACEDIM == 2
    det = aTransA[0][0] * aTransA[1][1] - aTransA[0][1] * aTransA[1][0];
    if (det < 1.e-15 && det > -1.e-15)
    {
      detZero = true;
    }
    else
    {
      invATransA[0][0] =  aTransA[1][1] / det;
      invATransA[0][1] = -aTransA[0][1] / det;
      invATransA[1][0] = -aTransA[1][0] / det;
      invATransA[1][1] =  aTransA[0][0] / det;
    }
#elif AMREX_SPACEDIM == 3
    det = aTransA[0][0] * ( aTransA[1][1] * aTransA[2][2]
                            - aTransA[1][2] * aTransA[2][1])
      + aTransA[0][1] * ( aTransA[1][2] * aTransA[2][0]
                          - aTransA[1][0] * aTransA[2][2])
      + aTransA[0][2] * ( aTransA[1][0] * aTransA[2][1]
                          - aTransA[1][1] * aTransA[2][0]);

    if (det < 1.e-15 && det > -1.e-15)
    {
      detZero = true;
    }
    else
    {
      invATransA[0][0] = ( aTransA[1][1] * aTransA[2][2]
                           - aTransA[1][2] * aTransA[2][1]) / det;
      invATransA[0][1] = ( aTransA[1][2] * aTransA[2][0]
                           - aTransA[1][0] * aTransA[2][2]) / det;
      invATransA[0][2] = ( aTransA[1][0] * aTransA[2][1]
                           - aTransA[1][1] * aTransA[2][0]) / det;
      invATransA[1][0] = ( aTransA[2][1] * aTransA[0][2]
                           - aTransA[2][2] * aTransA[0][1]) / det;
      invATransA[1][1] = ( aTransA[2][2] * aTransA[0][0]
                           - aTransA[2][0] * aTransA[0][2]) / det;
      invATransA[1][2] = ( aTransA[2][0] * aTransA[0][1]
                           - aTransA[2][1] * aTransA[0][0]) / det;
      invATransA[2][0] = ( aTransA[0][1] * aTransA[1][2]
                           - aTransA[0][2] * aTransA[1][1]) / det;
      invATransA[2][1] = ( aTransA[0][2] * aTransA[1][0]
                           - aTransA[0][0] * aTransA[1][2]) / det;
      invATransA[2][2] = ( aTransA[0][0] * aTransA[1][1]
                           - aTransA[0][1] * aTransA[1][0]) / det;
    }
#else
    THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_CH_SPACEDIM_IS_2_OR_3;
#endif

    //if (!detZero)
    {
      weightMatrix.resize(stenSize,RealVect::Zero);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        for (int isten = 0; isten < stenSize; isten++)
        {
          for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            weightMatrix[isten][idir] += invATransA[idir][jdir] * deltaX[isten][jdir];
          }
        }
      }
    }
  }
  void
  EBArith::calculateWeightingMatrixRed(RealVect           x00,
                                       Vector<RealVect>&  xpp,
                                       IntVect            dimm,
                                       Vector<RealVect>&  weightMatrix,
                                       bool&              deadRed)
  //CP: do the same thing for a reduced system, where some neighbors in the normal leastSquare stencil are covered
  //some dimensions might also have vanished. these need to be recorded outside
  //dimm[idir]==1, idir is valid, dimm[idir]==0, idir is missing

  {
    BL_PROFILE("EBArith::calculateWeightingMatrixRed");

    int stenSize = xpp.size();

    //now cast the problem to reduced dimension: make x0, xp
    int nr = 0;
    RealVect x0=RealVect::Zero;
    Vector<RealVect> xp(stenSize,RealVect::Zero);
    IntVect dirN;
    for (int idir = 0; idir< SpaceDim; idir++)
    {
      if (dimm[idir]>0)
      {
        nr++;
        x0[nr-1] = x00[idir];
        dirN[nr-1] = idir;
        for (int isten = 0; isten < stenSize; isten++)
        {
          xp[isten][nr-1]=xpp[isten][idir];
        }
      }
    }

    Vector<RealVect> deltaX = xp;
    for (int isten = 0; isten < stenSize; isten++)
    {
      deltaX[isten] -= x0;
    }

    Vector<RealVect> aTransA(SpaceDim,RealVect::Zero), invATransA(SpaceDim,RealVect::Zero);
    //CP: using the fact that nr <= SpaceDim

    for (int idir = 0; idir < nr; idir++)
    {
      for (int jdir = 0; jdir < nr; jdir++)
      {
        for (int isten = 0; isten < stenSize; isten++)
        {
          aTransA[idir][jdir] = aTransA[idir][jdir]
            + deltaX[isten][idir]*deltaX[isten][jdir];
        }
      }
    }

    Real det;
    if (nr == 1)
    {
      for (int isten = 0; isten< stenSize; isten++)
      {
        // this is worth more consideration when there is only one dimension
        // should all cells be included?
        invATransA[0][0] =  invATransA[0][0] + deltaX[isten][0]*deltaX[isten][0];
      }
      invATransA[0][0]=1/ invATransA[0][0];
      if ((invATransA[0][0]) == 0.0) deadRed = true;
    }
    else if (nr == 2)
    {
      det = aTransA[0][0] * aTransA[1][1] - aTransA[0][1] * aTransA[1][0];

      invATransA[0][0] =  aTransA[1][1] / det;
      invATransA[0][1] = -aTransA[0][1] / det;
      invATransA[1][0] = -aTransA[1][0] / det;
      invATransA[1][1] =  aTransA[0][0] / det;
      if ((det) == 0.0) deadRed = true;
    }
    else if (nr == 3)
    {
      det = aTransA[0][0] * ( aTransA[1][1] * aTransA[2][2]
                              - aTransA[1][2] * aTransA[2][1])
        + aTransA[0][1] * ( aTransA[1][2] * aTransA[2][0]
                            - aTransA[1][0] * aTransA[2][2])
        + aTransA[0][2] * ( aTransA[1][0] * aTransA[2][1]
                            - aTransA[1][1] * aTransA[2][0]);
      if (det< 0)
      {

      }
      else if (det == 0)
      {
        deadRed = true;
      }
      else
      {
        invATransA[0][0] = ( aTransA[1][1] * aTransA[2][2]
                             - aTransA[1][2] * aTransA[2][1]) / det;
        invATransA[0][1] = ( aTransA[1][2] * aTransA[2][0]
                             - aTransA[1][0] * aTransA[2][2]) / det;
        invATransA[0][2] = ( aTransA[1][0] * aTransA[2][1]
                             - aTransA[1][1] * aTransA[2][0]) / det;
        invATransA[1][0] = ( aTransA[2][1] * aTransA[0][2]
                             - aTransA[2][2] * aTransA[0][1]) / det;
        invATransA[1][1] = ( aTransA[2][2] * aTransA[0][0]
                             - aTransA[2][0] * aTransA[0][2]) / det;
        invATransA[1][2] = ( aTransA[2][0] * aTransA[0][1]
                             - aTransA[2][1] * aTransA[0][0]) / det;
        invATransA[2][0] = ( aTransA[0][1] * aTransA[1][2]
                             - aTransA[0][2] * aTransA[1][1]) / det;
        invATransA[2][1] = ( aTransA[0][2] * aTransA[1][0]
                             - aTransA[0][0] * aTransA[1][2]) / det;
        invATransA[2][2] = ( aTransA[0][0] * aTransA[1][1]
                             - aTransA[0][1] * aTransA[1][0]) / det;
      }
    }
    else
    {
      BL_ASSERT(nr<=3 && nr>0);
    }

    weightMatrix.resize(stenSize,RealVect::Zero);
    for (int idir = 0; idir < nr; idir++)
    {
      for (int isten = 0; isten < stenSize; isten++)
      {
        for (int jdir = 0; jdir < nr; jdir++)
        {
          weightMatrix[isten][dirN[idir]] += invATransA[idir][jdir] * deltaX[isten][jdir];
          // direction offset: set to original dimension
        }
      }
    }
  }


  /***************/
  bool
  EBArith::monotonePathVoFToCellVoF(VolIndex& a_vof2,
                                    const VolIndex& a_vof1,
                                    const IntVect& a_cell2,
                                    const EBISBox& a_ebisBox)
  {
    BL_PROFILE("EBArith::monotonePathVoFToCellVoF");
    IntVect diffVect = a_cell2 - a_vof1.gridIndex();
    int imaxdiff = 0;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (std::abs(diffVect[idir]) > imaxdiff)
        imaxdiff = std::abs(diffVect[idir]);
    }

    IntVect timesMoved = IntVect::TheZeroVector();
    IntVect pathSign   = IntVect::TheZeroVector();
    Vector<VolIndex> vofsStencil;
    getAllVoFsInMonotonePath(vofsStencil, timesMoved,
                             pathSign, a_vof1, a_ebisBox,
                             imaxdiff);

    bool found = isVoFHere(a_vof2,  vofsStencil, a_cell2);

    return found;
  }
  // This is the actual solution to the least squares problem with the inversion of matrices and solution in flux form
  void
  EBArith::
  getLeastSquaresGradSten(VoFStencil&     a_stencil,
                          Real&           a_weight,
                          const RealVect& a_normal,
                          const RealVect& a_centroid,
                          const IntVect&  a_quadrant,
                          const VolIndex& a_vof,
                          const EBISBox&  a_ebisBox,
                          const RealVect& a_dx,
                          const Box& a_domain,
                          int a_ivar)
  {
    BL_PROFILE("EBArith::getLeastSquaresGradSten");
    IntVect iv0 = a_vof.gridIndex();

#if   AMREX_SPACEDIM == 2
    int stenSize = 3;
#elif AMREX_SPACEDIM == 3
    int stenSize = 7;
#else
    THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_AMREX_SPACEDIM_IS_2_OR_3;
#endif
    //Box domainBox = a_domain;
    Vector<IntVect> ivSten(stenSize);

    ivSten[0] = iv0 + a_quadrant[0]*BASISV(0)                                                    ;
    ivSten[1] = iv0                           + a_quadrant[1]*BASISV(1)                          ;
    ivSten[2] = iv0 + a_quadrant[0]*BASISV(0) + a_quadrant[1]*BASISV(1)                          ;
#if AMREX_SPACEDIM == 3
    ivSten[3] = iv0                                                     + a_quadrant[2]*BASISV(2);
    ivSten[4] = iv0 + a_quadrant[0]*BASISV(0)                           + a_quadrant[2]*BASISV(2);
    ivSten[5] = iv0                           + a_quadrant[1]*BASISV(1) + a_quadrant[2]*BASISV(2);
    ivSten[6] = iv0 + a_quadrant[0]*BASISV(0) + a_quadrant[1]*BASISV(1) + a_quadrant[2]*BASISV(2);
#endif

    bool dropOrder = false;

    Vector<VolIndex> volSten(stenSize);
    for (int isten = 0; isten < stenSize; isten++)
    {
      //cp: it needs to be populated anyways
      if (a_ebisBox.getDomain().contains(ivSten[isten]))
      {
        if (a_ebisBox.getVoFs(ivSten[isten]).size() > 0)
        {
          volSten[isten] = a_ebisBox.getVoFs(ivSten[isten])[0];
        }
        else
        {
          dropOrder = true;
          volSten[isten] = VolIndex(IntVect(AMREX_D_DECL(0,0,0)),0);
          // break;
        }
      }
      else
      {
        volSten[isten] = VolIndex(IntVect(AMREX_D_DECL(0,0,0)),0);
        dropOrder = true;
        // break;
      }
    }

    std::vector<bool> monotonePath(stenSize);
    int nConn = stenSize;

    //restrictive because if one of the cells in the quadrant is not there then drop order
    for (int isten = 0; isten < stenSize; isten++)
    {
      monotonePath[isten] = EBArith::monotonePathVoFToCellVoF(volSten[isten],
                                                              a_vof,
                                                              ivSten[isten],
                                                              a_ebisBox);
      if (!monotonePath[isten])
      {
        dropOrder = true;
        nConn--;
      }
    }

    if (!dropOrder)
    {
      RealVect x0;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        x0[idir] = a_dx[idir] * (0.5 + a_centroid[idir] + iv0[idir]);
      }

      Vector<RealVect> xp(stenSize);
      for (int isten = 0; isten < stenSize; isten++)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          xp[isten][idir] = a_dx[idir] * (0.5 + ivSten[isten][idir]);
        }
      }

      Vector<RealVect> invATransAdeltaX(stenSize,RealVect::Zero);
      bool detZero = false;
      EBArith::calculateWeightingMatrix(x0, xp, invATransAdeltaX, detZero);

      a_stencil.clear();
      a_weight = 0.0;

      //if (!detZero)
      {
        for (int isten = 0; isten < stenSize; isten++)
        {
          Real dphidnWeight = 0.0;
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            dphidnWeight -= invATransAdeltaX[isten][idir] * a_normal[idir];
          }
              
          a_stencil.add(volSten[isten],dphidnWeight, a_ivar);
          a_weight -= dphidnWeight;
        }
      }
    }
    else
    {
      bool deadCell = true;
      if (nConn < 1)
      {
        deadCell = true;
      }
      else
      {
        //CP changed
        //collect all potentially usable points for the stencil
        RealVect x0;
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          x0[idir] = a_dx[idir] * (0.5 + a_centroid[idir] + iv0[idir]);
        }

        Vector<RealVect> xp;
        Vector<int> volStenIdx;
        int ns = 0;
        for (int isten = 0; isten < stenSize; isten++)
        {
          if (monotonePath[isten])
          {
            xp.push_back(RealVect::Zero);
            ns++;
            volStenIdx.push_back(isten);
            for (int idir = 0; idir < SpaceDim; idir++)
            {
              xp[ns-1][idir] = a_dx[idir] * (0.5 + ivSten[isten][idir]);
            }
          }
        }

        //determine which dimension to choose
        IntVect dimm = IntVect::TheZeroVector();
        IntVect mainDir = IntVect::TheZeroVector();
        //int diagDir[SpaceDim][2*SpaceDim - 1];
        int minDiag;
        int nDiagN = 2*(SpaceDim-1)-1; // the number of diagonal neighbors in one direction. Is this calculation correct? It is for 2D and 3D
        // mainDir and diagDir both contain neighbors' indices into monotonePath
        // main direction is its side neighbor, diagDir is the corner neighbors
        // how do I get diagDir? by looking at my 3D visualization box
        // by setting minDiag, we control, in the case of mainDir missing, how many corner neighbors must exist to make this dimension "valid". To always require the side neighbor, set this > 3 in 3D and > 1 in 2D

#if   AMREX_SPACEDIM == 2
        int diagDir[2][1] =
          {
            {
              2
            },
            {
              2
            }
          };
        mainDir   = IntVect(0,1);
        minDiag= 2; //minDiag = 2 this will always require the face neighbor
#elif AMREX_SPACEDIM == 3
        mainDir   = IntVect(0,1,3);
        int diagDir[3][3] =
          {
            {
              2,4,6
            },
            {
              2,5,6
            },
            {4,5,6
            }
          };
        minDiag= 2;
#else
        THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_AMREX_SPACEDIM_IS_2_OR_3;
#endif
        for (int iDir = 0; iDir < SpaceDim; iDir++)
        {
          if (monotonePath[mainDir[iDir]])
          {
            int diagCount = 0;
            for (int jDiag = 0; jDiag< nDiagN; jDiag++)
            {
              if (monotonePath[diagDir[iDir][jDiag] ]) diagCount++;
            }
            if (diagCount >= 1) dimm[iDir]=1;
          }
          else
            // face neighbor covered, counting number of diagonal neighbors
          {
            int diagCount = 0;
            for (int jDiag = 0; jDiag< nDiagN; jDiag++)
            {
              if (monotonePath[diagDir[iDir][jDiag] ]) diagCount++;
            }
            if (diagCount >= minDiag) dimm[iDir]=1;
          }
        }
        int dimsum=0;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          dimsum+= dimm[idir];
        }
        if (dimsum<1)
        {
          deadCell = true;
        }
        else
        {
          //calculate weights, corresponding to all stencil points
          Vector<RealVect> invATransAdeltaX(ns,RealVect::Zero);
          EBArith::calculateWeightingMatrixRed(x0,xp,dimm,invATransAdeltaX, deadCell);

          a_stencil.clear();
          a_weight = 0.0;

          for (int isten = 0; isten < xp.size(); isten++)
          {
            Real dphidnWeight = 0.0;
            for (int idir = 0; idir < SpaceDim; idir++)
            {
              dphidnWeight -= invATransAdeltaX[isten][idir] * a_normal[idir];
            }

            a_stencil.add(volSten[volStenIdx[isten]],dphidnWeight, a_ivar);
            a_weight -= dphidnWeight;
          }
        }
      }
      if (deadCell)
      {
        a_stencil.clear();
        a_weight = 0.0;
      }
    }
  }
  void
  EBArith::dataRayCast(bool&               a_dropOrder,
                       Vector<VoFStencil>& a_pointStencils,
                       Vector<Real>&       a_distanceAlongLine,
                       const RealVect&     a_normal,
                       const RealVect&     a_bndryCentroid,
                       const VolIndex&     a_vof,
                       const EBISBox&      a_ebisBox,
                       const RealVect&     a_dx,
                       const IntVectSet&   a_cfivs,
                       int a_ivar,
                       int a_numPoints)
  {
    a_dropOrder = false;

    BL_PROFILE("EBArith::johanStencil");
    IntVect iv = a_vof.gridIndex();

    //line starts at xb = location in physical space
    //of boundary centroid
    RealVect xb;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      xb[idir] = a_dx[idir]*(Real(iv[idir]) + 0.5 + a_bndryCentroid[idir]);
    }

    // Find InterpDirection and hiLo
    Real nMax = 0.0;
    int nMaxDir = 0;

    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      if (std::abs(a_normal[idir]) > nMax)
      {
        nMax = std::abs(a_normal[idir]);
        nMaxDir = idir;
      }
    }
    //sometimes normals can be zero
    if (std::abs(nMax) < 1.0e-15)
    {
      a_dropOrder = true;
      return;
    }
    int hiLo = 1;
    if (a_normal[nMaxDir] < 0.0)
    {
      hiLo = -1;
    }

    // Find Intersection of planes and ray
    //and the cells in which they live
    Vector<RealVect> intersectLoc(a_numPoints);
    Vector<IntVect>   intersectIV(a_numPoints);
    //equation of line
    //y = y_b + (ny/nx)*(x-xb)
    //we know x because that is the location of the
    //next cell center over in the nMaxDir direction
    //hiLo in the stencil depends on where the intersection point is in
    //relation to the cell center of the intersecting cell
    a_distanceAlongLine.resize(a_numPoints);
    a_pointStencils.resize(a_numPoints);
    const Box region = a_ebisBox.getRegion();

    for (int iinter = 0; iinter < a_numPoints; iinter++)
    {
      intersectIV[iinter] = iv + (iinter+1)*hiLo*BASISV(nMaxDir);
      // check whether intersectIV occurs outside of iv[idir!=nMaxDir]
      for (int idir=0; idir<SpaceDim; idir++)
      {
        if (idir != nMaxDir)
        {
          // what direction are we looking in?
          int isign = (a_normal[idir]<0.0) ? -1 : 1;
          // how far do we go in idir as a result of moving
          //  (iinter+1) cells in nMaxDir?
          Real xDist = std::abs((Real(iinter+1)-hiLo*a_bndryCentroid[nMaxDir])
                           *(a_normal[idir]/a_normal[nMaxDir]));
          // how far from the boundary centroid to the idir cell-edge
          Real xEdgeDist = std::abs(Real(isign)*0.5 - a_bndryCentroid[idir]);
          if (xDist > xEdgeDist)
          { // we're outside of iv[idir]. calculate by how many cells
            //  and adjust intersectIV accordingly
            intersectIV[iinter][idir] += isign*int(1+floor(xDist-xEdgeDist));
          }
        }
      }

      if (!region.contains(intersectIV[iinter]))
      {
        a_dropOrder = true;
        return;
      }
      if (a_ebisBox.numVoFs(intersectIV[iinter]) != 1)
      {
        a_dropOrder = true;
        return;
      }
      VolIndex centerVoF(intersectIV[iinter], 0);

      //small cells as centers can be bad
      if (a_ebisBox.volFrac(centerVoF) < 1.0e-7)
      {
        a_dropOrder = true;
        return;
      }

      Real xMaxDir  = a_dx[nMaxDir]*(Real(intersectIV[iinter][nMaxDir]) + 0.5);
      intersectLoc[iinter][nMaxDir] = xMaxDir;
      a_distanceAlongLine[iinter] = (xMaxDir-xb[nMaxDir])*(xMaxDir-xb[nMaxDir]);
      RealVect extrapDist = RealVect::Zero; //the initialization is important
      //get location of intersection and the distance along the line
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (idir != nMaxDir)
        {
          Real normalRat = a_normal[idir]/a_normal[nMaxDir];
          Real distDir =  normalRat*(xMaxDir-xb[nMaxDir]);
          Real spaceLoc  =  xb[idir] + distDir;
          intersectLoc[iinter][idir] = spaceLoc;
          a_distanceAlongLine[iinter] += distDir*distDir;

          //the nMaxDir value is set with RealVect::Zero initialization
          //no extrapolation in the nmax dir
          Real ccLoc = a_dx[idir]*(Real(intersectIV[iinter][idir]) + 0.5);
          extrapDist[idir] = spaceLoc - ccLoc;
        }
      }

      a_distanceAlongLine[iinter] = sqrt(a_distanceAlongLine[iinter]);

      int order = getExtrapolationStencil(a_pointStencils[iinter], extrapDist,
                                          a_dx, centerVoF, a_ebisBox, 2, nMaxDir,
                                          NULL, a_ivar);

      //the returned order is the order of the derivs taken.  can tolerate 1 or 2
      if (order == 0)
      {
        a_dropOrder = true;
        return;
      }
    }//end loop over intersection points

    return;
  }
  /****/
  int
  EBArith::
  getExtrapolationStencil(VoFStencil&     a_stencil,
                          const RealVect& a_dist,
                          const RealVect& a_dx,
                          const VolIndex& a_startVoF,
                          const EBISBox&  a_ebisBox,
                          int a_orderOfPolynomial,
                          int a_noExtrapThisDir,
                          IntVectSet*    a_cfivsPtr,
                          int ivar)
  {
    BL_PROFILE("EBArith::getExtrapolationStencil");
    int order = 2;

    a_stencil.clear();

    //zeroth order Taylor series
    a_stencil.add(a_startVoF,1.0,ivar);

    //do taylor series extrapolation stencil
    //get stencils for derivatives
    int derivorder;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if ((idir != a_noExtrapThisDir) && (a_orderOfPolynomial > 0))
      {
        VoFStencil firstDSten;
        derivorder = EBArith::getFirstDerivStencil(firstDSten,
                                                   a_startVoF, a_ebisBox,
                                                   idir, a_dx[idir], a_cfivsPtr, ivar);
        order = std::min(order, derivorder);

        firstDSten *= a_dist[idir];

        a_stencil += firstDSten;

        if (a_orderOfPolynomial > 1)
        {
          VoFStencil secondDSten;
          derivorder = EBArith::getSecondDerivStencil(secondDSten,
                                                      a_startVoF, a_ebisBox,
                                                      idir, a_dx[idir], a_cfivsPtr, ivar);
          order = std::min(order, derivorder);
          secondDSten *= (0.5*a_dist[idir]*a_dist[idir]);
          a_stencil += secondDSten;

#if AMREX_SPACEDIM==3
          int dir1, dir2;
          if(idir == 0)
          {
            dir1 = 1;
            dir2 = 2;
          }
          else if(idir == 1)
          {
            dir1 = 0;
            dir2 = 2;
          }
          else
          {
            dir1 = 0;
            dir2 = 1;
          }

          if ((dir1 != a_noExtrapThisDir) && (dir2 != a_noExtrapThisDir))
          {
            VoFStencil mixedDSten;
            derivorder = EBArith::getMixedDerivStencil(mixedDSten,
                                                       a_startVoF, a_ebisBox,
                                                       dir1, dir2, a_dx[dir1], a_dx[dir2], a_cfivsPtr, ivar);
            order = std::min(order, derivorder);
            mixedDSten *= (a_dist[dir1]*a_dist[dir2]);
            a_stencil += mixedDSten;
          }
#endif
        }
      }
    }

#if AMREX_SPACEDIM==2
    if (a_orderOfPolynomial > 1)
    {
      int dir1 = 0;
      int dir2 = 1;
      if ((dir1 != a_noExtrapThisDir) && (dir2 != a_noExtrapThisDir))
      {
        VoFStencil mixedDSten;
        derivorder = EBArith::getMixedDerivStencil(mixedDSten,
                                                   a_startVoF, a_ebisBox,
                                                   dir1, dir2, a_dx[dir1], a_dx[dir2], a_cfivsPtr, ivar);
        order = std::min(order, derivorder);
        mixedDSten *= (a_dist[dir1]*a_dist[dir2]);
        a_stencil += mixedDSten;
      }
    }
#endif

    return order;
  }
  /****/
  void
  EBArith::johanStencil(bool&               a_dropOrder,
                        Vector<VoFStencil>& a_pointStencils,
                        Vector<Real>&       a_distanceAlongLine,
                        const VolIndex&     a_vof,
                        const EBISBox&      a_ebisBox,
                        const RealVect&     a_dx,
                        const IntVectSet&   a_cfivs,
                        int a_ivar)
  {
    //BL_PROFILE("EBArith::johanStencil");
    //IntVect iv = a_vof.gridIndex();

    // Do Johansen-Colella cast-a-ray algorithm
    RealVect bndryCentroid = a_ebisBox.bndryCentroid(a_vof);
    RealVect normal     = a_ebisBox.normal(a_vof);
    //splitting up stuff this way to facillitate multifluid
    //which can have multiple normals and boundary centroids per cell.
    johanStencil(a_dropOrder, a_pointStencils, a_distanceAlongLine,
                 normal, bndryCentroid,
                 a_vof, a_ebisBox, a_dx, a_cfivs, a_ivar);
  }
  void
  EBArith::johanStencil(bool&               a_dropOrder,
                        Vector<VoFStencil>& a_pointStencils,
                        Vector<Real>&       a_distanceAlongLine,
                        const RealVect&     a_normal,
                        const RealVect&     a_bndryCentroid,
                        const VolIndex&     a_vof,
                        const EBISBox&      a_ebisBox,
                        const RealVect&     a_dx,
                        const IntVectSet&   a_cfivs,
                        int a_ivar)
  {
    int numExtrapPoints = 2;
    EBArith::dataRayCast(a_dropOrder, a_pointStencils, a_distanceAlongLine, a_normal, a_bndryCentroid,
                         a_vof, a_ebisBox, a_dx, a_cfivs, a_ivar, numExtrapPoints);
  }
}
