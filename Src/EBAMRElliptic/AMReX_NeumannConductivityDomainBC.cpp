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

#include "AMReX_NeumannConductivityDomainBC.H"

namespace amrex
{
  void 
  NeumannConductivityDomainBC::
  fillPhiGhost(FArrayBox&     a_phi,
               const Box&     a_valid,
               const Box&     a_domain,
               Real           a_dx,
               bool           a_homogeneous)
  {
    Box grownBox = a_valid;
    grownBox.grow(1);

    for (int idir=0; idir<CH_SPACEDIM; ++idir)
    {
      for(SideIterator sit; sit.ok(); ++sit)
      {
        Box choppedBox = grownBox;
        choppedBox.grow(idir,-1);
        Box toRegion = adjCellBox(choppedBox, idir, sit(), 1);

        if(!a_domain.contains(toRegion))
        {
          for (BoxIterator bit(toRegion); bit.ok(); ++bit)
          {
            const IntVect& iv = bit();
            //fake vof just to get the location
            VolIndex vof(iv, 0);
            RealVect loc = EBArith::getVoFLocation(vof, a_dx, RealVect::Zero);
            int isign = sign(sit());
            IntVect ivneigh = iv - isign*BASISV(idir);
            Real value = bcvaluefunc(loc, idir, sit());
            if(a_homogeneous) value = 0;
            a_phi(iv, 0) = a_phi(ivneigh, 0)  + a_dx*value;
          }
        }
      } 
    }//end loop over directions
  }
  /*****/
  void
  NeumannConductivityDomainBC::
  getFaceFlux(Real&                 a_faceFlux,
              const VolIndex&       a_vof,
              const MFIter    &     a_mfi,
              const EBCellFAB&      a_phi,
              const int&            a_idir,
              const Side::LoHiSide& a_side,
              const bool&           a_useHomogeneous)
  {
    const EBISBox& ebisBox = a_phi.getEBISBox();

    Real totalMassFlux = 0.0;
    Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
    for (int i = 0; i < faces.size(); i++)
    {
      const RealVect centroid = ebisBox.centroid(faces[i]);
      Real thisFaceFlux;
      getFaceGradPhi(thisFaceFlux,faces[i],mfi, a_phi, a_idir, a_side, a_useHomogeneous);
      totalMassFlux += thisFaceFlux;
    }
    if(faces.size() > 1)
    {
      totalMassFlux /= faces.size();
    }
    a_faceFlux = totalMassFlux;
   
    Real bcoave = 0;
    Real areaTot = 0.0;
    for (int iface = 0; iface < faces.size(); iface++)
    {
      Real areaFrac = ebisBox.areaFrac(faces[iface]);
      areaTot += areaFrac;
      Real bcoFace  = (*m_bcoef)[a_dit][a_idir](faces[iface], 0);
      bcoave += areaFrac*bcoFace;
    }
    if (areaTot > 1.0e-8)
    {
      bcoave /= areaTot;
    }
    a_faceFlux *= bcoave;
  }

  /*****/
  void
  NeumannConductivityDomainBC::
  getFaceGradPhi(Real&                 a_faceFlux,
                 const FaceIndex&      a_face,
                 const MFIter    &     a_mfi,
                 const EBCellFAB&      a_phi,
                 const int&            a_idir,
                 const Side::LoHiSide& a_side,
                 const bool&           a_useHomogeneous)
  {
    const int iside = -sign(a_side);

    const EBISBox& ebisBox = a_phi.getEBISBox();

    Real flux = -1.e99;
    if (a_useHomogeneous)
    {
      flux = 0.0;
    }
    else if (m_isFunction)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point  = a_centroid;
      point *= a_dx;
      point[a_face.direction()] = 0.0;//make this not depend on whatever ebisbox is returning for centroid in the face direction.
      point += EBArith::getFaceLocation(a_face,a_dx,a_probLo);
      flux = bcvaluefunc(point, a_idir, a_side);
    }
    else
    {
      flux = m_value;
    }
    a_faceFlux = iside*flux;
  }

}
