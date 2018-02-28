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

#include "AMReX_EBLoHiCenter.H"


namespace amrex
{
            
  void
  EBLoHiCenAllDirs(Box         a_loBox[AMREX_SPACEDIM],
                   Box         a_hiBox[AMREX_SPACEDIM],
                   Box         a_ceBox[AMREX_SPACEDIM],
                   int         a_hasLo[AMREX_SPACEDIM],
                   int         a_hasHi[AMREX_SPACEDIM],
                   int       & a_nearAnyBoundary,
                   const Box & a_inBox,
                   const Box & a_domain)
  {
    Box entirebox; //never saw the use in this one
    a_nearAnyBoundary = 0;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      eblohicenter(a_loBox[idir], a_hasLo[idir],
                   a_hiBox[idir], a_hasHi[idir],
                   a_ceBox[idir], entirebox,
                   a_inBox, a_domain, idir);

      if((a_hasLo[idir]==1) || (a_hasHi[idir]==1))
      {
        a_nearAnyBoundary = 1;
      }
    }
  }


  void getLHCBoxes(Box& a_loBox,
                   Box& a_hiBox,
                   Box& a_centerBox,
                   int& a_hasLo,
                   int& a_hasHi,
                   const Box& a_inBox,
                   const Box& a_dProblem,
                   const int& a_dir)
  {
    Box inputBox= a_inBox;
    inputBox &= a_dProblem;

    Box interiorDomain = a_dProblem;
    interiorDomain.grow(a_dir, -1);
    //the centerd diff thing is over cells not bordering domain
    a_centerBox = inputBox;
    a_centerBox &=  interiorDomain;
    Box tmpBox;
            
    //check the high side thing
    tmpBox = inputBox;
    tmpBox.shift(a_dir,1);
    tmpBox &= a_dProblem;
    if (tmpBox != inputBox)
    {
      a_hasHi = 1;
      a_hiBox = adjCellHi(tmpBox, a_dir);
    }
    else
    {
      a_hasHi = 0;
    }
            
    //check the low side thing
    tmpBox = inputBox;
    tmpBox.shift(a_dir,-1);
    tmpBox &= a_dProblem;
    if (tmpBox != inputBox)
    {
      a_hasLo = 1;
      a_hiBox = adjCellLo(tmpBox, a_dir);
    }
    else
    {
      a_hasHi = 0;
    }
  }
            
  void eblohicenter(Box&                 a_loBox,
                    int&                 a_hasLo,
                    Box&                 a_hiBox,
                    int&                 a_hasHi,
                    Box&                 a_centerBox,
                    Box&                 a_entireBox,
                    const Box&           a_inBox,
                    const Box& a_domain,
                    const int&           a_dir)
  {
    // Make a copy of the input box which can be modified
    Box inBox = a_inBox;
    inBox &= a_domain;
            
    // The centered difference box is always one smaller in a_dir
    a_centerBox = inBox;
    a_centerBox.grow(a_dir,-1);
            
            
    // The union of all the output boxes start off equal to the center
    // difference portion on the input box (intersected with the domain)
    a_entireBox = a_centerBox;
            
    // See if this chops off the high side of the input box
    Box tmp = inBox;
    tmp.shift(a_dir,1);
    tmp &= a_domain;
    tmp.shift(a_dir,-1);
            
    // If so, set up the high, one-sided difference box, a_hiBox, and expand
    // the entire box to include it
    if (tmp != inBox)
    {
      a_hasHi = 1;
      a_hiBox = adjCellHi(tmp,a_dir);
      a_entireBox.growHi(a_dir,1);
    }
    else
    {
      a_hasHi = 0;
    }
            
    // See if this chops off the low side of the input box
    tmp = inBox;
    tmp.shift(a_dir,-1);
    tmp &= a_domain;
    tmp.shift(a_dir,1);
            
    // If so, set up the low, one-sided difference box, a_loBox, and expand
    // the entire box to include it
    if (tmp != inBox)
    {
      a_hasLo = 1;
      a_loBox = adjCellLo(tmp,a_dir);
      a_entireBox.growLo(a_dir,1);
    }
    else
    {
      a_hasLo = 0;
    }
            
    // Make some simple sanity checks
    BL_ASSERT(a_entireBox.contains(a_centerBox));
            
    if (a_hasLo == 1)
    {
      BL_ASSERT(a_entireBox.contains(a_loBox));
    }
            
    if (a_hasHi == 1)
    {
      BL_ASSERT(a_entireBox.contains(a_hiBox));
    }
  }
            
  // This function is used when in direction a_dir a 2 point stencil of cell-
  // centered data is being used to compute something at the cell face between
  // the cell centers of the stencil.  The data for the stencil is valid in
  // a_inBox.  It uses a_inBox to compute a box (face-centered in a_dir) where
  // the full stencil can be used, a_centerBox, and boxes (face-centered in
  // a_dir) where a 1 point stencil can be used, a_loBox and a_hiBox based on
  // the current problem domain, a_domain, and the stencil direction, a_dir.
  // The union of these 1 and 2 point stencel boxes is returned as a_entireBox
  // (face-centered in a_dir).  The 1 point stencil boxes are one wide, at most,
  // and if they have been defined then the corresponding flag, a_hasLo or
  // a_hasHi, is set to one, otherwise these flags are zero.  All output boxes
  // lie within the domain.
            
  void eblohicenterFace(Box&                 a_loBox,
                        int&                 a_hasLo,
                        Box&                 a_hiBox,
                        int&                 a_hasHi,
                        Box&                 a_centerBox,
                        Box&                 a_entireBox,
                        const Box&           a_inBox,
                        const Box& a_domain,
                        const int&           a_dir)
  {
    // Make a copy of the input box which can be modified
    Box inBox = a_inBox;
    inBox &= a_domain;
            
    // The centered difference box is always one smaller in a_dir
    a_centerBox = inBox;
    a_centerBox.surroundingNodes(a_dir);
    a_centerBox.grow(a_dir,-1);
            
    // The union of all the output boxes start off equal to the center
    // difference portion on the input box (intersected with the domain)
    a_entireBox = a_centerBox;
            
    // See if this chops off the high side of the input box
    Box tmp = inBox;
    tmp.shift(a_dir,1);
    tmp &= a_domain;
    tmp.shift(a_dir,-1);
            
    // If so, set up the high, one-sided difference box, a_hiBox, and expand
    // the entire box to include it
    if (tmp != inBox)
    {
      a_hasHi = 1;
      a_hiBox = adjCellHi(tmp,a_dir);
      a_hiBox.shiftHalf(a_dir,1);
      a_entireBox.growHi(a_dir,1);
    }
    else
    {
      a_hasHi = 0;
    }
            
    // See if this chops off the low side of the input box
    tmp = inBox;
    tmp.shift(a_dir,-1);
    tmp &= a_domain;
    tmp.shift(a_dir,1);
            
    // If so, set up the low, one-sided difference box, a_loBox, and expand
    // the entire box to include it
    if (tmp != inBox)
    {
      a_hasLo = 1;
      a_loBox = adjCellLo(tmp,a_dir);
      a_loBox.shiftHalf(a_dir,-1);
      a_entireBox.growLo(a_dir,1);
    }
    else
    {
      a_hasLo = 0;
    }
            
    // Make some simple sanity checks
    BL_ASSERT(a_entireBox.contains(a_centerBox));
            
    if (a_hasLo == 1)
    {
      BL_ASSERT(a_entireBox.contains(a_loBox));
    }
            
    if (a_hasHi == 1)
    {
      BL_ASSERT(a_entireBox.contains(a_hiBox));
    }
  }
}
            
