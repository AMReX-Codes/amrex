#ifndef _EBLOHICENTER_H___
#define _EBLOHICENTER_H___

#include "AMReX_Box.H"

namespace amrex
{

  ///
  /**
     divide a box into its centered-diff boxes, its
     one-sided   boxes and tell whether the one sided boxes are there.
     The one-sided boxes are one wide at most.
     All Boxes are intersected with the domain.
     InBox should be the size of the domain of computation
     (the size of the gradient box).
  */
  extern void getLHCBoxes(Box        & a_loBox,
                          Box        & a_hiBox,
                          Box        & a_centerBox,
                          int        & a_haslo,
                          int        & a_hashi,
                          const Box  & a_inBox,
                          const Box  & a_dProblem,
                          const int  & a_dir);
                 
  ///
  /**
     Divide a box, a_inBox, into a box where centered differences can be used,
     a_centerBox, and boxes where one-sided difference can be used, a_loBox and
     a_hiBox based on the current problem domain, a_domain, and the difference
     direction, a_dir.  The union of these computation boxes are returned as
     a_entireBox.  The one-sided difference boxes are one wide at most and if
     they have been defined then the corresponding flag, a_hasLo or a_hasHi, is
     set to one, otherwise it is zero.  All boxes lie within the domain,
     a_domain.
  */
  extern void eblohicenter(Box       & a_loBox,
                           int       & a_hasLo,
                           Box       & a_hiBox,
                           int       & a_hasHi,
                           Box       & a_centerBox,
                           Box       & a_entireBox,
                           const Box & a_inBox,
                           const Box & a_domain,
                           const int & a_dir);
                 


  ///
  extern void EBLoHiCenAllDirs(Box         a_loBox [AMREX_SPACEDIM],
                               Box         a_hiBox [AMREX_SPACEDIM],
                               Box         a_cenBox[AMREX_SPACEDIM],
                               int         a_hasLo [AMREX_SPACEDIM],
                               int         a_hasHi [AMREX_SPACEDIM],
                               int       & a_nearAnyBoundary,
                               const Box & a_inBox,
                               const Box & a_domain);
                 
  ///
  /**
     This function is used when in direction a_dir a 2 point stencil of cell-
     centered data is being used to compute something at the cell face between
     the cell centers of the stencil.  The data for the stencil is valid in
     a_inBox.  It uses a_inBox to compute a box (face-centered in a_dir) where
     the full stencil can be used, a_centerBox, and boxes (face-centered in
     a_dir) where a 1 point stencil can be used, a_loBox and a_hiBox based on
     the current problem domain, a_domain, and the stencil direction, a_dir.
     The union of these 1 and 2 point stencel boxes is returned as a_entireBox
     (face-centered in a_dir).  The 1 point stencil boxes are one wide, at most,
     and if they have been defined then the corresponding flag, a_hasLo or
     a_hasHi, is set to one, otherwise these flags are zero.  All output boxes
     lie within the domain.
  */
  extern void eblohicenterFace(Box       & a_loBox,
                               int       & a_hasLo,
                               Box       & a_hiBox,
                               int       & a_hasHi,
                               Box       & a_centerBox,
                               Box       & a_entireBox,
                               const Box & a_inBox,
                               const Box & a_domain,
                               const int & a_dir);
                 
} 
#endif
