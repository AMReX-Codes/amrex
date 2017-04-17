#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTG  7-31-2001

#ifndef _EBFLUXFAB_H_
#define _EBFLUXFAB_H_

#include "Box.H"
#include "Vector.H"
#include "EBFaceFAB.H"
#include "NamespaceHeader.H"

/// A EBFaceFAB-like container for edge-centered fluxes
/**
    This is a class to contain face-centered fluxes on a box.
 */
class EBFluxFAB
{

public:

 ///
  /**
     This stuff required by LevelData in parallel:
  */
  static int preAllocatable()
  {
    return 1; // symmetric messaging
  }

  /// default constructor
  EBFluxFAB();

  /// destructor
  ~EBFluxFAB();

  /// this function returns the EBFluxFAB to the undefined state
  void clear();

  ///
  EBFluxFAB& operator+=(const EBFluxFAB& a_ebfab)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        (*m_fluxes[idir]) += (*a_ebfab.m_fluxes[idir]);
      }
    return *this;
  }


  ///
  EBFluxFAB& operator-=(const EBFluxFAB& a_ebfab)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        (*m_fluxes[idir]) -= (*a_ebfab.m_fluxes[idir]);
      }
    return *this;
  }


  ///
  EBFluxFAB& operator*=(const EBFluxFAB& a_ebfab)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        (*m_fluxes[idir]) *= (*a_ebfab.m_fluxes[idir]);
      }
    return *this;
  }


  ///
  EBFluxFAB& operator/=(const EBFluxFAB& a_ebfab)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        (*m_fluxes[idir]) /= (*a_ebfab.m_fluxes[idir]);
      }
    return *this;
  }
  
  const Box& box() const
  {
    return getRegion();
  }

  void define(const EBFluxFAB& a_fb)
  {
    define(a_fb.getEBISBox(), a_fb.box(), a_fb.nComp());
  }

  //! Returns the maximum value for the given component in the FAB for
  //! the given direction.
  //! \param a_dir The direction in which the maximum value will be returned.
  //! \param a_comp The component whose maximum value will be returned.
  Real max(int a_dir, int a_comp = 0) const
  {
     return m_fluxes[a_dir]->max(a_comp);
  }

  //! Returns the minimum value for the given component in the FAB for
  //! the given direction.
  //! \param a_dir The direction in which the minimum value will be returned.
  //! \param a_comp The component whose minimum value will be returned.
  Real min(int a_dir, int a_comp = 0) const
  {
     return m_fluxes[a_dir]->min(a_comp);
  }

  /// {\bf access functions}

  /// number of components
  int nComp() const;

  /// returns cell-centered box which defines fluxBox
  const Box& getRegion() const;

  /// returns EBFaceFAB in direction dir
  EBFaceFAB& operator[] (const int dir);

  /// constant version
  const EBFaceFAB& operator[] (const int dir) const;

  void setVal(const Real& val);

  ///defines this from input and copies the values
  void clone(const EBFluxFAB& a_input);

  ///
  /** Copies from a subsection of one box into another.
      Assumes the boxes are both in the same index space, and that
      box R is completely contained in both the src and destination
      boxes.  This function equired by BoxLayoutData */
  void copy(const Box& Rfrom,
            const Interval& Cdest,
            const Box& Rto,
            const EBFluxFAB& src,
            const Interval& Csrc);

  ///
  bool isDefined() const;

  ///
  EBFluxFAB(const EBISBox& a_ebisBox, const Box& bx, int n);

  ///
  void define(const EBISBox& a_ebisBox, const Box& bx, int n);

  int size(const Box& R, const Interval& comps) const;

  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  void linearIn(void* buf, const Box& R, const Interval& comps);

  ///
  /**
     Can only be used for an undefined ebfluxfab.
     Turns off deletion of pointers.
  */
  void alias(Vector<EBFaceFAB*> a_inputFAB);

  /// Invalid but necessary for LevelData<EBFluxFAB> to compile
  EBFluxFAB(const Box& bx, int n)
  {
    MayDay::Error("invalid constructor called for EBFluxFAB");
  }
  EBISBox getEBISBox() const
  {
    return m_ebisBox;
  }

protected:
  void setDefaultValues();

  ///true if facefabs set by hand, false otherwise
  bool m_aliased;
  ///
  Box m_region;

  ///
  EBISBox m_ebisBox;

  ///
  int m_nComp;

  bool m_isDefined;

  /// CH_SPACEDIM EBFaceFABes which hold fluxes
  Tuple<EBFaceFAB*, CH_SPACEDIM> m_fluxes;


private:
  /// these are disallowed

  EBFluxFAB (const EBFluxFAB&)
  {
    MayDay::Error("bogus constructor");
  }

  void operator = (const EBFluxFAB&)
  {
    MayDay::Error("bogus constructor");
  }

};

#include "NamespaceFooter.H"
#endif
