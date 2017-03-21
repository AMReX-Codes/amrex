#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#ifndef _SLABSERVICE_H_
#define _SLABSERVICE_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "Box.H"
#include "EBISBox.H"
#include "GeometryService.H"
#include "NamespaceHeader.H"
///
/**
   This is a pure base class that EBIndexSpace uses
   to interface with geometry-generation classes.
 */
class SlabService: public GeometryService
{
public:
  ///
  SlabService(const Box& a_coveredRegion);

  virtual ~SlabService();

  ///
  /**
      Return true if every cell in region is regular at the
      refinement described by dx.
  */
  virtual bool isRegular(const Box&      a_region,
                         const ProblemDomain&      a_domain,
                         const RealVect& a_origin,
                         const Real&     a_dx) const;

  ///
  /**
      Return true if every cell in region is covered at the
      refinement described by dx.
  */
  virtual bool isCovered(const Box&      a_region,
                         const ProblemDomain&      a_domain,
                         const RealVect& a_origin,
                         const Real&     a_dx) const;

  ///
  /**
     irregGraph contains a complete list of irregular vofs. \\
     validRegion is the the region of the graph \\
     domain is the domain of compuation \\
     regIrregCovered = 1 for regular, 0 for irregular, -1 for covered.
     regIrregCovered must contain the valid region grown by 1. \\
     In chombospeak, \\
     ghostregion = (grow(validRegion, 1) & domain);
   */
  virtual void fillGraph(BaseFab<int>&        a_regIrregCovered,
                         Vector<IrregNode>&   a_nodes,
                         const Box&           a_validRegion,
                         const Box&           a_ghostRegion,
                         const ProblemDomain&           a_domain,
                         const RealVect&      a_origin,
                         const Real&          a_dx) const;


  virtual InOut InsideOutside(const Box&           a_region,
                              const ProblemDomain& a_domain,
                              const RealVect&      a_origin,
                              const Real&          a_dx) const ;

protected:
  Box m_coveredRegion;

};
#include "NamespaceFooter.H"
#endif
