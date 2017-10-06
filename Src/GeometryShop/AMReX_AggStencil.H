/*
 *      .o.       ooo        ooooo ooooooooo.             ooooooo  ooooo 
 *     .888.      `88.       .888' `888   `Y88.            `8888    d8'  
 *    .8"888.      888b     d'888   888   .d88'  .ooooo.     Y888..8P    
 *   .8' `888.     8 Y88. .P  888   888ooo88P'  d88' `88b     `8888'     
 *  .88ooo8888.    8  `888'   888   888`88b.    888ooo888    .8PY888.    
 * .8'     `888.   8    Y     888   888  `88b.  888    .o   d8'  `888b   
 *o88o     o8888o o8o        o888o o888o  o888o `Y8bod8P' o888o  o88888o 
 *
 */
#ifndef _AGGSTENCIL_H_
#define _AGGSTENCIL_H_

#include "AMReX_Stencils.H"
#include "AMReX_BaseIndex.H"
#include "AMReX_REAL.H"

namespace amrex
{
/// Aggregated stencil
/**
   srcData_t classes need the following functions
   int   dataType(BaseIndex);
   int     offset(BaseIndex);
   Real* dataPtr(int dataType, int ivar);

   dstData_t classes need the following functions
   int   dataType(BaseIndex);
   int     offset(BaseIndex);
   Real* dataPtr(int dataType, int ivar);

   sten_t classes need the following functions
   srcIndex_t index(int isten)
   Real       weight(int isten)
 */
template <class srcData_t, class dstData_t>
class AggStencil
{
public:
  ///
  /**
     a_dstIndex is where the answers will go.
     a_srcStencil is the stencil into the src data
     a_srcData is to provide offsets for source data (can be a dummy but must be the right size)
     a_dstData is to provide offsets for destin data (can be a dummy but must be the right size)
   */
  AggStencil(const Vector<std::shared_ptr<BaseIndex  > >& a_dstVoFs,
             const Vector<std::shared_ptr<BaseStencil> >& a_stencil,
             const srcData_t                                 & a_srcData,
             const dstData_t                                 & a_dstData);

  ///
  ~AggStencil()
  {
  }

  ///
  void apply(dstData_t       & a_lph,
             const srcData_t & a_phi,
             const int       & a_varDest,
             const bool      & a_incrementOnly) const;

  ///
  /**
     for applying a scalar operator to multiple operators.
     Sets the initial variable for destination to a_varDest
     Sets dataPtr for a_phi initial variable to a_varSrc
     Runs over ncomp components
   */
  void apply(dstData_t       & a_lph,
             const srcData_t & a_phi,
             const int       & a_src,
             const int       & a_dst,
             const int       & a_nco,
             const bool      & a_incrementOnly) const;

  ///
  void   cache(const dstData_t& a_lph) const;

  ///
  void uncache(      dstData_t& a_phi) const;

  struct
  {
    size_t offset;
    int  dataID;
  } typedef access_t;

  typedef Vector<std::pair<access_t, Real> > stencil_t;

protected:

  int m_destVar;
  Vector<stencil_t>   m_ebstencil;
  Vector<access_t>    m_dstAccess;
  mutable Vector< Vector<Real> > m_cacheDst;

private:
  /// disallowed operators.   Without code because Jeff says that is better.

  ///weak construction is bad.
  AggStencil();

  ///deep copy for this object would kill performance
  AggStencil& operator=(const AggStencil& stenin);

  ///no reason to allow this one either.
  AggStencil(const AggStencil& stenin);

};

}
#include "AMReX_AggStencilI.H"
#endif
