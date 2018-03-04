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


#ifndef AMREX_BOXITERATOR_H_
#define AMREX_BOXITERATOR_H_

#include <cstdlib>
#include "AMReX_BLassert.H"
#include "AMReX_Box.H"
#include "AMReX_REAL.H"
#include "AMReX_SPACE.H"
#include "AMReX_IntVect.H"

namespace amrex
{
  ///iterates through the IntVects of a Box
  /**
     BoxIterator iterates through the IntVects of a box.  The actual
     sqeuence of IntVects is implementation-specific.
     Typical usage:

     Box b;
     ...
     BoxIterator bit (b);
     for (bit.begin(); bit.ok(); ++bit)
     {
     IntVect iv = bit();
     (do operations involving iv)
     }
  */
  class BoxIterator
  {
  public:
    ///
    /**
       Default constructor.  This constructs an invalid iterator.
       The user must call define before using.
    */
    BoxIterator();

    ///
    /**
       Constructs a BoxIterator and associates it with a Box.
       Arguments:
       a_bx (not modified) the Box to iterate over.
    */
    BoxIterator(const Box& a_bx);

    void setBox(const Box& a_bx);

    ///
    /**
       Associates a Box with this BoxIterator.
       Arguments:
       a_bx (not modified) the Box to iterate over.
    */
    void define(const Box& a_bx);

    ///
    /**
       Copy constructor.
       Arguments:
       a_iterIn (not modified) the BoxIterator to copy.
    */
    BoxIterator(const BoxIterator& a_iterIn);

    ///
    ~BoxIterator ()
    {
    }

    ///
    /**
       Sets this BoxIterator to the first IntVect in its Box.  The
       definition of the "first" IntVect is
       implementation-dependent.
    */
    void begin();

    ///
    /**
       Sets this BoxIterator to the first IntVect in its Box.  The
       definition of the "first" IntVect is
       implementation-dependent.
    */
    void reset();

    ///
    /**
       Modifies this BoxIterator to set it to the next location in its
       Box.  The definition of the "next location" of a Box is
       implementation-dependent.
    */
    void operator ++ ();

    void next();

    ///
    /**
       Returns the value of the InVect for the current location of this
       BoxIterator.
    */
    const IntVect& operator () () const;

    ///
    /**
       Returns true if this BoxIterator's location is within its Box.
    */
    bool ok();

  protected:
    IntVect m_current;
    IntVect m_boxLo;
    IntVect m_boxHi;
  };

  inline
  BoxIterator::BoxIterator()
  {
    m_current = IntVect::TheUnitVector();
    m_boxLo   = IntVect::TheUnitVector();
    m_boxHi   = IntVect::TheZeroVector();
  }

  inline
  BoxIterator::BoxIterator(const Box& a_bx)
  {
    define(a_bx);
  }

  inline
  BoxIterator::BoxIterator(const BoxIterator& a_iterIn)
  {
    m_current = a_iterIn.m_current;
    m_boxLo =   a_iterIn.m_boxLo;
    m_boxHi =   a_iterIn.m_boxHi;
  }

  inline
  void BoxIterator::begin()
  {
    if (m_boxLo <= m_boxHi)  m_current = m_boxLo;
  }

  inline
  void BoxIterator::reset()
  {
    begin();
  }

  inline
  void BoxIterator::operator ++ ()
  {
    next();
  }

  inline
  void BoxIterator::next()
  {
    m_current[0]++;
#if AMREX_SPACEDIM >= 2
    if (m_current[0] > m_boxHi[0])
      {
        m_current[0] = m_boxLo[0];
        m_current[1]++;
#if AMREX_SPACEDIM >= 3
        if (m_current[1] > m_boxHi[1])
          {
            m_current[1] = m_boxLo[1];
            m_current[2]++;
          }
#endif
      }      
#endif
  }

  inline
  const IntVect& BoxIterator::operator () () const
  {
    BL_ASSERT(m_current <= m_boxHi);
    BL_ASSERT(m_current >= m_boxLo);
    return m_current;
  }

  inline
  bool BoxIterator::ok()
  {
    return (m_current  <= m_boxHi);
  }
}
#endif
