#ifndef AMREX_BOXITERATOR_H_
#define AMREX_BOXITERATOR_H_
#include <AMReX_Config.H>

#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_IntVect.H>

#include <cstdlib>

namespace amrex
{
  /**
     \brief iterates through the IntVects of a Box


     BoxIterator iterates through the IntVects of a box.  The actual
     sequence of IntVects is implementation-specific.
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
    BoxIterator () noexcept = default;

    ///
    /**
       Constructs a BoxIterator and associates it with a Box.
       Arguments:
       a_bx (not modified) the Box to iterate over.
    */
    explicit BoxIterator (const Box& a_bx) noexcept;

    void setBox (const Box& a_bx) noexcept;

    ///
    /**
       Associates a Box with this BoxIterator.
       Arguments:
       a_bx (not modified) the Box to iterate over.
    */
    void define (const Box& a_bx) noexcept;

    ///
    /**
       Sets this BoxIterator to the first IntVect in its Box.  The
       definition of the "first" IntVect is
       implementation-dependent.
    */
    void begin () noexcept;

    ///
    /**
       Sets this BoxIterator to the first IntVect in its Box.  The
       definition of the "first" IntVect is
       implementation-dependent.
    */
    void reset () noexcept;

    ///
    /**
       Modifies this BoxIterator to set it to the next location in its
       Box.  The definition of the "next location" of a Box is
       implementation-dependent.
    */
    void operator ++ () noexcept;

    void next () noexcept;

    ///
    /**
       Returns the value of the InVect for the current location of this
       BoxIterator.
    */
    [[nodiscard]] const IntVect& operator () () const noexcept;

    ///
    /**
       Returns true if this BoxIterator's location is within its Box.
    */
    [[nodiscard]] bool ok () noexcept;

  protected:
    IntVect m_current = IntVect::TheUnitVector();
    IntVect m_boxLo   = IntVect::TheUnitVector();
    IntVect m_boxHi   = IntVect::TheZeroVector();
  };

  inline
  BoxIterator::BoxIterator (const Box& a_bx) noexcept
  {
    define(a_bx);
  }

  inline
  void BoxIterator::begin () noexcept
  {
    if (m_boxLo <= m_boxHi)  m_current = m_boxLo;
  }

  inline
  void BoxIterator::reset () noexcept
  {
    begin();
  }

  inline
  void BoxIterator::operator ++ () noexcept
  {
    next();
  }

  inline
  void BoxIterator::next () noexcept
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
  const IntVect& BoxIterator::operator () () const noexcept
  {
    BL_ASSERT(m_current <= m_boxHi);
    BL_ASSERT(m_current >= m_boxLo);
    return m_current;
  }

  inline
  bool BoxIterator::ok () noexcept
  {
    return (m_current  <= m_boxHi);
  }
}
#endif
