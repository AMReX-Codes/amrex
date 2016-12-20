#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOXITERATOR_H_
#define _BOXITERATOR_H_

#include <cstdlib>

#include "Box.H"
#include "REAL.H"
#include "SPACE.H"
#include "IntVect.H"
#include "NamespaceHeader.H"

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
  m_current = IntVect::Unit;
  m_boxLo   = IntVect::Unit;
  m_boxHi   = IntVect::Zero;
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
#if CH_SPACEDIM >= 2
  if (m_current[0] > m_boxHi[0])
    {
      m_current[0] = m_boxLo[0];
      m_current[1]++;
#if CH_SPACEDIM >= 3
      if (m_current[1] > m_boxHi[1])
        {
          m_current[1] = m_boxLo[1];
          m_current[2]++;
#if CH_SPACEDIM >= 4
          if (m_current[2] > m_boxHi[2])
            {
              m_current[2] = m_boxLo[2];
              m_current[3]++;
#if CH_SPACEDIM >= 5
              if (m_current[3] > m_boxHi[3])
                {
                  m_current[3] = m_boxLo[3];
                  m_current[4]++;
#if CH_SPACEDIM == 6
                  if (m_current[4] > m_boxHi[4])
                    {
                      m_current[4] = m_boxLo[4];
                      m_current[5]++;
                    }
#elif CH_SPACEDIM > 6
                  SPACEDIM > 6 undefined!;
#endif
                }
#endif
            }
#endif
        }
#endif
    }
#endif
}

inline
const IntVect& BoxIterator::operator () () const
{
  CH_assert(m_current <= m_boxHi);
  CH_assert(m_current >= m_boxLo);
  return m_current;
}

inline
bool BoxIterator::ok()
{
  return (m_current  <= m_boxHi);
}

#include "NamespaceFooter.H"
#endif
