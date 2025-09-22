#ifndef AMREX_BOXITERATOR_H_
#define AMREX_BOXITERATOR_H_
#include <AMReX_Config.H>

#include <AMReX_BLassert.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_IntVect.H>

#include <cstdlib>

namespace amrex
{

template<int dim>
class BoxND;

/**
* \brief iterates through the IntVects of a Box
*
* BoxIteratorND iterates through the IntVects of a box
* in the same order as a series of nested for loops
* with the innermost loop iterating over the first dimension.
*
* Typical usage using a range-based for loop:
* \code{.cpp}
*   Box b;
*   ...
*   for (IntVect iv : b.iterator()) {
*       // do operations involving iv
*   }
* \endcode
*
* Using a traditional for loop:
* \code{.cpp}
*   Box b;
*   ...
*   BoxIterator bit (b);
*   for (bit.begin(); bit.ok(); ++bit)
*   {
*       IntVect iv = bit();
*       // do operations involving iv
*   }
* \endcode
*
* Note that while regular iteration of the range-based version is reasonably performant,
* neither version is compatible with omp parallel for or simd autovectorization.
*/
template<int dim>
class BoxIteratorND
{
public:

    /**
    * Default constructor. This constructs an invalid iterator.
    * The user must call define before using.
    */
    BoxIteratorND () noexcept = default;

    /**
    * Constructs a BoxIteratorND and associates it with a Box.
    * Arguments:
    * a_bx (not modified) the Box to iterate over.
    */
    explicit BoxIteratorND (const BoxND<dim>& a_bx) noexcept {
        define(a_bx);
    }

    void setBox (const BoxND<dim>& a_bx) noexcept {
        define(a_bx);
    }

    /**
    * Associates a Box with this BoxIteratorND.
    * Arguments:
    * a_bx (not modified) the Box to iterate over.
    */
    void define (const BoxND<dim>& a_bx) noexcept {
        if (a_bx.ok()) {
            m_current = a_bx.smallEnd();
            m_boxLo   = a_bx.smallEnd();
            m_boxHi   = a_bx.bigEnd();
        } else {
            m_current = IntVectND<dim>::TheUnitVector();
            m_boxLo   = IntVectND<dim>::TheUnitVector();
            m_boxHi   = IntVectND<dim>::TheZeroVector();
        }
    }

    /**
    * Sets this BoxIteratorND to the first IntVect in its Box.
    * This is the smallEnd vector of the Box.
    */
    BoxIteratorND begin () noexcept {
        m_current = m_boxLo;
        return *this;
    }

    /**
    * Get a BoxIteratorND with its vector set to the end of the Box.
    * Should only be used internally by a range-based for loop.
    */
    [[nodiscard]] BoxIteratorND end () const noexcept {
        BoxIteratorND ret = *this;
        // set current of end to the past-the-end IntVect
        ret.m_current = m_boxLo;
        if (m_boxLo.allLE(m_boxHi)) {
            ret.m_current[dim-1] = m_boxHi[dim-1] + 1;
        }
        return ret;
    }

    /**
    * Sets this BoxIteratorND to the first IntVect in its Box.
    * This is the smallEnd vector of the Box.
    */
    void reset () noexcept {
        begin();
    }

    /**
    * Modifies this BoxIteratorND to set it to the next location in its Box.
    */
    BoxIteratorND& operator++ () noexcept {
        next();
        return *this;
    }

    /**
    * Modifies this BoxIteratorND to set it to the next location in its Box.
    */
    void next () noexcept {
        for (int i=0; i<dim; ++i) {
            ++m_current[i];
            if ((i+1<dim) && m_current[i] > m_boxHi[i]) {
                m_current[i] = m_boxLo[i];
                continue;
            }
            break;
        }
    }

    /**
    * Returns the value of the InVect for the current location of this BoxIteratorND.
    */
    [[nodiscard]] const IntVectND<dim>& operator () () const noexcept {
        AMREX_ASSERT(m_current.allLE(m_boxHi));
        AMREX_ASSERT(m_current.allGE(m_boxLo));
        return m_current;
    }

    /**
    * Returns the value of the InVect for the current location of this BoxIteratorND.
    */
    [[nodiscard]] IntVectND<dim> operator* () const noexcept {
        AMREX_ASSERT(m_current.allLE(m_boxHi));
        AMREX_ASSERT(m_current.allGE(m_boxLo));
        // Returning by value gives cleaner assembly output.
        return m_current;
    }

    /**
    * Returns true if this BoxIteratorND's location is within its Box.
    */
    [[nodiscard]] bool ok () const noexcept {
        return m_current.allLE(m_boxHi);
    }

    /**
    * Special operator for range-based for loops.
    * Only checks inequality of the last dimension of the current position for better performance.
    * Should only be used internally by a range-based for loop.
    */
    [[nodiscard]] friend
    bool operator != (const BoxIteratorND& b1, const BoxIteratorND& b2) noexcept {
        // Only check one value for better performance.
        return b1.m_current[dim-1] != b2.m_current[dim-1];
    }

private:
    IntVectND<dim> m_current = IntVectND<dim>::TheUnitVector();
    IntVectND<dim> m_boxLo   = IntVectND<dim>::TheUnitVector();
    IntVectND<dim> m_boxHi   = IntVectND<dim>::TheZeroVector();
};

using BoxIterator = BoxIteratorND<AMREX_SPACEDIM>;

}
#endif
