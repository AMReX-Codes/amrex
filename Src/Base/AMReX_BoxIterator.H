#ifndef AMREX_BOXITERATOR_H_
#define AMREX_BOXITERATOR_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_BoxIterator.H
 * \brief Iterator over all IntVects in a Box.
 */

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
 * \class amrex::BoxIteratorND
 * \ingroup amrex_index_space
 * \brief Iterator over all IntVects in a Box.
 *
 * Iterates through the IntVects of a Box in the same order as a series of
 * nested for loops, with the innermost loop over the first dimension.
 *
 * Typical usage with a range-based for loop:
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
 *
 * \tparam dim Spatial dimension.
 */
template<int dim>
class BoxIteratorND
{
public:

    //! Default constructor. Constructs an invalid iterator; call define() before use.
    BoxIteratorND () noexcept = default;

    /**
     * \brief Construct a BoxIteratorND associated with Box \p a_bx.
     * \param a_bx Box to iterate over.
     */
    explicit BoxIteratorND (const BoxND<dim>& a_bx) noexcept {
        define(a_bx);
    }

    //! Associate this iterator with Box \p a_bx (equivalent to define()).
    void setBox (const BoxND<dim>& a_bx) noexcept {
        define(a_bx);
    }

    /**
     * \brief Associate this iterator with Box \p a_bx.
     * \param a_bx Box to iterate over.
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

    //! Reset to the first IntVect (smallEnd) and return \c *this for range-based for.
    BoxIteratorND begin () noexcept {
        m_current = m_boxLo;
        return *this;
    }

    //! Return a past-the-end sentinel for use by range-based for loops.
    [[nodiscard]] BoxIteratorND end () const noexcept {
        BoxIteratorND ret = *this;
        // set current of end to the past-the-end IntVect
        ret.m_current = m_boxLo;
        if (m_boxLo.allLE(m_boxHi)) {
            ret.m_current[dim-1] = m_boxHi[dim-1] + 1;
        }
        return ret;
    }

    //! Reset to the first IntVect (smallEnd) of the Box.
    void reset () noexcept {
        begin();
    }

    //! Advance to the next IntVect in the Box (pre-increment).
    BoxIteratorND& operator++ () noexcept {
        next();
        return *this;
    }

    //! Advance to the next IntVect in the Box.
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

    //! Return the current IntVect.
    [[nodiscard]] const IntVectND<dim>& operator () () const noexcept {
        AMREX_ASSERT(m_current.allLE(m_boxHi));
        AMREX_ASSERT(m_current.allGE(m_boxLo));
        return m_current;
    }

    //! Dereference: return the current IntVect by value.
    [[nodiscard]] IntVectND<dim> operator* () const noexcept {
        AMREX_ASSERT(m_current.allLE(m_boxHi));
        AMREX_ASSERT(m_current.allGE(m_boxLo));
        // Returning by value gives cleaner assembly output.
        return m_current;
    }

    //! Return true if the iterator is within the Box.
    [[nodiscard]] bool ok () const noexcept {
        return m_current.allLE(m_boxHi);
    }

    //! Inequality test used by range-based for loops (checks only the last dimension).
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

/// \ingroup amrex_index_space
/// Alias for BoxIteratorND<AMREX_SPACEDIM>.
using BoxIterator = BoxIteratorND<AMREX_SPACEDIM>;

}
#endif
