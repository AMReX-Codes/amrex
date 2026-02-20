#ifndef AMREX_STACK_H_
#define AMREX_STACK_H_

/**
 * \file AMReX_Stack.H
 * \brief A compile-time-sized stack container.
 */

namespace amrex {

/**
 * \class amrex::Stack
 * \ingroup amrex_linear_containers
 * \brief A fixed-capacity stack with compile-time maximum size.
 *
 * \tparam T Value type.
 * \tparam N Maximum number of elements.
 */
template <typename T, int N>
struct Stack
{
public:
    //! Push \p v onto the top of the stack.
    constexpr void push (T v) { m_data[m_size++] = v; }
    //! Remove the top element from the stack.
    constexpr void pop () { --m_size; }
    //! Return true if the stack contains no elements.
    [[nodiscard]] constexpr bool empty () const { return m_size == 0; }
    //! Return the number of elements currently in the stack.
    [[nodiscard]] constexpr int size () const { return m_size; }
    //! Return a const reference to the top element.
    [[nodiscard]] constexpr T const& top () const { return m_data[m_size-1]; }
    //! Return a mutable reference to the top element.
    [[nodiscard]] constexpr T      & top ()       { return m_data[m_size-1]; }
    //! Return the element at index \p i (0-based, bottom of stack first).
    [[nodiscard]] constexpr T operator[] (int i) const { return m_data[i]; }
private:
    T m_data[N];
    int m_size = 0;
};

}

#endif
