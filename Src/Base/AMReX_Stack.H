#ifndef AMREX_STACK_H_
#define AMREX_STACK_H_

namespace amrex {

template <typename T, int N>
struct Stack
{
public:
    constexpr void push (T v) { m_data[m_size++] = v; }
    constexpr void pop () { --m_size; }
    [[nodiscard]] constexpr bool empty () const { return m_size == 0; }
    [[nodiscard]] constexpr int size () const { return m_size; }
    [[nodiscard]] constexpr T const& top () const { return m_data[m_size-1]; }
    [[nodiscard]] constexpr T      & top ()       { return m_data[m_size-1]; }
    [[nodiscard]] constexpr T operator[] (int i) const { return m_data[i]; }
private:
    T m_data[N];
    int m_size = 0;
};

}

#endif
