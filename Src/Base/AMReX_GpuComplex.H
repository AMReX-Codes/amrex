#ifndef AMREX_GPUCOMPLEX_H_
#define AMREX_GPUCOMPLEX_H_

#include <cmath>
#include <iostream>
#include <AMReX_Math.H>

namespace amrex {

template <typename T> struct GpuComplex;

template <typename T> 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE    
T norm (const GpuComplex<T>& a_z) noexcept;

/**
 * \brief A host / device complex number type, because std::complex doesn't
 *  work in device code with Cuda yet.
 *
 *  Should be bit-wise compatible with std::complex.
 */            
template <typename T>
struct GpuComplex 
{
    using value_type = T;

    T m_real, m_imag;

    /**
     * \brief Construct a complex number given the real and imaginary part.
     */    
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE    
    constexpr GpuComplex (const T& a_r = T(), const T& a_i = T()) noexcept
        : m_real(a_r), m_imag(a_i) {}

    /**
     * \brief Return the real part.
     */        
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr T real () const noexcept { return m_real; }

    /**
     * \brief Return the imaginary part.
     */            
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr T imag () const noexcept { return m_imag; }

   /**
     * \brief Add a real number to this complex number.
     */                
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuComplex<T>& operator+= (const T& a_t) noexcept
    {
        m_real += a_t;
        return *this;
    }

   /**
     * \brief Subtract a real number from this complex number.
     */                    
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuComplex<T>& operator-= (const T& a_t) noexcept
    {
        m_real -= a_t;
        return *this;
    }

    /**
     * \brief Multiply this complex number by a real.
     */                
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuComplex<T>& operator*= (const T& a_t) noexcept
    {
        m_real *= a_t;
        m_imag *= a_t;
        return *this;
    }

    /**
     * \brief Divide this complex number by a real.
     */                    
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuComplex<T>& operator/= (const T& a_t) noexcept
    {
        m_real /= a_t;
        m_imag /= a_t;
        return *this;
    }

    /**
     * \brief Add another complex number to this one.
     */                    
    template <typename U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuComplex<T>& operator+= (const GpuComplex<U>& a_z) noexcept
    {
        m_real += a_z.real();
        m_imag += a_z.imag();
        return *this;
    }

    /**
     * \brief Subtract another complex number from this one.
     */                        
    template <typename U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuComplex<T>& operator-= (const GpuComplex<U>& a_z) noexcept
    {
        m_real -= a_z.real();
        m_imag -= a_z.imag();
        return *this;
    }

    /**
     * \brief Multiply this complex number by another one.
     */                            
    template <typename U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuComplex<T>& operator*= (const GpuComplex<U>& a_z) noexcept
    {
        const T a_r = m_real * a_z.real() - m_imag * a_z.imag();
        m_imag = m_real * a_z.imag() + m_imag * a_z.real();
        m_real = a_r;
        return *this;
    }

    /**
     * \brief Divide this complex number by another one.
     */                                
    template <typename U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuComplex<T>& operator/= (const GpuComplex<U>& a_z) noexcept
    {
        const T r =  m_real * a_z.real() + m_imag * a_z.imag();
        const T n = amrex::norm(a_z);
        m_imag = (m_imag * a_z.real() - m_real * a_z.imag()) / n;
        m_real = r / n;
        return *this;
    }

    /**
     * \brief Print this complex number to an ostream.
     */                                    
    template <typename U>
    friend std::ostream & operator<< (std::ostream &out, const GpuComplex<U>& c);
};

template <typename U>
std::ostream& operator<< (std::ostream& out, const GpuComplex<U>& c)
{
    out << c.m_real;
    out << " + " << c.m_imag << "i";
    return out;
} 

/**
 * \brief Identity operation on a complex number.
 */                                        
template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator+ (const GpuComplex<T>& a_x) { return a_x; }

/**
 * \brief Negate a complex number.
 */    
template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator- (const GpuComplex<T>& a_x) { return GpuComplex<T>(-a_x.real(), -a_x.imag()); }

/**
 * \brief Subtract two complex numbers.
 */        
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator- (const GpuComplex<T>& a_x, const GpuComplex<T>& a_y) noexcept
{
    GpuComplex<T> r = a_x;
    r -= a_y;
    return r;
}

/**
 * \brief Subtract a real number from a complex one.
 */            
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator- (const GpuComplex<T>& a_x, const T& a_y) noexcept
{
    GpuComplex<T> r = a_x;
    r -= a_y;
    return r;
}

/**
 * \brief Subtract a complex number from a real one.
 */            
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator- (const T& a_x, const GpuComplex<T>& a_y) noexcept
{
    GpuComplex<T> r = a_y;
    r -= a_x;
    return -r;
}

/**
 * \brief Add two complex numbers.
 */            
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator+ (const GpuComplex<T>& a_x, const GpuComplex<T>& a_y) noexcept
{
    GpuComplex<T> r = a_x;
    r += a_y;
    return r;
}

/**
 * \brief Add a real number to a complex one
 */            
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator+ (const GpuComplex<T>& a_x, const T& a_y) noexcept
{
    GpuComplex<T> r = a_x;
    r += a_y;
    return r;
}

/**
 * \brief Add a complex number to a real one
 */            
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator+ (const T& a_x, const GpuComplex<T>& a_y) noexcept
{
    GpuComplex<T> r = a_y;
    r += a_x;
    return r;
}

/**
 * \brief Multiply two complex numbers.
 */            
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator* (const GpuComplex<T>& a_x, const GpuComplex<T>& a_y) noexcept
{
    GpuComplex<T> r = a_x;
    r *= a_y;
    return r;
}

/**
 * \brief Multiply a complex number by a real one.
 */            
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator* (const GpuComplex<T>& a_x, const T& a_y) noexcept
{
    GpuComplex<T> r = a_x;
    r *= a_y;
    return r;
}

/**
 * \brief Multiply a real number by a complex one.
 */                
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator* (const T& a_x, const GpuComplex<T>& a_y) noexcept
{
    GpuComplex<T> r = a_y;
    r *= a_x;
    return r;
}

/**
 * \brief Divide a complex number by another one.
 */                
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator/ (const GpuComplex<T>& a_x, const GpuComplex<T>& a_y) noexcept
{
    GpuComplex<T> r = a_x;
    r /= a_y;
    return r;
}

/**
 * \brief Divide a complex number by a real.
 */                    
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator/ (const GpuComplex<T>& a_x, const T& a_y) noexcept
{
    GpuComplex<T> r = a_x;
    r /= a_y;
    return r;
}

/**
 * \brief Divide a real number by a complex one.
 */                    
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> operator/ (const T& a_x, const GpuComplex<T>& a_y) noexcept
{
    GpuComplex<T> r = a_x;
    r /= a_y;
    return r;
}

/**
 * \brief Return a complex number given its polar representation.
 */                        
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> polar (const T& a_r, const T& a_theta) noexcept
{
    return GpuComplex<T>(a_r * std::cos(a_theta), a_r * std::sin(a_theta));
}    

/**
 * \brief Complex expotential function.
 */                            
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> exp (const GpuComplex<T>& a_z) noexcept
{
    return polar<T>(std::exp(a_z.real()), a_z.imag());
}    

/**
 * \brief Return the norm (magnitude squared) of a complex number
 */                                
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T norm (const GpuComplex<T>& a_z) noexcept
{
    const T x = a_z.real();
    const T y = a_z.imag();
    return x * x + y * y;
}    

/**
 * \brief Return the absolute value of a complex number
 */                                    
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T abs (const GpuComplex<T>& a_z) noexcept
{
    T x = a_z.real();
    T y = a_z.imag();

    const T s = amrex::max(amrex::Math::abs(x), amrex::Math::abs(y));
    if (s == T())
        return s;
    x /= s;
    y /= s;
    return s * std::sqrt(x * x + y * y);
}

/**
 * \brief Return the square root of a complex number
 */                                        
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> sqrt (const GpuComplex<T>& a_z) noexcept
{
    T x = a_z.real();
    T y = a_z.imag();

    if (x == T())
    {
        T t = std::sqrt(amrex::Math::abs(y) / 2);
        return GpuComplex<T>(t, y < T() ? -t : t);
    }
    else
    {
        T t = std::sqrt(2 * (amrex::abs(a_z) + amrex::Math::abs(x)));
        T u = t / 2;
        return x > T()
            ? GpuComplex<T>(u, y / t)
            : GpuComplex<T>(amrex::Math::abs(y) / t, y < T() ? -u : u);
    }
}

/**
 * \brief Return the angle of a complex number's polar representation.
 */                                        
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T arg (const GpuComplex<T>& a_z) noexcept
{
    return std::atan2(a_z.imag(), a_z.real());
}

/**
 * \brief Complex natural logarithm function.
 */                                            
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> log (const GpuComplex<T>& a_z) noexcept
{
    return GpuComplex<T>(std::log(amrex::abs(a_z)), amrex::arg(a_z));
}

/**
 * \brief Raise a complex number to a (real) power.
 */                                                
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> pow (const GpuComplex<T>& a_z, const T& a_y) noexcept
{
    if (a_z.imag() == T() && a_z.real() > T())
        return std::pow(a_z.real(), a_y);
    
    GpuComplex<T> t = amrex::log(a_z);
    return amrex::polar<T>(std::exp(a_y * t.real()), a_y * t.imag());
}

namespace detail
{
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuComplex<T> complex_pow_unsigned (GpuComplex<T> a_z, unsigned a_n)
    {
        GpuComplex<T> y = a_n % 2 ? a_z : GpuComplex<T>(1);
        
        while (a_n >>= 1)
        {
            a_z *= a_z;
            if (a_n % 2)
                y *= a_z;
        }
        
        return y;
    }
}

/**
 * \brief Raise a complex number to an integer power.
 */                                                    
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<T> pow (const GpuComplex<T>& a_z, int a_n) noexcept
{
    return a_n < 0
                 ? GpuComplex<T>(1) / detail::complex_pow_unsigned(a_z, (unsigned)(-a_n))
                 : detail::complex_pow_unsigned(a_z, a_n);
}
    
}

#endif
