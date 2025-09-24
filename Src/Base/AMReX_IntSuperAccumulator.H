#ifndef AMREX_INT_SUPERACCUMULATOR_H_
#define AMREX_INT_SUPERACCUMULATOR_H_

#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Math.H>

#include <cmath>
#include <cstdint>

namespace amrex {

/**
 * \brief Integer superaccumulator backed by IArrayBox storage.
 *
 * The accumulator stores the exact sum of single-precision floating-point
 * numbers by representing the scaled integer \f$\sum_i x_i 2^{-e_\text{min}}\f$
 * in a radix-\f$2^{30}\f$ positional system. Each digit is stored as a 32-bit int in
 * an `IArrayBox` component. A client that needs M independent sums must
 * allocate an accumulator with `M * digits_per_entry` components and index the
 * storage using `component = entry * digits_per_entry + digit`. Both
 * `accumulate` and `finalize` operate in single precision; callers may cast to
 * wider types if needed. All summands are assumed non-negative; subnormal
 * inputs are treated as zero to reduce the number of required digits per entry.
 */
class IntSuperAccumulatorFab
{
public:
    static constexpr int digit_bits = 30;                 //!< Width of each digit in bits.
    static constexpr int digit_base = 1 << digit_bits;    //!< Base of each digit.
    static constexpr int min_exponent = -149;             //!< Smallest exponent of scaled mantissa (binary32 minimum normal).
    static constexpr int max_exponent =  104;             //!< Largest exponent of scaled mantissa (binary32 maximum normal).
    static constexpr int significand_bits = 24;           //!< Binary32 significand bits including hidden bit.
    static constexpr int extra_digits = 1;                //!< Safety margin for carry propagation.
    static constexpr int digits_for_shift = ((max_exponent - min_exponent) / digit_bits) + 1; // covers exponent window
    static constexpr int digits_for_value = 2;            //!< Worst-case mantissa occupies two digits at radix 2^30.
    static constexpr int digits_per_entry = digits_for_shift + digits_for_value + extra_digits; //!< 12 digits per entry.

    IntSuperAccumulatorFab () = default;

    explicit IntSuperAccumulatorFab (Arena* arena) noexcept
        : m_storage(arena) {}

    IntSuperAccumulatorFab (const Box& box, int nentries, Arena* arena=nullptr)
        : m_storage(box, nentries * digits_per_entry, arena), m_entries(nentries) {}

    [[nodiscard]] int entries () const noexcept { return m_entries; }

    void define (const Box& box, int nentries, Arena* arena=nullptr)
    {
        m_entries = nentries;
        m_storage.resize(box, nentries * digits_per_entry, arena);
    }

    [[nodiscard]] Array4<int> array () noexcept { return m_storage.array(); }
    [[nodiscard]] Array4<int const> const_array () const noexcept { return m_storage.const_array(); }

    void reset ()
    {
        if (amrex::Gpu::inLaunchRegion()) {
            m_storage.template setVal<RunOn::Device>(0);
        } else {
            m_storage.template setVal<RunOn::Host>(0);
        }
    }

    void reset (Box const& bx, int nentries = -1)
    {
        const int entries = (nentries < 0) ? m_entries : nentries;
        const int ncomp = entries * digits_per_entry;
        if (amrex::Gpu::inLaunchRegion()) {
            m_storage.template setVal<RunOn::Device>(0, bx, 0, ncomp);
        } else {
            m_storage.template setVal<RunOn::Host>(0, bx, 0, ncomp);
        }
    }

    /**
     * \brief Accumulate a floating-point value into the superaccumulator.
     *
     * The caller must provide the digit array and logical comp index for the
     * entry being updated. All operations are carried out using integers; the
     * routine does not allocate additional storage.
     */
    template <typename Array>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void
    accumulate (Array const& digits, int i, int j, int k, int entry, float value,
                int entry_digits = digits_per_entry) noexcept
    {
        if (value <= 0.0f) {
            if (value < 0.0f) {
#ifndef AMREX_USE_GPU
                amrex::Abort("IntSuperAccumulatorFab::accumulate: negative input not supported");
#endif
            }
            return;
        }

        if (!amrex::Math::isfinite(value)) {
#ifndef AMREX_USE_GPU
            amrex::Abort("IntSuperAccumulatorFab::accumulate: non-finite input");
#else
            return;
#endif
        }

        const auto bits = detail::decompose(value);
        if (bits.is_zero) { return; }

        const int component_base = entry * entry_digits;

        const int shift = bits.exponent - min_exponent;
        int digit_index = shift / digit_bits;
        const int bit_offset = shift % digit_bits;

        const unsigned long long magnitude = bits.magnitude;
        unsigned long long chunk_value = magnitude << bit_offset;

        const int component_limit = component_base + entry_digits;

        while (chunk_value != 0ULL) {
            const unsigned long long chunk = chunk_value & static_cast<unsigned long long>(digit_base - 1);
            chunk_value >>= digit_bits;
            if (component_base + digit_index >= component_limit) {
#ifndef AMREX_USE_GPU
                amrex::Abort("IntSuperAccumulatorFab: digit overflow");
#else
                break;
#endif
            }
            detail::add_to_digit(digits, i, j, k,
                                 component_base + digit_index,
                                 component_limit,
                                 chunk);
            ++digit_index;
        }
    }

    /**
     * \brief Finalize the accumulator entry and return a binary32 result.
     */
    template <typename Array>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static float
    finalize (Array const& digits, int i, int j, int k, int entry,
              int entry_digits = digits_per_entry) noexcept
    {
        const int component_base = entry * entry_digits;

        double result = 0.0;
        for (int d = entry_digits - 1; d >= 0; --d) {
            const int digit = digits(i,j,k, component_base + d);
            if (digit == 0) { continue; }
            const int exponent = min_exponent + d * digit_bits;
            result += std::ldexp(static_cast<double>(digit), exponent);
        }
        return static_cast<float>(result);
    }

    template <typename Array>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void
    clear_entry (Array const& digits, int i, int j, int k, int entry,
                 int entry_digits = digits_per_entry) noexcept
    {
        const int component_base = entry * entry_digits;
        for (int d = 0; d < entry_digits; ++d) {
            digits(i,j,k, component_base + d) = 0;
        }
    }

private:
    struct Decomposed
    {
        unsigned long long magnitude;
        int exponent;
        bool is_zero;
    };

    struct detail
    {
        template <typename Array>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void
        add_to_digit (Array const& digits, int i, int j, int k,
                      int component, int component_limit, unsigned long long value) noexcept
        {
            unsigned long long carry = value;
            int current_component = component;
            constexpr unsigned long long base = static_cast<unsigned long long>(digit_base);

            while (carry != 0ULL) {
                if (current_component >= component_limit) {
#ifndef AMREX_USE_GPU
                    amrex::Abort("IntSuperAccumulatorFab: carry overflow");
#else
                    break;
#endif
                }

#if defined(AMREX_USE_GPU) && AMREX_DEVICE_COMPILE
                // When we compile for the device we cannot rely on a carry-save digit layout, so
                // each thread performs an immediate normalization using GPU atomics.
                int* const digit_ptr = &(digits(i,j,k,current_component));

                // Break the carry into the current radix chunk and the carry to forward.
                const unsigned long long chunk = carry & (base - 1ULL);
                carry >>= digit_bits;

                // Atomically accumulate this chunk.  Multiple threads may land on the same digit,
                // which is safe because atomicAdd enforces a total order.
                if (chunk != 0ULL) {
                    amrex::Gpu::Atomic::Add(digit_ptr, static_cast<int>(chunk));
                }

                // Track how many bases we successfully removed so they can be forwarded upward.
                unsigned long long propagated = 0ULL;
                while (true) {
                    // Snapshot the digit atomically without changing it.  This avoids racing
                    // against other threads that may be normalizing the same slot.
                    const int current_value = amrex::Gpu::Atomic::Add(digit_ptr, 0);
                    if (current_value < digit_base) {
                        break;
                    }

                    // Attempt to drain one full base from the digit.  The return value is the
                    // pre-subtraction contents, so we know whether we actually claimed the carry.
                    const int before = amrex::Gpu::Atomic::Add(digit_ptr, -digit_base);
                    if (before >= digit_base) {
                        // We consumed one unit of overflow and should forward it to the next digit.
                        ++propagated;
                        continue;
                    } else {
                        // Another thread already normalized the digit between the peek and our
                        // subtract; restore the value and exit.
                        amrex::Gpu::Atomic::Add(digit_ptr, digit_base);
                        break;
                    }
                }

                // Fold any carries we successfully removed into the higher-order digits.
                carry += propagated;
#else
                const unsigned long long sum = static_cast<unsigned long long>(digits(i,j,k,current_component)) + carry;
                const unsigned long long remainder = sum & (base - 1ULL);
                const unsigned long long quotient = sum >> digit_bits;

                digits(i,j,k,current_component) = static_cast<int>(remainder);
                carry = quotient;
#endif
                ++current_component;
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static Decomposed
        decompose (float value) noexcept
        {
            Decomposed bits{};

            if (value <= 0.0f) {
                bits.is_zero = true;
                return bits;
            }

            union {
                float f;
                std::uint32_t u;
            } pun{};
            pun.f = value;
            const std::uint32_t u = pun.u;

            const std::uint32_t exponent_field = (u >> 23) & 0xffu;
            const std::uint32_t fraction = u & ((1u << 23) - 1u);

            if (exponent_field == 0xffu) {
#ifndef AMREX_USE_GPU
                amrex::Abort("IntSuperAccumulatorFab::decompose: NaN or Inf not supported (float)");
#else
                Decomposed tmp{}; tmp.is_zero = true; return tmp;
#endif
            } else if (exponent_field == 0u) {
                // Treat subnormals as zero; caller is expected to pre-round.
                bits.is_zero = true;
                return bits;
            } else {
                bits.magnitude = (1u << 23) | fraction;
                bits.exponent = static_cast<int>(exponent_field) - 127 - significand_bits + 1;
            }

            bits.is_zero = false;
            return bits;
        }
    };

    IArrayBox m_storage;
    int m_entries = 0;
};

} // namespace amrex

#endif
