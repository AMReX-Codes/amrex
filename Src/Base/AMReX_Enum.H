#ifndef AMREX_ENUM_H_
#define AMREX_ENUM_H_

/**
 * \file AMReX_Enum.H
 * \brief Enum reflection utilities and the AMREX_ENUM macro.
 */

#include <AMReX_String.H>

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

template <typename T>
using amrex_enum_traits = decltype(amrex_get_enum_traits(std::declval<T>()));

namespace amrex
{
/// \cond DOXYGEN_IGNORE
namespace detail
{
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::vector<std::pair<std::string,T>> getNewEnumNameValuePairs ()
    {
        auto tmp = amrex::split(std::string(ET::enum_names), ",");
        std::vector<std::pair<std::string,T>> r;
        r.reserve(tmp.size());
        int next_value = 0;
        for (auto const& item : tmp) {
            auto const& kv = amrex::split(item, "=");
            if (kv.size() == 1) {
                r.emplace_back(amrex::trim(kv[0]), static_cast<T>(next_value));
                ++next_value;
            } else if (kv.size() == 2) {
                auto const& value_string = amrex::trim(kv[1]);
                // For k = v, find if v exists as a key. If not, v must be an int.
                auto it = std::find_if(r.begin(), r.end(),
                                       [&] (auto const& x) -> bool
                                           { return x.first == value_string; });
                T this_value;
                if (it != r.end()) {
                    this_value = it->second;
                } else { // v is an int
                    // Note that if value_string is not a valid string for int,
                    // it will fail at compile time, because the compiler
                    // will not allow something like `enum FOO : int { x = 3% }`.
                    this_value =  static_cast<T>(std::stoi(value_string));
                }
                r.emplace_back(amrex::trim(kv[0]), this_value);
                next_value = static_cast<int>(this_value) + 1;
            } else {
                std::string error_msg("amrex::getEnumNameIntPairs: AMREX_ENUM(");
                error_msg.append(ET::class_name).append(", ").append(ET::enum_names)
                    .append(").");
                throw std::runtime_error(error_msg);
            }
        }
        return r;
    }
}  // namespace detail
/// \endcond

    /**
     * \ingroup amrex_utilities
     * \brief Return all (name, value) pairs for an AMREX_ENUM-declared enum type.
     *
     * The result is cached the first time it is called for a given type \p T.
     *
     * \tparam T An enum type declared with AMREX_ENUM.
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::vector<std::pair<std::string,T>> const& getEnumNameValuePairs ()
    {
        // cache enum value strings for the whole type T
        static auto r = detail::getNewEnumNameValuePairs<T>();
        return r;
    }

    /**
     * \ingroup amrex_utilities
     * \brief Convert a string to an enum value.
     *
     * Throws `std::runtime_error` if \p s does not match any enumerator name.
     *
     * Example:
     * ```C++
     * AMREX_ENUM(Model,
     *     linear,
     *     nonlinear
     * );
     * std::string const model_str = "nonlinear";
     * Model const model = amrex::getEnum<Model>(model_str);
     * assert(model == Model::nonlinear);
     * ```
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    T getEnum (std::string_view const& s)
    {
        auto const& kv = getEnumNameValuePairs<T>();
        auto it = std::find_if(kv.begin(), kv.end(),
                               [&] (auto const& x) -> bool
                                   { return x.first == s; });
        if (it != kv.end()) {
            return it->second;
        } else {
            std::string error_msg("amrex::getEnum: Unknown enum: ");
            error_msg.append(s).append(" in AMREX_ENUM(").append(ET::class_name)
                .append(", ").append(ET::enum_names).append(").");
            throw std::runtime_error(error_msg);
            return T();
        }
    }

    /**
     * \ingroup amrex_utilities
     * \brief Convert a string case-insensitively to an enum value.
     *
     * Same as getEnum<T>, but uses a case-insensitive match.
     *
     * Example:
     * ```C++
     * AMREX_ENUM(Model,
     *     linear,
     *     nonlinear
     * );
     * std::string const model_str = "NonLinear";
     * Model const model = amrex::getEnumCaseInsensitive<Model>(model_str);
     * assert(model == Model::nonlinear);
     * ```
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    T getEnumCaseInsensitive (std::string_view const& s)
    {
        auto const& kv = getEnumNameValuePairs<T>();
        std::string ls = amrex::toLower(std::string(s));
        auto it = std::find_if(kv.begin(), kv.end(),
                               [&] (auto const& x) -> bool
                                   { return amrex::toLower(x.first) == ls; });
        if (it != kv.end()) {
            return it->second;
        } else {
            std::string error_msg("amrex::getEnumCaseInsensitive: Unknown enum: ");
            error_msg.append(s).append(" in AMREX_ENUM(").append(ET::class_name)
                .append(", ").append(ET::enum_names).append(").");
            throw std::runtime_error(error_msg);
            return T();
        }
    }

    /**
     * \ingroup amrex_utilities
     * \brief Get the name string of an enum value.
     *
     * Throws `std::runtime_error` if \p v has no matching enumerator name.
     *
     * Example:
     * ```C++
     * AMREX_ENUM(Model,
     *     linear,
     *     nonlinear
     * );
     * Model model = Model::linear;
     * std::string model_str = amrex::getEnumNameString(model);
     * assert(model_str == "linear");
     * ```
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::string getEnumNameString (T const& v)
    {
        auto const& kv = getEnumNameValuePairs<T>();
        auto it = std::find_if(kv.begin(), kv.end(),
                               [&] (auto const& x) -> bool
                                   { return x.second == v; });
        if (it != kv.end()) {
            return it->first;
        } else {
            std::string error_msg("amrex::getEnum: Unknown enum value: ");
            error_msg.append(std::to_string(static_cast<int>(v)))
                .append(" in AMREX_ENUM(").append(ET::class_name).append(", ")
                .append(ET::enum_names).append(").");
            throw std::runtime_error(error_msg);
            return std::string();
        }
    }

    /**
     * \ingroup amrex_utilities
     * \brief Return a list of all enumerator name strings for an AMREX_ENUM-declared type.
     *
     * \tparam T An enum type declared with AMREX_ENUM.
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::vector<std::string> getEnumNameStrings ()
    {
        auto r = amrex::split(std::string(ET::enum_names), ",");
        for (auto& s : r) {
            auto found = s.find('=');
            if (found != std::string::npos) {
                s.erase(found);
            }
            s = amrex::trim(s);
        }
        return r;
    }

    /**
     * \ingroup amrex_utilities
     * \brief Return the class name string of an AMREX_ENUM-declared enum type.
     *
     * \tparam T An enum type declared with AMREX_ENUM.
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::string getEnumClassName ()
    {
        return std::string(ET::class_name);
    }

    /**
     * \ingroup amrex_utilities
     * \brief Return the underlying integer value of an enum enumerator.
     *
     * Useful when building bitmasks.
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    constexpr auto toUnderlying (T v) noexcept
    {
        return static_cast<std::underlying_type_t<T>>(v);
    }
}

/**
 * \ingroup amrex_utilities
 * \brief Declare a scoped enum with reflection support.
 *
 * Defines `enum class CLASS : int` and registers its enumerator names so that
 * amrex::getEnum, amrex::getEnumCaseInsensitive, amrex::getEnumNameString,
 * and related utilities can convert between strings and values at run time.
 *
 * Example:
 * ```cpp
 * AMREX_ENUM(Color, red, green, blue);
 * Color c = amrex::getEnum<Color>("green");
 * ```
 */
#define AMREX_ENUM(CLASS, ...) \
    enum class CLASS : int { __VA_ARGS__ }; \
    /** \cond DOXYGEN_IGNORE */ \
    struct CLASS##_EnumTraits { \
        using enum_class_t = CLASS; \
        static constexpr bool value = true; \
        static constexpr std::string_view class_name{#CLASS}; \
        static constexpr std::string_view enum_names{#__VA_ARGS__}; \
    }; \
    /** \endcond */ \
    CLASS##_EnumTraits amrex_get_enum_traits(CLASS)

#endif
