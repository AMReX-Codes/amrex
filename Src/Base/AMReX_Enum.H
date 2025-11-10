#ifndef AMREX_ENUM_H_
#define AMREX_ENUM_H_

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

    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::vector<std::pair<std::string,T>> const& getEnumNameValuePairs ()
    {
        // cache enum value strings for the whole type T
        static auto r = detail::getNewEnumNameValuePairs<T>();
        return r;
    }

    /** Convert a string to an enum value
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

    /** Convert a string case insensitive to an enum value
     *
     * Same as getEnum<T>, but case insensitive match to enum value.
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

    /** Get a string from an enum value
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

    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::string getEnumClassName ()
    {
        return std::string(ET::class_name);
    }

    /** Return the underlying (u)int of an enum value
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

#define AMREX_ENUM(CLASS, ...) \
    enum class CLASS : int { __VA_ARGS__ }; \
    struct CLASS##_EnumTraits { \
        using enum_class_t = CLASS; \
        static constexpr bool value = true; \
        static constexpr std::string_view class_name{#CLASS}; \
        static constexpr std::string_view enum_names{#__VA_ARGS__}; \
    }; \
    CLASS##_EnumTraits amrex_get_enum_traits(CLASS)

#endif
