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
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::vector<std::pair<std::string,T>> getEnumNameValuePairs ()
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
                auto it = std::find_if(r.begin(), r.end(),
                                       [&] (auto const& x) -> bool
                                           { return x.first == value_string; });
                auto this_value = it->second;
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
