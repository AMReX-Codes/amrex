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

namespace amrex {
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    T getEnum (std::string_view const& s)
    {
        auto pos = ET::enum_names.find(s);
        if (pos == std::string_view::npos) {
            std::string error_msg("amrex::getEnum: Unknown enum: ");
            error_msg.append(s).append(" in AMREX_ENUM(").append(ET::class_name)
                .append(", ").append(ET::enum_names).append(").");
            throw std::runtime_error(error_msg);
        }
        auto count = std::count(ET::enum_names.begin(),
                                ET::enum_names.begin()+pos, ',');
        return static_cast<T>(count);
    }

    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::string getEnumNameString (T const& v)
    {
        auto n = static_cast<int>(v);
        std::size_t pos = 0;
        for (int i = 0; i < n; ++i) {
            pos = ET::enum_names.find(',', pos);
            if (pos == std::string::npos) {
                std::string error_msg("amrex::getEnum: Unknown enum value: ");
                error_msg.append(std::to_string(n)).append(" in AMREX_ENUM(")
                    .append(ET::class_name).append(", ").append(ET::enum_names)
                    .append(").");
                throw std::runtime_error(error_msg);
            }
            ++pos;
        }
        auto pos2 = ET::enum_names.find(',', pos);
        return amrex::trim(std::string(ET::enum_names.substr(pos,pos2-pos)));
    }

    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    std::vector<std::string> getEnumNameStrings ()
    {
        return amrex::split(std::string(ET::enum_names), ", ");
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
