#ifndef WARPX_amrex_qed_table_parser_helper_functions_h_
#define WARPX_amrex_qed_table_parser_helper_functions_h_

//This file contains helper functions to parse a char* array
//into a lookup table

#include <vector>
#include <tuple>

namespace QedUtils{
    template <class T>
    std::tuple<bool, std::vector<T>, const char*>parse_raw_data_vec(
        const char* p_data, size_t how_many, const char* const p_last)
    {
        std::vector<T> res;
        if(p_data + sizeof(T)*how_many > p_last)
            return std::make_tuple(false, res, nullptr);

        auto r_data = reinterpret_cast<const T*>(p_data);

        res.assign(r_data, r_data + how_many);

        p_data += sizeof(T)*how_many;
        return std::make_tuple(true, res, p_data);
    }

    template <class T>
    std::tuple<bool, T, const char*> parse_raw_data(
        const char* p_data, const char* const p_last)
    {
        T res;
        if(p_data + sizeof(T) > p_last)
            return std::make_tuple(false, res, nullptr);

        auto r_data = reinterpret_cast<const T*>(p_data);

        res = *r_data;

        p_data += sizeof(T);
        return std::make_tuple(true, res, p_data);
    }
};

#endif //WARPX_amrex_qed_table_parser_helper_functions_h_
