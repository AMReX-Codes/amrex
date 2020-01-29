/* Copyright 2019 Luca Fedeli, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_amrex_qed_table_parser_helper_functions_h_
#define WARPX_amrex_qed_table_parser_helper_functions_h_

/**
 * This header contains helper functions to safely extract data
 * (e.g. integers, floating point numbers) from raw binary data
 * (i.e. a char*) and to convert arrays into raw binary data.
 */

#include <AMReX_Vector.H>
#include <tuple>

namespace QedUtils{
    /**
    * This function safely extracts an amrex::Vector<T> from raw binary data.
    * T must be a simple datatype (e.g. an int, a float, a double...).
    *
    * @param[in] p_data a pointer to the binary stream
    * @param[in] how_many how many T should be read from stream
    * @param[in] p_last a pointer to the last element of the char* array
    * @return {a tuple containing
    * 1) flag (which is false if p_last is exceeded)
    * 2) a Vector of T
    * 3) a pointer to a new location of the binary data (after having read how_many T)}
    */
    template <class T>
    std::tuple<bool, amrex::Vector<T>, const char*>parse_raw_data_vec(
        const char* p_data, size_t how_many, const char* const p_last)
    {
        amrex::Vector<T> res;
        if(p_data + sizeof(T)*how_many > p_last)
            return std::make_tuple(false, res, nullptr);

        auto r_data = reinterpret_cast<const T*>(p_data);

        res.assign(r_data, r_data + how_many);

        p_data += sizeof(T)*how_many;
        return std::make_tuple(true, res, p_data);
    }

    /**
    * This function safely extracts a T from raw binary data.
    * T must be a simple datatype (e.g. an int, a float, a double...).
    *
    * @param[in] p_data a pointer to the binary stream
    * @param[in] p_last a pointer to the last element of the char* array
    * @return {a tuple containing
    * 1) flag (which is false if p_last is exceeded)
    * 2) a T
    * 3) a pointer to a new location of the binary data (after having read 1 T)}
    */
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

    /**
    * This function converts a C-style array of T into
    * a Vector<char> (i.e. raw binary data) and adds it
    * to an existing Vector<char> passed by reference
    * @param[in] p_data a pointer to the beginning of the array
    * @param[in] how_many number of elements of type T in the array
    * @param[in,out] raw_data data will be appended to this vector
    */
    template <class T>
    void add_data_to_vector_char (
        const T* p_data, size_t how_many, amrex::Vector<char>& raw_data)
    {
        raw_data.insert(
            raw_data.end(),
            reinterpret_cast<const char*>(p_data),
            reinterpret_cast<const char*>(p_data) +
            sizeof(T)*how_many
        );
    }
};

#endif //WARPX_amrex_qed_table_parser_helper_functions_h_
