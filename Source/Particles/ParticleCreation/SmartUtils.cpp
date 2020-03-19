/* Copyright 2019-2020 Andrew Myers, Axel Huebl,
 * Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SmartUtils.H"

PolicyVec getPolicies (const NameMap& names) noexcept
{
    PolicyVec policies;
    policies.resize(names.size());
    for (const auto& kv : names)
    {
        policies[kv.second] = initialization_policies[kv.first];
    }
    return policies;
}

SmartCopyTag getSmartCopyTag (const NameMap& src, const NameMap& dst) noexcept
{
    SmartCopyTag tag;

    // we use the fact that maps are sorted
    auto i_src = src.begin();
    auto i_dst = dst.begin();
    while ( (i_src != src.end()) and (i_dst != dst.end()) )
    {
        if (i_src->first < i_dst->first)
        {
            // names are not the same and src is lower
            ++i_src;
        }
        else if (i_src->first > i_dst->first)
        {
            // names are not the same and dst is lower
            ++i_dst;
        }
        else
        {
            // name is in both...
            tag.common_names.push_back(i_src->first);
            tag.src_comps.push_back(i_src->second);
            tag.dst_comps.push_back(i_dst->second);
            ++i_src;
            ++i_dst;
        }
    }

    return tag;
}
