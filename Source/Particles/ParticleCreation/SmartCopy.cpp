#include "SmartCopy.H"

SmartCopyTag getSmartCopyTag (const NameMap& src, const NameMap& dst)
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
