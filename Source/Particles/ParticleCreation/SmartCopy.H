#ifndef SMART_COPY_H_
#define SMART_COPY_H_

#include <map>
#include <string>
#include <vector>

using NameMap = std::map<std::string, int>;

struct SmartCopyTag
{
    std::vector<std::string> common_names;
    std::vector<int> src_comps;
    std::vector<int> dst_comps;
    
    int size () const noexcept { return common_names.size(); }
};

SmartCopyTag getSmartCopyTag (const NameMap& src, const NameMap& dst);

#endif
