//BL_COPYRIGHT_NOTICE

#ifndef _RegType_H_
#define _RegType_H_

#include <aString.H>

enum RegType
{
    interior, exterior, periodic, refWall, inflow, outflow, REG_TYPE_SIZE
};

inline
RegType
RegTypeSet(const aString& ass)
{
    if ( ass == "interior" ) return interior ;
    if ( ass == "exterior" ) return exterior ;
    if ( ass == "periodic" ) return periodic ;
    if ( ass == "refWall"  ) return refWall  ;
    if ( ass == "inflow"   ) return inflow   ;
    if ( ass == "outflow"  ) return outflow  ;
    return REG_TYPE_SIZE;
}


#endif /*_RegType_H_*/
