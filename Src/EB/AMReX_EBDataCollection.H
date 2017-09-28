#ifndef AMREX_EB_DATA_COLLECTION_H_
#define AMREX_EB_DATA_COLLECTION_H_

#include <AMReX_Geometry.H>

namespace amrex {

class EBDataCollection
{
public:

    EBDataCollection (const Geometry& a_geom);

    ~EBDataCollection () {;}

    EBDataCollection (const EBDataCollection&) = delete;
    EBDataCollection (EBDataCollection&&) = delete;
    EBDataCollection& operator= (const EBDataCollection&) = delete;
    EBDataCollection& operator= (EBDataCollection&&) = delete;

private:

    Geometry m_geom;
};

}

#endif
