#ifndef AMREX_EB2_2D_C_H_
#define AMREX_EB2_2D_C_H_

namespace amrex { namespace EB2 {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
amrex_eb2_build_types (Box const& tbx, Box const& bxg2,
                       Array4<Real const> const& s,
                       Array4<EBCellFlag> const& cell,
                       Array4<Type_t> const& fx,
                       Array4<Type_t> const& fy)
{
    Box bx = amrex::getThreadBox(tbx, bxg2);
    auto lo = amrex::lbound(bx);
    auto hi = amrex::ubound(bx);
    const int k = 0;
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        if (    s(i,j  ,k) < 0.0 and s(i+1,j  ,k) < 0.0
            and s(i,j+1,k) < 0.0 and s(i+1,j+1,k) < 0.0)
        {
            cell(i,j,k).setRegular();
        }
        else if (s(i,j  ,k) >= 0.0 and s(i+1,j  ,k) >= 0.0
            and  s(i,j+1,k) >= 0.0 and s(i+1,j+1,k) >= 0.0)
        {
            cell(i,j,k).setCovered();
        }
        else
        {
            cell(i,j,k).setSingleValued();
        }
    }}

    // x-face
    bx = amrex::getThreadBox(tbx, amrex::surroundingNodes(bxg2,0));
    lo = amrex::lbound(bx);
    hi = amrex::ubound(bx);
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        if (s(i,j,k) < 0.0 and s(i,j+1,k) < 0.0) {
            fx(i,j,k) = Type::regular;
        } else if (s(i,j,k) >= 0.0 and s(i,j+1,k) >= 0.0) {
            fx(i,j,k) = Type::covered;
        } else {
            fx(i,j,k) = Type::irregular;
        }
    }}

    // y-face
    bx = amrex::getThreadBox(tbx, amrex::surroundingNodes(bxg2,1));
    lo = amrex::lbound(bx);
    hi = amrex::ubound(bx);
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        if (s(i,j,k) < 0.0 and s(i+1,j,k) < 0.0) {
            fy(i,j,k) = Type::regular;
        } else if (s(i,j,k) >= 0.0 and s(i+1,j,k) >= 0.0) {
            fy(i,j,k) = Type::covered;
        } else {
            fy(i,j,k) = Type::irregular;
        }
    }}
}

}}

#endif
