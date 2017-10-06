#include <cmath>
#include <MyAmr.H>

using namespace amrex;

// Note that tags is built on level lev grids coarsened by bf_lev[lev].
void
MyAmr::ManualTagsPlacement (int lev, TagBoxArray& tags,
                            const Vector<IntVect>& bf_lev)
{
    // let's refine along the line r = (0,0,0) + (1.0,0.9,1.1)*k
    std::array<Real,3> n {1.0, 0.9, 1.1};
    // make it a unit vector
    {
        Real tmp = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        for (auto& x: n) x /= tmp;
    }

    for (MFIter mfi(tags); mfi.isValid(); ++mfi)
    {
        TagBox& tag = tags[mfi];
        const Box& bx = tag.box();
        for (IntVect cell=bx.smallEnd(); cell <= bx.bigEnd(); bx.next(cell))
        {
            std::array<Real,3> p {cell[0]+0.5, cell[1]+0.5, cell[2]+0.5}; // cell center
            // Compute distance bwtween p and the line
            Real pdotn = p[0]*n[0] + p[1]*n[1] + p[2]*n[2];
            // ||p - (p dot n) n||^2
            Real d = 0.0;
            for (int dim = 0; dim < 3; ++dim) {
                Real tmp = p[dim] - pdotn * n[dim];
                d += tmp*tmp;
            }

            if (d <= 3.0) {
                tag(cell) = TagBox::SET;
            }
        }
    }
}

