// Test for detail::color_tags_by_dfab coloring logic
#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_TagParallelFor.H>
#include <AMReX_FabArray.H> // pulls in AMReX_FabArrayCommI.H and AMReX_FBI.H inside namespace amrex
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>

#include <vector>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        using T = Real;

        // Create two destination fabs and three source fabs
        Box fullbx(IntVect(AMREX_D_DECL(0,0,0)), IntVect(AMREX_D_DECL(15,15,15)));
        FArrayBox dA(fullbx, 1);
        FArrayBox dB(fullbx, 1);
        FArrayBox s1(fullbx, 1);
        FArrayBox s2(fullbx, 1);
        FArrayBox s3(fullbx, 1);

        // Construct boxes: b0 and b1 overlap; b2 is disjoint from both
        Box b0(IntVect(AMREX_D_DECL(0,0,0)),  IntVect(AMREX_D_DECL(4,4,4)));
        Box b1(IntVect(AMREX_D_DECL(2,2,2)),  IntVect(AMREX_D_DECL(6,6,6))); // overlaps b0
        Box b2(IntVect(AMREX_D_DECL(8,8,8)),  IntVect(AMREX_D_DECL(9,9,9))); // disjoint

        // For second dfab group, construct overlapping pair
        Box b3 = b0; // overlaps with b4
        Box b4(IntVect(AMREX_D_DECL(3,3,3)),  IntVect(AMREX_D_DECL(5,5,5)));

        Vector<Array4CopyTag<T,T>> tags;
        tags.reserve(5);

        // Intentionally push in non-sorted order to test internal sorting
        tags.push_back(Array4CopyTag<T,T>{dA.array(), s3.const_array(), b2, Dim3{0,0,0}}); // idx 0
        tags.push_back(Array4CopyTag<T,T>{dA.array(), s1.const_array(), b0, Dim3{0,0,0}}); // idx 1
        tags.push_back(Array4CopyTag<T,T>{dA.array(), s2.const_array(), b1, Dim3{0,0,0}}); // idx 2
        tags.push_back(Array4CopyTag<T,T>{dB.array(), s1.const_array(), b3, Dim3{0,0,0}}); // idx 3
        tags.push_back(Array4CopyTag<T,T>{dB.array(), s2.const_array(), b4, Dim3{0,0,0}}); // idx 4

        Vector<int> color;
        int max_color = -1;
        detail::color_tags_by_dfab(tags, color, max_color);

        int errors = 0;

        // On dA: idx 1 (b0) and idx 2 (b1) overlap => different colors
        if (color[1] == color[2]) {
            amrex::Print() << "Error: overlapping tags on same dfab have same color\n";
            ++errors;
        }

        // On dA: idx 1 (b0) and idx 0 (b2) do not overlap => should be allowed same color
        if ( (color[0] == color[2]) && (color[1] == color[2]) ) {
            amrex::Print() << "Error: disjoint tag reused conflicting color\n";
            ++errors;
        }

        // On dB: idx 3 (b3) and idx 4 (b4) overlap => different colors
        if (color[3] == color[4]) {
            amrex::Print() << "Error: overlapping tags on second dfab have same color\n";
            ++errors;
        }

        // Max color must be at least 1 because we have overlaps
        if (max_color < 1) {
            amrex::Print() << "Error: max_color < 1 despite overlapping tags\n";
            ++errors;
        }

        if (errors != 0) {
            amrex::Abort("color_tags_by_dfab test failed");
        }
    }
    amrex::Finalize();
    return 0;
}
