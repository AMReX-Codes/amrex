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
        amrex::detail::color_tags_by_dfab(tags, color, max_color);

        int errors = 0;

        // Deterministic expectations for this setup with greedy coloring:
        // dA group (indices 1:b0, 2:b1, 0:b2 sorted by box): colors [0,1,0]
        if (!(color[1] == 0 && color[2] == 1 && color[0] == 0)) {
            amrex::Print() << "Error: unexpected colors for dA: got {"
                           << color[1] << "," << color[2] << "," << color[0]
                           << "}, expected {0,1,0}\n";
            ++errors;
        }

        // dB group (indices 3:b3, 4:b4): colors [0,1]
        if (!(color[3] == 0 && color[4] == 1)) {
            amrex::Print() << "Error: unexpected colors for dB: got {"
                           << color[3] << "," << color[4] << "}, expected {0,1}\n";
            ++errors;
        }

        if (max_color != 1) {
            amrex::Print() << "Error: unexpected max_color " << max_color << ", expected 1\n";
            ++errors;
        }

        if (errors != 0) {
            amrex::Abort("color_tags_by_dfab test failed");
        }

        // Second test: shuffled insertion order should yield same expected colors.
        Vector<Array4CopyTag<T,T>> tags2;
        tags2.reserve(5);
        // New insertion order: b4, b1, b0, b3, b2
        int idx_b4 = tags2.size(); tags2.push_back(Array4CopyTag<T,T>{dB.array(), s2.const_array(), b4, Dim3{0,0,0}});
        int idx_b1 = tags2.size(); tags2.push_back(Array4CopyTag<T,T>{dA.array(), s2.const_array(), b1, Dim3{0,0,0}});
        int idx_b0 = tags2.size(); tags2.push_back(Array4CopyTag<T,T>{dA.array(), s1.const_array(), b0, Dim3{0,0,0}});
        int idx_b3 = tags2.size(); tags2.push_back(Array4CopyTag<T,T>{dB.array(), s1.const_array(), b3, Dim3{0,0,0}});
        int idx_b2 = tags2.size(); tags2.push_back(Array4CopyTag<T,T>{dA.array(), s3.const_array(), b2, Dim3{0,0,0}});

        Vector<int> color2;
        int max_color2 = -1;
        amrex::detail::color_tags_by_dfab(tags2, color2, max_color2);

        int errors2 = 0;
        // Expect same mapping: b0->0, b1->1, b2->0; b3->0, b4->1
        if (!(color2[idx_b0] == 0 && color2[idx_b1] == 1 && color2[idx_b2] == 0)) {
            amrex::Print() << "Error: unexpected colors for dA (shuffled): got {"
                           << color2[idx_b0] << "," << color2[idx_b1] << "," << color2[idx_b2]
                           << "}, expected {0,1,0}\n";
            ++errors2;
        }
        if (!(color2[idx_b3] == 0 && color2[idx_b4] == 1)) {
            amrex::Print() << "Error: unexpected colors for dB (shuffled): got {"
                           << color2[idx_b3] << "," << color2[idx_b4] << "}, expected {0,1}\n";
            ++errors2;
        }
        if (max_color2 != 1) {
            amrex::Print() << "Error: unexpected max_color (shuffled) " << max_color2 << ", expected 1\n";
            ++errors2;
        }
        if (errors2 != 0) {
            amrex::Abort("color_tags_by_dfab shuffled-order test failed");
        }

        // Third test: three-box chain overlap on the same dfab (dA): c0 overlaps c1, c1 overlaps c2, c0 does not overlap c2.
        Box c0(IntVect(AMREX_D_DECL( 0, 0, 0)), IntVect(AMREX_D_DECL(3,3,3)));
        Box c1(IntVect(AMREX_D_DECL( 2, 2, 2)), IntVect(AMREX_D_DECL(5,5,5))); // overlaps c0
        Box c2(IntVect(AMREX_D_DECL( 4, 4, 4)), IntVect(AMREX_D_DECL(7,7,7))); // overlaps c1, disjoint from c0

        Vector<Array4CopyTag<T,T>> tags3;
        tags3.reserve(3);
        int id_c2 = tags3.size(); tags3.push_back(Array4CopyTag<T,T>{dA.array(), s3.const_array(), c2, Dim3{0,0,0}});
        int id_c0 = tags3.size(); tags3.push_back(Array4CopyTag<T,T>{dA.array(), s1.const_array(), c0, Dim3{0,0,0}});
        int id_c1 = tags3.size(); tags3.push_back(Array4CopyTag<T,T>{dA.array(), s2.const_array(), c1, Dim3{0,0,0}});

        Vector<int> color3;
        int max_color3 = -1;
        amrex::detail::color_tags_by_dfab(tags3, color3, max_color3);

        int errors3 = 0;
        // Expected greedy coloring after sort by box: c0->0, c1->1, c2->0
        if (!(color3[id_c0] == 0 && color3[id_c1] == 1 && color3[id_c2] == 0)) {
            amrex::Print() << "Error: unexpected colors for 3-chain: got {"
                           << color3[id_c0] << "," << color3[id_c1] << "," << color3[id_c2]
                           << "}, expected {0,1,0}\n";
            ++errors3;
        }
        if (max_color3 != 1) {
            amrex::Print() << "Error: unexpected max_color (3-chain) " << max_color3 << ", expected 1\n";
            ++errors3;
        }
        if (errors3 != 0) {
            amrex::Abort("color_tags_by_dfab 3-chain test failed");
        }
    }
    amrex::Finalize();
    return 0;
}
