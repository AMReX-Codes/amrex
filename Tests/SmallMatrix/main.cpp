#include <AMReX.H>
#include <AMReX_Math.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    static_assert(Order::C == Order::RowMajor &&
                  Order::F == Order::ColumnMajor);

    amrex::Initialize(argc, argv);
    {
        SmallMatrix<Real,3,4> m34{};
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 3; ++i) {
                AMREX_ALWAYS_ASSERT(m34(i,j) == 0.0_rt);
            }
        }
    }
    {
        SmallVector<Real,3> cv{};
        SmallRowVector<Real,3> rv{};
        SmallVector<int,3> cv2{1,2,3};
        SmallRowVector<int,3> rv2{0,10,20};
        SmallVector<int,5> cv3{0,1,2};
        for (int j = 0; j < 3; ++j) {
            AMREX_ALWAYS_ASSERT(cv(j) == 0.0_rt &&
                                rv(j) == 0.0_rt &&
                                cv2(j) == j+1 &&
                                rv2(j) == j*10 &&
                                cv3(j) == j);
        }
        AMREX_ALWAYS_ASSERT(cv3(3) == 0 && cv3(4) == 0);
    }
    {
        SmallMatrix<int,3,4> m34{{0,3,6,9},
                                 {1,4,7,10},
                                 {2,5,8,11}};
        int v = 0;
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 3; ++i) {
                AMREX_ALWAYS_ASSERT(m34(i,j) == v++);
            }
        }
        std::cout << m34;
    }
    {
        SmallMatrix<int,3,4,Order::C> m34{{0,1,2,3},
                                          {4,5,6,7},
                                          {8,9,10,11}};
        int v = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 4; ++j) {
                AMREX_ALWAYS_ASSERT(m34(i,j) == v++);
            }
        }
    }
    {
        auto v3 = SmallVector<double,3>::Zero();
        v3[0] = 1.;
        v3(1) = 2.;
        v3[2] = 3.;
        auto m33 = SmallMatrix<double,3,3>::Identity();
        auto r = m33*v3;
        AMREX_ALWAYS_ASSERT(almostEqual(r[0],v3[0]) &&
                            almostEqual(r[1],v3[1]) &&
                            almostEqual(r[2],v3[2]));
    }
    {
        SmallMatrix<int,4,3,Order::C> A{{1, 0, 1},
                                        {2, 1, 1},
                                        {0, 1, 1},
                                        {1, 1, 2}};
        SmallMatrix<int,3,3,Order::C> B{{1, 2, 1},
                                        {2, 3, 1},
                                        {4, 2, 2}};
        SmallMatrix<int,3,1,Order::C> C{10, 8, 6};
        auto ABC = A*B*C;
        AMREX_ALWAYS_ASSERT(ABC(0,0) == 100 &&
                            ABC(1,0) == 182 &&
                            ABC(2,0) == 118 &&
                            ABC(3,0) == 218);
    }
    {
        SmallMatrix<int,3,4,Order::F> A{{1, 2, 0, 1},
                                        {0, 1, 1, 1},
                                        {1, 1, 1, 2}};
        SmallMatrix<int,3,3,Order::F> B{{1, 2, 4},
                                        {2, 3, 2},
                                        {1, 1, 2}};
        SmallMatrix<int,1,3,Order::F> C{10, 8, 6};
        auto ABC = A.transpose()*B.transposeInPlace()*C.transpose();
        AMREX_ALWAYS_ASSERT(ABC(0,0) == 100 &&
                            ABC(1,0) == 182 &&
                            ABC(2,0) == 118 &&
                            ABC(3,0) == 218);
    }
    {
        SmallMatrix<int, 3, 4> m;
        m.setVal(2);
        using M = decltype(m);
        AMREX_ALWAYS_ASSERT(m.product() == Math::powi<M::row_size*M::column_size>(2));
        AMREX_ALWAYS_ASSERT(m.sum() == 2*m.row_size*m.column_size);
    }
    {
        SmallMatrix<double, 5, 5> m{{1.0, 3.4, 4.5, 5.6, 6.7},
                                    {1.3, 2.0, 4.5, 5.6, 6.7},
                                    {1.3, 1.0, 3.0, 5.6, 6.7},
                                    {1.3, 1.4, 4.5, 4.0, 6.7},
                                    {1.3, 1.0, 4.5, 5.6, 5.0}};
        AMREX_ALWAYS_ASSERT(m.trace() == double(1+2+3+4+5));
    }
    {
        SmallMatrix<int,2,3> a{{+1, +2, +3},
                               {+7, +8, +9}};
        SmallMatrix<int,2,3> b{{-1, -2, -3},
                               {-7, -8, -9}};
        auto c = a*2 + 2*b;
        for (auto const& x : c) {
            AMREX_ALWAYS_ASSERT(x == 0);
        }
    }
    {
        SmallMatrix<int,2,3> a{{+1, +2, +3},
                               {+7, +8, +9}};
        SmallMatrix<int,2,3> b{{-1, -2, -3},
                               {-7, -8, -9}};
        auto c = -a - b;
        for (auto const& x : c) {
            AMREX_ALWAYS_ASSERT(x == 0);
        }
    }
    {
        SmallMatrix<int,2,3> a{{+1, +2, +3},
                               {+7, +8, +9}};
        SmallMatrix<int,2,3> b;
        b.setVal(-1);
        AMREX_ALWAYS_ASSERT(a.dot(b) == -30);
    }
    amrex::Finalize();
}
