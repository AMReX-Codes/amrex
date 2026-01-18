#include <AMReX_Array4.H>
#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Random.H>
#include <AMReX_Reduce.H>

using namespace amrex;

template <typename T>
void test_array4 (Array4<T> const& a, T tot)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    Box box(IntVectShrink<AMREX_SPACEDIM,4>(a.begin),
            IntVectShrink<AMREX_SPACEDIM,4>(a.end-1));
    if (a.nComp() == 1) {
        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            auto v0 = a(i,j,k);
            auto v1 = a(i,j,k,0);
            auto v2 = a(IntVectND<4>(i,j,k,0));
            auto v3 = a(IntVect(AMREX_D_DECL(i,j,k)));
            auto v4 = a(Dim3{i,j,k});
            auto* p0 = a.ptr(i,j,k);
            auto* p1 = a.ptr(i,j,k,0);
            auto* p2 = a.ptr(IntVectND<4>(i,j,k,0));
            auto* p3 = a.ptr(IntVect(AMREX_D_DECL(i,j,k)));
            auto* p4 = a.ptr(Dim3{i,j,k});
            return (v0+v1+v2+v3+v4+*p0+*p1+*p2+*p3+*p4)/10;
        });
    } else {
        reduce_op.eval(box, a.nComp(), reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) -> ReduceTuple
        {
            auto v0 = a(i,j,k,n);
            auto v1 = a(IntVect(AMREX_D_DECL(i,j,k)),n);
            auto v2 = a(Dim3{i,j,k},n);
            auto* p0 = a.ptr(i,j,k,n);
            auto* p1 = a.ptr(IntVect(AMREX_D_DECL(i,j,k)),n);
            auto* p2 = a.ptr(Dim3{i,j,k},n);
            return (v0+v1+v2+*p0+*p1+*p2)/6;
        });
    }
    ReduceTuple hv = reduce_data.value(reduce_op);
    AMREX_ALWAYS_ASSERT(tot == amrex::get<0>(hv));
}

template <typename T, int N>
void test_comp (ArrayND<T,N,true> const& a, T tot)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    BoxND<N-1> box(IntVectShrink<N-1,N>(a.begin),
                   IntVectShrink<N-1,N>(a.end-1));
    reduce_op.eval(box, a.nComp(), reduce_data,
    [=] AMREX_GPU_DEVICE (IntVectND<N-1> const& iv, int n) -> ReduceTuple
    {
        auto iv_full = IntVectExpand<N,N-1>(iv,n);
        auto v0 = a(iv,n);
        auto v1 = a(iv_full);
        auto* p0 = a.ptr(iv,n);
        auto* p1 = a.ptr(iv_full);
        ArrayND<T const, N, true> b(a, n);
        T v2 = amrex::Apply([&] (auto&&... i) {
            return b(i...);
        }, iv);
        T const* p2 = amrex::Apply([&] (auto&&... i) {
            return b.ptr(i...,0);
        }, iv);
        return (v0+v1+v2+*p0+*p1+*p2)/6;
    });
    ReduceTuple hv = reduce_data.value(reduce_op);
    AMREX_ALWAYS_ASSERT(tot == amrex::get<0>(hv));
}

template <typename T, int N, bool C>
void test (ArrayND<T,N,C>& a)
{
    auto sz = a.size();
    auto* p = a.dataPtr();
    amrex::ParallelForRNG(sz, [=] AMREX_GPU_DEVICE (std::size_t i,
                                                    RandomEngine const& eng)
    {
        p[i] = amrex::Random_int(1000000, eng);
    });
    T tot = Reduce::Sum<T>(sz, p, 0);

    if constexpr (std::remove_reference_t<decltype(a)>::IsArray4_v) {
        test_array4(a, tot);
    }

    if constexpr (std::remove_reference_t<decltype(a)>::IsLastDimComponent_v) {
        test_comp(a, tot);
    }

    {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        BoxND<N> box(a.begin, a.end-1);
        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_DEVICE (IntVectND<N> const& iv) -> ReduceTuple
        {
            auto v0 = a(iv);
            auto* p0 = a.ptr(iv);
            auto v1 = amrex::Apply([&] (auto&&... i) {
                return a(i...);
            }, iv);
            auto* p1 = amrex::Apply([&] (auto&&... i) {
                return a.ptr(i...);
            }, iv);
            return (v0+v1+*p0+*p1)/4;
        });
        ReduceTuple hv = reduce_data.value(reduce_op);
        AMREX_ALWAYS_ASSERT(tot == amrex::get<0>(hv));
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        ArrayND<Real,5> a{};
        AMREX_ALWAYS_ASSERT(!a && !a.ok() && a.size() == 0);
    }

    {
        BoxND<1> box(IntVectND<1>(-1), IntVectND<1>(50));
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*sizeof(Long));
        ArrayND<Long, 1> a(p, box);
        test(a);
        The_Arena()->free(p);
    }

    {
        BoxND<2> box(IntVectND<2>(0,1), IntVectND<2>(100,99));
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*sizeof(Long));
        ArrayND<Long, 2> a(p, box);
        ArrayND<Long const, 2> b(a);
        test(a);
        AMREX_ALWAYS_ASSERT(a.size() == b.size() && a.dataPtr() == b.dataPtr());
        The_Arena()->free(p);
    }

    {
        BoxND<1> box(IntVectND<1>(1), IntVectND<1>(100));
        int ncomp = 3;
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*ncomp*sizeof(Long));
        ArrayND<Long, 2, true> a(p, box, ncomp);
        test(a);
        The_Arena()->free(p);
    }

    {
        BoxND<3> box(IntVectND<3>(-2), IntVectND<3>(100,-1,99));
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*sizeof(Long));
        ArrayND<Long, 3> a(p, box);
        test(a);
        The_Arena()->free(p);
    }

    {
        BoxND<2> box(IntVectND<2>(-2), IntVectND<2>(100));
        int ncomp = 1;
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*ncomp*sizeof(Long));
        ArrayND<Long, 3, true> a(p, box, ncomp);
        test(a);
        The_Arena()->free(p);
    }

    {
        BoxND<4> box(IntVectND<4>(100), IntVectND<4>(150));
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*sizeof(Long));
        ArrayND<Long, 4> a(p, box);
        test(a);
        static_assert(!std::remove_reference_t<decltype(a)>::IsArray4_v);
        The_Arena()->free(p);
    }

    {
        BoxND<3> box(IntVectND<3>(100), IntVectND<3>(150));
        int ncomp = 4;
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*ncomp*sizeof(Long));
        ArrayND<Long, 4, true> a(p, box, ncomp);
        static_assert(std::remove_reference_t<decltype(a)>::IsArray4_v);
        The_Arena()->free(p);
    }

    {
        Box box(IntVect(AMREX_D_DECL(10,12,13)), IntVect(AMREX_D_DECL(36,35,43)));
        int ncomp = 4;
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*ncomp*sizeof(Long));
        Array4<Long> a(p, box, ncomp);
        test(a);
        AMREX_ALWAYS_ASSERT(a.contains(box.smallEnd()) &&
                            !a.contains(box.bigEnd()+1));
        auto cell = amrex::begin(box);
        AMREX_ALWAYS_ASSERT(a.contains(cell) && a.contains(cell.x,cell.y,cell.z)
                            && !a.contains(-cell.x,cell.y,cell.z));
        static_assert(std::remove_reference_t<decltype(a)>::IsArray4_v);
        The_Arena()->free(p);
    }

    {
        Box box(IntVect(AMREX_D_DECL(10,12,13)), IntVect(AMREX_D_DECL(46,25,23)));
        int ncomp = 1;
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*ncomp*sizeof(Long));
        Array4<Long> a(p, amrex::begin(box), amrex::end(box), ncomp);
        test(a);
        static_assert(std::remove_reference_t<decltype(a)>::IsArray4_v);
        The_Arena()->free(p);
    }

    {
        BoxND<5> box(IntVectND<5>(-3), IntVectND<5>(10,11,12,13,5));
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*sizeof(Long));
        ArrayND a(p, box);
        test(a);
        The_Arena()->free(p);
    }

    {
        BoxND<4> box(IntVectND<4>(0), IntVectND<4>(10,11,3,13));
        int ncomp = 2;
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*2*sizeof(Long));
        ArrayND a(p, box, ncomp);
        test(a);
        The_Arena()->free(p);
    }

    {
        BoxND<6> box(IntVectND<6>(0), IntVectND<6>(4,5,1,5,1,3));
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*sizeof(Long));
        ArrayND a(p, box.smallEnd(), box.bigEnd()+1);
        test(a);
        AMREX_ALWAYS_ASSERT(Long(a.size()) == box.numPts());
        static_assert(!std::remove_reference_t<decltype(a)>::IsArray4_v &&
                      !std::remove_reference_t<decltype(a)>::IsLastDimComponent_v);
        The_Arena()->free(p);
    }

    {
        BoxND<5> box(IntVectND<5>(-1,0,1,2,3), IntVectND<5>(7));
        int ncomp = 2;
        auto* p = (Long*) The_Arena()->alloc(box.numPts()*ncomp*sizeof(Long));
        ArrayND a(p, box.smallEnd(), box.bigEnd()+1, 2);
        test(a);
        AMREX_ALWAYS_ASSERT(Long(a.size()) == box.numPts()*ncomp);
        static_assert(!std::remove_reference_t<decltype(a)>::IsArray4_v &&
                      std::remove_reference_t<decltype(a)>::IsLastDimComponent_v);
        The_Arena()->free(p);
    }

    amrex::Finalize();
}
