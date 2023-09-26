#ifndef AMREX_PAR_REDUCE_H_
#define AMREX_PAR_REDUCE_H_
#include <AMReX_Config.H>

#include <AMReX_Reduce.H>

namespace amrex {

/**
 * \brief Parallel reduce for MultiFab/FabArray.
 *
 * This performs reduction over a MultiFab's valid and specified ghost regions. For
 * example, the code below computes the minimum of the first MultiFab and the maximum of
 * the second MultiFab.
 \verbatim
     auto const& ma1 = mf1.const_arrays();
     auto const& ma2 = mf2.const_arrays();
     GpuTuple<Real,Real> mm = ParReduce(TypeList<ReduceOpMin,ReduceOpMax>{},
                                        TypeList<Real,Real>{},
                                        mf1, mf1.nGrowVect(),
     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
         -> GpuTuple<Real,Real>
     {
         return { ma1[box_no](i,j,k), ma2[box_no](i,j,k) };
     });
 \endverbatim
 *
 * \tparam Ops... reduce operators (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *                ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam Ts...  data types (e.g., Real, int, etc.)
 * \tparam FAB    MultiFab/FabArray type
 * \tparam F      callable type like a lambda function
 *
 * \param operation_list list of reduce operators
 * \param type_list      list of data types
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param nghost the number of ghost cells included in the iteration space
 * \param f      a callable object returning GpuTuple<Ts...>.  It takes four ints,
 *               where the first int is the local box index and the others are
 *               spatial indices for x, y, and z-directions.
 *
 * \return reduction result (GpuTuple<Ts...>)
 */
template <typename... Ops, typename... Ts, typename FAB, typename F,
          typename foo = std::enable_if_t<IsBaseFab<FAB>::value> >
typename ReduceData<Ts...>::Type
ParReduce (TypeList<Ops...> operation_list, TypeList<Ts...> type_list,
           FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    ReduceOps<Ops...> reduce_op;
    ReduceData<Ts...> reduce_data(reduce_op);
    reduce_op.eval(fa, nghost, reduce_data, std::forward<F>(f));
    return reduce_data.value(reduce_op);
}

/**
 * \brief Parallel reduce for MultiFab/FabArray.
 *
 * This performs reduction over a MultiFab's valid and specified ghost regions. For
 * example, the code below computes the sum of the processed data in a MultiFab.
 \verbatim
     auto const& ma = mf.const_arrays();
     Real ektot = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{},
                            mf, IntVect(0),
     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
         -> GpuTuple<Real>
     {
         auto rho = ma[box_no](i,j,k,0);
         auto mx = ma[box_no](i,j,k,1);
         auto my = ma[box_no](i,j,k,2);
         auto mz = ma[box_no](i,j,k,3);
         auto ek = (mx*mx+my*my+mz*mz)/(2.*rho);
         return { ek };
     });
 \endverbatim
 *
 * \tparam Op  Reduce operator (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *             ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam T   data type (e.g., Real, int, etc.)
 * \tparam FAB MultiFab/FabArray type
 * \tparam F   callable type like a lambda function
 *
 * \param operation_list a reduce operator stored in TypeList
 * \param type_list      a data type stored in TypeList
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param nghost the number of ghost cells included in the iteration space
 * \param f      a callable object returning GpuTuple<T>.  It takes four ints,
 *               where the first int is the local box index and the others are
 *               spatial indices for x, y, and z-directions.
 *
 * \return reduction result (T)
 */
template <typename Op, typename T, typename FAB, typename F,
          typename foo = std::enable_if_t<IsBaseFab<FAB>::value> >
T
ParReduce (TypeList<Op> operation_list, TypeList<T> type_list,
           FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    ReduceOps<Op> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    reduce_op.eval(fa, nghost, reduce_data, std::forward<F>(f));
    auto const& hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv);
}

/**
 * \brief Parallel reduce for MultiFab/FabArray.
 *
 * This performs reduction over a MultiFab's valid and specified ghost regions and
 * components.  For example, the code below computes the minimum of the first MultiFab
 * and the maximum of the second MultiFab.
 \verbatim
     auto const& ma1 = mf1.const_arrays();
     auto const& ma2 = mf2.const_arrays();
     GpuTuple<Real,Real> mm = ParReduce(TypeList<ReduceOpMin,ReduceOpMax>{},
                                        TypeList<Real,Real>{},
                                        mf1, mf1.nGrowVect(), mf1.nComp(),
     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
         -> GpuTuple<Real,Real>
     {
         return { ma1[box_no](i,j,k,n), ma2[box_no](i,j,k,n) };
     });
 \endverbatim
 *
 * \tparam Ops... reduce operators (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *                ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam Ts...  data types (e.g., Real, int, etc.)
 * \tparam FAB    MultiFab/FabArray type
 * \tparam F      callable type like a lambda function
 *
 * \param operation_list list of reduce operators
 * \param type_list      list of data types
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param nghost the number of ghost cells included in the iteration space
 * \param ncomp  the number of components in the iteration space
 * \param f      a callable object returning GpuTuple<Ts...>.  It takes five ints,
 *               where the first int is the local box index, the next three are
 *               spatial indices for x, y, and z-directions, and the last is for
 *               component.
 *
 * \return reduction result (GpuTuple<Ts...>)
 */
template <typename... Ops, typename... Ts, typename FAB, typename F,
          typename foo = std::enable_if_t<IsBaseFab<FAB>::value> >
typename ReduceData<Ts...>::Type
ParReduce (TypeList<Ops...> operation_list, TypeList<Ts...> type_list,
           FabArray<FAB> const& fa, IntVect const& nghost, int ncomp, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    ReduceOps<Ops...> reduce_op;
    ReduceData<Ts...> reduce_data(reduce_op);
    reduce_op.eval(fa, nghost, ncomp, reduce_data, std::forward<F>(f));
    return reduce_data.value(reduce_op);
}

/**
 * \brief Parallel reduce for MultiFab/FabArray.
 *
 * This performs reduction over a MultiFab's valid and specified ghost regions. For
 * example, the code below computes the sum of the data in a MultiFab.
 \verbatim
     auto const& ma = mf.const_arrays();
     Real ektot = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{},
                            mf, mf.nGrowVect(), mf.nComp(),
     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
         -> GpuTuple<Real>
     {
         return { ma[box_no](i,j,k,n) };
     });
 \endverbatim
 *
 * \tparam Op  Reduce operator (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *             ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam T   data type (e.g., Real, int, etc.)
 * \tparam FAB MultiFab/FabArray type
 * \tparam F   callable type like a lambda function
 *
 * \param operation_list a reduce operator stored in TypeList
 * \param type_list      a data type stored in TypeList
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param nghost the number of ghost cells included in the iteration space
 * \param ncomp  the number of components in the iteration space
 * \param f      a callable object returning GpuTuple<T>.  It takes four ints,
 *               where the first int is the local box index and the others are
 *               spatial indices for x, y, and z-directions.
 *
 * \return reduction result (T)
 */
template <typename Op, typename T, typename FAB, typename F,
          typename foo = std::enable_if_t<IsBaseFab<FAB>::value> >
T
ParReduce (TypeList<Op> operation_list, TypeList<T> type_list,
           FabArray<FAB> const& fa, IntVect const& nghost, int ncomp, F&& f)
{
    amrex::ignore_unused(operation_list,type_list);
    ReduceOps<Op> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    reduce_op.eval(fa, nghost, ncomp, reduce_data, std::forward<F>(f));
    auto const& hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv);
}

/**
 * \brief Parallel reduce for MultiFab/FabArray.
 *
 * This performs reduction over a MultiFab's valid region. For
 * example, the code below computes the minimum of the first MultiFab and the maximum of
 * the second MultiFab.
 \verbatim
     auto const& ma1 = mf1.const_arrays();
     auto const& ma2 = mf2.const_arrays();
     GpuTuple<Real,Real> mm = ParReduce(TypeList<ReduceOpMin,ReduceOpMax>{},
                                        TypeList<Real,Real>{}, mf1,
     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
         -> GpuTuple<Real,Real>
     {
         return { ma1[box_no](i,j,k), ma2[box_no](i,j,k) };
     });
 \endverbatim
 *
 * \tparam Ops... reduce operators (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *                ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam Ts...  data types (e.g., Real, int, etc.)
 * \tparam FAB    MultiFab/FabArray type
 * \tparam F      callable type like a lambda function
 *
 * \param operation_list list of reduce operators
 * \param type_list      list of data types
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param f      a callable object returning GpuTuple<Ts...>.  It takes four ints,
 *               where the first int is the local box index and the others are
 *               spatial indices for x, y, and z-directions.
 *
 * \return reduction result (GpuTuple<Ts...>)
 */
template <typename... Ops, typename... Ts, typename FAB, typename F,
          typename foo = std::enable_if_t<IsBaseFab<FAB>::value> >
typename ReduceData<Ts...>::Type
ParReduce (TypeList<Ops...> operation_list, TypeList<Ts...> type_list,
           FabArray<FAB> const& fa, F&& f)
{
    return ParReduce(operation_list, type_list, fa, IntVect(0), std::forward<F>(f));
}

/**
 * \brief Parallel reduce for MultiFab/FabArray.
 *
 * This performs reduction over a MultiFab's valid region. For
 * example, the code below computes the sum of the processed data in a MultiFab.
 \verbatim
     auto const& ma = mf.const_arrays();
     Real ektot = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mf,
     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
         -> GpuTuple<Real>
     {
         auto rho = ma[box_no](i,j,k,0);
         auto mx = ma[box_no](i,j,k,1);
         auto my = ma[box_no](i,j,k,2);
         auto mz = ma[box_no](i,j,k,3);
         auto ek = (mx*mx+my*my+mz*mz)/(2.*rho);
         return { ek };
     });
 \endverbatim
 *
 * \tparam Op  Reduce operator (e.g., ReduceOpSum, ReduceOpMin, ReduceOpMax,
 *             ReduceOpLogicalAnd, and ReduceOpLogicalOr)
 * \tparam T   data type (e.g., Real, int, etc.)
 * \tparam FAB MultiFab/FabArray type
 * \tparam F   callable type like a lambda function
 *
 * \param operation_list a reduce operator stored in TypeList
 * \param type_list      a data type stored in TypeList
 * \param fa     a MultiFab/FabArray object used to specify the iteration space
 * \param f      a callable object returning GpuTuple<T>.  It takes four ints,
 *               where the first int is the local box index and the others are
 *               spatial indices for x, y, and z-directions.
 *
 * \return reduction result (T)
 */
template <typename Op, typename T, typename FAB, typename F,
          typename foo = std::enable_if_t<IsBaseFab<FAB>::value> >
T
ParReduce (TypeList<Op> operation_list, TypeList<T> type_list,
           FabArray<FAB> const& fa, F&& f)
{
    return ParReduce(operation_list, type_list, fa, IntVect(0), std::forward<F>(f));
}

}
#endif
