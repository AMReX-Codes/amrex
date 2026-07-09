
#ifndef AMREX_BLASSERT_H_
#define AMREX_BLASSERT_H_
#include <AMReX_Config.H>

#include <AMReX.H>

/**
 * \file AMReX_BLassert.H
 * \brief Assertion macros used across AMReX for runtime consistency checks.
 *
 * In debug builds (or when `AMREX_USE_ASSERTION` is defined) the macros expand
 * to runtime checks that call `amrex::Assert` if they fail. In release builds
 * without assertion support, the non-\*_ALWAYS variants become no-ops.
 *
 * - `AMREX_ASSERT(expr)` evaluates \p expr and aborts on failure, emitting the
 *   expression, file, and line number.
 * - `AMREX_ASSERT_WITH_MESSAGE(expr,msg)` behaves the same but also forwards
 *   the user-provided \p msg to `amrex::Assert`, allowing custom diagnostics.
 * - The `_ALWAYS_` variants (`AMREX_ALWAYS_ASSERT`, etc.) remain active even in
 *   non-debug builds.
 * - `AMREX_GPU_ASSERT` provides the same semantics inside GPU code paths when
 *   `AMREX_USE_GPU` is defined and either `AMREX_DEBUG` or `AMREX_USE_ASSERTION`
 *   is enabled; otherwise it is a no-op.
 */

#ifdef AMREX_ASSERT
#undef AMREX_ASSERT
#endif

#ifdef BL_ASSERT
#undef BL_ASSERT
#endif

#if !defined(AMREX_DEBUG) && !defined(AMREX_USE_ASSERTION)

#define AMREX_ASSERT_WITH_MESSAGE(EX,MSG) ((void)0)
#define AMREX_ASSERT(EX) ((void)0)
#define    BL_ASSERT(EX) ((void)0)

#else

#define AMREX_ASSERT_WITH_MESSAGE(EX,MSG) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__ , MSG)
#define AMREX_ASSERT(EX) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__)
#define    BL_ASSERT(EX) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__)

#endif

#define AMREX_ALWAYS_ASSERT_WITH_MESSAGE(EX,MSG) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__ , MSG)
#define AMREX_ALWAYS_ASSERT(EX) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__)


#if !defined(AMREX_USE_GPU) || !defined(AMREX_DEBUG) && !defined(AMREX_USE_ASSERTION)
#define AMREX_GPU_ASSERT(EX) ((void)0)
#else
#define AMREX_GPU_ASSERT(EX) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__)
#endif

#endif /* AMREX_BLASSERT_H_ */
