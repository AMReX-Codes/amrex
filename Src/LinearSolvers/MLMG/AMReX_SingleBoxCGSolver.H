#ifndef AMREX_SINGLE_BOX_CG_SOLVER_H_
#define AMREX_SINGLE_BOX_CG_SOLVER_H_
#include <AMReX_Config.H>

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

#include <AMReX_BaseFab.H>
#include <AMReX_Print.H>

#include <iomanip>

namespace amrex
{

namespace detail {
    //! Everything bicgstab_solve copies back to the host, in one buffer.
    template <typename T>
    struct SingleBoxCGResult
    {
        int ret;     //!< return code, same meaning as in MLCGSolver
        int active;  //!< nonzero while the multi-kernel path keeps iterating
        int niters;  //!< number of iterations
        T rnorm0;    //!< initial residual norm
        T rnorm;     //!< final residual norm
    };
}

/**
 * \brief BiCGStab solver for a single Box that runs entirely on the GPU.
 *
 * This is currently available for CUDA and HIP only.  The return codes and
 * the verbose output follow MLCGSolver.
 *
 * \param box region for solve
 * \param x solution with initial value already set to zero
 * \param b RHS
 * \param lp operator
 * \param eps_rel relative tolerance
 * \param eps_abs absolute tolerance
 * \param maxiter max iteration
 * \param verbose verbosity, same meaning as in MLCGSolver
 * \param niters number of iterations performed (output)
 */
template <typename T, class LP>
int bicgstab_solve (Box const& box, BaseFab<T>& x, BaseFab<T> const& b, LP const& lp,
                    T eps_rel, T eps_abs, int maxiter, int verbose, int& niters)
{
    const int ncomp = x.nComp();

    enum Index {
        i_p, i_r, i_rh, i_v, i_t, nidx
    };

    BaseFab<T> fab(box, ncomp*nidx);
    auto const& a = fab.array();
    auto const& sol = x.array();
    auto const& rhs = b.const_array();

    // Single buffer for everything that goes back to the host.
    Gpu::Buffer<detail::SingleBoxCGResult<T>> res_buf(1);
    auto* resd = res_buf.data();

    Gpu::DeviceVector<T> reduce_buf;
    Gpu::DeviceVector<T> scalar_buf;

    constexpr int MT = AMREX_GPU_MAX_THREADS;

    const BoxIndexer indexer(box);

    bool const is_big = box.numPts() > (MT*4);

    if (is_big)
    {
        auto nblocks_long = (box.numPts() + MT - 1) / MT;
        AMREX_ALWAYS_ASSERT(nblocks_long <=
                            Long(std::numeric_limits<unsigned int>::max()/2));
        auto nblocks = static_cast<unsigned int>(nblocks_long);

        reduce_buf.resize(nblocks*2);
        enum ScalarIndex {
            s_rnorm0, s_rnorm, s_rho, s_rho_1, s_alpha, s_omega, s_beta, nscalar
        };
        scalar_buf.resize(nscalar);

        auto* reduce = reduce_buf.data();
        auto* scalar = scalar_buf.data();
        auto* active = &(resd->active);
        auto const npts = indexer.numPts();

        // Initialize r and rh, and write one rnorm contribution per block.
        amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
        {
            auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
            T rnorm = 0;
            if (ipt < npts) {
                auto const cell = indexer.intVect(ipt);
                for (int n = 0; n < ncomp; ++n) {
                    T v = rhs(cell,n);
                    lp.normalize(cell, n, v);
                    a(cell,i_r *ncomp+n) = v;
                    a(cell,i_rh*ncomp+n) = v;
                    rnorm = std::max(rnorm, std::abs(v));
                }
            }
            rnorm = Gpu::blockReduceMax<MT>(rnorm);
            if (threadIdx.x == 0) {
                reduce[blockIdx.x] = rnorm;
            }
        });

        // Finish the rnorm reduction with one block and initialize solver state.
        amrex::launch<MT>(1, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
        {
            T rnorm = 0;
            for (auto iblock = threadIdx.x; iblock < nblocks; iblock += blockDim.x) {
                rnorm = std::max(rnorm, reduce[iblock]);
            }
            rnorm = Gpu::blockReduceMax<MT>(rnorm);
            if (threadIdx.x == 0) {
                scalar[s_rnorm0] = rnorm;
                scalar[s_rnorm] = rnorm;
                scalar[s_rho] = 0;
                scalar[s_rho_1] = 0;
                scalar[s_alpha] = 0;
                scalar[s_omega] = 0;
                scalar[s_beta] = 0;
                resd->ret = 0;
                *active = !((rnorm == 0) || (rnorm < eps_abs));
            }
        });

        int iter = 1;
        for (; iter <= maxiter; ++iter) {
            // rho = dot(rh,r), with one partial sum per block.
            amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
                T rho = 0;
                if (ipt < npts) {
                    auto const cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        rho += lp.xdoty(cell, n, a(cell,i_rh*ncomp+n),
                                       a(cell,i_r*ncomp+n));
                    }
                }
                rho = Gpu::blockReduceSum<MT>(rho);
                if (threadIdx.x == 0) {
                    reduce[blockIdx.x] = rho;
                }
            });

            amrex::launch<MT>(1, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                T rho = 0;
                for (auto iblock = threadIdx.x; iblock < nblocks; iblock += blockDim.x) {
                    rho += reduce[iblock];
                }
                rho = Gpu::blockReduceSum<MT>(rho);
                if (threadIdx.x == 0) {
                    scalar[s_rho] = rho;
                    if (rho == 0) {
                        resd->ret = 1;
                        *active = 0;
                    } else if (iter > 1) {
                        scalar[s_beta] = (rho/scalar[s_rho_1])
                            * (scalar[s_alpha]/scalar[s_omega]);
                    }
                }
            });

            // Update p.  This is separate from applying the operator so that
            // all p values are globally visible before lp.apply reads them.
            amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
                if (ipt < npts) {
                    auto const cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        if (iter == 1) {
                            a(cell,i_p*ncomp+n) = a(cell,i_r*ncomp+n);
                        } else {
                            a(cell,i_p*ncomp+n) = a(cell,i_r*ncomp+n)
                                + scalar[s_beta] * (a(cell,i_p*ncomp+n)
                                                   - scalar[s_omega]*a(cell,i_v*ncomp+n));
                        }
                    }
                }
            });

            // Keep applying the operator separate from the reduction. Fusing them
            // increases register pressure for complex operators such as EB with
            // variable coefficients.
            // v = A p.
            amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
                if (ipt < npts) {
                    auto const cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        a(cell,i_v*ncomp+n) = lp.apply(cell, n, a, i_p*ncomp+n);
                        lp.normalize(cell, n, a(cell,i_v*ncomp+n));
                    }
                }
            });

            // rhTv = dot(rh,v).
            amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
                T rhTv = 0;
                if (ipt < npts) {
                    auto const cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        rhTv += lp.xdoty(cell, n, a(cell,i_rh*ncomp+n),
                                        a(cell,i_v*ncomp+n));
                    }
                }
                rhTv = Gpu::blockReduceSum<MT>(rhTv);
                if (threadIdx.x == 0) {
                    reduce[blockIdx.x] = rhTv;
                }
            });

            amrex::launch<MT>(1, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                T rhTv = 0;
                for (auto iblock = threadIdx.x; iblock < nblocks; iblock += blockDim.x) {
                    rhTv += reduce[iblock];
                }
                rhTv = Gpu::blockReduceSum<MT>(rhTv);
                if (threadIdx.x == 0) {
                    if (rhTv != 0) {
                        scalar[s_alpha] = scalar[s_rho] / rhTv;
                    } else {
                        resd->ret = 2;
                        *active = 0;
                    }
                }
            });

            // x += alpha*p; r -= alpha*v; rnorm = norminf(r).
            amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
                T rnorm = 0;
                if (ipt < npts) {
                    auto const cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        sol(cell,n) += scalar[s_alpha]*a(cell,i_p*ncomp+n);
                        a(cell,i_r*ncomp+n) -= scalar[s_alpha]*a(cell,i_v*ncomp+n);
                        rnorm = std::max(rnorm, std::abs(a(cell,i_r*ncomp+n)));
                    }
                }
                rnorm = Gpu::blockReduceMax<MT>(rnorm);
                if (threadIdx.x == 0) {
                    reduce[blockIdx.x] = rnorm;
                }
            });

            amrex::launch<MT>(1, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                T rnorm = 0;
                for (auto iblock = threadIdx.x; iblock < nblocks; iblock += blockDim.x) {
                    rnorm = std::max(rnorm, reduce[iblock]);
                }
                rnorm = Gpu::blockReduceMax<MT>(rnorm);
                if (threadIdx.x == 0) {
                    scalar[s_rnorm] = rnorm;
                    if (rnorm < eps_rel*scalar[s_rnorm0] || rnorm < eps_abs) {
                        *active = 0;
                    }
                }
            });

            // Keep applying the operator separate from the reduction for the same
            // register-pressure reason as above.
            // t = A r.
            amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
                if (ipt < npts) {
                    auto const cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        a(cell,i_t*ncomp+n) = lp.apply(cell, n, a, i_r*ncomp+n);
                        lp.normalize(cell, n, a(cell,i_t*ncomp+n));
                    }
                }
            });

            // t0 = dot(t,t) and t1 = dot(t,r).
            amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
                T t0 = 0;
                T t1 = 0;
                if (ipt < npts) {
                    auto const cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        t0 += lp.xdoty(cell, n, a(cell,i_t*ncomp+n),
                                      a(cell,i_t*ncomp+n));
                        t1 += lp.xdoty(cell, n, a(cell,i_t*ncomp+n),
                                      a(cell,i_r*ncomp+n));
                    }
                }
                t0 = Gpu::blockReduceSum<MT>(t0);
                t1 = Gpu::blockReduceSum<MT>(t1);
                if (threadIdx.x == 0) {
                    reduce[blockIdx.x] = t0;
                    reduce[nblocks+blockIdx.x] = t1;
                }
            });

            amrex::launch<MT>(1, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                T t0 = 0;
                T t1 = 0;
                for (auto iblock = threadIdx.x; iblock < nblocks; iblock += blockDim.x) {
                    t0 += reduce[iblock];
                    t1 += reduce[nblocks+iblock];
                }
                t0 = Gpu::blockReduceSum<MT>(t0);
                t1 = Gpu::blockReduceSum<MT>(t1);
                if (threadIdx.x == 0) {
                    if (t0 != 0) {
                        scalar[s_omega] = t1/t0;
                    } else {
                        resd->ret = 3;
                        *active = 0;
                    }
                }
            });

            // x += omega*r; r -= omega*t; rnorm = norminf(r).
            amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
                T rnorm = 0;
                if (ipt < npts) {
                    auto const cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        sol(cell,n) += scalar[s_omega]*a(cell,i_r*ncomp+n);
                        a(cell,i_r*ncomp+n) -= scalar[s_omega]*a(cell,i_t*ncomp+n);
                        rnorm = std::max(rnorm, std::abs(a(cell,i_r*ncomp+n)));
                    }
                }
                rnorm = Gpu::blockReduceMax<MT>(rnorm);
                if (threadIdx.x == 0) {
                    reduce[blockIdx.x] = rnorm;
                }
            });

            amrex::launch<MT>(1, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
            {
                if (*active == 0) { return; }
                T rnorm = 0;
                for (auto iblock = threadIdx.x; iblock < nblocks; iblock += blockDim.x) {
                    rnorm = std::max(rnorm, reduce[iblock]);
                }
                rnorm = Gpu::blockReduceMax<MT>(rnorm);
                if (threadIdx.x == 0) {
                    scalar[s_rnorm] = rnorm;
                    if (rnorm < eps_rel*scalar[s_rnorm0] || rnorm < eps_abs) {
                        *active = 0;
                    } else if (scalar[s_omega] == 0) {
                        resd->ret = 4;
                        *active = 0;
                    } else {
                        scalar[s_rho_1] = scalar[s_rho];
                    }
                }
            });

            if (res_buf.copyToHost()->active == 0) {
                break;
            }
        }

        // Match the return and solution-reset behavior of the one-block path.
        amrex::launch<MT>(1, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
        {
            if (threadIdx.x == 0) {
                T const rnorm = scalar[s_rnorm];
                T const rnorm0 = scalar[s_rnorm0];
                int ret = resd->ret;
                bool keep_solution = (rnorm0 == 0 || rnorm0 < eps_abs);
                if (!keep_solution) {
                    if (*active != 0 && ret == 0 &&
                        rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
                    {
                        ret = 8;
                    }
                    keep_solution = ((ret == 0 || ret == 8) && rnorm < rnorm0);
                    if (keep_solution && ret == 8) { ret = 9; }
                }
                resd->ret = ret;
                *active = keep_solution;
                resd->rnorm0 = rnorm0;
                resd->rnorm = rnorm;
            }
        });

        amrex::launch<MT>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
        {
            if (*active == 0) {
                auto const ipt = std::uint64_t(blockIdx.x)*blockDim.x + threadIdx.x;
                if (ipt < npts) {
                    auto const cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        sol(cell,n) = 0;
                    }
                }
            }
        });

        // Like MLCGSolver, iter is maxiter+1 if the loop was exhausted.  The
        // no-iteration case is fixed up below once the norms are on the host.
        niters = iter;
    }
    else
    {
        amrex::launch<MT>(1, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
        {
            auto npts = static_cast<unsigned int>(indexer.numPts());

            T rnorm = 0;
            for (auto ipt = threadIdx.x; ipt < npts; ipt += blockDim.x) {
                auto cell = indexer.intVect(ipt);
                for (int n = 0; n < ncomp; ++n) {
                    T v = rhs(cell,n);
                    lp.normalize(cell, n, v);
                    a(cell,i_r *ncomp+n) = v;
                    a(cell,i_rh*ncomp+n) = v;
                    rnorm = std::max(rnorm, std::abs(v));
                }
            }

            rnorm = Gpu::blockReduceMax<MT>(rnorm);

            __shared__ T shared_scalar[2];

            if (threadIdx.x == 0) { shared_scalar[0] = rnorm; }
            __syncthreads();
            T const rnorm0 = shared_scalar[0];

            if (rnorm0 == 0 || rnorm0 < eps_abs) {
                if (threadIdx.x == 0) {
                    resd->ret = 0;
                    resd->niters = 0;
                    resd->rnorm0 = rnorm0;
                    resd->rnorm = rnorm0;
                }
                return;
            }

            int ret = 0;
            T rho_1 = 0, alpha = 0, omega = 0;
            int iter = 1;
            for (; iter <= maxiter; ++iter) {
                T rho = 0;
                for (auto ipt = threadIdx.x; ipt < npts; ipt += blockDim.x) {
                    auto cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        rho += lp.xdoty(cell, n, a(cell,i_rh*ncomp+n), a(cell,i_r*ncomp+n));
                    }
                }
                rho = Gpu::blockReduceSum<MT>(rho);

                if (threadIdx.x == 0) { shared_scalar[0] = rho; }
                __syncthreads();
                rho = shared_scalar[0];

                if (rho == 0) {
                    ret = 1;
                    break;
                }

                if (iter == 1) {
                    for (auto ipt = threadIdx.x; ipt < npts; ipt += blockDim.x) {
                        auto cell = indexer.intVect(ipt);
                        for (int n = 0; n < ncomp; ++n) {
                            a(cell,i_p*ncomp+n) = a(cell,i_r*ncomp+n);
                        }
                    }
                } else {
                    T const beta = (rho/rho_1)*(alpha/omega);
                    for (auto ipt = threadIdx.x; ipt < npts; ipt += blockDim.x) {
                        auto cell = indexer.intVect(ipt);
                        for (int n = 0; n < ncomp; ++n) {
                            a(cell,i_p*ncomp+n) = a(cell,i_r*ncomp+n)
                                + beta * (a(cell,i_p*ncomp+n) - omega*a(cell,i_v*ncomp+n));
                        }
                    }
                }

                __syncthreads(); // needed for lp.apply

                T rhTv = 0;
                for (auto ipt = threadIdx.x; ipt < npts; ipt += blockDim.x) {
                    auto cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        a(cell,i_v*ncomp+n) = lp.apply(cell, n, a, i_p*ncomp+n);
                        lp.normalize(cell, n, a(cell,i_v*ncomp+n));
                        rhTv += lp.xdoty(cell, n, a(cell,i_rh*ncomp+n), a(cell,i_v*ncomp+n));
                    }
                }
                rhTv = Gpu::blockReduceSum<MT>(rhTv);

                if (threadIdx.x == 0) { shared_scalar[0] = rhTv; }
                __syncthreads();
                rhTv = shared_scalar[0];

                if (rhTv != 0) {
                    alpha = rho / rhTv;
                } else {
                    ret = 2;
                    break;
                }

                rnorm = 0;
                for (auto ipt = threadIdx.x; ipt < npts; ipt += blockDim.x) {
                    auto cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        sol(cell,n) += alpha * a(cell,i_p*ncomp+n);
                        a(cell,i_r*ncomp+n) -= alpha * a(cell,i_v*ncomp+n);
                        rnorm = std::max(rnorm, std::abs(a(cell,i_r*ncomp+n)));
                    }
                }
                rnorm = Gpu::blockReduceMax<MT>(rnorm);

                if (threadIdx.x == 0) { shared_scalar[0] = rnorm; }
                __syncthreads();
                rnorm = shared_scalar[0];

                if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

                T t0 = 0, t1 = 0;
                for (auto ipt = threadIdx.x; ipt < npts; ipt += blockDim.x) {
                    auto cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        a(cell,i_t*ncomp+n) = lp.apply(cell, n, a, i_r*ncomp+n);
                        lp.normalize(cell, n, a(cell,i_t*ncomp+n));
                        t0 += lp.xdoty(cell, n, a(cell,i_t*ncomp+n), a(cell,i_t*ncomp+n));
                        t1 += lp.xdoty(cell, n, a(cell,i_t*ncomp+n), a(cell,i_r*ncomp+n));
                    }
                }
                t0 = Gpu::blockReduceSum<MT>(t0);
                t1 = Gpu::blockReduceSum<MT>(t1);

                if (threadIdx.x == 0) {
                    shared_scalar[0] = t0;
                    shared_scalar[1] = t1;
                }
                __syncthreads();
                t0 = shared_scalar[0];
                t1 = shared_scalar[1];

                if (t0 != 0) {
                    omega = t1 / t0;
                } else {
                    ret = 3;
                    break;
                }

                rnorm = 0;
                for (auto ipt = threadIdx.x; ipt < npts; ipt += blockDim.x) {
                    auto cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        sol(cell,n) += omega * a(cell,i_r*ncomp+n);
                        a(cell,i_r*ncomp+n) -= omega * a(cell,i_t*ncomp+n);
                        rnorm = std::max(rnorm, std::abs(a(cell,i_r*ncomp+n)));
                    }
                }
                rnorm = Gpu::blockReduceMax<MT>(rnorm);

                if (threadIdx.x == 0) { shared_scalar[0] = rnorm; }
                __syncthreads();
                rnorm = shared_scalar[0];

                if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

                if ( omega == 0 ) {
                    ret = 4;
                    break;
                }
                rho_1 = rho;
            }

            if ( ret == 0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs) {
                ret = 8;
            }

            if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) ) {
                if (ret == 8) { ret = 9; }
            } else {
                for (auto ipt = threadIdx.x; ipt < npts; ipt += blockDim.x) {
                    auto cell = indexer.intVect(ipt);
                    for (int n = 0; n < ncomp; ++n) {
                        sol(cell,n) = 0;
                    }
                }
            }

            if (threadIdx.x == 0) {
                resd->ret = ret;
                resd->niters = iter;
                resd->rnorm0 = rnorm0;
                resd->rnorm = rnorm;
            }
        });
    }

    auto const* resh = res_buf.copyToHost(); // this has stream sync.
    int const ret = resh->ret;
    T const rnorm0 = resh->rnorm0;
    T const rnorm = resh->rnorm;
    if (is_big) {
        if (rnorm0 == 0 || rnorm0 < eps_abs) { niters = 0; }
    } else {
        niters = resh->niters;
    }

    // Host-side output only, mirroring MLCGSolver.  Only the rank that owns
    // the box gets here, so AllPrint is used instead of Print.
    if (verbose > 0) {
        amrex::AllPrint() << "SingleBoxCGSolver_BiCGStab: Initial error (error0) =        "
                          << rnorm0 << '\n';
        if (niters == 0) {
            amrex::AllPrint() << "SingleBoxCGSolver_BiCGStab: niter = 0, rnorm = " << rnorm
                              << ", eps_abs = " << eps_abs << '\n';
        } else {
            amrex::AllPrint() << "SingleBoxCGSolver_BiCGStab: Final: Iteration "
                              << std::setw(4) << niters
                              << " rel. err. " << rnorm/rnorm0 << '\n';
            if (ret == 8 || ret == 9) {
                amrex::Warning("SingleBoxCGSolver_BiCGStab:: failed to converge!");
            }
        }
    }

    return ret;
}

}

#endif
#endif
