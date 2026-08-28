#ifndef AMREX_SINGLE_BOX_CG_SOLVER_H_
#define AMREX_SINGLE_BOX_CG_SOLVER_H_
#include <AMReX_Config.H>

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

#include <AMReX_BaseFab.H>

namespace amrex
{

/**
 * \brief BiCGStab solver
 *
 * \param box region for solve
 * \param x solution with initial value already set to zero
 * \param b RHS
 * \param lp operator
 * \param eps_rel relative tolerance
 * \param eps_abs absolute tolerance
 * \param maxiter max iteration
 */
template <typename T, class LP>
int bicgstab_solve (Box const& box, BaseFab<T>& x, BaseFab<T> const& b, LP const& lp,
                    T eps_rel, T eps_abs, int maxiter)
{
    const int ncomp = x.nComp();

    enum Index {
        i_p, i_r, i_rh, i_v, i_t, nidx
    };

    BaseFab<T> fab(box, ncomp*nidx);
    auto const& a = fab.array();
    auto const& sol = x.array();
    auto const& rhs = const_cast<BaseFab<T>&>(b).array();

    Gpu::Buffer<int> ret_buf(2);
    auto* retd = ret_buf.data();

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
        auto* active = retd+1;
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
                retd[0] = 0;
                *active = !((rnorm == 0) || (rnorm < eps_abs));
            }
        });

        for (int iter = 1; iter <= maxiter; ++iter) {
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
                        retd[0] = 1;
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
                        retd[0] = 2;
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
                        retd[0] = 3;
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
                        retd[0] = 4;
                        *active = 0;
                    } else {
                        scalar[s_rho_1] = scalar[s_rho];
                    }
                }
            });

            if (ret_buf.copyToHost()[1] == 0) {
                break;
            }
        }

        // Match the return and solution-reset behavior of the one-block path.
        amrex::launch<MT>(1, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE ()
        {
            if (threadIdx.x == 0) {
                T const rnorm = scalar[s_rnorm];
                T const rnorm0 = scalar[s_rnorm0];
                int ret = retd[0];
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
                retd[0] = ret;
                *active = keep_solution;
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
                if (threadIdx.x == 0) { *retd = 0; }
                return;
            }

            int ret = 0;
            T rho_1 = 0, alpha = 0, omega = 0;
            for (int iter = 1; iter <= maxiter; ++iter) {
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

            if (threadIdx.x == 0) { *retd = ret; }
        });
    }

    auto* reth = ret_buf.copyToHost(); // this has stream sync.
    return *reth;
}

}

#endif
#endif
