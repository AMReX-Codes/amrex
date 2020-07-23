#include <AMReX_FilCC_C.H>

namespace amrex {

void fab_filcc (Box const& bx, Array4<Real> const& qn, int ncomp,
                Box const& domain, Real const* /*dx*/, Real const* /*xlo*/,
                BCRec const* bcn)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const Box qbx(qn);
    const auto qlo = amrex::lbound(qbx);
    const auto qhi = amrex::ubound(qbx);
    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

    const int is = std::max(qlo.x, domlo.x);
    const int ie = std::min(qhi.x, domhi.x);
    const int ilo = domlo.x;
    const int ihi = domhi.x;

#if AMREX_SPACEDIM >= 2
    const int js = std::max(qlo.y, domlo.y);
    const int je = std::min(qhi.y, domhi.y);
    const int jlo = domlo.y;
    const int jhi = domhi.y;
#endif

#if AMREX_SPACEDIM == 3
    const int ks = std::max(qlo.z, domlo.z);
    const int ke = std::min(qhi.z, domhi.z);
    const int klo = domlo.z;
    const int khi = domhi.z;
#endif

    for (int n = 0; n < ncomp; ++n)
    {
        Array4<Real> q(qn,n);
        BCRec const& bc = bcn[n];

        if (lo.x < ilo) {
           const int imin = lo.x;
           const int imax = ilo-1;
           if (bc.lo(0) == BCType::ext_dir) {
               // Do nothing.
           } else if (bc.lo(0) == BCType::foextrap) {
               for (int k = lo.z; k <= hi.z; ++k) {
               for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = imin; i <= imax; ++i) {
                   q(i,j,k) = q(ilo,j,k);
               }}}
           } else if (bc.lo(0) == BCType::hoextrap) {
               for (int k = lo.z; k <= hi.z; ++k) {
               for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = imin; i <= imax; ++i) {
                   if (i < ilo - 1) {
                       q(i,j,k) = q(ilo,j,k);
                   } else if (i == ilo - 1) {
                       if (ilo+2 <= ie) {
                           q(i,j,k) = (1./8.) * (15.*q(ilo,j,k) - 10.*q(ilo+1,j,k) + 3.*q(ilo+2,j,k));
                       } else {
                           q(i,j,k) = 0.5 * (3.*q(ilo,j,k) - q(ilo+1,j,k));
                       }
                   }
               }}}
           } else if (bc.lo(0) == BCType::hoextrapcc) {
               for (int k = lo.z; k <= hi.z; ++k) {
               for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = imin; i <= imax; ++i) {
                   q(i,j,k) = 2.*q(ilo,j,k) - q(ilo+1,j,k);
               }}}
           } else if (bc.lo(0) == BCType::reflect_even) {
               for (int k = lo.z; k <= hi.z; ++k) {
               for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = imin; i <= imax; ++i) {
                   q(i,j,k) = q(ilo+(ilo-i)-1,j,k);
               }}}
           } else if (bc.lo(0) == BCType::reflect_odd) {
               for (int k = lo.z; k <= hi.z; ++k) {
               for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = imin; i <= imax; ++i) {
                   q(i,j,k) = -q(ilo+(ilo-i)-1,j,k);
               }}}
           }
        }

        if (hi.x > ihi) {
            const int imin = ihi+1;
            const int imax = hi.x;

            if (bc.hi(0) == BCType::ext_dir) {
                // Do nothing.
            } else if (bc.hi(0) == BCType::foextrap) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = imin; i <= imax; ++i) {
                    q(i,j,k) = q(ihi,j,k);
                }}}
            } else if (bc.hi(0) == BCType::hoextrap) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = imin; i <= imax; ++i) {
                    if (i > ihi + 1) {
                        q(i,j,k) = q(ihi,j,k);
                    } else if (i == ihi + 1) {
                        if (ihi-2 >= is) {
                            q(i,j,k) = (1./8.) * (15.*q(ihi,j,k) - 10.*q(ihi-1,j,k) + 3.*q(ihi-2,j,k));
                        } else {
                            q(i,j,k) = 0.5 * (3.*q(ihi,j,k) - q(ihi-1,j,k));
                        }
                    }
                }}}
            } else if (bc.hi(0) == BCType::hoextrapcc) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = imin; i <= imax; ++i) {
                    q(i,j,k) = 2.*q(ihi,j,k) - q(ihi-1,j,k);
                }}}
            } else if (bc.hi(0) == BCType::reflect_even) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = imin; i <= imax; ++i) {
                    q(i,j,k) = q(ihi-(i-ihi)+1,j,k);
                }}}
            } else if (bc.hi(0) == BCType::reflect_odd) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = imin; i <= imax; ++i) {
                    q(i,j,k) = -q(ihi-(i-ihi)+1,j,k);
                }}}
            }
        }

#if AMREX_SPACEDIM >= 2

        if (lo.y < jlo) {
            const int jmin = lo.y;
            const int jmax = jlo-1;
            if (bc.lo(1) == BCType::ext_dir) {
                // Do nothing.
            } else if (bc.lo(1) == BCType::foextrap) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = q(i,jlo,k);
                }}}
            } else if (bc.lo(1) == BCType::hoextrap) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (j < jlo - 1) {
                        q(i,j,k) = q(i,jlo,k);
                    } else if (j == jlo - 1) {
                        if (jlo+2 <= je) {
                            q(i,j,k) = (1./8.) * (15.*q(i,jlo,k) - 10.*q(i,jlo+1,k) + 3.*q(i,jlo+2,k));
                        } else {
                            q(i,j,k) = 0.5 * (3.*q(i,jlo,k) - q(i,jlo+1,k));
                        }
                    }
                }}}
            } else if (bc.lo(1) == BCType::hoextrapcc) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = 2.*q(i,jlo,k) - q(i,jlo+1,k);
                }}}
            } else if (bc.lo(1) == BCType::reflect_even) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = q(i,jlo+(jlo-j)-1,k);
                }}}
            } else if (bc.lo(1) == BCType::reflect_odd) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = -q(i,jlo+(jlo-j)-1,k);
                }}}
            }
        }

        if (hi.y > jhi) {
            const int jmin = jhi+1;
            const int jmax = hi.y;
            if (bc.hi(1) == BCType::ext_dir) {
                // Do nothing.
            } else if (bc.hi(1) == BCType::foextrap) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = q(i,jhi,k);
                }}}
            } else if (bc.hi(1) == BCType::hoextrap) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (j > jhi + 1) {
                        q(i,j,k) = q(i,jhi,k);
                    } else if (j == jhi + 1) {
                        if (jhi-2 >= js) {
                            q(i,j,k) = (1./8.) * (15.*q(i,jhi,k) - 10.*q(i,jhi-1,k) + 3.*q(i,jhi-2,k));
                        } else {
                            q(i,j,k) = 0.5 * (3.*q(i,jhi,k) - q(i,jhi-1,k));
                        }
                    }
                }}}
            } else if (bc.hi(1) == BCType::hoextrapcc) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = 2.*q(i,jhi,k) - q(i,jhi-1,k);
                }}}
            } else if (bc.hi(1) == BCType::reflect_even) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = q(i,jhi-(j-jhi)+1,k);
                }}}
            } else if (bc.hi(1) == BCType::reflect_odd) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = -q(i,jhi-(j-jhi)+1,k);
                }}}
            }
        }
#endif

#if AMREX_SPACEDIM == 3

        if (lo.z < klo) {
            const int kmin = lo.z;
            const int kmax = klo-1;
            if (bc.lo(2) == BCType::ext_dir) {
                // Do nothing.
            } else if (bc.lo(2) == BCType::foextrap) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = q(i,j,klo);
                }}}
            } else if (bc.lo(2) == BCType::hoextrap) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (k < klo - 1) {
                        q(i,j,k) = q(i,j,klo);
                    } else if (k == klo - 1) {
                        if (klo+2 <= ke) {
                            q(i,j,k) = (1./8.) * (15.*q(i,j,klo) - 10.*q(i,j,klo+1) + 3.*q(i,j,klo+2));
                        } else {
                            q(i,j,k) = 0.5 * (3.*q(i,j,klo) - q(i,j,klo+1));
                        }
                    }
                }}}
            } else if (bc.lo(2) == BCType::hoextrapcc) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = 2.*q(i,j,klo) - q(i,j,klo+1);
                }}}
            } else if (bc.lo(2) == BCType::reflect_even) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = q(i,j,klo+(klo-k)-1);
                }}}
            } else if (bc.lo(2) == BCType::reflect_odd) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = -q(i,j,klo+(klo-k)-1);
                }}}
            }
        }

        if (hi.z > khi) {
            const int kmin = khi+1;
            const int kmax = hi.z;
            if (bc.hi(2) == BCType::ext_dir) {
                // Do nothing.
            } else if (bc.hi(2) == BCType::foextrap) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = q(i,j,khi);
                }}}
            } else if (bc.hi(2) == BCType::hoextrap) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (k > khi + 1) {
                        q(i,j,k) = q(i,j,khi);
                    } else if (k == khi + 1) {
                        if (khi-2 >= ks) {
                            q(i,j,k) = (1./8.) * (15.*q(i,j,khi) - 10.*q(i,j,khi-1) + 3.*q(i,j,khi-2));
                        } else {
                            q(i,j,k) = 0.5 * (3.*q(i,j,khi) - q(i,j,khi-1));
                        }
                    }
                }}}
            } else if (bc.hi(2) == BCType::hoextrapcc) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = 2.*q(i,j,khi) - q(i,j,khi-1);
                }}}
            } else if (bc.hi(2) == BCType::reflect_even) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = q(i,j,khi-(k-khi)+1);
                }}}
            } else if (bc.hi(2) == BCType::reflect_odd) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    q(i,j,k) = -q(i,j,khi-(k-khi)+1);
                }}}
            }
        }
#endif

#if AMREX_SPACEDIM >= 2
        // Now take care of the higher contributions

        //
        // First correct the i-j edges and all corners
        //
        if (bc.lo(0) == BCType::hoextrap and bc.lo(1) == BCType::hoextrap) {
            if (lo.x < ilo and lo.y < jlo) {
                const int imin = lo.x;
                const int imax = std::min(hi.x,ilo-1);
                const int jmin = lo.y;
                const int jmax = std::min(hi.y,jlo-1);
                const int i = ilo-1;
                const int j = jlo-1;
                if (i>=imin and i<=imax and j>=jmin and j<=jmax) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        if (jlo+2 <= je) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(ilo-1,jlo,k) - 10.*q(ilo-1,jlo+1,k) + 3.*q(ilo-1,jlo+2,k));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(ilo-1,jlo,k) - q(ilo-1,jlo+1,k));
                        }

                        if (ilo+2 <= ie) {
                            q(i,j,k) = q(ilo-1,jlo-1,k) +
                                0.5 * (1./8.) * (15.*q(ilo,jlo-1,k) - 10.*q(ilo+1,jlo-1,k) + 3.*q(ilo+2,jlo-1,k));
                        } else {
                            q(i,j,k) = q(ilo-1,jlo-1,k) + 0.5 * 0.5 * (3.*q(ilo,jlo-1,k) - q(ilo+1,jlo-1,k));
                        }

#if AMREX_SPACEDIM == 3

                        if (k == klo-1 and bc.lo(2) == BCType::hoextrap) {
                            if (klo+2 <= ke) {
                                q(i,j,k) = (1./8.) * ( (15.*q(ilo-1,jlo-1,klo) - 10.*q(ilo-1,jlo-1,klo+1) +
                                                          3.*q(ilo-1,jlo-1,klo+2)) );
                            } else {
                                q(i,j,k) = 0.5 * (3.*q(ilo-1,jlo-1,klo) - q(ilo-1,jlo-1,klo+1));
                            }
                        }

                        if (k == khi+1 and bc.hi(2) == BCType::hoextrap) {
                            if (khi-2 >= ks) {
                                q(i,j,k) = (1./8.) * ( (15.*q(ilo-1,jlo-1,khi) - 10.*q(ilo-1,jlo-1,khi-1) +
                                                          3.*q(ilo-1,jlo-1,khi-2)) );
                            } else {
                                q(i,j,k) = 0.5 * (3.*q(ilo-1,jlo-1,khi) - q(ilo-1,jlo-1,khi-1));
                            }
                        }
#endif
                    }
                }
            }
        }

        //
        // *********************************************************************
        //

        if (bc.lo(0) == BCType::hoextrap and bc.hi(1) == BCType::hoextrap) {
            if (lo.x < ilo and hi.y > jhi) {
                const int imin = lo.x;
                const int imax = std::min(hi.x,ilo-1);
                const int jmin = std::max(lo.y,jhi+1);
                const int jmax = hi.y;
                const int i = ilo-1;
                const int j = jhi+1;
                if (i>=imin and i<=imax and j>=jmin and j<=jmax) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        if (jhi-2 >= js) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(ilo-1,jhi,k) - 10.*q(ilo-1,jhi-1,k) + 3.*q(ilo-1,jhi-2,k));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(ilo-1,jhi,k) - q(ilo-1,jhi-1,k));
                        }

                        if (ilo+2 <= ie) {
                            q(i,j,k) = q(ilo-1,jhi+1,k) +
                                0.5 * (1./8.) * (15.*q(ilo,jhi+1,k) - 10.*q(ilo+1,jhi+1,k) + 3.*q(ilo+2,jhi+1,k));
                        } else {
                            q(i,j,k) = q(ilo-1,jhi+1,k) + 0.5 * 0.5 * (3.*q(ilo,jhi+1,k) - q(ilo+1,jhi+1,k));
                   }

#if (AMREX_SPACEDIM == 3)
                        if (k == klo-1 and bc.lo(2) == BCType::hoextrap) {
                            if (klo+2 <= ke) {
                                q(i,j,k) = (1./8.) * ( (15.*q(ilo-1,jhi+1,klo) - 10.*q(ilo-1,jhi+1,klo+1) +
                                                          3.*q(ilo-1,jhi+1,klo+2)) );
                            } else {
                                q(i,j,k) = 0.5 * (3.*q(ilo-1,jhi+1,klo) - q(ilo-1,jhi+1,klo+1));
                            }
                        }

                        if (k == khi+1 and bc.hi(2) == BCType::hoextrap) {
                            if (khi-2 >= ks) {
                                q(i,j,k) = (1./8.) * ( (15.*q(ilo-1,jhi+1,khi) - 10.*q(ilo-1,jhi+1,khi-1) +
                                                          3.*q(ilo-1,jhi+1,khi-2)) );
                            } else {
                                q(i,j,k) = 0.5 * (3.*q(ilo-1,jhi+1,khi) - q(ilo-1,jhi+1,khi-1));
                            }
                        }
#endif
                    }
                }
            }
        }

        //
        // **********************************************************************
        //

        if (bc.hi(0) == BCType::hoextrap and bc.lo(1) == BCType::hoextrap) {
            if (hi.x > ihi and lo.y < jlo) {
                const int imin = std::max(lo.x,ihi+1);
                const int imax = hi.x;
                const int jmin = lo.y;
                const int jmax = std::min(hi.y,jlo-1);
                const int i = ihi+1;
                const int j = jlo-1;
                if (i>=imin and i<=imax and j>=jmin and j<=jmax) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        if (jlo+2 <= je) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(ihi+1,jlo,k) - 10.*q(ihi+1,jlo+1,k) + 3.*q(ihi+1,jlo+2,k));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(ihi+1,jlo,k) - q(ihi+1,jlo+1,k));
                        }

                        if (ihi-2 >= is) {
                            q(i,j,k) = q(ihi+1,jlo-1,k) +
                                0.5 * (1./8.) * (15.*q(ihi,jlo-1,k) - 10.*q(ihi-1,jlo-1,k) + 3.*q(ihi-2,jlo-1,k));
                        } else {
                            q(i,j,k) = q(ihi+1,jlo-1,k) + 0.5 * 0.5 * (3.*q(ihi,jlo-1,k) - q(ihi-1,jlo-1,k));
                        }

#if (AMREX_SPACEDIM == 3)
                        if (k == klo-1 and bc.lo(2) == BCType::hoextrap) {
                            if (klo+2 <= ke) {
                                q(i,j,k) = (1./8.) * (15.*q(ihi+1,jlo-1,klo) - 10.*q(ihi+1,jlo-1,klo+1) + 3.*q(ihi+1,jlo-1,klo+2));
                            } else {
                                q(i,j,k) = 0.5 * (3.*q(ihi+1,jlo-1,klo) - q(ihi+1,jlo-1,klo+1));
                            }
                        }

                        if (k == khi+1 and bc.hi(2) == BCType::hoextrap) {
                            if (khi-2 >= ks) {
                                q(i,j,k) = (1./8.) * (15.*q(ihi+1,jlo-1,khi) - 10.*q(ihi+1,jlo-1,khi-1) + 3.*q(ihi+1,jlo-1,khi-2));
                            } else {
                                q(i,j,k) = 0.5 * (3.*q(ihi+1,jlo-1,khi) - q(ihi+1,jlo-1,khi-1));
                            }
                        }
#endif
                    }
                }
            }
        }

        //
        // **********************************************************************
        //

        if (bc.hi(0) == BCType::hoextrap and bc.hi(1) == BCType::hoextrap) {
            if (hi.x > ihi and hi.y > jhi) {
                const int imin = std::max(lo.x,ihi+1);
                const int imax = hi.x;
                const int jmin = std::max(lo.y,jhi+1);
                const int jmax = hi.y;
                const int i = ihi+1;
                const int j = jhi+1;
                if (i>=imin and i<=imax and j>=jmin and j<=jmax) {
                    for (int k = lo.z; k <= hi.z; ++k) {
                        if (jhi-2 >= js) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(ihi+1,jhi,k) - 10.*q(ihi+1,jhi-1,k) + 3.*q(ihi+1,jhi-2,k));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(ihi+1,jhi,k) - q(ihi+1,jhi-1,k));
                        }

                        if (ihi-2 >= is) {
                            q(i,j,k) = q(ihi+1,jhi+1,k) +
                                0.5 * (1./8.) * (15.*q(ihi,jhi+1,k) - 10.*q(ihi-1,jhi+1,k) + 3.*q(ihi-2,jhi+1,k));
                        } else {
                            q(i,j,k) = q(ihi+1,jhi+1,k) + 0.5 * 0.5 * (3.*q(ihi,jhi+1,k) - q(ihi-1,jhi+1,k));
                        }

#if (AMREX_SPACEDIM == 3)
                        if (k == klo-1 and bc.lo(2) == BCType::hoextrap) {
                            if (klo+2 <= ke) {
                                q(i,j,k) = (1./8.) * (15.*q(ihi+1,jhi+1,klo) - 10.*q(ihi+1,jhi+1,klo+1) + 3.*q(ihi+1,jhi+1,klo+2));
                            } else {
                                q(i,j,k) = 0.5 * (3.*q(ihi+1,jhi+1,klo) - q(ihi+1,jhi+1,klo+1));
                            }
                        }

                        if (k == khi+1 and bc.hi(2) == BCType::hoextrap) {
                            if (khi-2 >= ks) {
                                q(i,j,k) = (1./8.) * (15.*q(ihi+1,jhi+1,khi) - 10.*q(ihi+1,jhi+1,khi-1) + 3.*q(ihi+1,jhi+1,khi-2));
                            } else {
                                q(i,j,k) = 0.5 * (3.*q(ihi+1,jhi+1,khi) - q(ihi+1,jhi+1,khi-1));
                            }
                        }
#endif
                    }
                }
            }
        }
#endif

#if AMREX_SPACEDIM == 3
        //
        // Next correct the i-k edges
        //

        if (bc.lo(0) == BCType::hoextrap and bc.lo(2) == BCType::hoextrap) {
            if (lo.x < ilo and lo.z < klo) {
                const int imin = lo.x;
                const int imax = std::min(hi.x,ilo-1);
                const int kmin = lo.z;
                const int kmax = std::min(hi.z,klo-1);
                const int i = ilo-1;
                const int k = klo-1;
                if (i>=imin and i<=imax and k>=kmin and k<=kmax) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        if (klo+2 <= ke) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(ilo-1,j,klo) - 10.*q(ilo-1,j,klo+1) + 3.*q(ilo-1,j,klo+2));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(ilo-1,j,klo) - q(ilo-1,j,klo+1));
                        }

                        if (ilo+2 <= ie) {
                            q(i,j,k) = q(ilo-1,j,klo-1) +
                                0.5 * (1./8.) * (15.*q(ilo,j,klo-1) - 10.*q(ilo+1,j,klo-1) + 3.*q(ilo+2,j,klo-1));
                        } else {
                            q(i,j,k) = q(ilo-1,j,klo-1) + 0.5 * 0.5 * (3.*q(ilo,j,klo-1) - q(ilo+1,j,klo-1));
                        }
                    }
                }
            }
        }

        //
        // **********************************************************************
        //

        if (bc.lo(0) == BCType::hoextrap and bc.hi(2) == BCType::hoextrap) {
            if (lo.x < ilo and hi.z > khi) {
                const int imin = lo.x;
                const int imax = std::min(hi.x,ilo-1);
                const int kmin = std::max(lo.z,khi+1);
                const int kmax = hi.z;
                const int i = ilo-1;
                const int k = khi+1;
                if (i>=imin and i<=imax and k>=kmin and k<=kmax) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        if (khi-2 >= ks) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(ilo-1,j,khi) - 10.*q(ilo-1,j,khi-1) + 3.*q(ilo-1,j,khi-2));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(ilo-1,j,khi) - q(ilo-1,j,khi-1));
                        }

                        if (ilo+2 <= ie) {
                            q(i,j,k) = q(ilo-1,j,khi+1) +
                                0.5 * (1./8.) * (15.*q(ilo,j,khi+1) - 10.*q(ilo+1,j,khi+1) + 3.*q(ilo+2,j,khi+1));
                        } else {
                            q(i,j,k) = q(ilo-1,j,khi+1) + 0.5 * 0.5 * (3.*q(ilo,j,khi+1) - q(ilo+1,j,khi+1));
                        }
                    }
                }
            }
        }

        //
        // **********************************************************************
        //

        if (bc.hi(0) == BCType::hoextrap and bc.lo(2) == BCType::hoextrap) {
            if (hi.x > ihi and lo.z < klo) {
                const int imin = std::max(lo.x,ihi+1);
                const int imax = hi.x;
                const int kmin = lo.z;
                const int kmax = std::min(hi.z,klo-1);
                const int i = ihi+1;
                const int k = klo-1;
                if (i>=imin and i<=imax and k>=kmin and k<=kmax) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        if (klo+2 <= ke) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(ihi+1,j,klo) - 10.*q(ihi+1,j,klo+1) + 3.*q(ihi+1,j,klo+2));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(ihi+1,j,klo) - q(ihi+1,j,klo+1));
                        }

                        if (ihi-2 >= is) {
                            q(i,j,k) = q(ihi+1,j,klo-1) +
                                0.5 * (1./8.) * (15.*q(ihi,j,klo-1) - 10.*q(ihi-1,j,klo-1) + 3.*q(ihi-2,j,klo-1));
                        } else {
                            q(i,j,k) = q(ihi+1,j,klo-1) + 0.5 * 0.5 * (3.*q(ihi,j,klo-1) - q(ihi-1,j,klo-1));
                        }
                    }
                }
            }
        }

        //
        // **********************************************************************
        //

        if (bc.hi(0) == BCType::hoextrap and bc.hi(2) == BCType::hoextrap) {
            if (hi.x > ihi and hi.z > khi) {
                const int imin = std::max(lo.x,ihi+1);
                const int imax = hi.x;
                const int kmin = std::max(lo.z,khi+1);
                const int kmax = hi.z;
                const int i = ihi+1;
                const int k = khi+1;
                if (i>=imin and i<=imax and k>=kmin and k<=kmax) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        if (khi-2 >= ks) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(ihi+1,j,khi) - 10.*q(ihi+1,j,khi-1) + 3.*q(ihi+1,j,khi-2));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(ihi+1,j,khi) - q(ihi+1,j,khi-1));
                        }

                        if (ihi-2 >= is) {
                            q(i,j,k) = q(ihi+1,j,khi+1) +
                                0.5 * (1./8.) * (15.*q(ihi,j,khi+1) - 10.*q(ihi-1,j,khi+1) + 3.*q(ihi-2,j,khi+1));
                        } else {
                            q(i,j,k) = q(ihi+1,j,khi+1) + 0.5 * 0.5 * (3.*q(ihi,j,khi+1) - q(ihi-1,j,khi+1));
                        }
                    }
                }
            }
        }

        //
        // Next correct the j-k edges
        //

        if (bc.lo(1) == BCType::hoextrap and bc.lo(2) == BCType::hoextrap) {
            if (lo.y < jlo and lo.z < klo) {
                const int jmin = lo.y;
                const int jmax = std::min(hi.y,jlo-1);
                const int kmin = lo.z;
                const int kmax = std::min(hi.z,klo-1);
                const int j = jlo-1;
                const int k = klo-1;
                if (j>=jmin and j<=jmax and k>=kmin and k<=kmax) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        if (klo+2 <= ke) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(i,jlo-1,klo) - 10.*q(i,jlo-1,klo+1) + 3.*q(i,jlo-1,klo+2));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(i,jlo-1,klo) - q(i,jlo-1,klo+1));
                        }

                        if (jlo+2 <= je) {
                            q(i,j,k) = q(i,jlo-1,klo-1) +
                                0.5 * (1./8.) * (15.*q(i,jlo,klo-1) - 10.*q(i,jlo+1,klo-1) + 3.*q(i,jlo+2,klo-1));
                        } else {
                            q(i,j,k) = q(i,jlo-1,klo-1) + 0.5 * 0.5 * (3.*q(i,jlo,klo-1) - q(i,jlo+1,klo-1));
                        }
                    }
                }
            }
        }

        //
        // **********************************************************************
        //

        if (bc.lo(1) == BCType::hoextrap and bc.hi(2) == BCType::hoextrap) {
            if (lo.y < jlo and hi.z > khi) {
                const int jmin = lo.y;
                const int jmax = std::min(hi.y,jlo-1);
                const int kmin = std::max(lo.z,khi+1);
                const int kmax = hi.z;
                const int j = jlo-1;
                const int k = khi+1;
                if (j>=jmin and j<=jmax and k>=kmin and k<=kmax) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        if (khi-2 >= ks) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(i,jlo-1,khi) - 10.*q(i,jlo-1,khi-1) + 3.*q(i,jlo-1,khi-2));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(i,jlo-1,khi) - q(i,jlo-1,khi-1));
                        }

                        if (jlo+2 <= je) {
                            q(i,j,k) = q(i,jlo-1,khi+1) +
                                0.5 * (1./8.) * (15.*q(i,jlo,khi+1) - 10.*q(i,jlo+1,khi+1) + 3.*q(i,jlo+2,khi+1));
                        } else {
                            q(i,j,k) = q(i,jlo-1,khi+1) + 0.5 * 0.5 * (3.*q(i,jlo,khi+1) - q(i,jlo+1,khi+1));
                        }
                    }
                }
            }
        }

        //
        // **********************************************************************
        //

        if (bc.hi(1) == BCType::hoextrap and bc.lo(2) == BCType::hoextrap) {
            if (hi.y > jhi and lo.z < klo) {
                const int jmin = std::max(lo.y,jhi+1);
                const int jmax = hi.y;
                const int kmin = lo.z;
                const int kmax = std::min(hi.z,klo-1);
                const int j = jhi+1;
                const int k = klo-1;
                if (j>=jmin and j<=jmax and k>=kmin and k<=kmax) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        if (klo+2 <= ke) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(i,jhi+1,klo) - 10.*q(i,jhi+1,klo+1) + 3.*q(i,jhi+1,klo+2));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(i,jhi+1,klo) - q(i,jhi+1,klo+1));
                        }

                        if (jhi-2 >= js) {
                            q(i,j,k) = q(i,jhi+1,klo-1) +
                                0.5 * (1./8.) * (15.*q(i,jhi,klo-1) - 10.*q(i,jhi-1,klo-1) + 3.*q(i,jhi-2,klo-1));
                        } else {
                            q(i,j,k) = q(i,jhi+1,klo-1) + 0.5 * 0.5 * (3.*q(i,jhi,klo-1) - q(i,jhi-1,klo-1));
                        }
                    }
                }
            }
        }

        //
        // **********************************************************************
        //

        if (bc.hi(1) == BCType::hoextrap and bc.hi(2) == BCType::hoextrap) {
            if (hi.y > jhi and hi.z > khi) {
                const int jmin = std::max(lo.y,jhi+1);
                const int jmax = hi.y;
                const int kmin = std::max(lo.z,khi+1);
                const int kmax = hi.z;
                const int j = jhi+1;
                const int k = khi+1;
                if (j>=jmin and j<=jmax and k>=kmin and k<=kmax) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        if (khi-2 >= ks) {
                            q(i,j,k) = 0.5 * (1./8.) * (15.*q(i,jhi+1,khi) - 10.*q(i,jhi+1,khi-1) + 3.*q(i,jhi+1,khi-2));
                        } else {
                            q(i,j,k) = 0.5 * 0.5 * (3.*q(i,jhi+1,khi) - q(i,jhi+1,khi-1));
                        }

                        if (jhi-2 >= js) {
                            q(i,j,k) = q(i,jhi+1,khi+1) +
                                0.5 * (1./8.) * (15.*q(i,jhi,khi+1) - 10.*q(i,jhi-1,khi+1) + 3.*q(i,jhi-2,khi+1));
                        } else {
                            q(i,j,k) = q(i,jhi+1,khi+1) + 0.5 * 0.5 * (3.*q(i,jhi,khi+1) - q(i,jhi-1,khi+1));
                        }
                    }
                }
            }
        }
#endif
    }
}

}
