#include <WarpX.H>
#include <PhysicalParticleContainer.H>
#include <AMReX_Particles.H>

using namespace amrex;

/* \brief
*/
void
PhysicalParticleContainer::PartitionParticlesInBuffers(
    long& nfine_current, long& nfine_gather, long const np,
    WarpXParIter& pti, int const lev,
    iMultiFab const* current_masks,
    iMultiFab const* gather_masks,
    RealVector& uxp, RealVector& uyp, RealVector& uzp, RealVector& wp)
{
    BL_PROFILE("PPC::Evolve::partition");

    std::vector<bool> inexflag;
    inexflag.resize(np);

    auto& aos = pti.GetArrayOfStructs();

    // We need to partition the large buffer first
    iMultiFab const* bmasks =
        (WarpX::n_field_gather_buffer >= WarpX::n_current_deposition_buffer) ?
        gather_masks : current_masks;
    int i = 0;
    const auto& msk = (*bmasks)[pti];
    for (const auto& p : aos) {
        const IntVect& iv = Index(p, lev);
        inexflag[i++] = msk(iv);
    }

    Vector<long> pid;
    pid.resize(np);
    std::iota(pid.begin(), pid.end(), 0L);
    auto sep = std::stable_partition(pid.begin(), pid.end(),
                                     [&inexflag](long id) { return inexflag[id]; });

    if (WarpX::n_current_deposition_buffer == WarpX::n_field_gather_buffer) {
        nfine_current = nfine_gather = std::distance(pid.begin(), sep);
    } else if (sep != pid.end()) {
        int n_buf;
        if (bmasks == gather_masks) {
            nfine_gather = std::distance(pid.begin(), sep);
            bmasks = current_masks;
            n_buf = WarpX::n_current_deposition_buffer;
        } else {
            nfine_current = std::distance(pid.begin(), sep);
            bmasks = gather_masks;
            n_buf = WarpX::n_field_gather_buffer;
        }
        if (n_buf > 0)
        {
            const auto& msk2 = (*bmasks)[pti];
            for (auto it = sep; it != pid.end(); ++it) {
                const long id = *it;
                const IntVect& iv = Index(aos[id], lev);
                inexflag[id] = msk2(iv);
            }

            auto sep2 = std::stable_partition(sep, pid.end(),
                                              [&inexflag](long id) {return inexflag[id];});
            if (bmasks == gather_masks) {
                nfine_gather = std::distance(pid.begin(), sep2);
            } else {
                nfine_current = std::distance(pid.begin(), sep2);
            }
        }
    }

    // only deposit / gather to coarsest grid
    if (m_deposit_on_main_grid && lev > 0) {
        nfine_current = 0;
    }
    if (m_gather_from_main_grid && lev > 0) {
        nfine_gather = 0;
    }

    if (nfine_current != np || nfine_gather != np)
    {
        RealVector tmp;
        ParticleVector particle_tmp;

        particle_tmp.resize(np);
        for (long ip = 0; ip < np; ++ip) {
            particle_tmp[ip] = aos[pid[ip]];
        }
        std::swap(aos(), particle_tmp);

        tmp.resize(np);

        for (long ip = 0; ip < np; ++ip) {
            tmp[ip] = wp[pid[ip]];
        }
        std::swap(wp, tmp);

        for (long ip = 0; ip < np; ++ip) {
            tmp[ip] = uxp[pid[ip]];
        }
        std::swap(uxp, tmp);

        for (long ip = 0; ip < np; ++ip) {
            tmp[ip] = uyp[pid[ip]];
        }
        std::swap(uyp, tmp);

        for (long ip = 0; ip < np; ++ip) {
            tmp[ip] = uzp[pid[ip]];
        }
        std::swap(uzp, tmp);
    }

}
