#include <AMReX_GpuSequencer.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_BLassert.H>

#include <algorithm>

namespace amrex::Gpu {

void
Sequencer::fork (Fan& fan, ForkInfo const& info)
{
    AMREX_ASSERT(!fan.isOpen());
    auto& be = *m_backend;
    const int npool = std::min(be.numPoolStreams(), max_lanes);
    const int origin = be.current();
    const int width = std::max(1, std::min({info.nitems, info.max_width, npool}));

    fan.m_origin = origin;
    fan.m_width = width;
    // A fan that is not ordered runs on the pool streams and ends on stream
    // 0, unless it started on a stream that is not in the pool: that one
    // stays current until whoever set it removes it.
    fan.m_ordered = info.ordered || origin >= npool;
    // A thread that shares the loop joins on the host: the stream it would
    // join into is its own, not the one the work after the team is queued on.
    fan.m_join = (info.shared && info.join == StreamEdge::Gpu) ? StreamEdge::Host : info.join;

    if (width == 1 && fan.m_ordered) { return; } // The only lane is the origin.

    // Count the lanes before ordering them, so that a publish() on another
    // thread from now on includes them.
    if (width > 1) {
        for (int s = 0; s < width; ++s) { ++m_lane_users[s]; }
    }

    if (info.fork == StreamEdge::Host) {
        be.wait(origin);
    } else if (info.fork == StreamEdge::Gpu) {
        if (info.shared) {
            be.waitAll();
        } else {
            std::array<int, max_lanes> lanes{};
            for (int s = 0; s < width; ++s) { lanes[s] = s; }
            be.order(lanes.data(), width, origin);
        }
    }
}

void
Sequencer::lane (Fan const& fan, int i)
{
    AMREX_ASSERT(fan.isOpen());
    m_backend->setCurrent((fan.m_width > 1 || !fan.m_ordered) ? i % fan.m_width : fan.m_origin);
}

void
Sequencer::join (Fan& fan)
{
    if (!fan.isOpen()) { return; }
    auto& be = *m_backend;
    const int origin = fan.m_origin;
    const int width = std::exchange(fan.m_width, 0);

    if (width == 1 && fan.m_ordered) {
        if (fan.m_join == StreamEdge::Host) { be.wait(origin); }
    } else {
        for (int s = 0; s < width; ++s) {
            if (fan.m_join == StreamEdge::Host) {
                be.wait(s);
            } else if (fan.m_join == StreamEdge::Gpu && s != origin) {
                // Take the frees before the origin is ordered after the lane:
                // a free deferred to the lane later (by another thread) can
                // belong to work the origin does not wait for, so it stays
                // with the lane.
                auto frees = be.takeFrees(s);
                be.order(&origin, 1, s);
                be.giveFrees(origin, std::move(frees));
            }
        }
        if (width > 1) {
            for (int s = 0; s < width; ++s) { --m_lane_users[s]; }
        }
    }

    be.setCurrent(fan.m_ordered ? origin : 0);
}

void
Sequencer::publish (StreamEdge edge)
{
    auto& be = *m_backend;
    const int npool = std::min(be.numPoolStreams(), max_lanes);
    const int cur = be.current();
    std::array<int, max_lanes> lanes{};
    int n = 0;
    for (int s = 0; s < npool; ++s) {
        if (s != cur && m_lane_users[s].load(std::memory_order_relaxed) > 0) {
            lanes[n++] = s;
        }
    }
    if (n == 0 || edge == StreamEdge::None) { return; }
    if (edge == StreamEdge::Host) {
        be.wait(cur);
    } else {
        be.order(lanes.data(), n, cur);
    }
}

void
Sequencer::await ()
{
    m_backend->wait(m_backend->current());
}

int
Sequencer::openLanes () const noexcept
{
    for (int s = max_lanes; s > 0; --s) {
        if (m_lane_users[s-1].load(std::memory_order_relaxed) > 0) { return s; }
    }
    return 0;
}

void
await () noexcept
{
    Device::streamSynchronize();
}

void
publish () noexcept
{
#ifdef AMREX_USE_GPU
    // In a plain no-sync region the caller orders the streams, as always.
    theSequencer().publish(inStreamOrderedRegion() ? StreamEdge::Gpu
                           : (inNoSyncRegion() ? StreamEdge::None : StreamEdge::Host));
#endif
}

}
