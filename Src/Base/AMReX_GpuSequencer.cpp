#include <AMReX_GpuSequencer.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_BLassert.H>

#include <algorithm>

namespace amrex::Gpu {

namespace detail {

void
StreamSequencer::fork (LoopStreams& streams, ForkInfo const& info)
{
    AMREX_ASSERT(!streams.isOpen());
    auto& be = *m_backend;
    const int npool = std::min(be.numPoolStreams(), max_lanes);
    const int origin = be.current();
    const int width = std::max(1, std::min({info.nitems, info.max_width, npool}));

    streams.m_origin = origin;
    streams.m_width = width;
    // A loop that is not ordered runs on the pool streams and ends on stream
    // 0, unless it started on a stream that is not in the pool: that one
    // stays current until whoever set it removes it.
    streams.m_ordered = info.ordered || origin >= npool;
    // A thread that shares the loop joins on the host: the stream it would
    // join into is its own, not the one the work after the team is queued on.
    streams.m_join = (info.shared && info.join == StreamEdge::Gpu) ? StreamEdge::Host : info.join;

    if (width == 1 && streams.m_ordered) { return; } // The only lane is the origin.

    // Count the lanes before ordering them, so that an
    // orderLoopStreamsAfterCurrent() on another thread from now on includes
    // them.
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
StreamSequencer::lane (LoopStreams const& streams, int i)
{
    AMREX_ASSERT(streams.isOpen());
    m_backend->setCurrent((streams.m_width > 1 || !streams.m_ordered) ? i % streams.m_width : streams.m_origin);
}

void
StreamSequencer::join (LoopStreams& streams)
{
    if (!streams.isOpen()) { return; }
    auto& be = *m_backend;
    const int origin = streams.m_origin;
    const int width = std::exchange(streams.m_width, 0);

    if (width == 1 && streams.m_ordered) {
        if (streams.m_join == StreamEdge::Host) { be.wait(origin); }
    } else {
        for (int s = 0; s < width; ++s) {
            if (streams.m_join == StreamEdge::Host) {
                be.wait(s);
            } else if (streams.m_join == StreamEdge::Gpu && s != origin) {
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

    be.setCurrent(streams.m_ordered ? origin : 0);
}

void
StreamSequencer::orderLoopStreamsAfterCurrent (StreamEdge edge)
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
StreamSequencer::await ()
{
    m_backend->wait(m_backend->current());
}

int
StreamSequencer::openLanes () const noexcept
{
    for (int s = max_lanes; s > 0; --s) {
        if (m_lane_users[s-1].load(std::memory_order_relaxed) > 0) { return s; }
    }
    return 0;
}

} // namespace detail

void
orderLoopStreamsAfterCurrent () noexcept
{
#ifdef AMREX_USE_GPU
    // In a plain no-sync region the caller orders the streams, as always.
    using detail::StreamEdge;
    detail::theStreamSequencer().orderLoopStreamsAfterCurrent(
        inStreamOrderedRegion() ? StreamEdge::Gpu
        : (inNoSyncRegion() ? StreamEdge::None : StreamEdge::Host));
#endif
}

}
