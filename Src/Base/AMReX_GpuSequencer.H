#ifndef AMREX_GPU_SEQUENCER_H_
#define AMREX_GPU_SEQUENCER_H_
#include <AMReX_Config.H>

#include <AMReX_GpuControl.H>
#include <AMReX_Vector.H>

#include <array>
#include <atomic>
#include <utility>

namespace amrex {
    class Arena;
}

namespace amrex::Gpu {

// The machinery that keeps GPU loops in program order. Applications use
// Gpu::StreamOrderedRegion, Gpu::orderLoopStreamsAfterCurrent and
// Gpu::streamSynchronize instead.
namespace detail {

/**
 * \brief What Gpu::detail::StreamSequencer needs from the GPU runtime.
 *
 * Streams are numbered from 0. Streams 0 .. numPoolStreams()-1 are the ones
 * loops may run their iterations on; larger numbers are other streams work
 * can be queued on (e.g., a stream set with Gpu::setExternalGpuStream).
 */
class StreamSequencerBackend
{
public:
    using FreeList = Vector<std::pair<Arena*, void*>>;

    StreamSequencerBackend () = default;
    virtual ~StreamSequencerBackend () = default;
    StreamSequencerBackend (StreamSequencerBackend const&) = delete;
    StreamSequencerBackend (StreamSequencerBackend&&) = delete;
    StreamSequencerBackend& operator= (StreamSequencerBackend const&) = delete;
    StreamSequencerBackend& operator= (StreamSequencerBackend&&) = delete;

    //! Number of streams loops may run their iterations on.
    [[nodiscard]] virtual int numPoolStreams () const = 0;

    //! The stream the calling thread queues its work on.
    [[nodiscard]] virtual int current () const = 0;

    //! Makes \p stream the calling thread's current stream.
    virtual void setCurrent (int stream) = 0;

    /**
     * Makes the work queued after this call on each of the \p n streams
     * \p waiters run after all work queued so far on \p waited. The host does
     * not have to wait. Streams equal to \p waited are skipped.
     */
    virtual void order (int const* waiters, int n, int waited) = 0;

    //! The host waits for all work queued on \p stream. Memory whose free
    //! was deferred to the stream is then released.
    virtual void wait (int stream) = 0;

    //! The host waits for all work queued on the device.
    virtual void waitAll () = 0;

    //! Removes and returns the frees deferred to \p stream.
    [[nodiscard]] virtual FreeList takeFrees (int stream) = 0;

    //! Defers \p frees to \p stream.
    virtual void giveFrees (int stream, FreeList&& frees) = 0;
};

/**
 * \brief How the work on two streams is put in order.
 */
enum class StreamEdge {
    Host, //!< The host waits for the earlier work.
    Gpu,  //!< The later stream waits on the GPU; the host does not wait.
    None  //!< The work is not put in order. The caller is responsible for it.
};

/**
 * \brief The streams of a loop whose iterations may run on different
 * streams.
 *
 * StreamSequencer::fork opens them and StreamSequencer::join closes them.
 * While they are open, StreamSequencer::lane selects the stream of an
 * iteration. A LoopStreams that is moved from is no longer open.
 */
class LoopStreams
{
public:
    LoopStreams () = default;
    ~LoopStreams () = default;
    LoopStreams (LoopStreams const&) = delete;
    LoopStreams& operator= (LoopStreams const&) = delete;
    LoopStreams& operator= (LoopStreams&&) = delete;
    LoopStreams (LoopStreams&& rhs) noexcept
        : m_origin(rhs.m_origin), m_width(std::exchange(rhs.m_width, 0)),
          m_join(rhs.m_join), m_ordered(rhs.m_ordered)
    {}

    //! Whether the fan is open (forked and not yet joined).
    [[nodiscard]] bool isOpen () const noexcept { return m_width > 0; }

    //! Number of streams the iterations run on (0 if the fan is not open).
    [[nodiscard]] int width () const noexcept { return m_width; }

    //! The stream that was current when the fan was opened.
    [[nodiscard]] int origin () const noexcept { return m_origin; }

private:
    friend class StreamSequencer;
    int m_origin = 0;
    int m_width = 0;
    StreamEdge m_join = StreamEdge::Host;
    bool m_ordered = true;
};

/**
 * \brief Keeps GPU work on several streams in program order.
 *
 * Work queued on one stream runs in the order it was queued. Work on
 * different streams has no order unless something establishes one. The
 * sequencer is the only place that does so, and the only place that changes
 * the current stream of a loop. It maintains one invariant:
 *
 *   When the work queued so far on a thread's current stream has completed,
 *   all GPU work that precedes this point in program order has completed,
 *   except for the other iterations of the loops that are still open.
 *
 * Therefore waiting for the current stream (Gpu::streamSynchronize()) is
 * always enough before the host reads results, whatever the number of
 * streams and whether or not loops synchronize the host.
 *
 * The transitions are
 *  - fork:    opens the LoopStreams of a loop. Its lanes are ordered after
 *             the origin stream.
 *  - lane:    makes the lane of iteration i the current stream.
 *  - join:    closes the LoopStreams of a loop. The origin stream is ordered
 *             after its lanes and is the current stream again.
 *  - orderLoopStreamsAfterCurrent:
 *             orders the lanes of all open loops after the current stream,
 *             for work that other iterations depend on.
 *  - await:   the host waits for the current stream.
 *
 * Whether fork and join order the streams on the host (the host waits) or on
 * the GPU (the later stream waits) is a property of the loop. It changes the
 * cost, not the order.
 *
 * The invariant holds as long as every fork and join uses StreamEdge::Host
 * or StreamEdge::Gpu. A loop outside a Gpu::StreamOrderedRegion that does not
 * synchronize the host (MFItInfo::DisableDeviceSync, Gpu::NoSyncRegion) uses
 * StreamEdge::None, as it always has: its streams are not ordered, and the
 * caller is responsible for the order of its work.
 *
 * All functions are thread-safe. Every thread has its own current stream,
 * which is kept by the backend.
 */
class StreamSequencer
{
public:
    static constexpr int max_lanes = 64;

    explicit StreamSequencer (StreamSequencerBackend& backend) : m_backend(&backend) {}

    /**
     * \brief Options of a fork.
     */
    struct ForkInfo
    {
        int nitems = 0;     //!< number of iterations
        int max_width = 1;  //!< largest number of streams the loop may use
        StreamEdge fork = StreamEdge::Host; //!< how lanes follow the origin
        StreamEdge join = StreamEdge::Host; //!< how the origin follows the lanes
        /**
         * The calling thread shares the loop with other threads that each
         * fork and join for themselves (an OpenMP team). Such a thread's
         * current stream need not be the one the work before the team was
         * queued on, so a StreamEdge::Gpu is replaced by waiting on the host:
         * for all work on the device at the fork, and for the lanes at the
         * join.
         */
        bool shared = false;
        /**
         * Whether the loop is in a Gpu::StreamOrderedRegion. An ordered fan
         * with one lane stays on the origin, and join makes the origin the
         * current stream again. A fan that is not ordered always uses pool
         * streams 0 .. width-1 as its lanes (unless the origin is not a pool
         * stream), and join makes stream 0 the current stream.
         */
        bool ordered = true;
    };

    /**
     * \brief Opens \p fan from the current stream (the origin).
     *
     * The fan has min(nitems, max_width, numPoolStreams()) lanes, at least
     * one. With one lane, the lane of an ordered fan is the origin itself
     * and nothing has to be ordered. Otherwise the lanes are the pool streams
     * 0 .. width-1 and are ordered after the work queued so far on the origin
     * (unless the fork edge is StreamEdge::None). The current stream does not
     * change until lane() is called.
     */
    void fork (LoopStreams& fan, ForkInfo const& info);

    //! Makes the lane of iteration \p i of the open \p fan the current stream.
    void lane (LoopStreams const& fan, int i);

    /**
     * \brief Closes \p fan and makes its origin the current stream (stream
     * 0 for a fan that is not ordered).
     *
     * The origin is ordered after the lanes. With StreamEdge::Host, the host
     * waits for every lane (also for a fan with one lane). With
     * StreamEdge::Gpu, frees deferred to the lanes are handed to the origin.
     * With StreamEdge::None, nothing is ordered. Does nothing if the fan is
     * not open.
     */
    void join (LoopStreams& fan);

    /**
     * \brief Orders the lanes of all open fans (of all threads) after the
     * work queued so far on the current stream.
     *
     * For work queued in one iteration that other iterations depend on, such
     * as filling a cache that is shared by the iterations. Does nothing if no
     * open fan has more than one lane, or if \p edge is StreamEdge::None.
     */
    void orderLoopStreamsAfterCurrent (StreamEdge edge);

    //! The host waits for the current stream.
    void await ();

    //! Largest number of lanes of the open fans (of all threads) that have
    //! more than one lane, or 0 if there is none.
    [[nodiscard]] int openLanes () const noexcept;

    [[nodiscard]] StreamSequencerBackend& backend () noexcept { return *m_backend; }

private:
    StreamSequencerBackend* m_backend;
    //! m_lane_users[s]: number of open fans with more than one lane that
    //! have pool stream s as a lane.
    std::array<std::atomic<int>, max_lanes> m_lane_users{};
};

#ifdef AMREX_USE_GPU

//! The sequencer of AMReX's GPU streams.
[[nodiscard]] StreamSequencer& theStreamSequencer () noexcept;

#endif

} // namespace detail

/**
 * \brief Makes the GPU work queued so far on the current stream visible to
 * all iterations of the enclosing loops.
 *
 * Iterations of an MFIter or Gpu::StreamIter loop may run on different
 * streams, without an order between them. Call this after queuing work in an
 * iteration that the other iterations depend on (e.g., after filling a
 * cache they share). Outside a no-sync region the host waits for the work.
 * In a Gpu::StreamOrderedRegion the other streams wait for it on the GPU.
 * In a plain Gpu::NoSyncRegion nothing is ordered. Does nothing outside such
 * loops.
 */
void orderLoopStreamsAfterCurrent () noexcept;

}

#endif
