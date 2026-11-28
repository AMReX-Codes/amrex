// Tests Gpu::detail::StreamSequencer on a mock backend, so it also runs in
// CPU builds.
//
// The mock records every queued operation and every order between streams,
// which makes the order the GPU would be allowed to use an explicit graph.
// Random programs of nested loops are run through the sequencer, and three
// properties are checked on the graph (no timing is involved):
//
//  order:    every operation is ordered after all operations that precede it
//            in program order, except those in another iteration of the same
//            loop. A published operation precedes everything queued later.
//  await:    when await() returns, every operation that precedes it in
//            program order is complete.
//  lifetime: memory whose free was deferred is released only when every
//            operation that uses it is complete.
//
// Backends that break one of the functions the sequencer relies on must be
// detected, which shows that the checks can fail.
//
// The properties are checked for loops in a Gpu::StreamOrderedRegion (any mix
// of host and GPU edges, nested), and for loops outside it that synchronize
// the host (not nested). Loops outside it that do not synchronize are not
// ordered, which the checks must detect too.

#include <AMReX.H>
#include <AMReX_GpuSequencer.H>
#include <AMReX_Print.H>

#include <deque>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace amrex;

namespace {

constexpr int npool = 4;
constexpr int nstreams = npool + 1; // one stream outside the pool

enum class Sabotage { None, DropOrder, DropWait, FreeOnTake };

//! Ordered: loops in a Gpu::StreamOrderedRegion. Default: loops outside it
//! that synchronize the host. NoSync: loops outside it that do not.
enum class Mode { Ordered, Default, NoSync };

struct Violation { std::string what; };

class MockBackend final
    : public Gpu::detail::StreamSequencerBackend
{
public:
    explicit MockBackend (Sabotage a_sabotage) : m_sabotage(a_sabotage) {}

    [[nodiscard]] int numPoolStreams () const override { return npool; }
    [[nodiscard]] int current () const override { return m_current; }
    void setCurrent (int stream) override { m_current = stream; }

    void order (int const* waiters, int n, int waited) override {
        if (m_sabotage == Sabotage::DropOrder) { return; }
        const int event = m_tail[waited];
        for (int i = 0; i < n; ++i) {
            if (waiters[i] != waited && event >= 0) { node(waiters[i], event); }
        }
    }

    void wait (int stream) override {
        if (m_sabotage == Sabotage::DropWait) { return; }
        m_barrier = node(stream);
        release(stream);
    }

    void waitAll () override {
        if (m_sabotage == Sabotage::DropWait) { return; }
        int n = node(0);
        for (int s = 1; s < nstreams; ++s) {
            if (m_tail[s] >= 0) { m_deps[n].push_back(m_tail[s]); }
        }
        m_barrier = n;
    }

    [[nodiscard]] FreeList takeFrees (int stream) override {
        if (m_sabotage == Sabotage::FreeOnTake) { release(stream); }
        return std::exchange(m_frees[stream], {});
    }

    void giveFrees (int stream, FreeList&& frees) override {
        const FreeList given(std::move(frees));
        m_frees[stream].insert(m_frees[stream].end(), given.begin(), given.end());
    }

    // What the test does to the device.

    //! Queues an operation on the current stream.
    int launch () { return node(m_current); }

    void deferFree (int mem) {
        // The address of the memory's entry stands for the memory.
        m_frees[m_current].emplace_back(nullptr, &m_ids[mem]);
    }

    int newMemory () {
        m_users.emplace_back();
        m_ids.push_back(static_cast<int>(m_users.size())-1);
        return m_ids.back();
    }
    void use (int mem, int op) { m_users[mem].push_back(op); }

    //! Operations that are complete when \p n is (including \p n).
    [[nodiscard]] std::vector<char> before (int n) const {
        std::vector<char> seen(m_deps.size(), 0);
        std::vector<int> stack;
        if (n >= 0) { stack.push_back(n); }
        while (!stack.empty()) {
            const int x = stack.back();
            stack.pop_back();
            if (seen[x]) { continue; }
            seen[x] = 1;
            for (int d : m_deps[x]) { stack.push_back(d); }
        }
        return seen;
    }

    //! Operations the host knows to be complete.
    [[nodiscard]] std::vector<char> complete () const { return before(m_barrier); }

    void finalize () {
        for (int s = 0; s < nstreams; ++s) {
            m_barrier = node(s);
            release(s);
        }
    }

private:
    int node (int stream, int extra = -1) {
        const int n = static_cast<int>(m_deps.size());
        m_deps.emplace_back();
        for (int d : {m_tail[stream], m_barrier, extra}) {
            if (d >= 0) { m_deps[n].push_back(d); }
        }
        m_tail[stream] = n;
        return n;
    }

    void release (int stream) {
        const auto done = (m_sabotage == Sabotage::FreeOnTake) ? complete()
                                                               : before(m_tail[stream]);
        for (auto const& [arena, p] : m_frees[stream]) {
            const int mem = *static_cast<int const*>(p);
            for (int op : m_users[mem]) {
                if (!done[op]) {
                    throw Violation{"lifetime: memory released while an operation that uses it is pending"};
                }
            }
        }
        m_frees[stream].clear();
    }

    Sabotage m_sabotage;
    int m_current = 0;
    int m_barrier = -1;
    std::vector<std::vector<int>> m_deps;
    std::vector<int> m_tail = std::vector<int>(nstreams, -1);
    std::vector<FreeList> m_frees = std::vector<FreeList>(nstreams);
    std::vector<std::vector<int>> m_users;
    std::deque<int> m_ids; //!< m_ids[mem] == mem; a deque keeps the addresses
};

//! A random program of nested loops, run through the sequencer.
class Program
{
public:
    Program (unsigned seed, Sabotage sabotage, Mode mode)
        : m_rng(seed), m_backend(sabotage), m_sequencer(m_backend), m_mode(mode) {}

    void run () {
        if (m_mode == Mode::Ordered && pick(4) == 0) {
            m_backend.setCurrent(npool); // start outside the pool
        }
        block(0, {});
        await();
        m_backend.finalize();
        if (m_sequencer.openLanes() != 0) { throw Violation{"lanes left open"}; }
    }

private:
    using Path = std::vector<std::pair<int,int>>; // (loop, iteration), outermost first

    struct Op { int node; Path path; bool published; };

    int pick (int n) { return std::uniform_int_distribution<int>(0, n-1)(m_rng); }

    //! Whether an earlier operation at \p path has to precede the current point.
    [[nodiscard]] bool precedes (Path const& path) const {
        for (std::size_t i = 0; i < path.size() && i < m_path.size(); ++i) {
            if (path[i].first != m_path[i].first) { return true; }
            if (path[i].second != m_path[i].second) { return false; }
        }
        return true;
    }

    void check (std::vector<char> const& done, const char* what) const {
        for (auto const& op : m_ops) {
            if ((op.published || precedes(op.path)) && !done[op.node]) {
                throw Violation{what};
            }
        }
    }

    void launch (std::vector<int> const& mems, bool published) {
        const int n = m_backend.launch();
        check(m_backend.before(n), "order: operation not ordered after an earlier one");
        m_ops.push_back({n, m_path, published});
        for (int m : mems) { m_backend.use(m, n); }
        if (published) {
            m_sequencer.orderLoopStreamsAfterCurrent(
                m_mode == Mode::Ordered ? (pick(2) ? Gpu::detail::StreamEdge::Gpu
                                                   : Gpu::detail::StreamEdge::Host)
                : (m_mode == Mode::Default ? Gpu::detail::StreamEdge::Host
                                           : Gpu::detail::StreamEdge::None));
        }
    }

    void await () {
        m_sequencer.await();
        check(m_backend.complete(), "await: host continues before an earlier operation is complete");
    }

    void loop (int depth, std::vector<int> const& mems) {
        const int niter = pick(7);
        Gpu::detail::LoopStreams streams;
        Gpu::detail::StreamSequencer::ForkInfo info{
            .nitems = niter,
            .max_width = std::vector<int>{1,2,4,8}[pick(4)],
            .fork = pick(3) ? Gpu::detail::StreamEdge::Gpu : Gpu::detail::StreamEdge::Host,
            .join = pick(3) ? Gpu::detail::StreamEdge::Gpu : Gpu::detail::StreamEdge::Host,
            .shared = (pick(8) == 0)};
        if (m_mode != Mode::Ordered) {
            // Like an MFIter outside a Gpu::StreamOrderedRegion: it waits on
            // the host before the loop only if it uses more than one stream.
            const bool multi = std::min(niter, info.max_width) > 1;
            info.ordered = false;
            info.shared = false;
            info.fork = (m_mode == Mode::Default && multi) ? Gpu::detail::StreamEdge::Host
                                                           : Gpu::detail::StreamEdge::None;
            info.join = (m_mode == Mode::Default) ? Gpu::detail::StreamEdge::Host
                                                  : Gpu::detail::StreamEdge::None;
        }
        const int id = ++m_nloops;
        const Path outer = m_path;
        const int after = (m_mode == Mode::Ordered || niter == 0) ? m_backend.current() : 0;
        for (int i = 0; i < niter; ++i) {
            // Like MFIter, the loop opens at its first iteration.
            if (!streams.isOpen()) { m_sequencer.fork(streams, info); }
            m_sequencer.lane(streams, i);
            m_path = outer;
            m_path.emplace_back(id, i);
            block(depth+1, mems);
        }
        m_path = outer;
        m_sequencer.join(streams);
        if (m_backend.current() != after) {
            throw Violation{"join: wrong current stream after the loop"};
        }
    }

    //! Random statements. Memory allocated in a block is freed in it.
    void block (int depth, std::vector<int> const& outer_mems) {
        std::vector<int> mine;
        const int nstatements = 1 + pick(5);
        for (int k = 0; k < nstatements; ++k) {
            std::vector<int> live = outer_mems;
            live.insert(live.end(), mine.begin(), mine.end());
            const int r = pick(100);
            if (r < 40) {
                std::vector<int> used;
                for (int m : live) { if (pick(3) == 0) { used.push_back(m); } }
                launch(used, false);
            } else if (r < 50) {
                mine.push_back(m_backend.newMemory());
            } else if (r < 58) {
                launch({}, true);
            } else if (r < 66) {
                await();
            } else if (r < 72 && !mine.empty()) {
                const int m = mine.back();
                mine.pop_back();
                launch({m}, false);
                m_backend.deferFree(m);
            } else if (depth < (m_mode == Mode::Ordered ? 3 : 1)) {
                loop(depth, live);
            }
        }
        for (int m : mine) { m_backend.deferFree(m); }
    }

    std::mt19937 m_rng;
    MockBackend m_backend;
    Gpu::detail::StreamSequencer m_sequencer;
    Mode m_mode;
    std::vector<Op> m_ops;
    Path m_path;
    int m_nloops = 0;
};

int count_violations (int nprograms, Sabotage sabotage, Mode mode = Mode::Ordered,
                      std::string* first = nullptr)
{
    int nviolations = 0;
    for (int seed = 0; seed < nprograms; ++seed) {
        try {
            Program(static_cast<unsigned>(seed), sabotage, mode).run();
        } catch (Violation const& v) {
            if (nviolations++ == 0 && first) {
                *first = "seed " + std::to_string(seed) + ": " + v.what;
            }
        }
    }
    return nviolations;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int nfail = 0;
    {
        constexpr int nprograms = 5000;
        const std::pair<Mode, const char*> kept[] = {
            {Mode::Ordered, "in a StreamOrderedRegion"},
            {Mode::Default, "outside it, loops synchronize"}};
        for (auto const& [mode, name] : kept) {
            std::string first;
            const int n = count_violations(nprograms, Sabotage::None, mode, &first);
            amrex::Print() << "StreamSequencer (" << name << "): " << n << " of " << nprograms
                           << " random programs violate an invariant\n";
            if (n > 0) {
                amrex::Print() << "  " << first << "\n";
                ++nfail;
            }
        }
        {
            // Not ordered by design: the checks must see that.
            const int n = count_violations(nprograms/10, Sabotage::None, Mode::NoSync);
            amrex::Print() << "StreamSequencer (outside it, loops do not synchronize): unordered in "
                           << n << " of " << nprograms/10 << " programs\n";
            if (n == 0) { ++nfail; }
        }

        // The checks must be able to fail.
        const std::pair<Sabotage, const char*> broken[] = {
            {Sabotage::DropOrder, "order() does nothing"},
            {Sabotage::DropWait, "wait() does nothing"},
            {Sabotage::FreeOnTake, "frees released when handed to the origin"}};
        for (auto const& [sabotage, name] : broken) {
            const int nb = count_violations(nprograms/10, sabotage);
            amrex::Print() << "Broken backend (" << name << "): detected in " << nb
                           << " of " << nprograms/10 << " programs\n";
            if (nb == 0) { ++nfail; }
        }
    }
    if (nfail == 0) {
        amrex::Print() << "PASS\n";
    } else {
        amrex::Print() << "FAIL\n";
    }
    amrex::Finalize();
    return nfail;
}
