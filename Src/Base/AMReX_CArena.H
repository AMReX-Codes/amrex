#ifndef BL_CARENA_H
#define BL_CARENA_H
#include <AMReX_Config.H>

/**
 * \file AMReX_CArena.H
 * \brief Coalescing first-fit dynamic memory arena.
 */

#include <AMReX_Arena.H>

#include <atomic>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

namespace amrex {

/**
 * \class amrex::CArena
 * \ingroup amrex_memory
 * \brief Coalescing first-fit dynamic memory arena.
 *
 * Allocates large hunks from the heap and sub-allocates from them on
 * demand. Neighboring free blocks are merged on each \c free() call
 * to reduce fragmentation.
 */
class CArena
    :
    public Arena
{
public:
    /**
     * \brief Construct a coalescing memory manager.
     * \param hunk_size Minimum size of heap hunks to allocate.
     *        Defaults to \c DefaultHunkSize when 0.
     * \param info      Arena configuration flags.
     */
    CArena (std::size_t hunk_size = 0, ArenaInfo info = ArenaInfo());

    CArena (const CArena& rhs) = delete;
    CArena (CArena&& rhs) = delete;
    CArena& operator= (const CArena& rhs) = delete;
    CArena& operator= (CArena&& rhs) = delete;

    //! The destructor.
    ~CArena () override;

    //! Allocate some memory.
    [[nodiscard]] void* alloc (std::size_t nbytes) final;

    /**
     * \brief Allocate memory in-place if possible
     *
     * This function tries to allocate in-place by extending the capacity of
     * an existing pointer. If it's possible to extend the given pointer to
     * the minimum required size without performing a new allocation, this
     * function will return the original pointer back along with the new
     * size (which may be in the range of `[szmin,szmax]`); otherwise a new
     * allocation is performed to allocate the maximum requested size.
     *
     * Note that this function does NOT free the existing allocation, even
     * when a new memory allocation is performed. So it's the caller's
     * responsibility to free the original pointer when appropriate. A new
     * allocation can be detected by comparing the returned pointer with the
     * input pointer.
     *
     * \param pt existing pointer
     * \param szmin minimum required size in bytes
     * \param szmax maximum requested size in bytes
     * \return pair of pointer and allocation size in bytes
     */
    [[nodiscard]] std::pair<void*,std::size_t>
    alloc_in_place (void* pt, std::size_t szmin, std::size_t szmax) final;

    /**
     * \brief Shrink allocation size in-place
     *
     * This function shrinks the existing allocation of the given pointer to
     * a new size. It's a runtime error, if the requested new size is larger
     * than the original size. In the special case of pt == nullptr or
     * new_size == 0, this does not do anything other than returning nullptr.
     *
     * \param pt existing pointer
     * \param new_size new size after shrinking
     * \return either the original pointer or `nullptr` if the requested new
     *         size is zero.
     */
    [[nodiscard]] void*
    shrink_in_place (void* pt, std::size_t new_size) final;

    /**
    * \brief Free up allocated memory.  Merge neighboring free memory chunks
    * into largest possible chunk.
    *
    * \param vp Pointer previously returned by alloc() or alloc_in_place().
    */
    void free (void* vp) final;

    /**
    * \brief Free unused memory back to the system.  Return value is the
    * amount of memory freed.
    */
    std::size_t freeUnused () final;

    /**
     * \brief Return the amount of memory that can be released back to the heap.
     */
    [[nodiscard]] std::size_t freeableMemory () const;

    //! Return the largest free memory block already held by this CArena.
    [[nodiscard]] std::size_t largestFreeBlock () const final;

    /**
     * \brief Does the device have enough free memory for allocating this
     * much memory?  For CPU builds, this always return true.  This is not a
     * const function because it may attempt to release memory back to the
     * system.
     *
     * \param sz Number of bytes requested.
     */
    [[nodiscard]] bool hasFreeDeviceMemory (std::size_t sz) final;

    //! The current amount of heap space used by the CArena object. Thread-safe.
    std::size_t heap_space_used () const noexcept;

    //! Return the total amount of memory given out via alloc. Thread-safe.
    std::size_t heap_space_actually_used () const noexcept;

    /**
     * \brief Return the size of the allocation referenced by \p p (0 if unknown).
     *
     * \param p Pointer previously returned by alloc().
     */
    std::size_t sizeOf (void* p) const noexcept;

    /**
     * \brief Print memory usage information of this arena
     *
     * This function reports memory usage statistics for this `CArena`.  It
     * contains MPI collective calls and therefore must be called by all
     * processes. The output is printed on the I/O process only.
     *
     * Depending on the value of `print_max_usage`, the function prints
     * either the current memory usage (when false) or the maximum memory
     * usage recoreded (when true).
     *
     * \param name name of the arena used in printing.
     * \param print_max_usage whether to print the maximum usage or the
     *        current usage.
     */
    void PrintUsage (std::string const& name, bool print_max_usage) const;

    /**
     * \brief Print memory usage information of this arena to a given output
     * stream
     *
     * This function reports memory usage statistics this `CArena`. Unlike
     * `PrintUsage`, this does not contain MPI calls and therefore may be
     * called independently on individual processes.
     *
     * \param os output stream
     * \param name name of the arena used in printing.
     * \param space indentation prepended to printed lines
     */
    void PrintUsage (std::ostream& os, std::string const& name, std::string const& space) const;

    /// Reset the maximum usage counter.
    void ResetMaxUsageCounter () final {
        m_max_actually_used = 0;
    }

    //! The default memory hunk size to grab from the heap.
    constexpr static std::size_t DefaultHunkSize = 1024*1024*8;

protected:

    void* alloc_protected (std::size_t nbytes);

    std::size_t freeUnused_protected () final;

    //! The nodes in our free list and block list.
    class Node
    {
    public:
        Node (void* a_block, void* a_owner, std::size_t a_size, MemStat* a_stat=nullptr) noexcept
            :
            m_block(a_block), m_owner(a_owner), m_size(a_size), m_stat(a_stat) {}

        Node (const Node& rhs) = default;
        Node (Node&& rhs) = default;
        Node& operator= (const Node& rhs) = default;
        Node& operator= (Node&& rhs) = default;

        //! The "less-than" operator.
        bool operator< (const Node& rhs) const noexcept
        {
            return std::less<>{}(m_block, rhs.m_block);
        }

        //! The equality operator.
        bool operator== (const Node& rhs) const noexcept
        {
            return m_block == rhs.m_block;
        }

        //! The block address.
        [[nodiscard]] void* block () const noexcept { return m_block; }

        //! Set block address.
        void block (void* blk) noexcept { m_block = blk; }

        //! The size of the memory block.
        [[nodiscard]] std::size_t size () const noexcept { return m_size; }

        //! Set size.
        void size (std::size_t sz) noexcept { m_size = sz; }

        [[nodiscard]] void* owner () const noexcept { return m_owner; }

        [[nodiscard]] bool coalescable (const Node& rhs) const noexcept {
            return m_owner == rhs.m_owner;
        }

        //! Get the MemStat object of the function where this block was allocated.
        [[nodiscard]] MemStat* mem_stat () const { return m_stat; }

        //! Set MemStat
        void mem_stat (MemStat* a_stat) noexcept { m_stat = a_stat; }

        struct hash {
            std::size_t operator() (const Node& n) const noexcept {
                return std::hash<void*>{}(n.m_block);
            }
        };

    private:
        //! The block of memory we reference.
        void* m_block;
        //! The starting address of the original allocation
        void* m_owner;
        //! The size of the block we represent.
        std::size_t m_size;
        //! Used for profiling if this Node represents a user allocated block of memory.
        MemStat* m_stat;
    };


    class NodeTreeWithSize {
    public:

        NodeTreeWithSize () {
            m_nodes.reserve(1024);
        }

        int getFirstFit (std::size_t size) const {
            int curr = m_root;
            if (curr == -1 || m_nodes[curr].m_max_subtree_size < size) {
                return -1;
            }
            while (curr != -1) {
                int left = m_nodes[curr].m_left;
                if (left != -1 && m_nodes[left].m_max_subtree_size >= size) {
                    curr = left;
                } else if (m_nodes[curr].m_node.size() >= size) {
                    break;
                } else {
                    curr = m_nodes[curr].m_right;
                }
            }
            return curr;
        }

        std::size_t getLargestFreeBlock () const {
            if (m_root == -1) {
                return 0;
            }
            return m_nodes[m_root].m_max_subtree_size;
        }

        const Node& viewNode (int idx) const {
            return m_nodes.at(idx).m_node;
        }

        Node extractNode (int idx) {
            Node node = m_nodes.at(idx).m_node;
            m_root = remove(m_root, node);
            return node;
        }

        void insertNode (const Node& node) {
            m_root = insert(m_root, node);
        }

        int find (const Node& node) const {
            int curr = m_root;
            while (curr != -1) {
                if (node < m_nodes[curr].m_node) {
                    curr = m_nodes[curr].m_left;
                } else if (m_nodes[curr].m_node < node) {
                    curr = m_nodes[curr].m_right;
                } else {
                    return curr;
                }
            }
            return curr;
        }

        int findPredecessor (const Node& node) const {
            int curr = m_root;
            int pred = -1;
            while (curr != -1) {
                if (node < m_nodes[curr].m_node) {
                    curr = m_nodes[curr].m_left;
                } else {
                    pred = curr;
                    curr = m_nodes[curr].m_right;
                }
            }
            return pred;
        }

        int findSuccessor (const Node& node) const {
            int curr = m_root;
            int succ = -1;
            while (curr != -1) {
                if (m_nodes[curr].m_node < node) {
                    curr = m_nodes[curr].m_right;
                } else {
                    succ = curr;
                    curr = m_nodes[curr].m_left;
                }
            }
            return succ;
        }

        int size () const {
            return static_cast<int>(m_nodes.size())
                - static_cast<int>(m_free_node_indices.size());
        }

        template<class F>
        void iterate (const F& f) const {
            iterate_imp(m_root, f);
        }

    private:

        struct NodeWithSize {
            Node m_node;
            //! The max subtree block size lets us find the first-fit block quickly
            std::size_t m_max_subtree_size;
            int m_height;
            int m_left;
            int m_right;
        };

        std::vector<NodeWithSize> m_nodes;
        std::vector<int> m_free_node_indices;
        int m_root = -1;

        int createNode (const Node& node) {
            int idx;
            if (!m_free_node_indices.empty()) {
                idx = m_free_node_indices.back();
                m_free_node_indices.pop_back();
                m_nodes[idx] = {node, node.size(), 1, -1, -1};
            } else {
                idx = static_cast<int>(m_nodes.size());
                m_nodes.push_back(NodeWithSize{node, node.size(), 1, -1, -1});
            }
            return idx;
        }

        void destroyNode (int idx) {
            m_free_node_indices.push_back(idx);
        }

        std::size_t getMax (int idx) const {
            return idx != -1 ? m_nodes[idx].m_max_subtree_size : 0;
        }

        int getHeight (int idx) const {
            return idx != -1 ? m_nodes[idx].m_height : 0;
        }

        int getBalance (int idx) const {
            return idx != -1 ?
                getHeight(m_nodes[idx].m_left) - getHeight(m_nodes[idx].m_right) : 0;
        }

        void updateMetadata (int idx) {
            if (idx == -1) return;
            int l = m_nodes[idx].m_left;
            int r = m_nodes[idx].m_right;
            m_nodes[idx].m_height = 1 + std::max(getHeight(l), getHeight(r));
            m_nodes[idx].m_max_subtree_size = std::max({
                m_nodes[idx].m_node.size(),
                getMax(l),
                getMax(r)
            });
        }

        int rotateRight (int y) {
            int x = m_nodes[y].m_left;
            int T2 = m_nodes[x].m_right;

            m_nodes[x].m_right = y;
            m_nodes[y].m_left = T2;

            updateMetadata(y);
            updateMetadata(x);
            return x;
        }

        int rotateLeft (int x) {
            int y = m_nodes[x].m_right;
            int T2 = m_nodes[y].m_left;

            m_nodes[y].m_left = x;
            m_nodes[x].m_right = T2;

            updateMetadata(x);
            updateMetadata(y);
            return y;
        }

        int balanceNode (int idx) {
            updateMetadata(idx);
            int balance = getBalance(idx);
            if (balance > 1) {
                if (getBalance(m_nodes[idx].m_left) < 0) {
                    m_nodes[idx].m_left = rotateLeft(m_nodes[idx].m_left);
                }
                return rotateRight(idx);
            }
            if (balance < -1) {
                if (getBalance(m_nodes[idx].m_right) > 0) {
                    m_nodes[idx].m_right = rotateRight(m_nodes[idx].m_right);
                }
                return rotateLeft(idx);
            }
            return idx;
        }

        int insert (int idx, const Node& node) {
            if (idx == -1) return createNode(node);
            if (node < m_nodes[idx].m_node) {
                m_nodes[idx].m_left = insert(m_nodes[idx].m_left, node);
            } else {
                m_nodes[idx].m_right = insert(m_nodes[idx].m_right, node);
            }
            return balanceNode(idx);
        }

        int findMin (int idx) const {
            while (m_nodes[idx].m_left != -1) {
                idx = m_nodes[idx].m_left;
            }
            return idx;
        }

        int remove (int idx, const Node& node) {
            if (idx == -1) return -1;
            if (node < m_nodes[idx].m_node) {
                m_nodes[idx].m_left = remove(m_nodes[idx].m_left, node);
            } else if (m_nodes[idx].m_node < node) {
                m_nodes[idx].m_right = remove(m_nodes[idx].m_right, node);
            } else {
                if (m_nodes[idx].m_left == -1 || m_nodes[idx].m_right == -1) {
                    int temp = (m_nodes[idx].m_left != -1) ?
                        m_nodes[idx].m_left : m_nodes[idx].m_right;
                    destroyNode(idx);
                    return temp;
                } else {
                    int temp = findMin(m_nodes[idx].m_right);
                    m_nodes[idx].m_node = m_nodes[temp].m_node;
                    m_nodes[idx].m_max_subtree_size = m_nodes[temp].m_max_subtree_size;
                    m_nodes[idx].m_right = remove(m_nodes[idx].m_right, m_nodes[temp].m_node);
                }
            }
            return balanceNode(idx);
        }

        template<class F>
        void iterate_imp (int idx, const F& f) const {
            if (idx != -1) {
                iterate_imp(m_nodes[idx].m_left, f);
                f(m_nodes[idx].m_node);
                iterate_imp(m_nodes[idx].m_right, f);
            }
        }
    };

    //! The list of blocks allocated.
    std::vector<std::pair<void*,std::size_t> > m_alloc;

    /**
    * \brief The type of our freelist and blocklist.
    * We use a set sorted from lo to hi memory addresses.
    */
    using NL = NodeTreeWithSize;

    /**
    * \brief The free list of allocated but not currently used blocks.
    * Maintained in lo to hi memory sorted order.
    */
    NL m_freelist;

    /**
    * \brief The list of busy blocks.
    * A block is either on the freelist or on the blocklist, but not on both.
    */
//    NL m_busylist;
    std::unordered_set<Node, Node::hash> m_busylist;
    //! The minimal size of hunks to request from system
    std::size_t m_hunk;
    //! The amount of heap space currently allocated.
    //! Atomic, so it can be read without carena_mutex.
    std::atomic<std::size_t> m_used{0};
    //! The max amount of heap space currently allocated.
    std::size_t m_max_used{0};
    //! The amount of memory given out via alloc().
    //! Atomic, so it can be read without carena_mutex.
    std::atomic<std::size_t> m_actually_used{0};
    //! The max amount of memory given out via alloc().
    std::size_t m_max_actually_used{0};

    mutable std::mutex carena_mutex;

    //! Print a summary of the arena's free and busy lists to \p os.
    friend std::ostream& operator<< (std::ostream& os, const CArena& arena);
};

}

#endif /*BL_CARENA_H*/
