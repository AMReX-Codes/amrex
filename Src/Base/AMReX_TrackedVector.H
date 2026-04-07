#ifndef AMREX_TRACKED_VECTOR_H
#define AMREX_TRACKED_VECTOR_H

#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_GpuContainers.H>

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>


namespace amrex::Gpu
{
    /** Base for element data with dynamic host vectors mirrored lazily to GPU.
     *
     * Provides lifetimes independent of AMReX initialize/finalize cycles,
     * synchronization tracking, and data release APIs.
     *
     * This object is primarily for input handling, allowing to initialize
     * data even before AMReX was initialized and enabling workflows crossing
     * AMReX init/finalize cycles. GPU memory will always be bound to an
     * AMReX session, but CPU memory can be allocated and live arbitrarily long.
     *
     * For AMReX CPU builds, the host and device members point to the *same* memory
     * and the status will always be up to date.
     *
     * Usage contract:
     *   - Device data can only be allocated after AMReX initialize and before finalize.
     *   - You can call `release_gpu()` anytime, but we will call it during AMReX finalize
     *     to invalidate device().
     *   - Always access data via host/device[_const](). Do not cache references/pointers to
     *     host/device memory managed by this object, or you run the risk of stale memory access.
     *   - Read-only accessors (`host_const`, `device_const`) may trigger a synchronous
     *     copy when the opposite side is dirty. They are logically const but not
     *     side-effect-free.
     *   - This class is **not** thread-safe. All accesses (including `_const` variants)
     *     must be externally synchronized if used from multiple threads.
     */
    template <class T>
    struct TrackedVector
    {
        static_assert(std::is_trivially_copyable<T>(), "TrackedVector can only hold trivially copyable types");
        using value_type      = T;
        using size_type       = std::size_t;

        enum class Status {
            up_to_date,    //!< host and device data are in sync
            device_dirty,  //!< host data needs an update
            host_dirty     //!< device data needs an update
        };

        TrackedVector () = default;

        explicit TrackedVector (size_type a_size)
            : m_data(std::make_shared<Data>(std::vector<T>(a_size))) {}

        TrackedVector (size_type a_size, value_type const & a_value)
            : m_data(std::make_shared<Data>(std::vector<T>(a_size, a_value))) {}

        TrackedVector (std::initializer_list<T> a_initializer_list)
            : m_data(std::make_shared<Data>(std::vector<T>(a_initializer_list))) {}

        TrackedVector (std::vector<T> a_vector)
            : m_data(std::make_shared<Data>(std::move(a_vector))) {}

        TrackedVector (TrackedVector const & a_vector)
        {
            copy_from(a_vector);
        }

        /** Swap data of this and a_vector */
        TrackedVector (TrackedVector && a_vector) noexcept
        {
            std::swap(m_data, a_vector.m_data);
        }

        /** Copy data from a_vector */
        TrackedVector& operator= (TrackedVector const & a_vector) {
            if (this != &a_vector) {
                copy_from(a_vector);
            }
            return *this;
        }

        /** Swap the data of this and a_vector */
        TrackedVector& operator= (TrackedVector && a_vector) noexcept {
            if (this != &a_vector) {
                std::swap(m_data, a_vector.m_data);
            }
            return *this;
        }

        ~TrackedVector () = default;

        [[nodiscard]] Status status () const { return m_data->status; }

        /** Return writable host data
         *
         * If device data is newer, syncs device-to-host first (GPU builds).
         * Marks host as dirty so the next device access will sync back.
         */
        [[nodiscard]] std::vector<T> &
        host () {
#ifdef AMREX_USE_GPU
            if (m_data->status == Status::device_dirty) {
                to_host();
            }
            m_data->status = Status::host_dirty;
#endif
            return m_data->host;
        }

        /** Return read-only host data
         *
         * If device data is newer, syncs device-to-host first (GPU builds).
         */
        [[nodiscard]] std::vector<T> const &
        host_const () const {
#ifdef AMREX_USE_GPU
            if (m_data->status == Status::device_dirty) {
                to_host();
            }
#endif
            return m_data->host;
        }

#ifdef AMREX_USE_GPU
        /** Return writable device data
         *
         * If host data is newer, syncs host-to-device first.
         * Marks device as dirty so the next host access will sync back.
         */
        [[nodiscard]] amrex::Gpu::NonManagedDeviceVector<T> &
        device () {
            require_amrex("device");
            register_finalize();
            if (m_data->status == Status::host_dirty) {
                to_device();
            }
            m_data->status = Status::device_dirty;
            return m_data->device;
        }

        /** Return read-only device data
         *
         * If host data is newer, syncs host-to-device first.
         */
        [[nodiscard]] amrex::Gpu::NonManagedDeviceVector<T> const &
        device_const () const {
            require_amrex("device_const");
            register_finalize();
            if (m_data->status == Status::host_dirty) {
                to_device();
            }
            return m_data->device;
        }
#else
        /** Return writable device (==host) data, up to date by definition */
        [[nodiscard]] std::vector<T> &
        device () { return m_data->host; }

        /** Return read-only device (==host) data */
        [[nodiscard]] std::vector<T> const &
        device_const () const { return m_data->host; }
#endif

        /** Release GPU memory
         *
         * If the device side was modified (device_dirty), syncs device data
         * back to the host first so no updates are lost.
         * Host data preserved for reuse until destructor.
         * This enables use outside of and across AMReX init/finalize cycles.
         */
        void release_gpu ()
        {
#ifdef AMREX_USE_GPU
            Data::drain_device(*m_data);
#endif
        }

    private:

        struct Data {
            Status status = Status::up_to_date;
            std::vector<T> host;
#ifdef AMREX_USE_GPU
            amrex::Gpu::NonManagedDeviceVector<T> device;
            bool finalize_registered = false;
#endif
            Data () = default;

            /** Host owns authoritative data */
            explicit Data (std::vector<T> h)
                : host(std::move(h))
            {
#ifdef AMREX_USE_GPU
                status = Status::host_dirty;
#endif
            }

#ifdef AMREX_USE_GPU
            /** Sync device data back to host if dirty, then free device memory
             *
             * Leaves with the host side marked dirty, so a .device/.device_const
             * call in a future AMReX session resyncs to device.
             * A dirty marked device side remains authoritative even when it was
             * shrunk to size zero through the mutable device vector returned by device().
             *
             * Assumes the GPU context is still alive (called from ExecOnFinalize,
             * which runs before GPU device teardown in amrex::Finalize).
             */
            static void drain_device (Data& d) {
                if (d.status == Status::device_dirty) {
                    amrex::Gpu::streamSynchronize();
                    d.host.resize(d.device.size());
                    if (!d.device.empty()) {
                        amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                            d.device.begin(), d.device.end(), d.host.begin());
                    }
                }
                d.device.clear();
                d.device.shrink_to_fit();
                d.status = Status::host_dirty;
            }
#endif
        };

        /** Check amrex::Initialized or throw */
        static void require_amrex (char const * func) {
            if (!amrex::Initialized()) {
                throw std::runtime_error(
                    std::string("TrackedVector::") + func +
                    " called outside of AMReX initialize/finalize");
            }
        }

        /** Helper to copy Data and register conditionally */
        void copy_from (TrackedVector const & a_vector) {
#ifdef AMREX_USE_GPU
            // Outside an AMReX session, host data is the only
            // usable copy.  A device_dirty source has stale host
            // data and no accessible device, so the copy would
            // silently lose updates.
            if (!amrex::Initialized()
                && a_vector.m_data->status == Status::device_dirty)
            {
                throw std::runtime_error(
                    "TrackedVector::copy_from: source is device_dirty "
                    "outside an AMReX session (host data is stale)");
            }
            bool const was_finalize_registered = m_data->finalize_registered;
#endif
            *m_data = *a_vector.m_data;
#ifdef AMREX_USE_GPU
            m_data->finalize_registered = was_finalize_registered;
            if (amrex::Initialized()) {
                // Don't copy stale device data — it will be resynced on next access.
                if (m_data->status == Status::host_dirty) {
                    m_data->device.clear();
                    m_data->device.shrink_to_fit();
                }
                if (!was_finalize_registered) {
                    register_finalize();
                }
            }
#endif
        }

        /** Synchronize host data to device */
        void to_device () const {
#ifdef AMREX_USE_GPU
            require_amrex("to_device");
            auto const size = m_data->host.size();
            if (size > 0U) {
                m_data->device.resize(size);
                amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                    m_data->host.begin(), m_data->host.end(), m_data->device.begin());
            } else {
                m_data->device.clear();
                m_data->device.shrink_to_fit();
            }
            m_data->status = Status::up_to_date;
#endif
        }

        /** Synchronize device data to host */
        void to_host () const {
#ifdef AMREX_USE_GPU
            require_amrex("to_host");
            amrex::Gpu::streamSynchronize();
            m_data->host.resize(m_data->device.size());
            if (!m_data->device.empty()) {
                amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                    m_data->device.begin(), m_data->device.end(), m_data->host.begin());
            }
            m_data->status = Status::up_to_date;
#endif
        }

        /** Register with AMReX to free device memory on Finalize()
         *
         * This function is idempotent and drains device memory to
         * host first if dirty.
         */
        void register_finalize () const {
#ifdef AMREX_USE_GPU
            require_amrex("register_finalize");
            if (m_data->finalize_registered) { return; }
            m_data->finalize_registered = true;

            std::weak_ptr<Data> weak_data = m_data;

            amrex::ExecOnFinalize([weak_data]() {
                auto data = weak_data.lock();
                if (!data) { return; }
                Data::drain_device(*data);
                data->finalize_registered = false;
            });
#endif
        }

        // mutable: lazy-sync accessors (host_const, device_const) need to
        // update status and trigger copies through the shared_ptr even on
        // const TrackedVector references.
        mutable std::shared_ptr<Data> m_data = std::make_shared<Data>();
    };

}

#endif
