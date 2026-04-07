#ifndef AMREX_TRACKED_VECTOR_H
#define AMREX_TRACKED_VECTOR_H

#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_GpuContainers.H>

#include <memory>
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
        {
            m_data->host.resize(a_size);
#ifdef AMREX_USE_GPU
            m_data->status = Status::host_dirty;
#endif
        }

        TrackedVector (size_type a_size, value_type const & a_value)
        {
            m_data->host.assign(a_size, a_value);
#ifdef AMREX_USE_GPU
            m_data->status = Status::host_dirty;
#endif
        }

        TrackedVector (std::initializer_list<T> a_initializer_list)
        {
            m_data->host = a_initializer_list;
#ifdef AMREX_USE_GPU
            m_data->status = Status::host_dirty;
#endif
        }

        TrackedVector (std::vector<T> a_vector)
        {
            m_data->host = std::move(a_vector);
#ifdef AMREX_USE_GPU
            m_data->status = Status::host_dirty;
#endif
        }

        TrackedVector (TrackedVector const & a_vector)
            : m_data(std::make_shared<Data>())
        {
            m_data->status = a_vector.m_data->status;
            m_data->host = a_vector.m_data->host;
#ifdef AMREX_USE_GPU
            m_data->device = a_vector.m_data->device;
            if (amrex::Initialized()) {
                register_finalize();
            }
#endif
        }

        /** Swap data of this and a_vector */
        TrackedVector (TrackedVector && a_vector) noexcept
        {
            std::swap(m_data, a_vector.m_data);
        }

        TrackedVector& operator= (TrackedVector const & a_vector) {
            if (this != &a_vector) {
                m_data->status = a_vector.m_data->status;
                m_data->host = a_vector.m_data->host;
#ifdef AMREX_USE_GPU
                m_data->device = a_vector.m_data->device;
                if (amrex::Initialized()) {
                    register_finalize();
                }
#endif
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
        [[nodiscard]] amrex::Gpu::DeviceVector<T> &
        device () {
            if (!amrex::Initialized()) {
                throw std::runtime_error("TrackedVector::device() called before AMReX initialize/after AMReX finalize");
            }
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
        [[nodiscard]] amrex::Gpu::DeviceVector<T> const &
        device_const () const {
            if (!amrex::Initialized()) {
                throw std::runtime_error("TrackedVector::device_const() called before AMReX initialize/after AMReX finalize");
            }
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
            if (m_data->status == Status::device_dirty) {
                to_host();
            }
            m_data->device.clear();
            m_data->device.shrink_to_fit();
            m_data->status = Status::host_dirty;
#endif
        }

    private:

        struct Data {
            Status status = Status::up_to_date;
            std::vector<T> host;
#ifdef AMREX_USE_GPU
            amrex::Gpu::DeviceVector<T> device;
            bool finalize_registered = false;
#endif
        };

        /** Synchronize host data to device */
        void to_device () const {
#ifdef AMREX_USE_GPU
            if (!amrex::Initialized()) {
                throw std::runtime_error("TrackedVector::to_device() called outside of AMReX initialize/finalize");
            }
            auto const size = m_data->host.size();
            if (size > 0U) {
                m_data->device.resize(size);
                amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                    m_data->host.begin(), m_data->host.end(), m_data->device.begin());
            } else {
                m_data->device.clear();
            }
            m_data->status = Status::up_to_date;
#endif
        }

        /** Synchronize device data to host */
        void to_host () const {
#ifdef AMREX_USE_GPU
            if (!amrex::Initialized()) {
                throw std::runtime_error("TrackedVector::to_host() called outside of AMReX initialize/finalize");
            }
            m_data->host.resize(m_data->device.size());
            amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                m_data->device.begin(), m_data->device.end(), m_data->host.begin());
            m_data->status = Status::up_to_date;
#endif
        }

        void register_finalize () const {
#ifdef AMREX_USE_GPU
            if (!amrex::Initialized()) {
                throw std::runtime_error("TrackedVector::register_finalize() called outside of AMReX initialize/finalize");
            }

            if (m_data->finalize_registered) { return; }
            m_data->finalize_registered = true;

            std::weak_ptr<Data> weak_data = m_data;

            amrex::ExecOnFinalize([weak_data]() {
                auto data = weak_data.lock();
                if (!data) { return; }

                if (data->status == Status::device_dirty
                    && !data->device.empty())
                {
                    data->host.resize(data->device.size());
                    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                        data->device.begin(), data->device.end(), data->host.begin());
                }
                data->device.clear();
                data->device.shrink_to_fit();
                data->status = Status::host_dirty;
                data->finalize_registered = false;
            });
#endif
        }

        mutable std::shared_ptr<Data> m_data = std::make_shared<Data>();
    };

}

#endif
