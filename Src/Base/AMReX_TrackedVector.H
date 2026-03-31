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

    private:
        void register_finalize () {
#ifdef AMREX_USE_GPU
            std::weak_ptr<amrex::Gpu::DeviceVector<T>> weak_device = m_device;
            std::weak_ptr<Status> weak_status = m_status;
            amrex::ExecOnFinalize([weak_device, weak_status]() {
                // see: release_gpu()
                if (auto self = weak_device.lock()) {
                    self->clear();
                    self->shrink_to_fit();
                }
                if (auto self = weak_status.lock()) {
                    *self = Status::host_dirty;
                }
            });
#endif
        }
    public:

        constexpr TrackedVector () {
            register_finalize();
        }

        explicit TrackedVector (size_type a_size)
            : m_host(a_size)
        {
            register_finalize();
            host_dirty();
        }

        TrackedVector (size_type a_size, value_type const & a_value)
            : m_host(a_size, a_value)
        {
            register_finalize();
            host_dirty();
        }

        TrackedVector (std::initializer_list<T> a_initializer_list)
            : m_host(a_initializer_list)
        {
            register_finalize();
            host_dirty();
        }

        TrackedVector (std::vector<T> a_vector)
            : m_host(std::move(a_vector))
        {
            register_finalize();
            host_dirty();
        }

        TrackedVector (TrackedVector const & a_vector)
            : m_host(a_vector.m_host)
        {
            *m_status = *a_vector.m_status;
#ifdef AMREX_USE_GPU
            *m_device = *a_vector.m_device;
#endif

            register_finalize();
        }

        /** Swap the empty data of this and a_vector */
        TrackedVector (TrackedVector && a_vector) noexcept
        {
            std::swap(m_status, a_vector.m_status);
            std::swap(m_host, a_vector.m_host);
#ifdef AMREX_USE_GPU
            std::swap(m_device, a_vector.m_device);
#endif

            // We inherit a_vector's original registration.
            // But we need to register a_vector's new data now:
            //   the shared ptrs that we owned briefly up to here.
            // That way, a_vector can be either reused/assigned values
            // or safely destructed.
            a_vector.register_finalize();
        }

        TrackedVector& operator= (TrackedVector const & a_vector) {
            if (this != &a_vector) {
                *m_status = *a_vector.m_status;
                m_host = a_vector.m_host;
#ifdef AMREX_USE_GPU
                *m_device = *a_vector.m_device;
#endif
            }
            return *this;
        }

        /** Swap the data of this and a_vector */
        TrackedVector& operator= (TrackedVector && a_vector) noexcept {
            if (this != &a_vector) {
                std::swap(m_host, a_vector.m_host);
                std::swap(m_status, a_vector.m_status);
#ifdef AMREX_USE_GPU
                std::swap(m_device, a_vector.m_device);
#endif
            }
            return *this;
        }

        ~TrackedVector () = default;

        enum class Status {
            up_to_date,    //<! host and device data are in sync
            device_dirty,  //!< host data needs an update
            host_dirty     //!< device data needs an update
        };

        [[nodiscard]] Status status () const { return *m_status; }

        /** Return writable host data, mark host as dirty */
        [[nodiscard]] std::vector<T> &
        host () {
            host_dirty();
            return m_host;
        }

        /** Return read-only host data */
        [[nodiscard]] std::vector<T> const &
        host_const () const { return m_host; }

#ifdef AMREX_USE_GPU
        /** Return writable device data, mark device as dirty */
        [[nodiscard]] amrex::Gpu::DeviceVector<T> &
        device () {
            if (!amrex::Initialized()) {
                throw std::runtime_error("TrackedVector::device() called before AMReX initialize/after AMReX finalize");
            }
            device_dirty();
            return *m_device;
        }

        /** Return read-only device data */
        [[nodiscard]] amrex::Gpu::DeviceVector<T> const &
        device_const () const {
            if (!amrex::Initialized()) {
                throw std::runtime_error("TrackedVector::device() called before AMReX initialize/after AMReX finalize");
            }
            return *m_device;
        }
#else
        /** Return writable device (==host) data, up to date by definition */
        [[nodiscard]] std::vector<T> &
        device () { return m_host; }

        /** Return read-only device (==host) data */
        [[nodiscard]] std::vector<T> const &
        device_const () const { return m_host; }
#endif

        /** Release GPU memory
         *
         * Host data preserved for reuse until destructor.
         * This enables use outside of and across AMReX init/finalize cycles.
         */
        void release_gpu ()
        {
#ifdef AMREX_USE_GPU
            device().clear();
            device().shrink_to_fit();
            *m_status = Status::host_dirty;
#endif
        }

        /** Mark host data as modified; next sync will copy to GPU. */
        void host_dirty () {
#ifdef AMREX_USE_GPU
            *m_status = Status::host_dirty;
#endif
        }

        /** Mark device data as modified; next sync will copy to CPU. */
        void device_dirty () {
#ifdef AMREX_USE_GPU
            *m_status = Status::device_dirty;
#endif
        }

        /** Conditionally synchronize host data to device
         *
         * If either the host or device side are dirty, this will
         * copy the host data over the device data.
         *
         * @param[in] force  If true, force a copy even if the status is up to date.
         */
        void to_device (bool force=false) {
#ifdef AMREX_USE_GPU
            if (status() != Status::up_to_date || force) {
                if (!amrex::Initialized()) {
                    throw std::runtime_error("TrackedVector::to_device() called outside of AMReX initialize/finalize");
                }
                auto const size = m_host.size();
                if (size > 0U) {
                    m_device->resize(size);
                    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                        m_host.begin(), m_host.end(), m_device->begin());
                } else {
                    m_device->clear();
                }
            }
            *m_status = Status::up_to_date;
#else
            amrex::ignore_unused(force);
#endif
        }

        /** Conditionally synchronize device data to host
         *
         * If either the host or device side are dirty, this will
         * copy the device data over the host data.
         *
         * @param[in] force  If true, force a copy even if the status is up to date.
         */
        void to_host (bool force=false) {
#ifdef AMREX_USE_GPU
            if (status() != Status::up_to_date || force) {
                if (!amrex::Initialized()) {
                    throw std::runtime_error("TrackedVector::to_host() called outside of AMReX initialize/finalize");
                }
                m_host.resize(m_device->size());
                amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                    m_device->begin(), m_device->end(), m_host.begin());
            }
            *m_status = Status::up_to_date;
#else
            amrex::ignore_unused(force);
#endif
        }

        /** Ensure equivalent host-device data (sync copy)
         *
         * This performs a conditional synchronous copy from whatever
         * host()/device() side was accessed last.
         */
        void ensure_same ()
        {
#ifdef AMREX_USE_GPU
            if (status() == Status::host_dirty) {
                to_device(false);
            } else if (status() == Status::device_dirty) {
                to_host(false);
            }
#endif
        }

    private:
        std::shared_ptr<Status> m_status = std::make_shared<Status>();
        std::vector<T> m_host;
#ifdef AMREX_USE_GPU
        std::shared_ptr<amrex::Gpu::DeviceVector<T>> m_device = std::make_shared<amrex::Gpu::DeviceVector<T>>();
#endif
    };

}

#endif
