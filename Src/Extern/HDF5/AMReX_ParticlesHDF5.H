#ifndef AMREX_PARTICLES_HDF5_H_
#define AMREX_PARTICLES_HDF5_H_

/**
 * \file AMReX_ParticlesHDF5.H
 *
 * HDF5-specific ParticleContainer hooks (checkpointing, plotfile output, and
 * supporting helpers shared by all particle types).
 */

public:
    /**
     * \brief Write a contiguous chunk of real particle data to an HDF5 file.
     *
     * \param data A pointer to the start of the buffer to write
     * \param size The number of elements to write
     * \param fid The HDF5 file ID
     */
    void WriteParticleRealDataHDF5 (void* data, size_t size, hid_t fid) const;

    /**
     * \brief Read a contiguous chunk of real particle data from HDF5 file.
     *
     * \param data A pointer to the start of the buffer into which to read
     * \param size The number of elements to read
     * \param fid The istream from which to read the data
     */
    void ReadParticleRealDataHDF5 (void* data, size_t size, hid_t fid);

    /**
     * \brief Writes a particle checkpoint to file, suitable for restarting.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     */
    void CheckpointHDF5 (const std::string& dir, const std::string& name,
                         const std::string& compression = "None@0") const;

    /**
     * \brief Writes a particle checkpoint to file, suitable for restarting.
     *
     * \param dir              Directory that owns the per-type subdirectory.
     * \param name             Particle type subdirectory (e.g., "Tracer").
     * \param is_checkpoint    Reserved for parity with ASCII checkpoints (ignored because HDF5 always writes id/cpu).
     * \param real_comp_names  Optional labels for each real component.
     * \param int_comp_names   Optional labels for each integer component.
     * \param compression      Compression scheme, e.g., `"ZFP_ACCURACY@0.001"`.
     */
    void CheckpointHDF5 (const std::string& dir, const std::string& name, bool is_checkpoint,
                         const Vector<std::string>& real_comp_names = Vector<std::string>(),
                         const Vector<std::string>& int_comp_names = Vector<std::string>(),
                         const std::string& compression = "None@0") const;

     /**
      * \brief Writes particle data to disk in the HDF5 format.
      *
      * \tparam F function type
      *
      * \param dir The base directory into which to write (i.e. "plt00000")
      * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
      * \param write_real_comp for each real component, whether or not we include that component in the file
      * \param write_int_comp for each integer component, whether or not we include that component in the file
      * \param real_comp_names for each real component, a name to label the data with
      * \param int_comp_names for each integer component, a name to label the data with
      * \param f callable that returns whether a given particle should be written or not
      * \param is_checkpoint Whether the particle id and cpu are included in the file.
      * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
      */
    template <class F>
    void WriteHDF5ParticleData (const std::string& dir,
                                const std::string& name,
                                const Vector<int>& write_real_comp,
                                const Vector<int>& write_int_comp,
                                const Vector<std::string>& real_comp_names,
                                const Vector<std::string>&  int_comp_names,
                                const std::string& compression,
                                F&& f, bool is_checkpoint=false) const;


    /**
     * \brief Prepare metadata used by deferred HDF5 checkpoint writes.
     *
     * Tallies particles/NextID values per level so CheckpointPostHDF5() can
     * finish the header without revisiting particle data.
     */
    void CheckpointPreHDF5 ();

    /**
     * \brief Finalize a checkpoint header handed off by CheckpointPreHDF5().
     *
     * Aggregates distributed offsets/counts, appends per-level bookkeeping,
     * and removes empty per-file fragments when requested.
     */
    void CheckpointPostHDF5 ();

    /**
     * \brief Restart from checkpoint
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     */
    void RestartHDF5 (const std::string& dir, const std::string& file);

    /**
     * \brief Older version, for backwards compatibility
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param is_checkpoint Whether the particle id and cpu are included in the file.
     */
    void RestartHDF5 (const std::string& dir, const std::string& file, bool is_checkpoint);

    /**
     * \brief This version of WritePlotFile writes all components and assigns component names
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     */
    void WritePlotFileHDF5 (const std::string& dir, const std::string& name, const std::string& compression) const;

    /**
     * \brief This version of WritePlotFile writes all components and assigns component names
     *
     * This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     * \param f callable that returns whether or not to write each particle
     */
    template <class F, std::enable_if_t<!std::is_same_v<F, Vector<std::string>>>* = nullptr>
    void WritePlotFileHDF5 (const std::string& dir, const std::string& name, const std::string& compression, F&& f) const;

    /**
     * \brief This version of WritePlotFile writes all components and allows the user to specify the names of the components.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     */
    void WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                            const Vector<std::string>& real_comp_names,
                            const Vector<std::string>&  int_comp_names,
                            const std::string& compression) const;

    /**
     * \brief This version of WritePlotFile writes all components and allows the user to specify the names of the components.
     *
     * This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     * \param f callable that returns whether or not to write each particle
     */
    template <class F>
    void WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                            const Vector<std::string>& real_comp_names,
                            const Vector<std::string>&  int_comp_names,
                            const std::string& compression, F&& f) const;

    /**
     * \brief This version of WritePlotFile writes all components and allows the user to specify
     * the names of the components.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     */
    void WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                            const Vector<std::string>& real_comp_names,
                            const std::string& compression) const;

    /**
     * \brief This version of WritePlotFile writes all components and allows the user to specify
     * the names of the components.
     *
     * This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     * \param f callable that returns whether or not to write each particle
     */
    template <class F, std::enable_if_t<!std::is_same_v<F, Vector<std::string>>>* = nullptr>
    void WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                            const Vector<std::string>& real_comp_names,
                            const std::string& compression, F&& f) const;

    /**
     * \brief This version of WritePlotFile assigns component names, but allows the user to pass
     * in a vector of ints that toggle on / off the writing of specific components.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     */
    void WritePlotFileHDF5 (const std::string& dir,
                            const std::string& name,
                            const Vector<int>& write_real_comp,
                            const Vector<int>& write_int_comp,
                            const std::string& compression) const;

    /**
     * \brief This version of WritePlotFile assigns component names, but allows the user to pass
     * in a vector of ints that toggle on / off the writing of specific components.
     *
     * This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     * \param f callable that returns whether or not to write each particle
     */
    template <class F>
    void WritePlotFileHDF5 (const std::string& dir,
                            const std::string& name,
                            const Vector<int>& write_real_comp,
                            const Vector<int>& write_int_comp,
                            const std::string& compression, F&& f) const;

    /**
     * \brief This is the most general version of WritePlotFile, which takes component
     * names and flags for whether to write each variable as components. Note that
     * the user should pass in vectors containing names of all the components, whether
     * they are written or not.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     */
    void WritePlotFileHDF5 (const std::string& dir,
                            const std::string& name,
                            const Vector<int>& write_real_comp,
                            const Vector<int>& write_int_comp,
                            const Vector<std::string>& real_comp_names,
                            const Vector<std::string>& int_comp_names,
                            const std::string& compression) const;

    /**
     * \brief This is the most general version of WritePlotFile, which takes component
     * names and flags for whether to write each variable as components. Note that
     * the user should pass in vectors containing names of all the components, whether
     * they are written or not.
     *
     *  This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param name The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     * \param f callable that returns whether or not to write each particle
     */
    template <class F>
    void WritePlotFileHDF5 (const std::string& dir,
                            const std::string& name,
                            const Vector<int>& write_real_comp,
                            const Vector<int>& write_int_comp,
                            const Vector<std::string>& real_comp_names,
                            const Vector<std::string>&  int_comp_names,
                            const std::string& compression,
                            F&& f) const;

    /**
     * \brief Plotfile-friendly alias for CheckpointPreHDF5().
     */
    void WritePlotFilePreHDF5 ();

    /**
     * \brief Plotfile-friendly alias for CheckpointPostHDF5().
     */
    void WritePlotFilePostHDF5 ();

    /**
     * \brief Write all particles from one AMR level into the given HDF5 group.
     *
     * \param level            AMR level whose particles are written.
     * \param grp              Open HDF5 group corresponding to `level_%d`.
     * \param which            Placeholder vector kept for API compatibility (unused by the HDF5 backend).
     * \param count            Per-grid particle counts updated in place.
     * \param where            Placeholder vector kept for API compatibility (unused by the HDF5 backend).
     * \param write_real_comp  Flags enabling each real component.
     * \param write_int_comp   Flags enabling each integer component.
     * \param compression      Compression string (`ZFP_*`, `SZ`, etc.).
     * \param particle_io_flags Per-grid/tile masks returned by the filter functor.
     * \param is_checkpoint    Reserved; the HDF5 writer always preserves id/cpu regardless of this flag.
     */
    void
    WriteParticlesHDF5 (int level, hid_t grp,
                        Vector<int>& which, Vector<int>& count, Vector<Long>& where,
                        const Vector<int>& write_real_comp, const Vector<int>& write_int_comp,
                        const std::string& compression,
                        const Vector<std::map<std::pair<int, int>,IntVector>>& particle_io_flags, bool is_checkpoint) const;

protected:

/**
 * \brief Read particle records from HDF5 datasets into a level/tile buffer.
 *
 * \tparam RTYPE  Floating-point storage type stored in the file.
 * \param offset  Starting element offset inside the dataset.
 * \param count   Number of entries to load.
 * \param grd     Grid index being populated.
 * \param lev     AMR level index being populated.
 * \param int_dset HDF5 dataset handle for integer components.
 * \param real_dset HDF5 dataset handle for real components.
 * \param finest_level_in_file Finest level encoded in the file (for sanity checks).
 * \param convert_ids True when reading `Version_Two_Dot_One` checkpoints whose 64-bit id/cpu pairs were split into two 32-bit words and must be reassembled.
 */
template <class RTYPE>
void ReadParticlesHDF5 (hsize_t offset, hsize_t count, int grd, int lev, hid_t int_dset, hid_t real_dset, int finest_level_in_file, bool convert_ids);

#endif
