#ifndef AMREX_PARTICLES_HDF5_H_
#define AMREX_PARTICLES_HDF5_H_

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
     * \param fid The HDF5 file ID
     * \param compression compression parameter (i.e. "ZFP_ACCURACY@0.001", "SZ@sz.config")
     */
    void CheckpointHDF5 (const std::string& dir, const std::string& name,
                         const std::string& compression = "None@0") const;

    /**
     * \brief Writes a particle checkpoint to file, suitable for restarting.
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


    void CheckpointPreHDF5 ();

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
     * \param f callable that returns whether or not to write each particle
     */
    template <class F, typename std::enable_if<!std::is_same<F, Vector<std::string>>::value>::type* = nullptr>
    void WritePlotFileHDF5 (const std::string& dir, const std::string& name, const std::string& compression, F&& f) const;

    /**
     * \brief This version of WritePlotFile writes all components and allows the user to specify the names of the components.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
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
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
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
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
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
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param real_comp_names for each real component, a name to label the data with
     * \param f callable that returns whether or not to write each particle
     */
    template <class F, typename std::enable_if<!std::is_same<F, Vector<std::string>>::value>::type* = nullptr>
    void WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                            const Vector<std::string>& real_comp_names,
                            const std::string& compression, F&& f) const;

    /**
     * \brief This version of WritePlotFile assigns component names, but allows the user to pass
     * in a vector of ints that toggle on / off the writing of specfic components.
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     */
    void WritePlotFileHDF5 (const std::string& dir,
                            const std::string& name,
                            const Vector<int>& write_real_comp,
                            const Vector<int>& write_int_comp,
                            const std::string& compression) const;

    /**
     * \brief This version of WritePlotFile assigns component names, but allows the user to pass
     * in a vector of ints that toggle on / off the writing of specfic components.
     *
     * This version also lets you pass in a functor to toggle whether each particle gets output.
     *
     * \tparam F function type
     *
     * \param dir The base directory into which to write (i.e. "plt00000")
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
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
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
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
     * \param file The name of the sub-directory for this particle type (i.e. "Tracer")
     * \param write_real_comp for each real component, whether to include that comp in the file
     * \param write_int_comp for each integer component, whether to include that comp in the file
     * \param real_comp_names for each real component, a name to label the data with
     * \param int_comp_names for each integer component, a name to label the data with
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

    void WritePlotFilePreHDF5 ();

    void WritePlotFilePostHDF5 ();

    void
    WriteParticlesHDF5 (int level, hid_t grp,
                        Vector<int>& which, Vector<int>& count, Vector<Long>& where,
                        const Vector<int>& write_real_comp, const Vector<int>& write_int_comp,
                        const std::string& compression,
                        const Vector<std::map<std::pair<int, int>,IntVector>>& particle_io_flags, bool is_checkpoint) const;

protected:

template <class RTYPE>
void ReadParticlesHDF5 (hsize_t offset, hsize_t count, int grd, int lev, hid_t int_dset, hid_t real_dset, int finest_level_in_file, bool convert_ids);

#endif
