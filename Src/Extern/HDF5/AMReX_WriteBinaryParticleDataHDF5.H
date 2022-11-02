#ifndef AMREX_WRITE_BINARY_PARTICLE_DATA_HDF5_H
#define AMREX_WRITE_BINARY_PARTICLE_DATA_HDF5_H

static int CreateWriteHDF5Attr(hid_t loc, const char *name, hsize_t n, void *data, hid_t dtype)
{
    herr_t ret;
    hid_t attr, attr_space;
    hsize_t dims = n;

    attr_space = H5Screate_simple(1, &dims, NULL);

    attr = H5Acreate(loc, name, dtype, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        amrex::Print() << " Error with H5Acreate [" << name << "]\n";
        return -1;
    }

    ret  = H5Awrite(attr, dtype, (void*)data);
    if (ret < 0) {
        amrex::Print() << " Error with H5Awrite [" << name << "]\n";
        return -1;
    }
    H5Sclose(attr_space);
    H5Aclose(attr);
    return 1;
}

static int CreateWriteHDF5AttrString(hid_t loc, const char *name, const char* str)
{
    hid_t attr, atype, space;
    herr_t ret;

    AMREX_ASSERT(name);
    AMREX_ASSERT(str);

    space = H5Screate(H5S_SCALAR);
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(str)+1);
    H5Tset_strpad(atype,H5T_STR_NULLTERM);
    attr = H5Acreate(loc, name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        amrex::Print() << " Error with H5Acreate [" << name << "]\n";
        return -1;
    }

    ret = H5Awrite(attr, atype, str);
    if (ret < 0) {
        amrex::Print() << " Error with H5Awrite [" << name << "]\n";
        return -1;
    }

    H5Tclose(atype);
    H5Sclose(space);
    H5Aclose(attr);

    return 1;
}

static int ReadHDF5Attr(hid_t loc, const char *name, void *data, hid_t dtype)
{
    herr_t ret;
    hid_t attr;

    attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr < 0) {
        amrex::Print() << " Error with H5Aopen [" << name << "]\n";
        return -1;
    }

    ret  = H5Aread(attr, dtype, data);
    if (ret < 0) {
        amrex::Print() << " Error with H5Aread [" << name << "]\n";
        return -1;
    }
    H5Aclose(attr);
    return 1;
}

#ifdef BL_USE_MPI
static void SetHDF5fapl(hid_t fapl, MPI_Comm comm)
#else
static void SetHDF5fapl(hid_t fapl)
#endif
{
#ifdef BL_USE_MPI
    H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL);

    // Alignment and metadata block size
    int alignment = 16 * 1024 * 1024;
    int blocksize =  4 * 1024 * 1024;
    H5Pset_alignment(fapl, alignment, alignment);
    H5Pset_meta_block_size(fapl, blocksize);

    // Collective metadata ops
    H5Pset_coll_metadata_write(fapl, true);
    H5Pset_all_coll_metadata_ops(fapl, true);

    // Defer cache flush
    H5AC_cache_config_t cache_config;
    cache_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    H5Pget_mdc_config(fapl, &cache_config);
    cache_config.set_initial_size = 1;
    cache_config.initial_size = 16 * 1024 * 1024;
    cache_config.evictions_enabled = 0;
    cache_config.incr_mode = H5C_incr__off;
    cache_config.flash_incr_mode = H5C_flash_incr__off;
    cache_config.decr_mode = H5C_decr__off;
    H5Pset_mdc_config (fapl, &cache_config);
#else
    H5Pset_fapl_sec2(fapl);
#endif

}
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
void WriteHDF5ParticleDataSync (PC const& pc,
                                  const std::string& dir, const std::string& name,
                                  const Vector<int>& write_real_comp,
                                  const Vector<int>& write_int_comp,
                                  const Vector<std::string>& real_comp_names,
                                  const Vector<std::string>& int_comp_names,
                                  const std::string& compression,
                                  F&& f, bool is_checkpoint)
{
    BL_PROFILE("WriteHDF5ParticleDataSync()");
    AMREX_ASSERT(pc.OK());

    AMREX_ASSERT(sizeof(typename PC::ParticleType::RealType) == 4 ||
                 sizeof(typename PC::ParticleType::RealType) == 8);

    constexpr int NStructReal = PC::NStructReal;
    constexpr int NStructInt  = PC::NStructInt;

    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();

    AMREX_ALWAYS_ASSERT(real_comp_names.size() == pc.NumRealComps() + NStructReal);
    AMREX_ALWAYS_ASSERT( int_comp_names.size() == pc.NumIntComps() + NStructInt);

#ifdef AMREX_USE_HDF5_ASYNC
    // For HDF5 async VOL, block and wait previous tasks have all completed
    if (es_par_g != 0)
        async_vol_es_wait_particle();
    else {
        ExecOnFinalize(async_vol_es_wait_close_particle);
        es_par_g = H5EScreate();
    }
#endif

    std::string pdir = dir;
    if ( ! pdir.empty() && pdir[pdir.size()-1] != '/') pdir += '/';
    pdir += name;

    if ( ! pc.GetLevelDirectoriesCreated()) {
        if (ParallelDescriptor::IOProcessor())
        {
            if ( ! amrex::UtilCreateDirectory(pdir, 0755))
            {
                amrex::CreateDirectoryFailed(pdir);
            }
        }
        ParallelDescriptor::Barrier();
    }

    /* std::ofstream HdrFile; */
    hid_t fid, grp, fapl, comp_dtype;
    int status;

    Long nparticles = 0;
    Long maxnextid;

    comp_dtype = H5Tcreate (H5T_COMPOUND, 2 * AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_j", 3 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_k", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 3 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_j", 4 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_k", 5 * sizeof(int), H5T_NATIVE_INT);
    }

    // evaluate f for every particle to determine which ones to output
    Vector<std::map<std::pair<int, int>, typename PC::IntVector > >
        particle_io_flags(pc.GetParticles().size());
    for (int lev = 0; lev < pc.GetParticles().size();  lev++)
    {
        const auto& pmap = pc.GetParticles(lev);
        for (const auto& kv : pmap)
        {
            auto& flags = particle_io_flags[lev][kv.first];
            particle_detail::fillFlags(flags, kv.second, std::forward<F>(f));
        }
    }

    Gpu::streamSynchronize();

    if(pc.GetUsePrePost())
    {
        nparticles = pc.GetNParticlesPrePost();
        maxnextid  = pc.GetMaxNextIDPrePost();
    }
    else
    {
        nparticles = particle_detail::countFlags(particle_io_flags, pc);
        maxnextid  = PC::ParticleType::NextID();
        ParallelDescriptor::ReduceLongSum(nparticles, IOProcNumber);
        PC::ParticleType::NextID(maxnextid);
        ParallelDescriptor::ReduceLongMax(maxnextid, IOProcNumber);
    }

    std::string HDF5FileName = pdir;
    if ( ! HDF5FileName.empty() && HDF5FileName[HDF5FileName.size()-1] != '/')
        HDF5FileName += '/';

    HDF5FileName += name;
    HDF5FileName += ".h5";
    pc.HdrFileNamePrePost = HDF5FileName;

    BL_PROFILE_VAR("H5FileCreate", h5fg);

    if (ParallelDescriptor::IOProcessor())
    {
        int set_stripe = 0;
        char setstripe[1024];
        int stripe_count = 128;
        int stripe_size = 1;
        char *stripe_count_str = getenv("HDF5_STRIPE_COUNT");
        char *stripe_size_str  = getenv("HDF5_STRIPE_SIZE");
        if (stripe_count_str) {
            stripe_count = atoi(stripe_count_str);
            set_stripe = 1;
        }
        if (stripe_size_str) {
            stripe_size = atoi(stripe_size_str);
            set_stripe = 1;
        }
        if (set_stripe == 1) {
            snprintf(setstripe, sizeof setstripe, "lfs setstripe -c %d -S %dm %s", stripe_count, stripe_size, pdir.c_str());
            std::cout << "Setting stripe parameters for HDF5 output: " << setstripe << std::endl;
            amrex::ignore_unused(std::system(setstripe));
        }

        fid = H5Fcreate(HDF5FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (fid < 0) amrex::FileOpenFailed(HDF5FileName.c_str());

        //
        // First thing written is our Checkpoint/Restart version string.
        // We append "_single" or "_double" to the version string indicating
        // whether we're using "float" or "double" floating point data in the
        // particles so that we can Restart from the checkpoint files.
        //
        std::string versionName = is_checkpoint ? PC::CheckpointVersion() : PC::PlotfileVersion();
        if (sizeof(typename PC::ParticleType::RealType) == 4)
            versionName += "_single";
        else
            versionName += "_double";

        CreateWriteHDF5AttrString(fid, "version_name", versionName.c_str());

        int num_output_real = 0;
        for (int i = 0; i < pc.NumRealComps() + NStructReal; ++i)
            if (write_real_comp[i]) ++num_output_real;

        int num_output_int = 0;
        for (int i = 0; i < pc.NumIntComps() + NStructInt; ++i)
            if (write_int_comp[i]) ++num_output_int;

        // AMREX_SPACEDIM and N for sanity checking.
        int ndim = AMREX_SPACEDIM;
        grp = H5Gcreate(fid, "Chombo_global", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        CreateWriteHDF5Attr(grp, "SpaceDim", 1, &ndim, H5T_NATIVE_INT);
        H5Gclose(grp);

        // The number of extra real parameters
        int ncomp = num_output_real + num_output_int;
        CreateWriteHDF5Attr(fid, "num_component", 1, &ncomp, H5T_NATIVE_INT);
        CreateWriteHDF5Attr(fid, "num_component_real", 1, &num_output_real, H5T_NATIVE_INT);

        char comp_name[128];
        // Real component names
        for (int i = 0; i < NStructReal + pc.NumRealComps(); ++i ) {
            if (write_real_comp[i]) {
                /* HdrFile << real_comp_names[i] << '\n'; */
                snprintf(comp_name, sizeof comp_name, "real_component_%d", i);
                CreateWriteHDF5AttrString(fid, comp_name, real_comp_names[i].c_str());
            }
        }

        // The number of extra int parameters
        CreateWriteHDF5Attr(fid, "num_component_int", 1, &num_output_int, H5T_NATIVE_INT);

        // int component names
        for (int i = 0; i < NStructInt + pc.NumIntComps(); ++i ) {
            if (write_int_comp[i]) {
                snprintf(comp_name, sizeof comp_name, "int_component_%d", i);
                CreateWriteHDF5AttrString(fid, comp_name, int_comp_names[i].c_str());
            }
        }

        // The total number of particles.
        CreateWriteHDF5Attr(fid, "nparticles", 1, &nparticles, H5T_NATIVE_LLONG);

        // The value of nextid that we need to restore on restart.
        CreateWriteHDF5Attr(fid, "maxnextid", 1, &maxnextid, H5T_NATIVE_LLONG);

        // Then the finest level of the AMR hierarchy.
        int finest_level = pc.finestLevel();
        CreateWriteHDF5Attr(fid, "finest_level", 1, &finest_level, H5T_NATIVE_INT);

        char level_name[128];
        int ngrids;

        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);

        for (int lev = 0; lev <= finest_level; ++lev) {
            snprintf(level_name, sizeof level_name, "level_%d", lev);

            grp = H5Gcreate(fid, level_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (grp < 0) {
                std::cout << "H5Gcreate [" << level_name << "] failed!" << std::endl;
                continue;
            }

            // Then the number of grids at each level.
            ngrids = pc.ParticleBoxArray(lev).size();
            CreateWriteHDF5Attr(grp, "ngrids", 1, &ngrids, H5T_NATIVE_INT);

            /* int ratio = 1; */
            /* CreateWriteHDF5Attr(grp, "ref_ratio", 1, &ratio); */

            int mfs_size = 2 * AMREX_SPACEDIM;
            hsize_t mfs_dim = (hsize_t)ngrids;

            hid_t mfs_dset_space = H5Screate_simple(1, &mfs_dim, NULL);
            hid_t mfs_dset= H5Dcreate(grp, "boxes", comp_dtype, mfs_dset_space, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

            Vector<int> vbox(ngrids * mfs_size);
            for(int j = 0; j < pc.ParticleBoxArray(lev).size(); ++j) {
                for(int i = 0; i < AMREX_SPACEDIM; ++i) {
                    vbox[(j * mfs_size) + i] = pc.ParticleBoxArray(lev)[j].smallEnd(i);
                    vbox[(j * mfs_size) + i + AMREX_SPACEDIM] = pc.ParticleBoxArray(lev)[j].bigEnd(i);
                }
            }

            status = H5Dwrite(mfs_dset, comp_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vbox[0]));
            if (status < 0) {
                std::string msg("ParticleContainer::WriteHDF5ParticleDataSync(): unable to write boxes dataset");
                amrex::Abort(msg.c_str());
            }

            H5Sclose(mfs_dset_space);
            H5Dclose(mfs_dset);

            H5Gclose(grp);
        }
        H5Pclose(dcpl_id);
        H5Fclose(fid);
    }

    ParallelDescriptor::Barrier();
    BL_PROFILE_VAR_STOP(h5fg);

    // We want to write the data out in parallel.
    // We'll allow up to nOutFiles active writers at a time.
    int nOutFiles(-1);

    ParmParse pp("particles");
    pp.queryAdd("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) nOutFiles = NProcs;
    /* nOutFiles = std::max(1, std::min(nOutFiles,NProcs)); */
    pc.nOutFilesPrePost = nOutFiles;

    fapl = H5Pcreate (H5P_FILE_ACCESS);
#ifdef BL_USE_MPI
    SetHDF5fapl(fapl, ParallelDescriptor::Communicator());
#else
    SetHDF5fapl(fapl);
#endif

#ifdef AMREX_USE_HDF5_ASYNC
    fid = H5Fopen_async(HDF5FileName.c_str(), H5F_ACC_RDWR, fapl, es_par_g);
#else
    fid = H5Fopen(HDF5FileName.c_str(), H5F_ACC_RDWR, fapl);
#endif
    if (fid < 0)
        FileOpenFailed(HDF5FileName.c_str());

    char level_name[64];
    for (int lev = 0; lev <= pc.finestLevel(); lev++)
    {
        snprintf(level_name, sizeof level_name, "level_%d", lev);
#ifdef AMREX_USE_HDF5_ASYNC
        grp = H5Gopen_async(fid, level_name, H5P_DEFAULT, es_par_g);
#else
        grp = H5Gopen(fid, level_name, H5P_DEFAULT);
#endif

        bool gotsome;
        if(pc.usePrePost)
            gotsome = (pc.nParticlesAtLevelPrePost[lev] > 0);
        else
            gotsome = (pc.NumberOfParticlesAtLevel(lev) > 0);

        MFInfo info;
        info.SetAlloc(false);
        MultiFab state(pc.ParticleBoxArray(lev),
                       pc.ParticleDistributionMap(lev),
                       1,0,info);

        // We eventually want to write out the file name and the offset
        // into that file into which each grid of particles is written.
        Vector<int>  which(state.size(),0);
        Vector<int > count(state.size(),0);
        Vector<Long> where(state.size(),0);

        if(pc.usePrePost)
            pc.filePrefixPrePost[lev] = HDF5FileName;

        if (gotsome)
        {
            pc.WriteParticlesHDF5(lev, grp, which, count, where,
                                  write_real_comp, write_int_comp,
                                  compression,
                                  particle_io_flags, is_checkpoint);

            if(pc.usePrePost) {
                pc.whichPrePost[lev] = which;
                pc.countPrePost[lev] = count;
                pc.wherePrePost[lev] = where;
            } else {
                ParallelDescriptor::ReduceIntSum (which.dataPtr(), which.size(), IOProcNumber);
                ParallelDescriptor::ReduceIntSum (count.dataPtr(), count.size(), IOProcNumber);
                ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProcNumber);
            }
        }

#ifdef AMREX_USE_HDF5_ASYNC
        H5Gclose_async(grp, es_par_g);
#else
        H5Gclose(grp);
#endif
    } // end for (lev...)

    H5Tclose(comp_dtype);
    H5Pclose(fapl);

#ifdef AMREX_USE_HDF5_ASYNC
    H5Fclose_async(fid, es_par_g);
#else
    H5Fclose(fid);
#endif

    return;
} // End WriteHDF5ParticleDataSync

#endif /*AMREX_WRITE_BINARY_PARTICLE_DATA_HDF5_H*/
