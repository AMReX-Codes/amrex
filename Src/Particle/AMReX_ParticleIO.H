#ifndef AMREX_PARTICLEIO_H
#define AMREX_PARTICLEIO_H

#ifdef AMREX_USE_HDF5
#include <hdf5.h>
#endif

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WriteParticleRealData (void* data, size_t size,
                         std::ostream& os, const RealDescriptor& rd) const
{
    if (sizeof(typename ParticleType::RealType) == 4) {
        writeFloatData((float*) data, size, os, ParticleRealDescriptor);
    } 
    else if (sizeof(typename ParticleType::RealType) == 8) {
        writeDoubleData((double*) data, size, os, ParticleRealDescriptor);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::ReadParticleRealData (void* data, size_t size,
                         std::istream& is, const RealDescriptor& rd)
{
    if (sizeof(typename ParticleType::RealType) == 4) {
        readFloatData((float*) data, size, is, ParticleRealDescriptor);
    } 
    else if (sizeof(typename ParticleType::RealType) == 8) {
        readDoubleData((double*) data, size, is, ParticleRealDescriptor);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Checkpoint (const std::string& dir,
              const std::string& name, bool is_checkpoint,
              const Vector<std::string>& real_comp_names,
              const Vector<std::string>& int_comp_names) const
{
    Vector<int> write_real_comp;
    Vector<std::string> tmp_real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        if (real_comp_names.size() == 0)
        {
            std::stringstream ss;
            ss << "real_comp" << i;
            tmp_real_comp_names.push_back(ss.str());
        }
        else
        {
            tmp_real_comp_names.push_back(real_comp_names[i]);
        }
    }
    
    Vector<int> write_int_comp;
    Vector<std::string> tmp_int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        if (int_comp_names.size() == 0)
        {
            std::stringstream ss;
            ss << "int_comp" << i;
            tmp_int_comp_names.push_back(ss.str());
        }
        else
        {
            tmp_int_comp_names.push_back(int_comp_names[i]);
        }
    }

    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            tmp_real_comp_names, tmp_int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p) -> int
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Checkpoint (const std::string& dir, const std::string& name) const
{
    Vector<int> write_real_comp;
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }
    
    Vector<int> write_int_comp;
    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }
    
    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p) -> int
                            {
                                return p.id() > 0;
                            });
}

#ifdef AMREX_USE_HDF5
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::CheckpointHDF5 (const std::string& dir,
              const std::string& name, bool is_checkpoint,
              const Vector<std::string>& real_comp_names,
              const Vector<std::string>& int_comp_names) const
{
    Vector<int> write_real_comp;
    Vector<std::string> tmp_real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        if (real_comp_names.size() == 0)
        {
            std::stringstream ss;
            ss << "real_comp" << i;
            tmp_real_comp_names.push_back(ss.str());
        }
        else
        {
            tmp_real_comp_names.push_back(real_comp_names[i]);
        }
    }
    
    Vector<int> write_int_comp;
    Vector<std::string> tmp_int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        if (int_comp_names.size() == 0)
        {
            std::stringstream ss;
            ss << "int_comp" << i;
            tmp_int_comp_names.push_back(ss.str());
        }
        else
        {
            tmp_int_comp_names.push_back(int_comp_names[i]);
        }
    }

    WriteHDF5ParticleData(dir, name, write_real_comp, write_int_comp,
                            tmp_real_comp_names, tmp_int_comp_names);
    
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::CheckpointHDF5 (const std::string& dir, const std::string& name) const
{
    Vector<int> write_real_comp;
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }
    
    Vector<int> write_int_comp;
    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }
    
    WriteHDF5ParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names);
}

static int CreateWriteHDF5AttrDouble(hid_t loc, const char *name, hsize_t n, const double *data)
{
    herr_t ret;
    hid_t attr, attr_space;
    hsize_t dims = n;

    attr_space = H5Screate_simple(1, &dims, NULL);

    attr = H5Acreate(loc, name, H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        amrex::Print() << " Error with H5Acreate [" << name << "]\n";
        return -1;
    }

    ret  = H5Awrite(attr, H5T_NATIVE_DOUBLE, (void*)data);
    if (ret < 0) {
        amrex::Print() << " Error with H5Awrite [" << name << "]\n";
        return -1;
    }
    H5Sclose(attr_space);
    H5Aclose(attr);
    return 1;
}

static int CreateWriteHDF5AttrLong(hid_t loc, const char *name, hsize_t n, const Long *data)
{
    herr_t ret;
    hid_t attr, attr_space;
    hsize_t dims = n;

    attr_space = H5Screate_simple(1, &dims, NULL);

    attr = H5Acreate(loc, name, H5T_NATIVE_LLONG, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        amrex::Print() << " Error with H5Acreate [" << name << "]\n";
        return -1;
    }

    ret  = H5Awrite(attr, H5T_NATIVE_LLONG, (void*)data);
    if (ret < 0) {
        amrex::Print() << " Error with H5Awrite [" << name << "]\n";
        return -1;
    }
    H5Sclose(attr_space);
    H5Aclose(attr);
    return 1;
}


static int CreateWriteHDF5AttrInt(hid_t loc, const char *name, hsize_t n, const int *data)
{
    herr_t ret;
    hid_t attr, attr_space;
    hsize_t dims = n;

    attr_space = H5Screate_simple(1, &dims, NULL);

    attr = H5Acreate(loc, name, H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        amrex::Print() << " Error with H5Acreate [" << name << "]\n";
        return -1;
    }

    ret  = H5Awrite(attr, H5T_NATIVE_INT, (void*)data);
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

static int CreateWriteDsetDouble(hid_t loc, const char *name, hsize_t n, const double *data)
{
    herr_t ret;
    hid_t dset, dset_space;
    hsize_t dims = n;

    dset_space = H5Screate_simple(1, &dims, NULL);

    dset = H5Dcreate(loc, name, H5T_NATIVE_DOUBLE, dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        amrex::Print() << " Error with H5Dcreate [" << name << "]\n";
        return -1;
    }

    ret  = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)data);
    if (ret < 0) {
        amrex::Print() << " Error with H5Dwrite [" << name << "]\n";
        return -1;
    }
    H5Sclose(dset_space);
    H5Aclose(dset);
    return 1;
}

static int ReadHDF5AttrDouble(hid_t loc, const char *name, void *data)
{
    herr_t ret;
    hid_t attr;

    attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr < 0) {
        amrex::Print() << " Error with H5Aopen [" << name << "]\n";
        return -1;
    }

    ret  = H5Aread(attr, H5T_NATIVE_DOUBLE, data);
    if (ret < 0) {
        amrex::Print() << " Error with H5Aread [" << name << "]\n";
        return -1;
    }
    H5Aclose(attr);
    return 1;
}

static int ReadHDF5AttrInt(hid_t loc, const char *name, void *data)
{
    herr_t ret;
    hid_t attr;

    attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr < 0) {
        amrex::Print() << " Error with H5Aopen [" << name << "]\n";
        return -1;
    }

    ret  = H5Aread(attr, H5T_NATIVE_INT, data);
    if (ret < 0) {
        amrex::Print() << " Error with H5Aread [" << name << "]\n";
        return -1;
    }
    H5Aclose(attr);
    return 1;
}

static int ReadHDF5AttrLong(hid_t loc, const char *name, void *data)
{
    herr_t ret;
    hid_t attr;

    attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr < 0) {
        amrex::Print() << " Error with H5Aopen [" << name << "]\n";
        return -1;
    }

    ret  = H5Aread(attr, H5T_NATIVE_LLONG, data);
    if (ret < 0) {
        amrex::Print() << " Error with H5Aread [" << name << "]\n";
        return -1;
    }
    H5Aclose(attr);
    return 1;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WriteHDF5ParticleData (const std::string& dir, const std::string& name,
                           const Vector<int>& write_real_comp,
                           const Vector<int>& write_int_comp,
                           const Vector<std::string>& real_comp_names,
                           const Vector<std::string>& int_comp_names) const
{
    BL_PROFILE("ParticleContainer::WriteHDF5ParticleData()");
    BL_ASSERT(OK());
    
    BL_ASSERT(sizeof(typename ParticleType::RealType) == 4 ||
              sizeof(typename ParticleType::RealType) == 8);
    
    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = amrex::second();
    
    AMREX_ALWAYS_ASSERT(real_comp_names.size() == NumRealComps() + NStructReal);
    AMREX_ALWAYS_ASSERT( int_comp_names.size() == NumIntComps() + NStructInt);

    std::string pdir = dir;
    if ( not pdir.empty() and pdir[pdir.size()-1] != '/') pdir += '/';
    
    if ( ! levelDirectoriesCreated)
    {
        if (ParallelDescriptor::IOProcessor()) {
            if ( ! amrex::UtilCreateDirectory(pdir, 0755)) 
                amrex::CreateDirectoryFailed(pdir);

            int set_stripe = 0;
            char setstripe[1024];
            int stripe_count = 128;
            int stripe_size = 32;
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
                sprintf(setstripe, "lfs setstripe -c %d -S %dm %s", stripe_count, stripe_size, pdir.c_str());
                std::cout << "Setting stripe parameters for HDF5 output: " << setstripe << std::endl; 
                amrex::ignore_unused(std::system(setstripe));
            }
        }
        ParallelDescriptor::Barrier();
    }

    Long nparticles = 0;
    int maxnextid;
    
    if(usePrePost)
    {
        nparticles = nparticlesPrePost;
        maxnextid  = maxnextidPrePost;
    } else
    {
        nparticles = 0;
        maxnextid  = ParticleType::NextID();
        
        for (int lev = 0; lev < m_particles.size();  lev++)
        {
            const auto& pmap = m_particles[lev];
            for (const auto& kv : pmap)
            {
                const auto& aos = kv.second.GetArrayOfStructs();
                for (int k = 0; k < aos.size(); ++k) 
                {
                    // Only count (and checkpoint) valid particles.
                    const ParticleType& p = aos[k];
                    if (p.m_idata.id > 0) nparticles++;
                }
            }
        }
        
        ParallelDescriptor::ReduceLongSum(nparticles, IOProcNumber);
        ParticleType::NextID(maxnextid);
        ParallelDescriptor::ReduceIntMax(maxnextid, IOProcNumber);
    }

    hid_t fapl, fid, grp;
    int status;

    hid_t comp_dtype = H5Tcreate (H5T_COMPOUND, 2 * AMREX_SPACEDIM * sizeof(int));
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


    std::string HDF5FileName = pdir;
    if ( ! HDF5FileName.empty() && HDF5FileName[HDF5FileName.size()-1] != '/')
        HDF5FileName += '/';
    
    HDF5FileName += name;
    if (std::string::npos == HDF5FileName.find("h5") && std::string::npos == HDF5FileName.find("hdf5")) 
        HDF5FileName += ".h5";
    HdrFileNamePrePost = HDF5FileName;
	
    if (ParallelDescriptor::IOProcessor()) {
   
        // Have only one rank to create and write metadata (header)
        fapl = H5Pcreate (H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(fapl, MPI_COMM_SELF, MPI_INFO_NULL);

        // Create the HDF5 file
        fid = H5Fcreate(HDF5FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        if (fid < 0) 
            amrex::FileOpenFailed(HDF5FileName.c_str());

        //
        // First thing written is our Checkpoint/Restart version string.
        // We append "_single" or "_double" to the version string indicating
        // whether we're using "float" or "double" floating point data in the
        // particles so that we can Restart from the checkpoint files.
        //
        std::string versionName = ParticleType::Version();
        if (sizeof(typename ParticleType::RealType) == 4)
            versionName += "_single";
        else
            versionName += "_double";

        CreateWriteHDF5AttrString(fid, "version_name", versionName.c_str());

        int num_output_real = 0;
        for (int i = 0; i < NumRealComps() + NStructReal; ++i)
            if (write_real_comp[i]) ++num_output_real;
        
        int num_output_int = 0;
        for (int i = 0; i < NumIntComps() + NStructInt; ++i)
            if (write_int_comp[i]) ++num_output_int;
        
        // AMREX_SPACEDIM and N for sanity checking.
        int ndim = AMREX_SPACEDIM;
        grp = H5Gcreate(fid, "Chombo_global", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        CreateWriteHDF5AttrInt(grp, "SpaceDim", 1, &ndim);
        H5Gclose(grp);
	
	// The number of extra real parameters
        int ncomp = num_output_real + num_output_int;
        CreateWriteHDF5AttrInt(fid, "num_component", 1, &ncomp);
        CreateWriteHDF5AttrInt(fid, "num_component_real", 1, &num_output_real);
        
        char comp_name[128];
        // Real component names
        for (int i = 0; i < NStructReal + NumRealComps(); ++i ) {
            if (write_real_comp[i]) {
                /* HdrFile << real_comp_names[i] << '\n'; */
                sprintf(comp_name, "real_component_%d", i);
                CreateWriteHDF5AttrString(fid, comp_name, real_comp_names[i].c_str());
            }
        }
        
	// The number of extra int parameters
        CreateWriteHDF5AttrInt(fid, "num_component_int", 1, &num_output_int);
        
        // int component names
        for (int i = 0; i < NStructInt + NumIntComps(); ++i ) {
            if (write_int_comp[i]) {
                /* HdrFile << int_comp_names[i] << '\n'; */
                sprintf(comp_name, "int_component_%d", i);
                CreateWriteHDF5AttrString(fid, comp_name, int_comp_names[i].c_str());
            }
        }


        // The total number of particles.
        CreateWriteHDF5AttrLong(fid, "nparticles", 1, &nparticles);

        // The value of nextid that we need to restore on restart.
        CreateWriteHDF5AttrInt(fid, "maxnextid", 1, &maxnextid);

        // Then the finest level of the AMR hierarchy.
        int finest_level = finestLevel();
        CreateWriteHDF5AttrInt(fid, "finest_level", 1, &finest_level);


        char level_name[128];
        int ngrids;
        for (int lev = 0; lev <= finest_level; ++lev) {
            sprintf(level_name, "level_%d", lev);

            grp = H5Gcreate(fid, level_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (grp < 0) {
                std::cout << "H5Gcreate [" << level_name << "] failed!" << std::endl;
                continue;
            }

            // Then the number of grids at each level.
            ngrids = ParticleBoxArray(lev).size();
            CreateWriteHDF5AttrInt(grp, "ngrids", 1, &ngrids);

            /* int ratio = 1; */
            /* CreateWriteHDF5AttrInt(grp, "ref_ratio", 1, &ratio); */

            int mfs_size = 2 * AMREX_SPACEDIM;
            hsize_t mfs_dim = (hsize_t)ngrids;
        
            hid_t mfs_dset_space = H5Screate_simple(1, &mfs_dim, NULL);
            hid_t mfs_dset= H5Dcreate(grp, "boxes", comp_dtype, mfs_dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            Vector<int> vbox(ngrids * mfs_size);
            for(int j = 0; j < ParticleBoxArray(lev).size(); ++j) {
                for(int i = 0; i < AMREX_SPACEDIM; ++i) {
                    vbox[(j * mfs_size) + i] = ParticleBoxArray(lev)[j].smallEnd(i);
                    vbox[(j * mfs_size) + i + AMREX_SPACEDIM] = ParticleBoxArray(lev)[j].bigEnd(i);
                }
            }

            status = H5Dwrite(mfs_dset, comp_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vbox[0]));
            if (status < 0) {
                std::string msg("ParticleContainer::CheckpointHDF5(): unable to write mfs dataset");
                amrex::Abort(msg.c_str());
            }

            H5Sclose(mfs_dset_space);
            H5Dclose(mfs_dset);

            H5Gclose(grp);
        }

        H5Pclose(fapl);
        H5Fclose(fid);
    }

    ParallelDescriptor::Barrier();

    fapl = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl,  ParallelDescriptor::Communicator(), MPI_INFO_NULL);

    // All process open the file
    fid = H5Fopen(HDF5FileName.c_str(), H5F_ACC_RDWR, fapl);
    if (fid < 0) 
        FileOpenFailed(HDF5FileName.c_str());

    char level_name[32];
    for (int lev = 0; lev <= finestLevel(); lev++)
    {
        sprintf(level_name, "level_%d", lev);
        grp = H5Gopen(fid, level_name, H5P_DEFAULT);

        bool gotsome;
	if(usePrePost)
            gotsome = (nParticlesAtLevelPrePost[lev] > 0);
        else
            gotsome = (NumberOfParticlesAtLevel(lev) > 0);


	MFInfo info;
	info.SetAlloc(false);
	MultiFab state(ParticleBoxArray(lev), ParticleDistributionMap(lev), 1, 0, info);

        // We eventually want to write out the file name and the offset
        // into that file into which each grid of particles is written.
        Vector<int>  which(state.size(),0);
        Vector<int > count(state.size(),0);
        Vector<Long> where(state.size(),0);
	

	bool groupSets(false), setBuf(true);
        
        if (gotsome)
	{
            WriteParticlesHDF5(grp, lev, count, where);
            
	    if(usePrePost) {
                whichPrePost[lev] = which;
                countPrePost[lev] = count;
                wherePrePost[lev] = where;
	    } else {
                ParallelDescriptor::ReduceIntSum (which.dataPtr(), which.size(), IOProcNumber);
                ParallelDescriptor::ReduceIntSum (count.dataPtr(), count.size(), IOProcNumber);
                ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProcNumber);
	    }
        }

        H5Gclose(grp);
    }            // ---- end for(lev...)

    H5Tclose(comp_dtype);
    H5Pclose(fapl);
    H5Fclose(fid);
   
    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;
        
        ParallelDescriptor::ReduceRealMax(stoptime, IOProcNumber);
        
        amrex::Print() << "ParticleContainer::CheckpointHDF5() time: " << stoptime << '\n';
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WriteParticlesHDF5 ( hid_t grp, int lev, Vector<int>& count, Vector<Long>& where) const
{
    BL_PROFILE("ParticleContainer::WriteParticlesHDF5()");

    std::map<int, Vector<int> > tile_map;

    for (const auto& kv : m_particles[lev])
    {
        const int grid = kv.first.first;
        const int tile = kv.first.second;
        tile_map[grid].push_back(tile);

        // Only write out valid particles.
        int cnt = 0;	
	for (int k = 0; k < kv.second.GetArrayOfStructs().size(); ++k)
	{
	    const ParticleType& p = kv.second.GetArrayOfStructs()[k];
  	    if (p.m_idata.id > 0) {
                cnt++;
	    }	    
	}

        count[grid] += cnt;
    }
	
    MFInfo info;
    info.SetAlloc(false);
    MultiFab state(ParticleBoxArray(lev), ParticleDistributionMap(lev), 1, 0, info);

    int            my_mfi_cnt = 0;
    ULong          my_mfi_int_total_size = 0, my_mfi_real_total_size = 0, int_size, real_size;
    Vector<int>    all_mfi_cnt(ParallelDescriptor::NProcs());
    Vector<int>    my_mfi_real_size;
    Vector<int>    my_mfi_int_size;
    Vector<int>    my_nparticles;
    Vector<ULong>  all_mfi_real_total_size(ParallelDescriptor::NProcs());
    Vector<ULong>  all_mfi_int_total_size(ParallelDescriptor::NProcs());
    hid_t          offset_id, offset_space, real_mem_space, real_dset_space, real_dset_id;
    hid_t          int_mem_space, int_dset_id, int_dset_space;
    hsize_t        total_mfi = 0, total_real_size = 0, total_int_size = 0, real_file_offset = 0, int_file_offset = 0;
    hsize_t        my_int_offset, my_int_count, my_real_offset, my_real_count;

    // Get the size for each mf so we know the amount of data from each rank
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const int grid = mfi.index();
        if (count[grid] == 0) 
            continue;

        int_size  = count[grid] * (2 + NStructInt + NArrayInt);
        my_mfi_int_size.push_back(int_size);
        my_nparticles.push_back(count[grid]);
        my_mfi_int_total_size += int_size;

        real_size = count[grid] * (AMREX_SPACEDIM + NStructReal + NArrayReal);
        my_mfi_real_size.push_back(real_size);
        my_mfi_real_total_size += real_size;
        my_mfi_cnt++;
    }

    // Collect the number of mf and total size of mf from each rank
    MPI_Allgather(&my_mfi_cnt, 1, MPI_INT, &(all_mfi_cnt[0]), 1, MPI_INT, ParallelDescriptor::Communicator());
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++) 
        total_mfi += all_mfi_cnt[i];

    // Create the int data 
    MPI_Allgather(&my_mfi_int_total_size, 1, Mpi_typemap<ULong>::type(), 
                  &(all_mfi_int_total_size[0]), 1, Mpi_typemap<ULong>::type(), ParallelDescriptor::Communicator());
    
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++) 
        total_int_size += all_mfi_int_total_size[i]; 

    int_dset_space = H5Screate_simple(1, &total_int_size, NULL); 
    int_dset_id  = H5Dcreate(grp, "data:datatype=0", H5T_NATIVE_INT, int_dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(int_dset_space);

    int_file_offset = 0;
    for (int i = 0; i < ParallelDescriptor::MyProc(); i++) 
        int_file_offset += all_mfi_int_total_size[i];
    my_int_offset = int_file_offset;
    my_int_count  = 0;
    

    // Create the real data 
    MPI_Allgather(&my_mfi_real_total_size, 1, Mpi_typemap<ULong>::type(), 
                  &(all_mfi_real_total_size[0]), 1, Mpi_typemap<ULong>::type(), ParallelDescriptor::Communicator());

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++) 
        total_real_size += all_mfi_real_total_size[i];

    real_dset_space = H5Screate_simple(1, &total_real_size, NULL); 
    if (sizeof(typename ParticleType::RealType) == 4) 
        real_dset_id  = H5Dcreate(grp, "data:datatype=1", H5T_NATIVE_FLOAT, real_dset_space, 
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    else 
        real_dset_id  = H5Dcreate(grp, "data:datatype=1", H5T_NATIVE_DOUBLE, real_dset_space, 
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(real_dset_space);

    real_file_offset = 0;
    for (int i = 0; i < ParallelDescriptor::MyProc(); i++) 
        real_file_offset += all_mfi_real_total_size[i];
    my_real_offset = real_file_offset;
    my_real_count  = 0;

    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const int grid = mfi.index();
        
        if (count[grid] == 0) 
            continue;
        
        // First write out the integer data in binary.
        const int      iChunkSize = 2 + NStructInt + NArrayInt;
        Vector<int> istuff(count[grid]*iChunkSize);
        int* iptr = istuff.dataPtr();

        for (unsigned i = 0; i < tile_map[grid].size(); i++) {
            const auto& pbox = m_particles[lev].at(std::make_pair(grid, tile_map[grid][i]));
            for (int pindex = 0; pindex < pbox.GetArrayOfStructs().size(); ++pindex) {
                const ParticleType& p = pbox.GetArrayOfStructs()[pindex];
                if (p.m_idata.id > 0) {
                    for (int j = 0; j < 2 + NStructInt; j++) {
                        iptr[j] = p.m_idata.arr[j];
                    }
                    iptr += 2 + NStructInt;
                    const auto& soa  = pbox.GetStructOfArrays();
                    for (int j = 0; j < NArrayInt; j++) {
                        iptr[j] = soa.GetIntData(j)[pindex];
                    }
                    iptr += NArrayInt;
                }
            }
        }

        /* writeIntData(istuff.dataPtr(), istuff.size(), ofs); */
        my_int_count = istuff.size();
        int_mem_space = H5Screate_simple(1, &my_int_count, NULL); 
        /* std::cout << "Rank " << ParallelDescriptor::MyProc() << ": my_int_offset = " << */ 
        /*         my_int_offset << ", my_int_count = " << my_int_count << ", total_int_size = " << total_int_size << '\n'; */
        int_dset_space = H5Screate_simple(1, &total_int_size, NULL); 
        H5Sselect_hyperslab (int_dset_space, H5S_SELECT_SET, &my_int_offset, NULL, &my_int_count, NULL);
        H5Dwrite(int_dset_id, H5T_NATIVE_INT, int_mem_space, int_dset_space, H5P_DEFAULT, istuff.dataPtr());
        H5Sclose(int_mem_space);
        H5Sclose(int_dset_space);

        my_int_offset += my_int_count;

        
        // Write the Real data in binary.
        const int      rChunkSize = AMREX_SPACEDIM + NStructReal + NArrayReal;
        Vector<typename ParticleType::RealType> rstuff(count[grid]*rChunkSize);
        typename ParticleType::RealType* rptr = rstuff.dataPtr();
        
        for (unsigned i = 0; i < tile_map[grid].size(); i++) {
            const auto& pbox = m_particles[lev].at(std::make_pair(grid, tile_map[grid][i]));
            for (int pindex = 0; pindex < pbox.GetArrayOfStructs().size(); ++pindex) {
                const ParticleType& p = pbox.GetArrayOfStructs()[pindex];
                if (p.m_idata.id > 0) {
                    for (int j = 0; j < AMREX_SPACEDIM + NStructReal; j++) {
                        rptr[j] = p.m_rdata.arr[j];
                    }
                    rptr += AMREX_SPACEDIM + NStructReal;
                    const auto& soa  = pbox.GetStructOfArrays();
                    for (int j = 0; j < NArrayReal; j++) {
                        rptr[j] = (typename ParticleType::RealType) soa.GetRealData(j)[pindex];
                    }
                    rptr += NArrayReal;
                }
            }
        }

        my_real_count = rstuff.size();
        real_mem_space = H5Screate_simple(1, &my_real_count, NULL); 
        /* std::cout << "Rank " << ParallelDescriptor::MyProc() << ": my_real_offset = " << */ 
        /*         my_real_offset << ", my_real_count = " << my_real_count << ", total_real_size = " << total_real_size << '\n'; */
        real_dset_space = H5Screate_simple(1, &total_real_size, NULL); 
        H5Sselect_hyperslab (real_dset_space, H5S_SELECT_SET, &my_real_offset, NULL, &my_real_count, NULL);
        H5Dwrite(real_dset_id, H5T_NATIVE_DOUBLE, real_mem_space, real_dset_space, H5P_DEFAULT, rstuff.dataPtr());
        H5Sclose(real_mem_space);
        H5Sclose(real_dset_space);

        my_real_offset += my_real_count;
        /* WriteParticleRealData(rstuff.dataPtr(), rstuff.size(), ofs, ParticleRealDescriptor); */
    } // End for mfi(state)

    H5Dclose(real_dset_id);

    H5Dclose(int_dset_id);

    // Create and write the size dataset
    offset_space = H5Screate_simple(1, &total_mfi, NULL); 
    offset_id  = H5Dcreate(grp, "nparticles_grid", H5T_NATIVE_INT, offset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    my_int_offset = 0;
    for (int i = 0; i < ParallelDescriptor::MyProc(); i++) 
        my_int_offset += all_mfi_cnt[i];
    my_int_count  = my_mfi_cnt;
    int_mem_space = H5Screate_simple(1, &my_int_count, NULL); 
    /* std::cout << "Rank " << ParallelDescriptor::MyProc() << ": my_int_offset = " << */ 
    /*         my_int_offset << ", my_int_count = " << my_int_count << ", total_mfi = " << total_mfi << '\n'; */
    H5Sselect_hyperslab (offset_space, H5S_SELECT_SET, &my_int_offset, NULL, &my_int_count, NULL);

    H5Dwrite(offset_id, H5T_NATIVE_INT, int_mem_space, offset_space, H5P_DEFAULT, &(my_nparticles[0]));

    H5Sclose(int_mem_space);
    H5Dclose(offset_id);
    H5Sclose(offset_space);

} // End WriteParticlesHDF5

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::RestartHDF5 (const std::string& dir, const std::string& file, bool is_checkpoint)
{
    Restart(dir, file);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::RestartHDF5 (const std::string& dir, const std::string& file)
{
    BL_PROFILE("ParticleContainer::RestartHDF5()");
    AMREX_ASSERT(!dir.empty());
    AMREX_ASSERT(!file.empty());
    
    const Real strttime = amrex::second();
    
    /* int DATA_Digits_Read(5); */
    /* ParmParse pp("particles"); */
    /* pp.query("datadigits_read",DATA_Digits_Read); */

    std::string fullname = dir;
    if (!fullname.empty() && fullname[fullname.size()-1] != '/')
        fullname += '/';
    fullname += file;

    if (std::string::npos == fullname.find("h5") && std::string::npos == fullname.find("hdf5")) 
        fullname += ".h5";

    hid_t fid, dset, grp, fapl, attr, atype, dspace, int_dset, real_dset;
    int ret;

    fapl = H5Pcreate (H5P_FILE_ACCESS);
    /* H5Pset_fapl_mpio(fapl, MPI_COMM_SELF, MPI_INFO_NULL); */
    H5Pset_fapl_mpio(fapl, ParallelDescriptor::Communicator(), MPI_INFO_NULL);
    fid = H5Fopen(fullname.c_str(), H5F_ACC_RDONLY, fapl);
    if (fid < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to open file: ");
        msg += fullname;
        amrex::Abort(msg.c_str());
    }

    attr = H5Aopen(fid, "version_name", H5P_DEFAULT);
    if (attr < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to open version attribute ");
        amrex::Abort(msg.c_str());
    }

    hsize_t attr_len = H5Aget_storage_size(attr);

    std::string version;
    version.resize(attr_len+1);

    atype = H5Aget_type(attr);
    ret = H5Aread(attr, atype, &version[0]);
    H5Tclose(atype);

    H5Aclose(attr);

    AMREX_ASSERT(!version.empty());

    // What do our version strings mean?
    // "Version_One_Dot_Zero" -- hard-wired to write out in double precision.
    // "Version_One_Dot_One" -- can write out either as either single or double precision.
    // Appended to the latter version string are either "_single" or "_double" to
    // indicate how the particles were written.
    // "Version_Two_Dot_Zero" -- this is the AMReX particle file format
    std::string how;
    if (version.find("Version_One_Dot_Zero") != std::string::npos) {
        how = "double";
    }
    else if (version.find("Version_One_Dot_One")  != std::string::npos or
             version.find("Version_Two_Dot_Zero") != std::string::npos) {
        if (version.find("_single") != std::string::npos) {
            how = "single";
        }
        else if (version.find("_double") != std::string::npos) {
            how = "double";
        }
        else {
            std::string msg("ParticleContainer::Restart(): bad version string: ");
            msg += version;
            amrex::Error(version.c_str());
        }
    }
    else {
        std::string msg("ParticleContainer::Restart(): unknown version string: ");
        msg += version;
        amrex::Abort(msg.c_str());
    }

    std::string gname = "Chombo_global";
    std::string aname = "SpaceDim";
    int dm;
    grp = H5Gopen(fid, gname.c_str(), H5P_DEFAULT);
    if (grp < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to open group ");
        msg += gname;
        amrex::Abort(msg.c_str());
    }
    ret = ReadHDF5AttrInt(grp, aname.c_str(), &dm);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    H5Gclose(grp);

    aname = "num_component_real";
    int nr;
    ret = ReadHDF5AttrInt(fid, aname.c_str(), &nr);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    if (nr != NStructReal + NumRealComps())
        amrex::Abort("ParticleContainer::Restart(): nr != NStructReal + NumRealComps()");
    
    aname = "num_component_int";
    int ni;
    ret = ReadHDF5AttrInt(fid, aname.c_str(), &ni);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    if (ni != NStructInt + NumIntComps())
        amrex::Abort("ParticleContainer::Restart(): ni != NStructInt");
    
    aname = "nparticles";
    Long nparticles;
    ret = ReadHDF5AttrLong(fid, aname.c_str(), &nparticles);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    AMREX_ASSERT(nparticles >= 0);
    
    aname = "maxnextid";
    int maxnextid;
    ret = ReadHDF5AttrInt(fid, aname.c_str(), &maxnextid);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    AMREX_ASSERT(maxnextid > 0);
    ParticleType::NextID(maxnextid);
    
    aname = "finest_level";
    int finest_level_in_file;
    ret = ReadHDF5AttrInt(fid, aname.c_str(), &finest_level_in_file);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    AMREX_ASSERT(finest_level_in_file >= 0);
    
    // Determine whether this is a dual-grid restart or not.
    Vector<BoxArray> particle_box_arrays(finest_level_in_file + 1);
    bool dual_grid = false;

    hid_t comp_dtype = H5Tcreate (H5T_COMPOUND, 2 * AMREX_SPACEDIM * sizeof(int));
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

    for (int lev = 0; lev <= finest_level_in_file; lev++)
    {
        if (lev > finestLevel())
        {
            dual_grid = true;
            break;
        }
 
        gname = "level_" + std::to_string(lev);
        grp = H5Gopen(fid, gname.c_str(), H5P_DEFAULT);
        if (grp < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open group ");
            msg += gname;
            amrex::Abort(msg.c_str());
        }
 
        dset = H5Dopen(grp, "boxes", H5P_DEFAULT);
        dspace = H5Dget_space(dset);
        hsize_t ngrid;
        H5Sget_simple_extent_dims(dspace, &ngrid, NULL);

        Vector<int> boxes(ngrid*AMREX_SPACEDIM*2);
        ret = H5Dread(dset, comp_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(boxes[0]));
        if (ret < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to read nparticles_grid dataset");
            amrex::Abort(msg.c_str());
        }
        H5Sclose(dspace);
        H5Dclose(dset);
        H5Gclose(grp);
      
        /* particle_box_arrays[lev].readFrom(phdr_file); */
        particle_box_arrays[lev] = ParticleBoxArray(lev);
        for (int i = 0; i < ngrid; i++) {
            Box tmp ({AMREX_D_DECL(boxes[i*AMREX_SPACEDIM], boxes[i*AMREX_SPACEDIM+1], boxes[i*AMREX_SPACEDIM+2])}, {AMREX_D_DECL(boxes[i*AMREX_SPACEDIM*2], boxes[i*AMREX_SPACEDIM*2+1], boxes[i*AMREX_SPACEDIM*2+2])});
            particle_box_arrays[lev][i] = tmp;
        }

        if (not particle_box_arrays[lev].CellEqual(ParticleBoxArray(lev))) 
            dual_grid = true;
    }
    if (dual_grid) {
        for (int lev = 0; lev <= finestLevel(); lev++) {
            SetParticleBoxArray(lev, particle_box_arrays[lev]);
            DistributionMapping pdm(particle_box_arrays[lev]);
            SetParticleDistributionMap(lev, pdm);
        }
    }
    
    Vector<int> ngrids(finest_level_in_file+1);
    for (int lev = 0; lev <= finest_level_in_file; lev++) {
        gname = "level_" + std::to_string(lev);
        grp = H5Gopen(fid, gname.c_str(), H5P_DEFAULT);
        if (grp < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open group ");
            msg += gname;
            amrex::Abort(msg.c_str());
        }

        aname = "ngrids";
        ret = ReadHDF5AttrInt(grp, aname.c_str(), &ngrids[lev]);
        if (ret < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
            msg += aname;
            amrex::Abort(msg.c_str());
        }

        AMREX_ASSERT(ngrids[lev] > 0);
        if (lev <= finestLevel()) {
            AMREX_ASSERT(ngrids[lev] == int(ParticleBoxArray(lev).size()));
        }

        H5Gclose(grp);
    }
    
    resizeData();
    
    if (finest_level_in_file > finestLevel()) {
        m_particles.resize(finest_level_in_file+1);
    }
    
    for (int lev = 0; lev <= finest_level_in_file; lev++) {

        gname = "level_" + std::to_string(lev);
        grp = H5Gopen(fid, gname.c_str(), H5P_DEFAULT);
        if (grp < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open group ");
            msg += gname;
            amrex::Abort(msg.c_str());
        }

        dset = H5Dopen(grp, "nparticles_grid", H5P_DEFAULT);
        Vector<int> count(ngrids[lev]);
        ret = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(count[0]));
        if (ret < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to read nparticles_grid dataset");
            amrex::Abort(msg.c_str());
        }
        H5Dclose(dset);

        Vector<hsize_t>  offset(ngrids[lev]);
        offset[0] = 0;
        for (int i = 1; i < ngrids[lev]; i++) {
            offset[i] = offset[i-1] + count[i-1];
        }

        int_dset  = H5Dopen(grp, "data:datatype=0", H5P_DEFAULT);
        if (int_dset < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open int dataset");
            amrex::Abort(msg.c_str());
        }
        real_dset = H5Dopen(grp, "data:datatype=1", H5P_DEFAULT);
        if (real_dset < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open int dataset");
            amrex::Abort(msg.c_str());
        }

        Vector<int> grids_to_read;
        if (lev <= finestLevel()) {
            for (MFIter mfi(*m_dummy_mf[lev]); mfi.isValid(); ++mfi) {
                grids_to_read.push_back(mfi.index());
            }
        } else {
            
            // we lost a level on restart. we still need to read in particles
            // on finer levels, and put them in the right place via Redistribute()
            
            const int rank = ParallelDescriptor::MyProc();
            const int NReaders = ParticleType::MaxReaders();
            if (rank >= NReaders) return;
            
            const int Navg = ngrids[lev] / NReaders;
            const int Nleft = ngrids[lev] - Navg * NReaders;
            
            int lo, hi;
            if (rank < Nleft) {
                lo = rank*(Navg + 1);
                hi = lo + Navg + 1;
            } 
            else {
                lo = rank * Navg + Nleft;
                hi = lo + Navg;
            }
            
            for (int i = lo; i < hi; ++i) {
                grids_to_read.push_back(i);
            }
        }
        
        for(int igrid = 0; igrid < static_cast<int>(grids_to_read.size()); ++igrid) {
            const int grid = grids_to_read[igrid];
            
            if (how == "single") {
                ReadParticlesHDF5<float>(offset[grid], count[grid], grid, lev, int_dset, real_dset, finest_level_in_file);
            }
            else if (how == "double") {
                ReadParticlesHDF5<double>(offset[grid], count[grid], grid, lev, int_dset, real_dset, finest_level_in_file);
            }
            else {
                std::string msg("ParticleContainer::Restart(): bad parameter: ");
                msg += how;
                amrex::Error(msg.c_str());
            }
        }

        H5Dclose(int_dset);
        H5Dclose(real_dset);
        H5Gclose(grp);
    } // end for level
    
    H5Fclose(fid);
    
    Redistribute();
    
    AMREX_ASSERT(OK());
    
    if (m_verbose > 1) {
        Real stoptime = amrex::second() - strttime;	
        ParallelDescriptor::ReduceRealMax(stoptime, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "ParticleContainer::Restart() time: " << stoptime << '\n';
    }
}

// Read a batch of particles from the checkpoint file
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class RTYPE>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::ReadParticlesHDF5 (hsize_t offset, hsize_t cnt, int grd, int lev, hid_t int_dset, hid_t real_dset, int finest_level_in_file)
{
    BL_PROFILE("ParticleContainer::ReadParticlesHDF5()");
    AMREX_ASSERT(cnt > 0);
    AMREX_ASSERT(lev < int(m_particles.size()));

    hid_t int_dspace, int_fspace, real_dspace, real_fspace;
    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    const int iChunkSize = 2 + NStructInt + NumIntComps();
    Vector<int> istuff(cnt*iChunkSize);
    int_fspace = H5Dget_space(int_dset);
    hsize_t int_cnt = cnt*iChunkSize;
    hsize_t int_offset = offset*iChunkSize;
    int_dspace = H5Screate_simple(1, &int_cnt, NULL); 
    H5Sselect_hyperslab (int_fspace, H5S_SELECT_SET, &int_offset, NULL, &int_cnt, NULL);
    H5Dread(int_dset, H5T_NATIVE_INT, int_dspace, int_fspace, H5P_DEFAULT, istuff.dataPtr());

    H5Sclose(int_fspace);
    H5Sclose(int_dspace);
    
    // Then the real data in binary.
    const int rChunkSize = AMREX_SPACEDIM + NStructReal + NumRealComps();
    Vector<RTYPE> rstuff(cnt*rChunkSize);
    real_fspace = H5Dget_space(real_dset);
    hsize_t real_cnt = cnt*rChunkSize;
    hsize_t real_offset = offset*rChunkSize;
    real_dspace = H5Screate_simple(1, &real_cnt, NULL); 
    H5Sselect_hyperslab (real_fspace, H5S_SELECT_SET, &real_offset, NULL, &real_cnt, NULL);
    H5Dread(real_dset, H5T_NATIVE_DOUBLE, real_dspace, real_fspace, H5P_DEFAULT, rstuff.dataPtr());

    H5Sclose(real_fspace);
    H5Sclose(real_dspace);
    
    // Now reassemble the particles.
    int*   iptr = istuff.dataPtr();
    RTYPE* rptr = rstuff.dataPtr();
    
    ParticleType p;
    ParticleLocData pld;

    Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
    host_particles.reserve(15);
    host_particles.resize(finest_level_in_file+1);

    Vector<std::map<std::pair<int, int>,
                    std::vector<Gpu::HostVector<Real> > > > host_real_attribs;
    host_real_attribs.reserve(15);
    host_real_attribs.resize(finest_level_in_file+1);

    Vector<std::map<std::pair<int, int>,
                    std::vector<Gpu::HostVector<int> > > > host_int_attribs;
    host_int_attribs.reserve(15);
    host_int_attribs.resize(finest_level_in_file+1);

    for (int i = 0; i < cnt; i++) {
        p.m_idata.id   = iptr[0];
        p.m_idata.cpu  = iptr[1];
        
        iptr += 2;
            
        for (int j = 0; j < NStructInt; j++)
        {
            p.m_idata.arr[2+j] = *iptr;
            ++iptr;
        }

        AMREX_ASSERT(p.m_idata.id > 0);
        
        AMREX_D_TERM(p.m_rdata.pos[0] = rptr[0];,
                     p.m_rdata.pos[1] = rptr[1];,
                     p.m_rdata.pos[2] = rptr[2];);
        
        rptr += AMREX_SPACEDIM;
        
        for (int j = 0; j < NStructReal; j++)
        {
            p.m_rdata.arr[AMREX_SPACEDIM+j] = *rptr;
            ++rptr;
        }

        locateParticle(p, pld, 0, finestLevel(), 0);
        
	std::pair<int, int> ind(grd, pld.m_tile);

        host_real_attribs[lev][ind].resize(NumRealComps());
        host_int_attribs[lev][ind].resize(NumIntComps());
        
	// add the struct
	host_particles[lev][ind].push_back(p);

	// add the real...
	for (int icomp = 0; icomp < NumRealComps(); icomp++) {
            host_real_attribs[lev][ind][icomp].push_back(*rptr);
            ++rptr;
	}
        
	// ... and int array data
	for (int icomp = 0; icomp < NumIntComps(); icomp++) {
            host_int_attribs[lev][ind][icomp].push_back(*iptr);
            ++iptr;
	}        
    }

    for (int host_lev = 0; host_lev < static_cast<int>(host_particles.size()); ++host_lev)
      {
	for (auto& kv : host_particles[host_lev]) {
	  auto grid = kv.first.first;
	  auto tile = kv.first.second;
	  const auto& src_tile = kv.second;
          
	  auto& dst_tile = DefineAndReturnParticleTile(host_lev, grid, tile);
	  auto old_size = dst_tile.GetArrayOfStructs().size();
	  auto new_size = old_size + src_tile.size();
	  dst_tile.resize(new_size);
                
	  Gpu::copy(Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
                    dst_tile.GetArrayOfStructs().begin() + old_size);
	  
	  for (int i = 0; i < NumRealComps(); ++i) {
              Gpu::copy(Gpu::hostToDevice,
                        host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                        host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                        dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
	  }
	  
	  for (int i = 0; i < NumIntComps(); ++i) {
              Gpu::copy(Gpu::hostToDevice,
                        host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                        host_int_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                        dst_tile.GetStructOfArrays().GetIntData(i).begin() + old_size);
	  }
	}
      }
    
    Gpu::streamSynchronize();
}

#endif

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir, const std::string& name) const
{
    Vector<int> write_real_comp;
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }
    
    Vector<int> write_int_comp;
    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }
    
    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names,
                 const Vector<std::string>& int_comp_names) const
{    
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());
    AMREX_ASSERT( int_comp_names.size() == NStructInt  + NumIntComps() );

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names) const
{    
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }
        
    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir,
                 const std::string& name,
                 const Vector<int>& write_real_comp,
                 const Vector<int>& write_int_comp) const
{    
    AMREX_ASSERT(write_real_comp.size() == NStructReal + NumRealComps());
    AMREX_ASSERT(write_int_comp.size()  == NStructInt  + NArrayInt );
    
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
WritePlotFile (const std::string& dir, const std::string& name,
               const Vector<int>& write_real_comp,
               const Vector<int>& write_int_comp,    
               const Vector<std::string>& real_comp_names,
               const Vector<std::string>&  int_comp_names) const
{
    BL_PROFILE("ParticleContainer::WritePlotFile()");
    
    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir, const std::string& name, F&& f) const
{
    Vector<int> write_real_comp;
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }
    
    Vector<int> write_int_comp;
    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }
    
    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names,
                 const Vector<std::string>& int_comp_names, F&& f) const
{    
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());
    AMREX_ASSERT( int_comp_names.size() == NStructInt  + NArrayInt );

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names, F&& f) const
{    
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }
        
    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir,
                 const std::string& name,
                 const Vector<int>& write_real_comp,
                 const Vector<int>& write_int_comp, F&& f) const
{    
    AMREX_ASSERT(write_real_comp.size() == NStructReal + NumRealComps());
    AMREX_ASSERT(write_int_comp.size()  == NStructInt  + NumIntComps() );
    
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
WritePlotFile (const std::string& dir, const std::string& name,
               const Vector<int>& write_real_comp,
               const Vector<int>& write_int_comp,    
               const Vector<std::string>& real_comp_names,
               const Vector<std::string>&  int_comp_names,
               F&& f) const
{
    BL_PROFILE("ParticleContainer::WritePlotFile()");
    
    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WriteBinaryParticleData (const std::string& dir, const std::string& name,
                           const Vector<int>& write_real_comp,
                           const Vector<int>& write_int_comp,
                           const Vector<std::string>& real_comp_names,
                           const Vector<std::string>& int_comp_names,
                           F&& f) const
{
    BL_PROFILE("ParticleContainer::WriteBinaryParticleData()");
    AMREX_ASSERT(OK());
    
    AMREX_ASSERT(sizeof(typename ParticleType::RealType) == 4 ||
                 sizeof(typename ParticleType::RealType) == 8);
    
    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = amrex::second();
    
    AMREX_ALWAYS_ASSERT(real_comp_names.size() == NumRealComps() + NStructReal);
    AMREX_ALWAYS_ASSERT( int_comp_names.size() == NumIntComps() + NStructInt);

    std::string pdir = dir;
    if ( not pdir.empty() and pdir[pdir.size()-1] != '/') pdir += '/';
    pdir += name;
    
    if ( ! levelDirectoriesCreated)
    {
        if (ParallelDescriptor::IOProcessor()) 
            if ( ! amrex::UtilCreateDirectory(pdir, 0755)) 
                amrex::CreateDirectoryFailed(pdir);
        ParallelDescriptor::Barrier();
    }
    
    std::ofstream HdrFile;
    
    Long nparticles = 0;
    int maxnextid;

    // evaluate f for every particle to determine which ones to output
    Vector<std::map<std::pair<int, int>, Gpu::DeviceVector<int> > > particle_io_flags(m_particles.size());		
    for (int lev = 0; lev < m_particles.size();  lev++)
    {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap)
        {
            const auto ptd = kv.second.getConstParticleTileData();
            const auto np = kv.second.numParticles();
            particle_io_flags[lev][kv.first].resize(np, 0);
            auto pflags = particle_io_flags[lev][kv.first].data();
            AMREX_HOST_DEVICE_FOR_1D( np, k, 
            {
                const auto p = ptd.getSuperParticle(k);
                pflags[k] = f(p);
            });
        }
    }

    Gpu::Device::synchronize();
    
    if(usePrePost)
    {
        nparticles = nparticlesPrePost;
        maxnextid  = maxnextidPrePost;
    }
    else
    {
        nparticles = 0;
        maxnextid  = ParticleType::NextID();

        for (int lev = 0; lev < m_particles.size();  lev++)
        {
            const auto& pmap = m_particles[lev];
            for (const auto& kv : pmap)
            {
                const auto& pflags = particle_io_flags[lev][kv.first];
                for (int k = 0; k < kv.second.numParticles(); ++k) 
                {
                    if (pflags[k]) nparticles++;
                }
            }
        }
        
        ParallelDescriptor::ReduceLongSum(nparticles, IOProcNumber);
        ParticleType::NextID(maxnextid);
        ParallelDescriptor::ReduceIntMax(maxnextid, IOProcNumber);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        std::string HdrFileName = pdir;
	
        if ( ! HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
            HdrFileName += '/';
        
        HdrFileName += "Header";
        HdrFileNamePrePost = HdrFileName;
	
        HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);
	
        if ( ! HdrFile.good()) amrex::FileOpenFailed(HdrFileName);

        //
        // First thing written is our Checkpoint/Restart version string.
        // We append "_single" or "_double" to the version string indicating
        // whether we're using "float" or "double" floating point data in the
        // particles so that we can Restart from the checkpoint files.
        //
        if (sizeof(typename ParticleType::RealType) == 4)
        {
            HdrFile << ParticleType::Version() << "_single" << '\n';
        }
        else
        {
            HdrFile << ParticleType::Version() << "_double" << '\n';
        }

        int num_output_real = 0;
        for (int i = 0; i < NumRealComps() + NStructReal; ++i)
            if (write_real_comp[i]) ++num_output_real;
        
        int num_output_int = 0;
        for (int i = 0; i < NumIntComps() + NStructInt; ++i)
            if (write_int_comp[i]) ++num_output_int;
        
        // AMREX_SPACEDIM and N for sanity checking.
        HdrFile << AMREX_SPACEDIM << '\n';
	
		// The number of extra real parameters
        HdrFile << num_output_real << '\n';
        
        // Real component names
        for (int i = 0; i < NStructReal + NumRealComps(); ++i )
            if (write_real_comp[i]) HdrFile << real_comp_names[i] << '\n';
        
		// The number of extra int parameters
        HdrFile << num_output_int << '\n';
        
        // int component names
        for (int i = 0; i < NStructInt + NumIntComps(); ++i )
            if (write_int_comp[i]) HdrFile << int_comp_names[i] << '\n';

        bool is_checkpoint = true; // legacy
        HdrFile << is_checkpoint << '\n';

        // The total number of particles.
        HdrFile << nparticles << '\n';

        // The value of nextid that we need to restore on restart.
        HdrFile << maxnextid << '\n';

        // Then the finest level of the AMR hierarchy.
        HdrFile << finestLevel() << '\n';

        // Then the number of grids at each level.
        for (int lev = 0; lev <= finestLevel(); lev++)
            HdrFile << ParticleBoxArray(lev).size() << '\n';
    }

    // We want to write the data out in parallel.
    // We'll allow up to nOutFiles active writers at a time.
    int nOutFiles(256);

    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) nOutFiles = NProcs;
    nOutFiles = std::max(1, std::min(nOutFiles,NProcs));
    nOutFilesPrePost = nOutFiles;

    for (int lev = 0; lev <= finestLevel(); lev++)
    {
        bool gotsome;
		if(usePrePost)
        {
            gotsome = (nParticlesAtLevelPrePost[lev] > 0);
		}
        else
        {
            gotsome = (NumberOfParticlesAtLevel(lev) > 0);
		}

        // We store the particles at each level in their own subdirectory.
        std::string LevelDir = pdir;
        
        if (gotsome)
        {
            if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') LevelDir += '/';
	    
            LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);
			
			if ( ! levelDirectoriesCreated) {
                if (ParallelDescriptor::IOProcessor()) 
                    if ( ! amrex::UtilCreateDirectory(LevelDir, 0755)) 
                        amrex::CreateDirectoryFailed(LevelDir);
                //
                // Force other processors to wait until directory is built.
                //
                ParallelDescriptor::Barrier();
			}
        }
	
        // Write out the header for each particle
        if (gotsome and ParallelDescriptor::IOProcessor()) {
            std::string HeaderFileName = LevelDir;
            HeaderFileName += "/Particle_H";
            std::ofstream ParticleHeader(HeaderFileName);
            
            ParticleBoxArray(lev).writeOn(ParticleHeader);
            ParticleHeader << '\n';
            
            ParticleHeader.flush();
            ParticleHeader.close();
        }
        
		MFInfo info;
		info.SetAlloc(false);
		MultiFab state(ParticleBoxArray(lev),
					   ParticleDistributionMap(lev),
					   1,0,info);

		// We eventually want to write out the file name and the offset
		// into that file into which each grid of particles is written.
		Vector<int>  which(state.size(),0);
		Vector<int > count(state.size(),0);
		Vector<Long> where(state.size(),0);
		
		std::string filePrefix(LevelDir);
		filePrefix += '/';
		filePrefix += ParticleType::DataPrefix();
		if(usePrePost) {
            filePrefixPrePost[lev] = filePrefix;
		}
		bool groupSets(false), setBuf(true);
        
		if (gotsome)
		{
			for(NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf); nfi.ReadyToWrite(); ++nfi)
			{
                            std::ofstream& myStream = (std::ofstream&) nfi.Stream();
                            WriteParticles(lev, myStream, nfi.FileNumber(), which, count, where,
                                           write_real_comp, write_int_comp, particle_io_flags);
			}
            
			if(usePrePost) {
				whichPrePost[lev] = which;
				countPrePost[lev] = count;
				wherePrePost[lev] = where;
			} else {
				ParallelDescriptor::ReduceIntSum (which.dataPtr(), which.size(), IOProcNumber);
				ParallelDescriptor::ReduceIntSum (count.dataPtr(), count.size(), IOProcNumber);
				ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProcNumber);
			}
		}
		
		if (ParallelDescriptor::IOProcessor())
		{
			if(usePrePost) {
				// ---- write to the header and unlink in CheckpointPost
			} else {
				for (int j = 0; j < state.size(); j++)
				{
					HdrFile << which[j] << ' ' << count[j] << ' ' << where[j] << '\n';
				}
				
				if (gotsome && doUnlink)
				{
					// Unlink any zero-length data files.
					Vector<Long> cnt(nOutFiles,0);
					
					for (int i = 0, N=count.size(); i < N; i++) {
						cnt[which[i]] += count[i];
					}
					
					for (int i = 0, N=cnt.size(); i < N; i++)
					{
						if (cnt[i] == 0)
						{
							std::string FullFileName = NFilesIter::FileName(i, filePrefix);
							amrex::UnlinkFile(FullFileName.c_str());
						}
					}
				}                
			}            
		}
    }
	
    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile.flush();
        HdrFile.close();
        if ( ! HdrFile.good())
        {
            amrex::Abort("ParticleContainer::Checkpoint(): problem writing HdrFile");
        }
    }
    
    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;        
        ParallelDescriptor::ReduceRealMax(stoptime, IOProcNumber);        
        amrex::Print() << "ParticleContainer::Checkpoint() time: " << stoptime << '\n';
    }
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::CheckpointPre ()
{
    if( ! usePrePost) {
        return;
    }
    
    BL_PROFILE("ParticleContainer::CheckpointPre()");
    
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    Long nparticles = 0;
    int  maxnextid  = ParticleType::NextID();
    
    for (int lev = 0; lev < m_particles.size();  lev++) {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
			for (int k = 0; k < aos.size(); ++k) {
				const ParticleType& p = aos[k];
                if (p.m_idata.id > 0) {
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
				}
            }
        }
    }
    ParallelDescriptor::ReduceLongSum(nparticles, IOProcNumber);

    ParticleType::NextID(maxnextid);
    ParallelDescriptor::ReduceIntMax(maxnextid, IOProcNumber);
    
    nparticlesPrePost = nparticles;
    maxnextidPrePost  = maxnextid;
    
    nParticlesAtLevelPrePost.clear();
    nParticlesAtLevelPrePost.resize(finestLevel() + 1, 0);
    for(int lev(0); lev <= finestLevel(); ++lev) {
        nParticlesAtLevelPrePost[lev] = NumberOfParticlesAtLevel(lev);
    }
    
    whichPrePost.clear();
    whichPrePost.resize(finestLevel() + 1);
    countPrePost.clear();
    countPrePost.resize(finestLevel() + 1);
    wherePrePost.clear();
    wherePrePost.resize(finestLevel() + 1);
    
    filePrefixPrePost.clear();
    filePrefixPrePost.resize(finestLevel() + 1);
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::CheckpointPost ()
{
    if( ! usePrePost) {
        return;
    }
    
    BL_PROFILE("ParticleContainer::CheckpointPost()");
    
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    std::ofstream HdrFile;
    HdrFile.open(HdrFileNamePrePost.c_str(), std::ios::out | std::ios::app);
    
    for(int lev(0); lev <= finestLevel(); ++lev) {
        ParallelDescriptor::ReduceIntSum (whichPrePost[lev].dataPtr(), whichPrePost[lev].size(), IOProcNumber);
        ParallelDescriptor::ReduceIntSum (countPrePost[lev].dataPtr(), countPrePost[lev].size(), IOProcNumber);
        ParallelDescriptor::ReduceLongSum(wherePrePost[lev].dataPtr(), wherePrePost[lev].size(), IOProcNumber);
        
        
        if(ParallelDescriptor::IOProcessor()) {
            for(int j(0); j < whichPrePost[lev].size(); ++j) {
                HdrFile << whichPrePost[lev][j] << ' ' << countPrePost[lev][j] << ' ' << wherePrePost[lev][j] << '\n';
            }
            
            const bool gotsome = (nParticlesAtLevelPrePost[lev] > 0);
            if(gotsome && doUnlink) {
//            BL_PROFILE_VAR("PC<NNNN>::Checkpoint:unlink", unlink_post);
                // Unlink any zero-length data files.
                Vector<Long> cnt(nOutFilesPrePost,0);
                
                for(int i(0), N = countPrePost[lev].size(); i < N; ++i) {
                    cnt[whichPrePost[lev][i]] += countPrePost[lev][i];
                }
                
                for(int i(0), N = cnt.size(); i < N; ++i) {
                    if(cnt[i] == 0) {
                        std::string FullFileName = NFilesIter::FileName(i, filePrefixPrePost[lev]);
                        amrex::UnlinkFile(FullFileName.c_str());
                    }
                }
            }
        }
    }
    
    if(ParallelDescriptor::IOProcessor()) {
        HdrFile.flush();
        HdrFile.close();
        if( ! HdrFile.good()) {
            amrex::Abort("ParticleContainer::CheckpointPost(): problem writing HdrFile");
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFilePre ()
{
    CheckpointPre();
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFilePost ()
{
    CheckpointPost();
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WriteParticles (int lev, std::ofstream& ofs, int fnum,
                  Vector<int>& which, Vector<int>& count, Vector<Long>& where,
                  const Vector<int>& write_real_comp,
                  const Vector<int>& write_int_comp,
                  const Vector<std::map<std::pair<int, int>, Gpu::DeviceVector<int>>>& particle_io_flags) const
{
    BL_PROFILE("ParticleContainer::WriteParticles()");

    // For a each grid, the tiles it contains
    std::map<int, Vector<int> > tile_map;

    for (const auto& kv : m_particles[lev])
    {
        const int grid = kv.first.first;
        const int tile = kv.first.second;
        tile_map[grid].push_back(tile);
        const auto& pflags = particle_io_flags[lev].at(kv.first);
		
        // Only write out valid particles.
        int cnt = 0;	
        for (int k = 0; k < kv.second.GetArrayOfStructs().size(); ++k)
        {
            if (pflags[k]) cnt++;
        }
		
        count[grid] += cnt;
    }
    
    MFInfo info;
    info.SetAlloc(false);
    MultiFab state(ParticleBoxArray(lev), ParticleDistributionMap(lev), 1,0,info);
    
    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const int grid = mfi.index();
        
        which[grid] = fnum;
        where[grid] = VisMF::FileOffset(ofs);
        
        if (count[grid] == 0) continue;
      
        // First write out the integer data in binary.
        int num_output_int = 0;
        for (int i = 0; i < NumIntComps() + NStructInt; ++i)
            if (write_int_comp[i]) ++num_output_int;
        
        const int iChunkSize = 2 + num_output_int;
        Vector<int> istuff(count[grid]*iChunkSize);
        int* iptr = istuff.dataPtr();
        
        for (unsigned i = 0; i < tile_map[grid].size(); i++) {
            auto ptile_index = std::make_pair(grid, tile_map[grid][i]);
            const auto& pbox = m_particles[lev].at(ptile_index);
            const auto& pflags = particle_io_flags[lev].at(ptile_index);
            for (int pindex = 0; pindex < pbox.GetArrayOfStructs().size(); ++pindex) {
                const auto& aos = pbox.GetArrayOfStructs();
                const auto& p = aos[pindex];
                if (pflags[pindex])
                {
                    // always write these
                    for (int j = 0; j < 2; j++) iptr[j] = p.m_idata.arr[j];
                    iptr += 2;
                    
                    // optionally write these
                    for (int j = 0; j < NStructInt; j++)
                    {
                        if (write_int_comp[j])
                        {
                            *iptr = p.m_idata.arr[2+j];
                            ++iptr;
                        }
                    }
                    
                    const auto& soa  = pbox.GetStructOfArrays();
                    for (int j = 0; j < NumIntComps(); j++)
                    {
                        if (write_int_comp[NStructInt+j])
                        {
                            *iptr = soa.GetIntData(j)[pindex];
                            ++iptr;
                        }
                    }
                }
            }
        }
                
        writeIntData(istuff.dataPtr(), istuff.size(), ofs);
        ofs.flush();  // Some systems require this flush() (probably due to a bug)
        
        // Write the Real data in binary.
        int num_output_real = 0;
        for (int i = 0; i < NumRealComps() + NStructReal; ++i)
            if (write_real_comp[i]) ++num_output_real;
        
        const int rChunkSize = AMREX_SPACEDIM + num_output_real;
        Vector<typename ParticleType::RealType> rstuff(count[grid]*rChunkSize);
        typename ParticleType::RealType* rptr = rstuff.dataPtr();
        
        for (unsigned i = 0; i < tile_map[grid].size(); i++) {
			auto ptile_index = std::make_pair(grid, tile_map[grid][i]);
            const auto& pbox = m_particles[lev].at(ptile_index);
			const auto& pflags = particle_io_flags[lev].at(ptile_index);
            for (int pindex = 0; pindex < pbox.GetArrayOfStructs().size(); ++pindex) {
				const auto& aos = pbox.GetArrayOfStructs();
				const auto& p = aos[pindex];				
                if (pflags[pindex])
                {
                    // always write these
                    for (int j = 0; j < AMREX_SPACEDIM; j++) rptr[j] = p.m_rdata.arr[j];
                    rptr += AMREX_SPACEDIM;
                    
                    // optionally write these
                    for (int j = 0; j < NStructReal; j++)
                    {
                        if (write_real_comp[j])
                        {
                            *rptr = p.m_rdata.arr[AMREX_SPACEDIM+j];
                            ++rptr;
                        }
                    }
                    
                    const auto& soa  = pbox.GetStructOfArrays();
                    for (int j = 0; j < NumRealComps(); j++)
                    {
                        if (write_real_comp[NStructReal+j])
                        {
                            *rptr = (typename ParticleType::RealType) soa.GetRealData(j)[pindex];
                            ++rptr;
                        }
                    }
                }
            }
        }
        
        WriteParticleRealData(rstuff.dataPtr(), rstuff.size(), ofs, ParticleRealDescriptor);
        ofs.flush();  // Some systems require this flush() (probably due to a bug)
    }
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Restart (const std::string& dir, const std::string& file, bool is_checkpoint)
{
    Restart(dir, file);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Restart (const std::string& dir, const std::string& file)
{
    BL_PROFILE("ParticleContainer::Restart()");
    AMREX_ASSERT(!dir.empty());
    AMREX_ASSERT(!file.empty());
    
    const Real strttime = amrex::second();
    
    int DATA_Digits_Read(5);
    ParmParse pp("particles");
    pp.query("datadigits_read",DATA_Digits_Read);

    std::string fullname = dir;
    if (!fullname.empty() && fullname[fullname.size()-1] != '/')
        fullname += '/';
    fullname += file;
    std::string HdrFileName = fullname;
    if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
        HdrFileName += '/';
    HdrFileName += "Header";
    
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(HdrFileName, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream HdrFile(fileCharPtrString, std::istringstream::in);
  
    std::string version;
    HdrFile >> version;
    AMREX_ASSERT(!version.empty());
    
    // What do our version strings mean?
    // "Version_One_Dot_Zero" -- hard-wired to write out in double precision.
    // "Version_One_Dot_One" -- can write out either as either single or double precision.
    // Appended to the latter version string are either "_single" or "_double" to
    // indicate how the particles were written.
    // "Version_Two_Dot_Zero" -- this is the AMReX particle file format
    std::string how;
    if (version.find("Version_One_Dot_Zero") != std::string::npos) {
        how = "double";
    }
    else if (version.find("Version_One_Dot_One")  != std::string::npos or
             version.find("Version_Two_Dot_Zero") != std::string::npos) {
        if (version.find("_single") != std::string::npos) {
            how = "single";
        }
        else if (version.find("_double") != std::string::npos) {
            how = "double";
        }
        else {
            std::string msg("ParticleContainer::Restart(): bad version string: ");
            msg += version;
            amrex::Error(version.c_str());
        }
    }
    else {
        std::string msg("ParticleContainer::Restart(): unknown version string: ");
        msg += version;
        amrex::Abort(msg.c_str());
    }
    
    int dm;
    HdrFile >> dm;
    if (dm != AMREX_SPACEDIM)
        amrex::Abort("ParticleContainer::Restart(): dm != AMREX_SPACEDIM");
    
    int nr;
    HdrFile >> nr;
    if (nr != NStructReal + NumRealComps())
        amrex::Abort("ParticleContainer::Restart(): nr != NStructReal + NumRealComps()");
    
    std::string comp_name;
    for (int i = 0; i < nr; ++i)
        HdrFile >> comp_name;
    
    int ni;
    HdrFile >> ni;
    if (ni != NStructInt + NumIntComps())
        amrex::Abort("ParticleContainer::Restart(): ni != NStructInt");
    
    for (int i = 0; i < ni; ++i)
        HdrFile >> comp_name;
    
    bool checkpoint;
    HdrFile >> checkpoint;
    
    Long nparticles;
    HdrFile >> nparticles;
    AMREX_ASSERT(nparticles >= 0);
    
    int maxnextid;
    HdrFile >> maxnextid;
    AMREX_ASSERT(maxnextid > 0);
    ParticleType::NextID(maxnextid);
    
    int finest_level_in_file;
    HdrFile >> finest_level_in_file;
    AMREX_ASSERT(finest_level_in_file >= 0);
    
    // Determine whether this is a dual-grid restart or not.
    Vector<BoxArray> particle_box_arrays(finest_level_in_file + 1);
    bool dual_grid = false;

    bool have_pheaders = false;
    for (int lev = 0; lev <= finest_level_in_file; lev++)
    {
        std::string phdr_name = fullname;
        phdr_name = amrex::Concatenate(phdr_name + "/Level_", lev, 1);
        phdr_name += "/Particle_H";

        if (amrex::FileExists(phdr_name)) {
            have_pheaders = true;
            break;
        }        
    }

    if (have_pheaders)
    {
        for (int lev = 0; lev <= finest_level_in_file; lev++)
        {
            std::string phdr_name = fullname;
            phdr_name = amrex::Concatenate(phdr_name + "/Level_", lev, 1);
            phdr_name += "/Particle_H";

            if (not amrex::FileExists(phdr_name)) continue;
            
            Vector<char> phdr_chars;
            ParallelDescriptor::ReadAndBcastFile(phdr_name, phdr_chars);
            std::string phdr_string(phdr_chars.dataPtr());
            std::istringstream phdr_file(phdr_string, std::istringstream::in);
            
            if (lev > finestLevel())
            {
                dual_grid = true;
                break;
            }
        
            particle_box_arrays[lev].readFrom(phdr_file);
            if (not particle_box_arrays[lev].CellEqual(ParticleBoxArray(lev))) dual_grid = true;
        }
    } else // if no particle box array information exists in the file, we assume a single grid restart
    {
        dual_grid = false;
    }

    if (dual_grid) {
        for (int lev = 0; lev <= finestLevel(); lev++) {
            SetParticleBoxArray(lev, particle_box_arrays[lev]);
            DistributionMapping pdm(particle_box_arrays[lev]);
            SetParticleDistributionMap(lev, pdm);
        }
    }
    
    Vector<int> ngrids(finest_level_in_file+1);
    for (int lev = 0; lev <= finest_level_in_file; lev++) {
        HdrFile >> ngrids[lev];
        AMREX_ASSERT(ngrids[lev] > 0);
        if (lev <= finestLevel()) {
            AMREX_ASSERT(ngrids[lev] == int(ParticleBoxArray(lev).size()));
        }
    }
    
    resizeData();
    
    if (finest_level_in_file > finestLevel()) {
        m_particles.resize(finest_level_in_file+1);
    }
    
    for (int lev = 0; lev <= finest_level_in_file; lev++) {
        Vector<int>  which(ngrids[lev]);
        Vector<int>  count(ngrids[lev]);
        Vector<Long> where(ngrids[lev]);
        for (int i = 0; i < ngrids[lev]; i++) {
            HdrFile >> which[i] >> count[i] >> where[i];
        }
        
        Vector<int> grids_to_read;
        if (lev <= finestLevel()) {
            for (MFIter mfi(*m_dummy_mf[lev]); mfi.isValid(); ++mfi) {
                grids_to_read.push_back(mfi.index());
            }
        } else {
            
            // we lost a level on restart. we still need to read in particles
            // on finer levels, and put them in the right place via Redistribute()
            
            const int rank = ParallelDescriptor::MyProc();
            const int NReaders = ParticleType::MaxReaders();
            if (rank >= NReaders) return;
            
            const int Navg = ngrids[lev] / NReaders;
            const int Nleft = ngrids[lev] - Navg * NReaders;
            
            int lo, hi;
            if (rank < Nleft) {
                lo = rank*(Navg + 1);
                hi = lo + Navg + 1;
            } 
            else {
                lo = rank * Navg + Nleft;
                hi = lo + Navg;
            }
            
            for (int i = lo; i < hi; ++i) {
                grids_to_read.push_back(i);
            }
        }
        
        for(int igrid = 0; igrid < static_cast<int>(grids_to_read.size()); ++igrid) {
            const int grid = grids_to_read[igrid];
            
            if (count[grid] <= 0) continue;
            
            // The file names in the header file are relative.
            std::string name = fullname;
            
            if (!name.empty() && name[name.size()-1] != '/')
                name += '/';
            
            name += "Level_";
            name += amrex::Concatenate("", lev, 1);
            name += '/';
            name += ParticleType::DataPrefix();
            name += amrex::Concatenate("", which[grid], DATA_Digits_Read);
            
            std::ifstream ParticleFile;
            
            ParticleFile.open(name.c_str(), std::ios::in | std::ios::binary);
            
            if (!ParticleFile.good())
                amrex::FileOpenFailed(name);
            
            ParticleFile.seekg(where[grid], std::ios::beg);
            
            if (how == "single") {
                ReadParticles<float>(count[grid], grid, lev, ParticleFile, finest_level_in_file);
            }
            else if (how == "double") {
                ReadParticles<double>(count[grid], grid, lev, ParticleFile, finest_level_in_file);
            }
            else {
                std::string msg("ParticleContainer::Restart(): bad parameter: ");
                msg += how;
                amrex::Error(msg.c_str());
            }
            
            ParticleFile.close();
            
            if (!ParticleFile.good())
                amrex::Abort("ParticleContainer::Restart(): problem reading particles");
        }
    }
    
    Redistribute();
    
    AMREX_ASSERT(OK());
    
    if (m_verbose > 1) {
        Real stoptime = amrex::second() - strttime;	
        ParallelDescriptor::ReduceRealMax(stoptime, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "ParticleContainer::Restart() time: " << stoptime << '\n';
    }
}

// Read a batch of particles from the checkpoint file
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class RTYPE>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::ReadParticles (int cnt, int grd, int lev, std::ifstream& ifs, int finest_level_in_file)
{
    BL_PROFILE("ParticleContainer::ReadParticles()");
    AMREX_ASSERT(cnt > 0);
    AMREX_ASSERT(lev < int(m_particles.size()));

    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    const int iChunkSize = 2 + NStructInt + NumIntComps();
    Vector<int> istuff(cnt*iChunkSize);
    readIntData(istuff.dataPtr(), istuff.size(), ifs, FPC::NativeIntDescriptor());
    
    // Then the real data in binary.
    const int rChunkSize = AMREX_SPACEDIM + NStructReal + NumRealComps();
    Vector<RTYPE> rstuff(cnt*rChunkSize);
    ReadParticleRealData(rstuff.dataPtr(), rstuff.size(), ifs, ParticleRealDescriptor);
    
    // Now reassemble the particles.
    int*   iptr = istuff.dataPtr();
    RTYPE* rptr = rstuff.dataPtr();
    
    ParticleType p;
    ParticleLocData pld;

    Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
    host_particles.reserve(15);
    host_particles.resize(finest_level_in_file+1);

    Vector<std::map<std::pair<int, int>,
                    std::vector<Gpu::HostVector<Real> > > > host_real_attribs;
    host_real_attribs.reserve(15);
    host_real_attribs.resize(finest_level_in_file+1);

    Vector<std::map<std::pair<int, int>,
                    std::vector<Gpu::HostVector<int> > > > host_int_attribs;
    host_int_attribs.reserve(15);
    host_int_attribs.resize(finest_level_in_file+1);

    for (int i = 0; i < cnt; i++) {
        p.m_idata.id   = iptr[0];
        p.m_idata.cpu  = iptr[1];
        
        iptr += 2;
            
        for (int j = 0; j < NStructInt; j++)
        {
            p.m_idata.arr[2+j] = *iptr;
            ++iptr;
        }

        AMREX_ASSERT(p.m_idata.id > 0);
        
        AMREX_D_TERM(p.m_rdata.pos[0] = rptr[0];,
                     p.m_rdata.pos[1] = rptr[1];,
                     p.m_rdata.pos[2] = rptr[2];);
        
        rptr += AMREX_SPACEDIM;
        
        for (int j = 0; j < NStructReal; j++)
        {
            p.m_rdata.arr[AMREX_SPACEDIM+j] = *rptr;
            ++rptr;
        }

        locateParticle(p, pld, 0, finestLevel(), 0);
        
	std::pair<int, int> ind(grd, pld.m_tile);

        host_real_attribs[lev][ind].resize(NumRealComps());
        host_int_attribs[lev][ind].resize(NumIntComps());
        
	// add the struct
	host_particles[lev][ind].push_back(p);

	// add the real...
	for (int icomp = 0; icomp < NumRealComps(); icomp++) {
            host_real_attribs[lev][ind][icomp].push_back(*rptr);
            ++rptr;
	}
        
	// ... and int array data
	for (int icomp = 0; icomp < NumIntComps(); icomp++) {
            host_int_attribs[lev][ind][icomp].push_back(*iptr);
            ++iptr;
	}        
    }

    for (int host_lev = 0; host_lev < static_cast<int>(host_particles.size()); ++host_lev)
      {
	for (auto& kv : host_particles[host_lev]) {
	  auto grid = kv.first.first;
	  auto tile = kv.first.second;
	  const auto& src_tile = kv.second;
          
	  auto& dst_tile = DefineAndReturnParticleTile(host_lev, grid, tile);
	  auto old_size = dst_tile.GetArrayOfStructs().size();
	  auto new_size = old_size + src_tile.size();
	  dst_tile.resize(new_size);
                
	  Gpu::copy(Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
                    dst_tile.GetArrayOfStructs().begin() + old_size);
	  
	  for (int i = 0; i < NumRealComps(); ++i) {
              Gpu::copy(Gpu::hostToDevice,
                        host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                        host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                        dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
	  }
	  
	  for (int i = 0; i < NumIntComps(); ++i) {
              Gpu::copy(Gpu::hostToDevice,
                        host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                        host_int_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                        dst_tile.GetStructOfArrays().GetIntData(i).begin() + old_size);
	  }
	}
      }
    
    Gpu::streamSynchronize();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::WriteAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer::WriteAsciiFile()");
    AMREX_ASSERT(!filename.empty());

    const Real strttime = amrex::second();
    //
    // Count # of valid particles.
    //
    Long nparticles = 0;

    for (int lev = 0; lev < m_particles.size();  lev++) {
        auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
	    auto np = aos.numParticles();
	    Gpu::HostVector<ParticleType> host_aos(np);
	    Gpu::copy(Gpu::deviceToHost, aos.begin(), aos.end(), host_aos.begin());
	    for (int k = 0; k < np; ++k) {
	        const ParticleType& p = host_aos[k];
                if (p.m_idata.id > 0)
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
            }
        }
    }
    
    //
    // And send count to I/O processor.
    //
    ParallelDescriptor::ReduceLongSum(nparticles,ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Have I/O processor open file and write out particle metadata.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out|std::ios::trunc);

        if (!File.good())
            amrex::FileOpenFailed(filename);

        File << nparticles  << '\n';
        File << NStructReal << '\n';
        File << NStructInt  << '\n';
        File << NumRealComps()  << '\n';
        File << NumIntComps()   << '\n';
            
        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("ParticleContainer::WriteAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
    {
        if (MyProc == proc)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            std::ofstream File;

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                amrex::FileOpenFailed(filename);

	    for (int lev = 0; lev < m_particles.size();  lev++) {
	      auto& pmap = m_particles[lev];
	      for (const auto& kv : pmap) {
                const auto& aos = kv.second.GetArrayOfStructs();
                const auto& soa = kv.second.GetStructOfArrays();

		auto np = aos.numParticles();
		Gpu::HostVector<ParticleType> host_aos(np);
		Gpu::copy(Gpu::deviceToHost, aos.begin(), aos.end(), host_aos.begin());

		for (int index = 0; index < np; ++index) {
		    const ParticleType* it = &host_aos[index];
		    if (it->m_idata.id > 0) {

                        // write out the particle struct first... 
                        AMREX_D_TERM(File << it->m_rdata.pos[0] << ' ',
                               << it->m_rdata.pos[1] << ' ',
                               << it->m_rdata.pos[2] << ' ');

                        for (int i = AMREX_SPACEDIM; i < AMREX_SPACEDIM + NStructReal; i++)
                            File << it->m_rdata.arr[i] << ' ';

                        File << it->m_idata.id  << ' ';
                        File << it->m_idata.cpu << ' ';
                        
                        for (int i = 2; i < 2 + NStructInt; i++)
                            File << it->m_idata.arr[i] << ' ';
		      
                        // then the particle attributes.
                        for (int i = 0; i < NumRealComps(); i++)
                            File << soa.GetRealData(i)[index] << ' ';
                        
                        for (int i = 0; i < NumIntComps(); i++)
                            File << soa.GetIntData(i)[index] << ' ';
                        
                        File << '\n';                                                    
                    }
                }
              }
            }
	    
            File.flush();
	    
            File.close();
            
            if (!File.good())
                amrex::Abort("ParticleContainer::WriteAsciiFile(): problem writing file");
	    
        }
	
        ParallelDescriptor::Barrier();
    }
    
    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;
        
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
        
        amrex::Print() << "ParticleContainer::WriteAsciiFile() time: " << stoptime << '\n';
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal,NStructInt,NArrayReal, NArrayInt>::WriteCoarsenedAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer::WriteCoarsenedAsciiFile()");
    AMREX_ASSERT(!filename.empty());

    const Real strttime = amrex::second();
 
    //
    // Count # of valid particles.
    //
    Long nparticles = 0;

    for (int lev = 0; lev < m_particles.size();  lev++) {
        auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
            for (int k = 0; k < aos.size(); ++k) {
	        const ParticleType& p = aos[k];
                if (p.m_idata.id > 0)
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
            }
        }
    }
 
    //
    // And send count to I/O processor.
    //
    ParallelDescriptor::ReduceLongSum(nparticles,ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Have I/O processor open file and write out particle count.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out|std::ios::trunc);

        if (!File.good())
            amrex::FileOpenFailed(filename);

        File << nparticles << '\n';
            
        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("ParticleContainer::WriteCoarsenedAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
    {
        if (MyProc == proc)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            std::ofstream File;

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                amrex::FileOpenFailed(filename);

	    for (int lev = 0; lev < m_particles.size();  lev++) {
	      auto& pmap = m_particles[lev];
	      for (auto& kv : pmap) {
                  auto& aos = kv.second.GetArrayOfStructs();
                  auto& soa = kv.second.GetStructOfArrays();
                  
                  int index = 0;
                  ParticleLocData pld;
                  for (int k = 0; k < aos.size(); ++k) {
					  ParticleType* it = &aos[k];
                      locateParticle(*it, pld, 0, finestLevel(), 0);
                      // Only keep particles in even cells
                      if (it->id() > 0 &&
                          (pld.m_cell[0])%2 == 0 && (pld.m_cell[1])%2 == 0 && (pld.m_cell[2])%2 == 0)
                      {
                          
                          // Only keep particles in even cells               
                          if (it->m_idata.id > 0) {
                              
                              File << it->m_idata.id  << ' ';
                              File << it->m_idata.cpu << ' ';
                              
                              AMREX_D_TERM(File << it->m_rdata.pos[0] << ' ',
                                     << it->m_rdata.pos[1] << ' ',
                                     << it->m_rdata.pos[2] << ' ');
                              
                              for (int i = 0; i < NumRealComps(); i++) {
                                  File << soa.GetRealData(i)[index] << ' ';
                              }
                              index++;
                              
                              for (int i = AMREX_SPACEDIM; i < AMREX_SPACEDIM + NStructReal; i++) {
                                  char ws = (i == AMREX_SPACEDIM + NStructReal - 1) ? '\n' : ' ';
                                  if (i == AMREX_SPACEDIM) {
                                      // Multiply mass by 8 since we are only taking 1/8 of the 
                                      // total particles and want to keep the mass in the domain the same.
                                      File << 8.0* it->m_rdata.arr[i] << ws;
                                  }
                                  else {
                                      File << it->m_rdata.arr[i] << ws;
                                  }
                              }
                          }
                      }
                  }
              }
            }

            File.flush();

            File.close();

            if (!File.good())
                amrex::Abort("ParticleContainer::WriteCoarsenedAsciiFile(): problem writing file");

        }

        ParallelDescriptor::Barrier();
    }

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "ParticleContainer::WriteCoarsenedAsciiFile() time: " << stoptime << '\n';
    }
}
#endif /*AMREX_PARTICLEIO_H*/
