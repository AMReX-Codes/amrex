#ifndef AMREX_WRITE_BINARY_PARTICLE_DATA_H
#define AMREX_WRITE_BINARY_PARTICLE_DATA_H
#include <AMReX_Config.H>

#include <AMReX_TypeTraits.H>
#include <AMReX_Particles.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_GpuDevice.H>

struct KeepValidFilter
{
    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int i) const noexcept
    {
        return (src.m_aos[i].id() > 0);
    }
};

namespace particle_detail {

template <typename ParticleReal>
std::size_t PSizeInFile (const Vector<int>& wrc, const Vector<int>& wic)
{
    std::size_t rsize = sizeof(ParticleReal)*std::accumulate(wrc.begin(), wrc.end(), 0);
    std::size_t isize = sizeof(int)*std::accumulate(wic.begin(), wic.end(), 0);
    return rsize + isize + AMREX_SPACEDIM*sizeof(ParticleReal) + 2*sizeof(int);
}

template <template <class, class> class Container,
          class Allocator,
          class PTile,
          class F>
typename std::enable_if<RunOnGpu<typename Container<int, Allocator>::allocator_type>::value>::type
fillFlags (Container<int, Allocator>& pflags, const PTile& ptile, F&& f)
{
    const auto ptd = ptile.getConstParticleTileData();
    const auto np = ptile.numParticles();
    pflags.resize(np, 0);
    auto flag_ptr = pflags.data();
    amrex::ParallelForRNG(np,
        [=] AMREX_GPU_DEVICE (int k, amrex::RandomEngine const& engine) noexcept
        {
            const auto p = ptd.getSuperParticle(k);
            flag_ptr[k] = particle_detail::call_f(f,p,engine);
        });
}

template <template <class, class> class Container,
          class Allocator,
          class PTile,
          class F>
typename std::enable_if<!RunOnGpu<typename Container<int, Allocator>::allocator_type>::value>::type
fillFlags (Container<int, Allocator>& pflags, const PTile& ptile, F&& f)
{
    const auto ptd = ptile.getConstParticleTileData();
    const auto np = ptile.numParticles();
    pflags.resize(np, 0);
    auto flag_ptr = pflags.data();
    for (int k = 0; k < np; ++k) {
        const auto p = ptd.getSuperParticle(k);
        flag_ptr[k] = particle_detail::call_f(f,p,RandomEngine{});
    }
}

template <template <class, class> class Container, class Allocator, class PC>
typename std::enable_if<RunOnGpu<typename Container<int, Allocator>::allocator_type>::value, amrex::Long>::type
countFlags (const Vector<std::map<std::pair<int,int>,Container<int,Allocator>>>& particle_io_flags, const PC& pc)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (int lev = 0; lev < pc.GetParticles().size();  lev++)
    {
        const auto& pmap = pc.GetParticles(lev);
        for (const auto& kv : pmap)
        {
            const auto& pflags = particle_io_flags[lev].at(kv.first);
            const auto flag_ptr = pflags.data();
            reduce_op.eval(pflags.size(), reduce_data,
                [=] AMREX_GPU_DEVICE (const int i) -> ReduceTuple
                {
                    return flag_ptr[i] ? 1 : 0;
                });
        }
    }
    ReduceTuple hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv);
}

template <template <class, class> class Container, class Allocator>
typename std::enable_if<RunOnGpu<typename Container<int, Allocator>::allocator_type>::value, int>::type
countFlags (const Container<int,Allocator>& pflags)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    const auto flag_ptr = pflags.data();
    reduce_op.eval(pflags.size(), reduce_data,
        [=] AMREX_GPU_DEVICE (const int i) -> ReduceTuple
        {
            return flag_ptr[i] ? 1 : 0;
        });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv);
}

template <template <class, class> class Container, class Allocator, class PC>
typename std::enable_if<!RunOnGpu<typename Container<int, Allocator>::allocator_type>::value, amrex::Long>::type
countFlags (const Vector<std::map<std::pair<int,int>,Container<int,Allocator>>>& particle_io_flags, const PC& pc)
{
    amrex::Long nparticles = 0;
    for (int lev = 0; lev < pc.GetParticles().size();  lev++)
    {
        const auto& pmap = pc.GetParticles(lev);
        for (const auto& kv : pmap)
        {
            const auto& pflags = particle_io_flags[lev].at(kv.first);
            for (int k = 0; k < kv.second.numParticles(); ++k)
            {
                if (pflags[k]) nparticles++;
            }
        }
    }
    return nparticles;
}

template <template <class, class> class Container, class Allocator>
typename std::enable_if<!RunOnGpu<typename Container<int, Allocator>::allocator_type>::value, int>::type
countFlags (const Container<int,Allocator>& pflags)
{
    int nparticles = 0;
    for (std::size_t k = 0; k < pflags.size(); ++k)
    {
        if (pflags[k]) nparticles++;
    }
    return nparticles;
}

template <class PC>
typename std::enable_if<RunOnGpu<typename PC::template AllocatorType<int>>::value>::type
packIOData (Vector<int>& idata, Vector<ParticleReal>& rdata, const PC& pc, int lev, int grid,
            const Vector<int>& write_real_comp, const Vector<int>& write_int_comp,
            const Vector<std::map<std::pair<int, int>, typename PC::IntVector>>& particle_io_flags,
            const Vector<int>& tiles, int np)
{
    int num_output_int = 0;
    for (int i = 0; i < pc.NumIntComps() + PC::NStructInt; ++i)
        if (write_int_comp[i]) ++num_output_int;

    const int iChunkSize = 2 + num_output_int;
    idata.resize(np*iChunkSize);

    int num_output_real = 0;
    for (int i = 0; i < pc.NumRealComps() + PC::NStructReal; ++i)
        if (write_real_comp[i]) ++num_output_real;

    const int rChunkSize = AMREX_SPACEDIM + num_output_real;
    rdata.resize(np*rChunkSize);

    typename PC::IntVector  idata_d(idata.size());
    typename PC::RealVector rdata_d(rdata.size());

    typename PC::IntVector write_int_comp_d(write_int_comp.size());
    typename PC::IntVector write_real_comp_d(write_real_comp.size());
    Gpu::copy(Gpu::hostToDevice, write_int_comp.begin(), write_int_comp.end(),
              write_int_comp_d.begin());
    Gpu::copy(Gpu::hostToDevice, write_real_comp.begin(), write_real_comp.end(),
              write_real_comp_d.begin());
    Gpu::Device::synchronize();

    const auto write_int_comp_d_ptr = write_int_comp_d.data();
    const auto write_real_comp_d_ptr = write_real_comp_d.data();

    std::size_t poffset = 0;
    for (int i = 0; i < tiles.size(); i++) {
        const auto& ptile = pc.ParticlesAt(lev, grid, tiles[i]);
        const auto& pflags = particle_io_flags[lev].at(std::make_pair(grid, tiles[i]));
        int np_tile = ptile.GetArrayOfStructs().numParticles();
        typename PC::IntVector offsets(np_tile);
        int num_copies = Scan::ExclusiveSum(np_tile, pflags.begin(), offsets.begin(), Scan::retSum);

        const auto flag_ptr = pflags.data();

        auto idata_d_ptr = idata_d.data();
        auto rdata_d_ptr = rdata_d.data();

        const auto ptd = ptile.getConstParticleTileData();
        amrex::ParallelFor(num_copies,
        [=] AMREX_GPU_DEVICE (int pindex) noexcept
        {
            // might be worth using shared memory here
            const auto p = ptd.getSuperParticle(pindex);

            if (flag_ptr[pindex]) {
                std::size_t iout_index = (pindex+poffset)*iChunkSize;
                idata_d_ptr[iout_index] = p.id();
                idata_d_ptr[iout_index+1] = p.cpu();
                iout_index += 2;

                std::size_t rout_index = (pindex+poffset)*rChunkSize;
                for (int j = 0; j < AMREX_SPACEDIM; j++) {
                  rdata_d_ptr[rout_index] = p.pos(j);
                  rout_index++;
                }

                for (int j = 0; j < PC::SuperParticleType::NInt; j++) {
                    if (write_int_comp_d_ptr[j]) {
                        idata_d_ptr[iout_index] = p.idata(j);
                        iout_index++;
                    }
                }

                for (int j = 0; j < PC::SuperParticleType::NReal; j++) {
                    if (write_real_comp_d_ptr[j]) {
                        rdata_d_ptr[rout_index] = p.rdata(j);
                        rout_index++;
                    }
                }
            }
        });

        poffset += num_copies;
    }

    Gpu::copy(Gpu::deviceToHost, idata_d.begin(), idata_d.end(), idata.begin());
    Gpu::copy(Gpu::deviceToHost, rdata_d.begin(), rdata_d.end(), rdata.begin());
    Gpu::Device::synchronize();
}

template <class PC>
typename std::enable_if<!RunOnGpu<typename PC::template AllocatorType<int>>::value>::type
packIOData (Vector<int>& idata, Vector<ParticleReal>& rdata, const PC& pc, int lev, int grid,
            const Vector<int>& write_real_comp, const Vector<int>& write_int_comp,
            const Vector<std::map<std::pair<int, int>, typename PC::IntVector>>& particle_io_flags,
            const Vector<int>& tiles, int np)
{
    int num_output_int = 0;
    for (int i = 0; i < pc.NumIntComps() + PC::NStructInt; ++i)
        if (write_int_comp[i]) ++num_output_int;

    const int iChunkSize = 2 + num_output_int;
    idata.resize(np*iChunkSize);

    int num_output_real = 0;
    for (int i = 0; i < pc.NumRealComps() + PC::NStructReal; ++i)
        if (write_real_comp[i]) ++num_output_real;

    const int rChunkSize = AMREX_SPACEDIM + num_output_real;
    rdata.resize(np*rChunkSize);

    int* iptr = idata.dataPtr();
    ParticleReal* rptr = rdata.dataPtr();
    for (unsigned i = 0; i < tiles.size(); i++) {
        const auto& ptile = pc.ParticlesAt(lev, grid, tiles[i]);
        const auto& pflags = particle_io_flags[lev].at(std::make_pair(grid, tiles[i]));
        for (int pindex = 0; pindex < ptile.GetArrayOfStructs().numParticles(); ++pindex) {
            const auto& aos = ptile.GetArrayOfStructs();
            const auto& p = aos[pindex];
            if (pflags[pindex]) {
                *iptr = p.id(); ++iptr;
                *iptr = p.cpu(); ++iptr;
                for (int j = 0; j < AMREX_SPACEDIM; j++) rptr[j] = p.pos(j);
                rptr += AMREX_SPACEDIM;

                for (int j = 0; j < PC::NStructInt; j++) {
                    if (write_int_comp[j]) {
                        *iptr = p.idata(j);
                        ++iptr;
                    }
                }

                const auto& soa  = ptile.GetStructOfArrays();
                for (int j = 0; j < pc.NumIntComps(); j++) {
                    if (write_int_comp[PC::NStructInt+j]) {
                        *iptr = soa.GetIntData(j)[pindex];
                        ++iptr;
                    }
                }

                for (int j = 0; j < PC::NStructReal; j++) {
                    if (write_real_comp[j]) {
                        *rptr = p.rdata(j);
                        ++rptr;
                    }
                }

                for (int j = 0; j < pc.NumRealComps(); j++) {
                    if (write_real_comp[PC::NStructReal+j]) {
                        *rptr = (ParticleReal) soa.GetRealData(j)[pindex];
                        ++rptr;
                    }
                }
            }
        }
    }
}
}

template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
void WriteBinaryParticleDataSync (PC const& pc,
                                  const std::string& dir, const std::string& name,
                                  const Vector<int>& write_real_comp,
                                  const Vector<int>& write_int_comp,
                                  const Vector<std::string>& real_comp_names,
                                  const Vector<std::string>& int_comp_names,
                                  F&& f)
{
    BL_PROFILE("WriteBinaryParticleData()");
    AMREX_ASSERT(pc.OK());

    AMREX_ASSERT(sizeof(typename PC::ParticleType::RealType) == 4 ||
                 sizeof(typename PC::ParticleType::RealType) == 8);

    constexpr int NStructReal = PC::NStructReal;
    constexpr int NStructInt  = PC::NStructInt;

    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();

    AMREX_ALWAYS_ASSERT(real_comp_names.size() == pc.NumRealComps() + NStructReal);
    AMREX_ALWAYS_ASSERT( int_comp_names.size() == pc.NumIntComps() + NStructInt);

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

    std::ofstream HdrFile;

    Long nparticles = 0;
    Long maxnextid;

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

    Gpu::Device::synchronize();

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

    if (ParallelDescriptor::IOProcessor())
    {
        std::string HdrFileName = pdir;

        if ( ! HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
            HdrFileName += '/';

        HdrFileName += "Header";
        pc.HdrFileNamePrePost = HdrFileName;

        HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);

        if ( ! HdrFile.good()) amrex::FileOpenFailed(HdrFileName);

        //
        // First thing written is our Checkpoint/Restart version string.
        // We append "_single" or "_double" to the version string indicating
        // whether we're using "float" or "double" floating point data in the
        // particles so that we can Restart from the checkpoint files.
        //
        if (sizeof(typename PC::ParticleType::RealType) == 4)
        {
            HdrFile << PC::Version() << "_single" << '\n';
        }
        else
        {
            HdrFile << PC::Version() << "_double" << '\n';
        }

        int num_output_real = 0;
        for (int i = 0; i < pc.NumRealComps() + NStructReal; ++i)
            if (write_real_comp[i]) ++num_output_real;

        int num_output_int = 0;
        for (int i = 0; i < pc.NumIntComps() + NStructInt; ++i)
            if (write_int_comp[i]) ++num_output_int;

        // AMREX_SPACEDIM and N for sanity checking.
        HdrFile << AMREX_SPACEDIM << '\n';

        // The number of extra real parameters
        HdrFile << num_output_real << '\n';

        // Real component names
        for (int i = 0; i < NStructReal + pc.NumRealComps(); ++i )
            if (write_real_comp[i]) HdrFile << real_comp_names[i] << '\n';

        // The number of extra int parameters
        HdrFile << num_output_int << '\n';

        // int component names
        for (int i = 0; i < NStructInt + pc.NumIntComps(); ++i )
            if (write_int_comp[i]) HdrFile << int_comp_names[i] << '\n';

        bool is_checkpoint = true; // legacy
        HdrFile << is_checkpoint << '\n';

        // The total number of particles.
        HdrFile << nparticles << '\n';

        // The value of nextid that we need to restore on restart.
        HdrFile << maxnextid << '\n';

        // Then the finest level of the AMR hierarchy.
        HdrFile << pc.finestLevel() << '\n';

        // Then the number of grids at each level.
        for (int lev = 0; lev <= pc.finestLevel(); lev++)
            HdrFile << pc.ParticleBoxArray(lev).size() << '\n';
    }

    // We want to write the data out in parallel.
    // We'll allow up to nOutFiles active writers at a time.
    int nOutFiles(256);

    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) nOutFiles = NProcs;
    nOutFiles = std::max(1, std::min(nOutFiles,NProcs));
    pc.nOutFilesPrePost = nOutFiles;

    for (int lev = 0; lev <= pc.finestLevel(); lev++)
    {
        bool gotsome;
        if(pc.usePrePost)
        {
            gotsome = (pc.nParticlesAtLevelPrePost[lev] > 0);
        }
        else
        {
            gotsome = (pc.NumberOfParticlesAtLevel(lev) > 0);
        }

        // We store the particles at each level in their own subdirectory.
        std::string LevelDir = pdir;

        if (gotsome)
        {
            if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') LevelDir += '/';

            LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);

            if ( ! pc.GetLevelDirectoriesCreated())
            {
                if (ParallelDescriptor::IOProcessor())
                    if ( ! amrex::UtilCreateDirectory(LevelDir, 0755))
                        amrex::CreateDirectoryFailed(LevelDir);
                ParallelDescriptor::Barrier();
            }
        }

        // Write out the header for each particle
        if (gotsome && ParallelDescriptor::IOProcessor()) {
            std::string HeaderFileName = LevelDir;
            HeaderFileName += "/Particle_H";
            std::ofstream ParticleHeader(HeaderFileName);

            pc.ParticleBoxArray(lev).writeOn(ParticleHeader);
            ParticleHeader << '\n';

            ParticleHeader.flush();
            ParticleHeader.close();
        }

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

        std::string filePrefix(LevelDir);
        filePrefix += '/';
        filePrefix += PC::DataPrefix();
        if(pc.usePrePost) {
            pc.filePrefixPrePost[lev] = filePrefix;
        }
        bool groupSets(false), setBuf(true);

        if (gotsome)
        {
            for(NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf); nfi.ReadyToWrite(); ++nfi)
            {
                std::ofstream& myStream = (std::ofstream&) nfi.Stream();
                pc.WriteParticles(lev, myStream, nfi.FileNumber(), which, count, where,
                                  write_real_comp, write_int_comp, particle_io_flags);
            }

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

        if (ParallelDescriptor::IOProcessor())
        {
            if(pc.GetUsePrePost()) {
                // ---- write to the header and unlink in CheckpointPost
            } else {
                for (int j = 0; j < state.size(); j++)
                {
                    HdrFile << which[j] << ' ' << count[j] << ' ' << where[j] << '\n';
                }

                if (gotsome && pc.doUnlink)
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
                            FileSystem::Remove(FullFileName);
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
}

template <class PC, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
void WriteBinaryParticleDataAsync (PC const& pc,
                                   const std::string& dir, const std::string& name,
                                   const Vector<int>& write_real_comp,
                                   const Vector<int>& write_int_comp,
                                   const Vector<std::string>& real_comp_names,
                                   const Vector<std::string>& int_comp_names)
{
    BL_PROFILE("WriteBinaryParticleDataAsync");
    AMREX_ASSERT(pc.OK());

    AMREX_ASSERT(sizeof(typename PC::ParticleType::RealType) == 4 ||
                 sizeof(typename PC::ParticleType::RealType) == 8);

    constexpr int NStructReal = PC::NStructReal;
    constexpr int NStructInt  = PC::NStructInt;
    constexpr int NArrayReal  = PC::NArrayReal;
    constexpr int NArrayInt   = PC::NArrayInt;

    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProcNumber = NProcs - 1;

    AMREX_ALWAYS_ASSERT(real_comp_names.size() == pc.NumRealComps() + NStructReal);
    AMREX_ALWAYS_ASSERT( int_comp_names.size() == pc.NumIntComps() + NStructInt);

    Vector<LayoutData<Long> > np_per_grid_local(pc.finestLevel()+1);
    for (int lev = 0; lev <= pc.finestLevel(); lev++)
    {
        np_per_grid_local[lev].define(pc.ParticleBoxArray(lev), pc.ParticleDistributionMap(lev));
        using ParIter = typename PC::ParConstIterType;
        for (ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            int gid = pti.index();
            const auto& ptile = pc.ParticlesAt(lev, pti);
            const auto& aos = ptile.GetArrayOfStructs();
            const auto pstruct = aos().dataPtr();
            const int np = ptile.numParticles();

            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            reduce_op.eval(np, reduce_data,
            [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
            {
                return (pstruct[i].id() > 0) ? 1 : 0;
            });

            int np_valid = amrex::get<0>(reduce_data.value(reduce_op));
            np_per_grid_local[lev][gid] += np_valid;
        }
    }

    Vector<Vector<Long> > np_per_grid_global(pc.finestLevel()+1);
    Long total_np = 0;
    Vector<Long> np_per_level(pc.finestLevel()+1);
    for (int lev = 0; lev <= pc.finestLevel(); lev++)
    {
        np_per_grid_global[lev].resize(np_per_grid_local[lev].size());
        ParallelDescriptor::GatherLayoutDataToVector(np_per_grid_local[lev],
                                                     np_per_grid_global[lev],
                                                     IOProcNumber);
        np_per_level[lev] = std::accumulate(np_per_grid_global[lev].begin(),
                                            np_per_grid_global[lev].end(), 0L);
        total_np += np_per_level[lev];
    }

    std::string pdir = dir;
    if ( ! pdir.empty() && pdir[pdir.size()-1] != '/') pdir += '/';
    pdir += name;

    if (MyProc == IOProcNumber)
    {
        if ( ! pc.GetLevelDirectoriesCreated())
        {
            if ( ! amrex::UtilCreateDirectory(pdir, 0755))
            {
                amrex::CreateDirectoryFailed(pdir);
            }
        }

        for (int lev = 0; lev <= pc.finestLevel(); lev++)
        {
            std::string LevelDir = pdir;
            bool gotsome = np_per_level[lev];

            if (gotsome)
            {
                if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') LevelDir += '/';

                LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);

                if ( ! pc.GetLevelDirectoriesCreated())
                {
                    if ( ! amrex::UtilCreateDirectory(LevelDir, 0755))
                    {
                        amrex::CreateDirectoryFailed(LevelDir);
                    }
                }

                std::string HeaderFileName = LevelDir;
                HeaderFileName += "/Particle_H";
                std::ofstream ParticleHeader(HeaderFileName);

                pc.ParticleBoxArray(lev).writeOn(ParticleHeader);
                ParticleHeader << '\n';

                ParticleHeader.flush();
                ParticleHeader.close();
            }
        }
    }
    ParallelDescriptor::Barrier();

    Long maxnextid = PC::ParticleType::NextID();
    ParallelDescriptor::ReduceLongMax(maxnextid, IOProcNumber);

    Vector<Long> np_on_rank(NProcs, 0L);
    std::size_t psize = particle_detail::PSizeInFile<ParticleReal>(write_real_comp, write_int_comp);
    Vector<int64_t> rank_start_offset(NProcs);
    if (MyProc == IOProcNumber)
    {
        for (int lev = 0; lev <= pc.finestLevel(); lev++)
        {
            for (int k = 0; k < pc.ParticleBoxArray(lev).size(); ++k)
            {
                int rank = pc.ParticleDistributionMap(lev)[k];
                np_on_rank[rank] += np_per_grid_global[lev][k];
            }
        }

        for (int ip = 0; ip < NProcs; ++ip)
        {
            auto info = AsyncOut::GetWriteInfo(ip);
            rank_start_offset[ip] = (info.ispot == 0) ? 0 : rank_start_offset[ip-1] + np_on_rank[ip-1]*psize;
        }
    }

    // make tmp particle tiles in pinned memory to write
    using PinnedPTile = ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt,
                                     PinnedArenaAllocator>;
    auto myptiles = std::make_shared<Vector<std::map<std::pair<int, int>,PinnedPTile> > >();
    myptiles->resize(pc.finestLevel()+1);
    for (int lev = 0; lev <= pc.finestLevel(); lev++)
    {
        for (MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto& new_ptile = (*myptiles)[lev][std::make_pair(mfi.index(),
                                                              mfi.LocalTileIndex())];

            if (np_per_grid_local[lev][mfi.index()] > 0)
            {
                const auto& ptile = pc.ParticlesAt(lev, mfi);
                new_ptile.resize(np_per_grid_local[lev][mfi.index()]);
                amrex::filterParticles(new_ptile, ptile, KeepValidFilter());
            }
        }
    }

    int finest_level = pc.finestLevel();
    Vector<BoxArray> bas;
    Vector<DistributionMapping> dms;
    for (int lev = 0; lev <= pc.finestLevel(); lev++)
    {
        bas.push_back(pc.ParticleBoxArray(lev));
        dms.push_back(pc.ParticleDistributionMap(lev));
    }

    int nrc = pc.NumRealComps();
    int nic = pc.NumIntComps();

    auto RD = pc.ParticleRealDescriptor;

    AsyncOut::Submit([=] ()
#if defined(__GNUC__) && (__GNUC__ == 8) && (__GNUC_MINOR__ == 1)
                     mutable // workaround for bug in gcc 8.1
#endif
    {
        if (MyProc == IOProcNumber)
        {
            std::string HdrFileName = pdir;
            std::ofstream HdrFile;

            if ( ! HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
                HdrFileName += '/';

            HdrFileName += "Header";

            HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);

            if ( ! HdrFile.good()) amrex::FileOpenFailed(HdrFileName);

            if (sizeof(typename PC::ParticleType::RealType) == 4)
            {
                HdrFile << PC::Version() << "_single" << '\n';
            }
            else
            {
                HdrFile << PC::Version() << "_double" << '\n';
            }

            int num_output_real = 0;
            for (int i = 0; i < nrc + NStructReal; ++i)
                if (write_real_comp[i]) ++num_output_real;

            int num_output_int = 0;
            for (int i = 0; i < nic + NStructInt; ++i)
                if (write_int_comp[i]) ++num_output_int;

            // AMREX_SPACEDIM and N for sanity checking.
            HdrFile << AMREX_SPACEDIM << '\n';

            // The number of extra real parameters
            HdrFile << num_output_real << '\n';

            // Real component names
            for (int i = 0; i < NStructReal + nrc; ++i )
                if (write_real_comp[i]) HdrFile << real_comp_names[i] << '\n';

            // The number of extra int parameters
            HdrFile << num_output_int << '\n';

            // int component names
            for (int i = 0; i < NStructInt + nic; ++i )
                if (write_int_comp[i]) HdrFile << int_comp_names[i] << '\n';

            bool is_checkpoint = true; // legacy
            HdrFile << is_checkpoint << '\n';

            // The total number of particles.
            HdrFile << total_np << '\n';

            // The value of nextid that we need to restore on restart.
            HdrFile << maxnextid << '\n';

            // Then the finest level of the AMR hierarchy.
            HdrFile << finest_level << '\n';

            // Then the number of grids at each level.
            for (int lev = 0; lev <= finest_level; lev++)
                HdrFile << dms[lev].size() << '\n';

            for (int lev = 0; lev <= finest_level; lev++)
            {
                Vector<int64_t> grid_offset(NProcs, 0);
                for (int k = 0; k < bas[lev].size(); ++k)
                {
                    int rank = dms[lev][k];
                    auto info = AsyncOut::GetWriteInfo(rank);
                    HdrFile << info.ifile << ' '
                            << np_per_grid_global[lev][k] << ' '
                            << grid_offset[rank] + rank_start_offset[rank] << '\n';
                    grid_offset[rank] += np_per_grid_global[lev][k]*psize;
                }
            }

            HdrFile.flush();
            HdrFile.close();
            if ( ! HdrFile.good())
            {
                amrex::Abort("ParticleContainer::Checkpoint(): problem writing HdrFile");
            }
        }

        AsyncOut::Wait();  // Wait for my turn

        for (int lev = 0; lev <= finest_level; lev++)
        {
            // For a each grid, the tiles it contains
            std::map<int, Vector<int> > tile_map;

            for (const auto& kv : (*myptiles)[lev])
            {
                const int grid = kv.first.first;
                const int tile = kv.first.second;
                tile_map[grid].push_back(tile);
            }

            std::string LevelDir = pdir;
            if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') LevelDir += '/';
            LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);
            std::string filePrefix(LevelDir);
            filePrefix += '/';
            filePrefix += PC::DataPrefix();
            auto info = AsyncOut::GetWriteInfo(MyProc);
            std::string file_name = amrex::Concatenate(filePrefix, info.ifile, 5);
            std::ofstream ofs;
            ofs.open(file_name.c_str(), (info.ispot == 0) ? (std::ios::binary | std::ios::trunc)
                     : (std::ios::binary | std::ios::app));

            for (int k = 0; k < bas[lev].size(); ++k)
            {
                int rank = dms[lev][k];
                if (rank != MyProc) continue;
                const int grid = k;
                if (np_per_grid_local[lev][grid] == 0) continue;

                // First write out the integer data in binary.
                int num_output_int = 0;
                for (int i = 0; i < nic + NStructInt; ++i)
                    if (write_int_comp[i]) ++num_output_int;

                const int iChunkSize = 2 + num_output_int;
                Vector<int> istuff(np_per_grid_local[lev][grid]*iChunkSize);
                int* iptr = istuff.dataPtr();

                for (unsigned i = 0; i < tile_map[grid].size(); i++) {
                    auto ptile_index = std::make_pair(grid, tile_map[grid][i]);
                    const auto& pbox = (*myptiles)[lev][ptile_index];
                    for (int pindex = 0;
                         pindex < pbox.GetArrayOfStructs().numParticles(); ++pindex)
                    {
                        const auto& aos = pbox.GetArrayOfStructs();
                        const auto& p = aos[pindex];

                        if (p.id() <= 0) continue;

                        // always write these
                        *iptr = p.id(); ++iptr;
                        *iptr = p.cpu(); ++iptr;

                        // optionally write these
                        for (int j = 0; j < NStructInt; j++)
                        {
                            if (write_int_comp[j])
                            {
                                *iptr = p.idata(j);
                                ++iptr;
                            }
                        }

                        const auto& soa  = pbox.GetStructOfArrays();
                        for (int j = 0; j < nic; j++)
                        {
                            if (write_int_comp[NStructInt+j])
                            {
                                *iptr = soa.GetIntData(j)[pindex];
                                ++iptr;
                            }
                        }
                    }
                }

                writeIntData(istuff.dataPtr(), istuff.size(), ofs);
                ofs.flush();  // Some systems require this flush() (probably due to a bug)

                // Write the Real data in binary.
                int num_output_real = 0;
                for (int i = 0; i < nrc + NStructReal; ++i)
                    if (write_real_comp[i]) ++num_output_real;

                const int rChunkSize = AMREX_SPACEDIM + num_output_real;
                Vector<typename PC::ParticleType::RealType> rstuff(np_per_grid_local[lev][grid]*rChunkSize);
                typename PC::ParticleType::RealType* rptr = rstuff.dataPtr();

                for (unsigned i = 0; i < tile_map[grid].size(); i++) {
                    auto ptile_index = std::make_pair(grid, tile_map[grid][i]);
                    const auto& pbox = (*myptiles)[lev][ptile_index];
                    for (int pindex = 0;
                         pindex < pbox.GetArrayOfStructs().numParticles(); ++pindex)
                    {
                        const auto& aos = pbox.GetArrayOfStructs();
                        const auto& p = aos[pindex];

                        if (p.id() <= 0) continue;

                        // always write these
                        for (int j = 0; j < AMREX_SPACEDIM; j++) rptr[j] = p.pos(j);
                        rptr += AMREX_SPACEDIM;

                        // optionally write these
                        for (int j = 0; j < NStructReal; j++)
                        {
                            if (write_real_comp[j])
                            {
                                *rptr = p.rdata(j);
                                ++rptr;
                            }
                        }

                        const auto& soa  = pbox.GetStructOfArrays();
                        for (int j = 0; j < nrc; j++)
                        {
                            if (write_real_comp[NStructReal+j])
                            {
                                *rptr = (typename PC::ParticleType::RealType) soa.GetRealData(j)[pindex];
                                ++rptr;
                            }
                        }
                    }
                }

                if (sizeof(typename PC::ParticleType::RealType) == 4) {
                    writeFloatData((float*) rstuff.dataPtr(), rstuff.size(), ofs, RD);
                }
                else if (sizeof(typename PC::ParticleType::RealType) == 8) {
                    writeDoubleData((double*) rstuff.dataPtr(), rstuff.size(), ofs, RD);
                }

                ofs.flush();  // Some systems require this flush() (probably due to a bug)
            }
        }
        AsyncOut::Notify();  // Notify others I am done
    });
}

#ifdef AMREX_USE_HDF5
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
                                  F&& f)
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
    hid_t fid, grp, fapl, dxpl, comp_dtype;
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

    Gpu::Device::synchronize();

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
            sprintf(setstripe, "lfs setstripe -c %d -S %dm %s", stripe_count, stripe_size, pdir.c_str());
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
        std::string versionName = Version();
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
                sprintf(comp_name, "real_component_%d", i);
                CreateWriteHDF5AttrString(fid, comp_name, real_comp_names[i].c_str());
            }
        }

        // The number of extra int parameters
        CreateWriteHDF5Attr(fid, "num_component_int", 1, &num_output_int, H5T_NATIVE_INT);

        // int component names
        for (int i = 0; i < NStructInt + pc.NumIntComps(); ++i ) {
            if (write_int_comp[i]) {
                sprintf(comp_name, "int_component_%d", i);
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
        for (int lev = 0; lev <= finest_level; ++lev) {
            sprintf(level_name, "level_%d", lev);

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
            hid_t mfs_dset= H5Dcreate(grp, "boxes", comp_dtype, mfs_dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

        H5Fclose(fid);
    }

    ParallelDescriptor::Barrier();
    BL_PROFILE_VAR_STOP(h5fg);

    // We want to write the data out in parallel.
    // We'll allow up to nOutFiles active writers at a time.
    int nOutFiles(-1);

    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
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
    // For HDF5 async VOL, block and wait previous tasks have all completed
    if (es_par_g != 0)
        async_vol_es_wait_particle();
    else {
        ExecOnFinalize(async_vol_es_wait_close_particle);
        es_par_g = H5EScreate();
    }
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
        sprintf(level_name, "level_%d", lev);
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

        bool groupSets(false), setBuf(true);

        if (gotsome)
        {
            pc.WriteParticlesHDF5(lev, grp, which, count, where,
                                  write_real_comp, write_int_comp, particle_io_flags);

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
}
#endif /*AMREX_USE_HDF5*/
#endif /*AMREX_WRITE_BINARY_PARTICLE_DATA_H*/
