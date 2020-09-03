#ifndef AMREX_WRITE_BINARY_PARTICLE_DATA_H
#define AMREX_WRITE_BINARY_PARTICLE_DATA_H

#include <AMReX_TypeTraits.H>
#include <AMReX_Particles.H>

// note - namespace

struct KeepValidFilter
{    
    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int i) const noexcept
    {
        return (src.m_aos[i].id() > 0);
    }
};

template <typename ParticleReal>
std::size_t PSizeInFile (const Vector<int>& wrc, const Vector<int>& wic)
{
    std::size_t rsize = sizeof(ParticleReal)*std::accumulate(wrc.begin(), wrc.end(), 0);
    std::size_t isize = sizeof(int)*std::accumulate(wic.begin(), wic.end(), 0);
    return rsize + isize + AMREX_SPACEDIM*sizeof(ParticleReal) + 2*sizeof(int);
}

template <class PC, class F, EnableIf_t<IsParticleContainer<PC>::value, int> foo = 0>
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
    if ( not pdir.empty() and pdir[pdir.size()-1] != '/') pdir += '/';
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
    Vector<std::map<std::pair<int, int>, Gpu::DeviceVector<int> > > particle_io_flags(pc.GetParticles().size());
    for (int lev = 0; lev < pc.GetParticles().size();  lev++)
    {
        const auto& pmap = pc.GetParticles(lev);
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

    if(pc.GetUsePrePost())
    {
        nparticles = pc.GetNParticlesPrePost();
        maxnextid  = pc.GetMaxNextIDPrePost();
    }
    else
    {
        nparticles = 0;
        maxnextid  = PC::ParticleType::NextID();

        for (int lev = 0; lev < pc.GetParticles().size();  lev++)
        {
            const auto& pmap = pc.GetParticles(lev);
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
            HdrFile << PC::ParticleType::Version() << "_single" << '\n';
        }
        else
        {
            HdrFile << PC::ParticleType::Version() << "_double" << '\n';
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
        if (gotsome and ParallelDescriptor::IOProcessor()) {
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
        filePrefix += PC::ParticleType::DataPrefix();
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

template <class PC, EnableIf_t<IsParticleContainer<PC>::value, int> foo = 0>
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
        for (MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int gid = mfi.index();
            const auto& ptile = pc.ParticlesAt(lev, mfi);
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

            int np_valid = amrex::get<0>(reduce_data.value());
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
    if ( not pdir.empty() and pdir[pdir.size()-1] != '/') pdir += '/';
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
    std::size_t psize = PSizeInFile<ParticleReal>(write_real_comp, write_int_comp);
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

            if (sizeof(typename PC::ParticleType) == 4)
            {
                HdrFile << PC::ParticleType::Version() << "_single" << '\n';
            }
            else
            {
                HdrFile << PC::ParticleType::Version() << "_double" << '\n';
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
            filePrefix += PC::ParticleType::DataPrefix();
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

#endif
