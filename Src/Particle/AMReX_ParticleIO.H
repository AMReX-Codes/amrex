#ifndef BL_PARTICLEIO_H
#define BL_PARTICLEIO_H

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
                            tmp_real_comp_names, tmp_int_comp_names);
    
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
                            real_comp_names, int_comp_names);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir, const std::string& name) const
{
    Vector<int> write_real_comp;
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NArrayReal; ++i )
    {
        write_real_comp.push_back(1);
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }
    
    Vector<int> write_int_comp;
    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NArrayInt; ++i )
    {
        write_int_comp.push_back(1);
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }
    
    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names,
                 const Vector<std::string>& int_comp_names) const
{    
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NArrayReal);
    AMREX_ASSERT( int_comp_names.size() == NStructInt  + NArrayInt );

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NArrayReal; ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NArrayInt; ++i) write_int_comp.push_back(1);

    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names);        
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names) const
{    
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NArrayReal);

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NArrayReal; ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NArrayInt; ++i) write_int_comp.push_back(1);

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NArrayInt; ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }
        
    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names);        
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFile (const std::string& dir,
                 const std::string& name,
                 const Vector<int>& write_real_comp,
                 const Vector<int>& write_int_comp) const
{    
    AMREX_ASSERT(write_real_comp.size() == NStructReal + NArrayReal);
    AMREX_ASSERT(write_int_comp.size()  == NStructInt  + NArrayInt );
    
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NArrayReal; ++i )
    {
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NArrayInt; ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names);
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
                            real_comp_names, int_comp_names);
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WriteBinaryParticleData (const std::string& dir, const std::string& name,
                           const Vector<int>& write_real_comp,
                           const Vector<int>& write_int_comp,
                           const Vector<std::string>& real_comp_names,
                           const Vector<std::string>& int_comp_names) const
{
    BL_PROFILE("ParticleContainer::WriteBinaryParticleData()");
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
    pdir += name;
    
    if ( ! levelDirectoriesCreated)
    {
        if (ParallelDescriptor::IOProcessor()) 
            if ( ! amrex::UtilCreateDirectory(pdir, 0755)) 
                amrex::CreateDirectoryFailed(pdir);
        ParallelDescriptor::Barrier();
    }
    
    std::ofstream HdrFile;
    
    long nparticles = 0;
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
        Vector<long> where(state.size(),0);
	
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
                //
                // Write out all the valid particles we own at the specified level.
                // Do it grid block by grid block remembering the seek offset
                // for the start of writing of each block of data.
                //
                WriteParticles(lev, myStream, nfi.FileNumber(), which, count, where,
                               write_real_comp, write_int_comp);
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
                    //
                    // We now write the which file, the particle count, and the
                    // file offset into which the data for each grid was written,
                    // to the header file.
                    //
                    HdrFile << which[j] << ' ' << count[j] << ' ' << where[j] << '\n';
                }
                
                if (gotsome && doUnlink)
                {
//		BL_PROFILE_VAR("PC<NNNN>::Checkpoint:unlink", unlink);
                    //
                    // Unlink any zero-length data files.
                    //
                    Vector<long> cnt(nOutFiles,0);
                    
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
    }            // ---- end for(lev...)

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
    long nparticles = 0;
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
                Vector<long> cnt(nOutFilesPrePost,0);
                
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
                  Vector<int>& which, Vector<int>& count, Vector<long>& where,
                  const Vector<int>& write_real_comp,
                  const Vector<int>& write_int_comp) const
{
    BL_PROFILE("ParticleContainer::WriteParticles()");

    // For a each grid, the tiles it contains
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
            const auto& pbox = m_particles[lev].at(std::make_pair(grid, tile_map[grid][i]));
            for (int pindex = 0; pindex < pbox.GetArrayOfStructs().size(); ++pindex) {
                const ParticleType& p = pbox.GetArrayOfStructs()[pindex];
                if (p.m_idata.id > 0)
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
            const auto& pbox = m_particles[lev].at(std::make_pair(grid, tile_map[grid][i]));
            for (int pindex = 0; pindex < pbox.GetArrayOfStructs().size(); ++pindex) {
                const ParticleType& p = pbox.GetArrayOfStructs()[pindex];
                if (p.m_idata.id > 0)
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
    BL_ASSERT(!dir.empty());
    BL_ASSERT(!file.empty());
    
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
    BL_ASSERT(!version.empty());
    
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
    
    long nparticles;
    HdrFile >> nparticles;
    BL_ASSERT(nparticles >= 0);
    
    int maxnextid;
    HdrFile >> maxnextid;
    BL_ASSERT(maxnextid > 0);
    ParticleType::NextID(maxnextid);
    
    int finest_level_in_file;
    HdrFile >> finest_level_in_file;
    BL_ASSERT(finest_level_in_file >= 0);
    
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
        BL_ASSERT(ngrids[lev] > 0);
        if (lev <= finestLevel()) {
            BL_ASSERT(ngrids[lev] == int(ParticleBoxArray(lev).size()));
        }
    }
    
    resizeData();
    
    if (finest_level_in_file > finestLevel()) {
        m_particles.resize(finest_level_in_file+1);
    }
    
    for (int lev = 0; lev <= finest_level_in_file; lev++) {
        Vector<int>  which(ngrids[lev]);
        Vector<int>  count(ngrids[lev]);
        Vector<long> where(ngrids[lev]);
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
    
    BL_ASSERT(OK());
    
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
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < int(m_particles.size()));

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
                    std::vector<Cuda::HostVector<Real> > > > host_real_attribs;
    host_real_attribs.reserve(15);
    host_real_attribs.resize(finest_level_in_file+1);

    Vector<std::map<std::pair<int, int>,
                    std::vector<Cuda::HostVector<int> > > > host_int_attribs;
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

        BL_ASSERT(p.m_idata.id > 0);
        
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
                
	  Cuda::thrust_copy(src_tile.begin(),
			    src_tile.end(),
			    dst_tile.GetArrayOfStructs().begin() + old_size);
	  
	  for (int i = 0; i < NumRealComps(); ++i) {
              Cuda::thrust_copy(host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                                host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                                dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
	  }
	  
	  for (int i = 0; i < NumIntComps(); ++i) {
              Cuda::thrust_copy(host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
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
    BL_ASSERT(!filename.empty());

    const Real strttime = amrex::second();
    //
    // Count # of valid particles.
    //
    long nparticles = 0;

    for (int lev = 0; lev < m_particles.size();  lev++) {
        auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
	    auto np = aos.numParticles();
	    Gpu::HostVector<ParticleType> host_aos(np);
	    Gpu::thrust_copy(aos.begin(), aos.end(), host_aos.begin());
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
		Gpu::thrust_copy(aos.begin(), aos.end(), host_aos.begin());

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
    BL_ASSERT(!filename.empty());

    const Real strttime = amrex::second();
 
    //
    // Count # of valid particles.
    //
    long nparticles = 0;

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
#endif /*BL_PARTICLEIO_H*/
