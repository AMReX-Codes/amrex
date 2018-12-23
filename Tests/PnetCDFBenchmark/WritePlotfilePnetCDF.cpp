/* This file is translated form Tests/HDF5Benchmark/WritePlotfileHDF5.cpp
 * A HDF5 variable maps to a NetCDF variable of the same size
 * A HDF attribute maps to a NetCDF attribute 
 * A HDF5 data space is mapped to a set of NetCDF dimensions representing every dimension of the data space
 * PnetCDF does not support composite datatype in attributes and variables. Luckily, all composite datatype 
 * used in this benchmark are sinmply repetiion of a single elementary datatype. We simply add one dimension 
 * to the variable to accomodate multiple values per cell.
 * Since NetCDF classic format does not have the concept of groups, groups are represented by appending the 
 * group name to the name of all member objects similar to a path. For example: variable X under group Y will 
 * be represented by a global variable called Y/X
 * PnetCDF have data and define mode. Objects can only be created in define mode and variable can only be 
 * accessed in data mode. We follow the flow of the HDF5 benchmark and switch mode as we need.
 */

#include <WritePlotfilePnetCDF.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <unistd.h>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

#ifdef BL_PNETCDF
#include <pnetcdf.h>
#endif

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, ncmpi_strerror(err));}}

using namespace amrex;

#ifdef BL_PNETCDF

namespace {
    
    /* Write std amp content as NetCDF attributes
     * m_int: IN: integer attributes 
     * m_real: IN: double attributes 
     * m_string: IN: string attributes 
     */
    int VWriteToLocation(int ncid, 
                            std::string prefix,
                            std::map<std::string, int>  &m_int,
                            std::map<std::string, Real> &m_real,
                            std::map<std::string, std::string> &m_string)
    {
        int err, status = NC_NOERR;
        std::string attname;
    
        // Write integer attribute
        for (std::map<std::string, int>::const_iterator p = m_int.begin(); p!= m_int.end(); ++p){
            attname = prefix + p->first;
            std::replace(attname.begin(), attname.end(), '/', '_'); 
            err = ncmpi_put_att_int(ncid, NC_GLOBAL, attname.c_str(), NC_INT, 1, &(p->second)); ERR
        }

        // Write real valued attribute
        for (std::map<std::string, double>::const_iterator p = m_real.begin(); p!= m_real.end(); ++p){
            double val;

            attname = std::string(prefix) + p->first;
            std::replace(attname.begin(), attname.end(), '/', '_'); 
            val = (double)(p->second);  // NetCDF has no specific datatype for Real, so we cast it as double
            err = ncmpi_put_att_double(ncid, NC_GLOBAL, attname.c_str(), NC_DOUBLE, 1, &(val)); ERR
        }

        // Write string attribute
        for (std::map<std::string, std::string>::const_iterator p = m_string.begin(); p!= m_string.end(); ++p){
            attname = std::string(prefix) + p->first;
            std::replace(attname.begin(), attname.end(), '/', '_'); 
            err = ncmpi_put_att_text(ncid, NC_GLOBAL, attname.c_str(), p->second.length(), p->second.c_str()); ERR
        }

        return status;
    }
}

/* Write NetCDF formated plot file
 */
void WriteMultiLevelPlotfilePNETCDF (const std::string &plotfilename,
                                  int nlevels,
				  const Vector<const MultiFab*> &mf,
				  const Vector<std::string> &varnames,
				  const Vector<Geometry> &geom,
				  Real time, Real dt,
				  const Vector<IntVect> &ref_ratio)
{
    BL_PROFILE("WriteMultiLevelPlotfilePNETCDF");
    
    int myProc(ParallelDescriptor::MyProc());
    int nProcs(ParallelDescriptor::NProcs());
    
    const int nComp = mf[0]->nComp();
    const int nGrow = 0;
    std::string filename(plotfilename + ".nc");

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime80(ParallelDescriptor::second());

    // ---- start PnetCDF part
    //herr_t  ret;
    int err, ret;
    IntVect iv1;
    Box b3(geom[0].Domain());

    int  b3int[2 * BL_SPACEDIM];
    for(int i(0); i < BL_SPACEDIM; ++i) {
        b3int[i] = b3.smallEnd(i);
        b3int[i + BL_SPACEDIM] = b3.bigEnd(i);
    }

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime741(ParallelDescriptor::second());
    double dPlotFileTime742(dPlotFileTime741 - dPlotFileTime80);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime742);
    if ( ParallelDescriptor::IOProcessor() )
    {
        //std::cout << "Write_H5M_time_1_0 = " << dPlotFileTime742 << "  seconds." << std::endl;
    }
    
#if BL_SPACEDIM == 1
    amrex::Abort("WriteMultiLevelPlotfilePNETCDF not implemented in 1d.");
#elif BL_SPACEDIM == 2
    amrex::Abort("WriteMultiLevelPlotfilePNETCDF not implemented in 2d.");
#endif
    
    // ASim@lbl.gov 6/15/2016
    double dPlotFileTime711(ParallelDescriptor::second());
    double dPlotFileTime712(dPlotFileTime711 - dPlotFileTime741);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime712);
    if ( ParallelDescriptor::IOProcessor() )
    {
        //std::cout << "Write_H5M_time_1_1 = " << dPlotFileTime712 << "  seconds." << std::endl;
    }
    
    std::string vGroupName = "_";
    int vFile;
    
    // Write one level at a time
    for (int level = 0; level < nlevels; ++level) {

    // Create the file only whne it is the first level, otherwise, open existing file
    if ( level == 0 ) {
        std::string filedescriptor("VanillaAMRFileType");

        err = ncmpi_create(ParallelDescriptor::Communicator(), filename.c_str(), NC_CLOBBER | NC_64BIT_DATA, MPI_INFO_NULL, &vFile); ERR

        std::string vGlobalGroupName = "Chombo_global";
        
        // ASim@lbl.gov 6/15/2017
        double dPlotFileTime713(ParallelDescriptor::second());
        double dPlotFileTime714(dPlotFileTime713 - dPlotFileTime711);
        ParallelDescriptor::ReduceRealMax(dPlotFileTime714);
        if(ParallelDescriptor::IOProcessor()) {
          std::cout << "Write_NCHEADER_time_1_2_0 = " << dPlotFileTime714 << "  seconds." << std::endl;
          std::cout << "#%$: ncmpi_create_time: " << dPlotFileTime714 << std::endl;
        }

        std::map<std::string, int>  vMInt;
        std::map<std::string, Real> vMReal;
        std::map<std::string, std::string> vMString;
        
        vMInt["SpaceDim"] = amrex::SpaceDim;
        vMReal["testReal"] = 0.0;
        vMString["testString"] = "vMString::testString";
        ret = VWriteToLocation(vFile, vGlobalGroupName, vMInt, vMReal, vMString);
        if(ret < 0) {
            std::cout << myProc << "**** Error 0:  ret = " << ret << std::endl;
        }
        
        // ASim@lbl.gov 6/15/2017
        double dPlotFileTime7131(ParallelDescriptor::second());
        double dPlotFileTime7141(dPlotFileTime7131 - dPlotFileTime713);
        ParallelDescriptor::ReduceRealMax(dPlotFileTime7141);
        if( ParallelDescriptor::IOProcessor() )
        {
            //std::cout << "Write_NCATT_time_1_2 = " << dPlotFileTime7141 << "  seconds." << std::endl;
        }
        
        vMInt.clear();
        vMReal.clear();
        
        vMString ["filetype"]    = filedescriptor;
        vMInt ["num_levels"]     = nlevels;
        vMInt ["num_components"] = nComp;
        for(int ivar(0); ivar < nComp; ++ivar) {
            char labelChSt[100];
            sprintf(labelChSt, "component_%d", ivar);
            std::string label(labelChSt);
            vMString[label] = varnames[ivar];
        }
        
        ret = VWriteToLocation(vFile, vGroupName, vMInt, vMReal, vMString);
        if(ret < 0) {
            std::cout << myProc << "**** Error 1:  ret = " << ret << std::endl;
        }
        
        // ASim@lbl.gov 6/15/2017
        double dPlotFileTime715(ParallelDescriptor::second());
        double dPlotFileTime716(dPlotFileTime715 - dPlotFileTime7131);
        ParallelDescriptor::ReduceRealMax(dPlotFileTime716);
        if( ParallelDescriptor::IOProcessor() )
        {
            std::cout << "Write_NCATT_time_1_3 = " << dPlotFileTime716 << "  seconds." << std::endl;
            std::cout << "#%$: file_init_att_time: " << dPlotFileTime716 + dPlotFileTime716 << std::endl;
        }
        
    } else {
        err = ncmpi_open(ParallelDescriptor::Communicator(), filename.c_str(), NC_WRITE, MPI_INFO_NULL, &vFile); ERR

        // ASim@lbl.gov 6/15/2017
        double dPlotFileTime713(ParallelDescriptor::second());
        double dPlotFileTime714(dPlotFileTime713 - dPlotFileTime711);
        ParallelDescriptor::ReduceRealMax(dPlotFileTime714);
        if( ParallelDescriptor::IOProcessor() )
        {
            std::cout << "Write_NCHEADER_time_1_4 = " << dPlotFileTime714 << "  seconds." << std::endl;
        }
    }
    
    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime71(ParallelDescriptor::second());
    double dPlotFileTime72(dPlotFileTime71 - dPlotFileTime80);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime72);
    if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Write_NCINIT_time_1 = " << dPlotFileTime72 << "  seconds." << std::endl;
        std::cout << "#%$: create_open_time: " << dPlotFileTime72 << std::endl;
    }
        
    char levelName[10];
    sprintf(levelName, "_level_%i", level);
    std::string gL(vGroupName + levelName);
    std::string gLDA(gL + "_data_attributes");
    err = ncmpi_put_att_int(vFile, NC_GLOBAL, (gLDA + "comps").c_str(), NC_INT, 1, &nComp); ERR

    IntVect giv;
    int gint[BL_SPACEDIM];
    for(int gi(0); gi < BL_SPACEDIM; ++gi) {
        giv[gi] = 0;
        gint[gi] = 0;
    }
    int gintvect_id[3];
    err = ncmpi_put_att_int(vFile, NC_GLOBAL, (gLDA + "_" + "ghost").c_str(), NC_INT, 3, gintvect_id); ERR
    
    const Real *a_dx = geom[level].CellSize();
    Real vData(dt);
    std::string vName("dt");
    
    err = ncmpi_put_att_double(vFile, NC_GLOBAL, (gL + "_" + vName).c_str(), NC_DOUBLE, 1, &vData); ERR
    err = ncmpi_put_att_double(vFile, NC_GLOBAL, (gL + "_" + "dx").c_str(), NC_DOUBLE, 1, &(a_dx[level])); ERR
    err = ncmpi_put_att_double(vFile, NC_GLOBAL, (gL + "_" + "time").c_str(), NC_DOUBLE, 1, &time); ERR
    err = ncmpi_put_att_int(vFile, NC_GLOBAL, (gL + "_" + "prob_domain").c_str(), NC_INT, 6, b3int); ERR

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime73(ParallelDescriptor::second());
    double dPlotFileTime74(dPlotFileTime73 - dPlotFileTime71);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime74);
    if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Write_NCATT_time_2 = " << dPlotFileTime74 << "  seconds." << std::endl;
        std::cout << "#%$: write_att_time: " << dPlotFileTime74 << std::endl;
    }
    
    // ---- "boxes" and "Processors" data
    Vector<int> procMap = mf[level]->DistributionMap().ProcessorMap();
    const BoxArray& grids = mf[level]->boxArray();
    int procdataset;
    int procdataspacedimids[1];
    int boxdataset;
    int boxdataspacedimids[2];
    int offsetdataset;
    int offsetdataspacedimids[1];
    std::string pdsname("Processors");
    std::string bdsname("boxes");
    std::string odsname("data:offsets=0");
    MPI_Offset flatdims[1], count[1], ocount[1];
    
    flatdims[0] = (MPI_Offset)grids.size();
    err = ncmpi_def_dim(vFile, "procdataspace_0", flatdims[0], procdataspacedimids); ERR
    err = ncmpi_def_var(vFile, (gL + "_" + pdsname).c_str(), NC_INT, 1, procdataspacedimids, &procdataset); ERR
 
    
    int bab3int[2 * BL_SPACEDIM];
    for(int i(0); i < BL_SPACEDIM; ++i) {
        bab3int[i] = 0;
        bab3int[i + BL_SPACEDIM] = 1;
    }
    int  boxSize(2 * BL_SPACEDIM);

    flatdims[0] = (MPI_Offset)grids.size();

    err = ncmpi_def_dim(vFile, "boxdataspace_0", flatdims[0], boxdataspacedimids); ERR
    err = ncmpi_def_dim(vFile, "boxdataspace_1", 6, boxdataspacedimids + 1); ERR
    err = ncmpi_def_var(vFile, (gL + "_" + bdsname).c_str(), NC_INT, 2, boxdataspacedimids, &boxdataset); ERR

    int iRefRatio(1);
    if(level < nlevels-1) {
        iRefRatio = ref_ratio[level][0];
    }
    err = ncmpi_put_att_int(vFile, NC_GLOBAL, (gL + "_" + "ref_ratio").c_str(), NC_INT, 1, &iRefRatio); ERR

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime75(ParallelDescriptor::second());
    double dPlotFileTime76(dPlotFileTime75 - dPlotFileTime73);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime76);
    if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Write_NCATT_time_3 = " << dPlotFileTime76 << "  seconds." << std::endl;
        std::cout << "#%$: def_meta_var_dim_time: " << dPlotFileTime76 << std::endl;
    }
    
    // ---- create a boxarray sorted by rank
    std::map<int, Vector<Box> > gridMap;
    for(int i(0); i < grids.size(); ++i) {
        int gridProc(procMap[i]);
        Vector<Box> &boxesAtProc = gridMap[gridProc];
        boxesAtProc.push_back(grids[i]);
    }
    BoxArray sortedGrids(grids.size());
    Vector<int> sortedProcs(grids.size());
    int bIndex(0);
    for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
        int proc = it->first;
        Vector<Box> &boxesAtProc = it->second;
        for(int ii(0); ii < boxesAtProc.size(); ++ii) {
            sortedGrids.set(bIndex, boxesAtProc[ii]);
            sortedProcs[bIndex] = proc;
            ++bIndex;
        }
    }
    
    MPI_Offset oflatdims[1];
    oflatdims[0] = (MPI_Offset)sortedGrids.size() + 1;
    err = ncmpi_def_dim(vFile, "offsetdataspace_0", oflatdims[0], offsetdataspacedimids); ERR
    err = ncmpi_def_var(vFile, (gL + "_" + odsname).c_str(), NC_INT64, 1, offsetdataspacedimids, &offsetdataset); ERR

    Vector<unsigned long long> offsets(sortedGrids.size() + 1);
    unsigned long long currentOffset(0L);
    ocount[0] = (MPI_Offset)sortedGrids.size() + 1;
    for(int b(0); b < sortedGrids.size(); ++b) {
        offsets[b] = currentOffset;
        currentOffset += sortedGrids[b].numPts() * nComp;
    }
    offsets[sortedGrids.size()] = currentOffset;
    
    Vector<unsigned long long> procOffsets(nProcs);
    int posCount(0);
    Vector<unsigned long long> procBufferSize(nProcs);
    unsigned long long totalOffset(0);
    for(auto it = gridMap.begin(); it != gridMap.end(); ++it) {
        int proc = it->first;
        Vector<Box> &boxesAtProc = it->second;
        BL_ASSERT(posCount == proc);
        procOffsets[posCount] = totalOffset;
        ++posCount;
        procBufferSize[proc] = 0L;
        for(int b(0); b < boxesAtProc.size(); ++b) {
            procBufferSize[proc] += boxesAtProc[b].numPts() * nComp;
        }
        totalOffset += procBufferSize[proc];
    }

    // Some metadata are stored in small variables, enter ddata mode to write variables
    err = ncmpi_enddef(vFile); ERR

    // Since only part of the processes write the variable, we use independent mode
    err = ncmpi_begin_indep_data(vFile); ERR
    
    if(ParallelDescriptor::IOProcessor()) {
        int vbCount(0);
        Vector<int> vbox(sortedGrids.size() * boxSize);
        Vector<int> pid(sortedGrids.size());
        count[0] = (MPI_Offset)sortedGrids.size();
        for(int b(0); b < sortedGrids.size(); ++b) {
            for(int i(0); i < BL_SPACEDIM; ++i) {
                vbox[(vbCount * boxSize) + i] = sortedGrids[b].smallEnd(i);
                vbox[(vbCount * boxSize) + i + BL_SPACEDIM] = sortedGrids[b].bigEnd(i);
            }
            ++vbCount;
            pid[b] = sortedProcs[b];
        }
        // ASim@lbl.gov CollectiveMetaData
	    // commented out 10/3/2017
        //ret = H5Pset_all_coll_metadata_ops(bmemdataspace, true);
        //ret = H5Pset_coll_metadata_write(bmemdataspace, true);
        //ret = H5Pset_all_coll_metadata_ops(pmemdataspace, true);
        //ret = H5Pset_coll_metadata_write(pmemdataspace, true);
        
        /*
        // ASim@lbl.gov 03/20/2017 for collective io setting: H5FD_MPIO_COLLECTIVE
        // Collective IO does not work here. H5FD_MPIO_INDEPENDENT = H5P_DEFAULT 
        // Is there a reason why this array needs to be written 
        // on this particular IOProcessor?
        // leaving along dxfer_template = H5P_DEFAULT;
        hid_t dxfer_template;
        dxfer_template = H5Pcreate(H5P_DATASET_XFER);
        // H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_INDEPENDENT);
        */
        
        if(vbox.size() > 0) {
            ret = ncmpi_put_var_ulonglong(vFile, offsetdataset, &(offsets[0])); ERR
            ret = ncmpi_put_var_int(vFile, boxdataset, &(vbox[0])); ERR
            ret = ncmpi_put_var_int(vFile, procdataset, &(pid[0])); ERR
        } else {
            /*
            hid_t dxfer_template;
            dxfer_template = H5P_DEFAULT;
            ret = H5Dwrite(offsetdataset, H5T_NATIVE_LLONG, omemdataspace, offsetdataspace,
                           dxfer_template, NULL);
            if(ret < 0) { std::cout << "_here 3:  ret = " << ret << std::endl; }
            ret = H5Dwrite(boxdataset, babox_id, bmemdataspace, boxdataspace,
                           dxfer_template, NULL);
            if(ret < 0) { std::cout << "_here 4:  ret = " << ret << std::endl; }
            ret = H5Dwrite(procdataset, H5T_NATIVE_INT, pmemdataspace, procdataspace,
                           dxfer_template, NULL);
            if(ret < 0) { std::cout << "_here 5:  ret = " << ret << std::endl; }
            */
        }
        
        /*
        // ASim@lbl.gov 03/20/2017 for closing collective io
        H5Pclose(dxfer_template);
        */
    }

    // Back to collective mode
    err = ncmpi_end_indep_data(vFile); ERR
    
    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime77(ParallelDescriptor::second());
    double dPlotFileTime78(dPlotFileTime77 - dPlotFileTime75);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime78);
    if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Write_NCVAR_time_4 = " << dPlotFileTime78 << "  seconds." << std::endl;
        std::cout << "#%$: write_meta_var_time: " << dPlotFileTime78 << std::endl;
    }
    
    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime85(ParallelDescriptor::second());
    double dPlotFileTime86(dPlotFileTime85 - dPlotFileTime80);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime86);
    if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Write_NCATT_time = " << dPlotFileTime86 << "  seconds." << std::endl;
        std::cout << "#%$: write_all_meta_time: " << dPlotFileTime86 << std::endl;
    }
    
    // Now it's time for the main data, we enter define mode to define variable for the main data
    err = ncmpi_redef(vFile); ERR

    {  // ---- data write
        BL_PROFILE_VAR("NCVarPut", h5dwd);

        MPI_Offset hs_procsize[1], hs_allprocsize[1], ch_offset[1];
        
        ch_offset[0]      = (MPI_Offset)procOffsets[myProc];          // ---- offset on this proc
        hs_procsize[0]    = (MPI_Offset)procBufferSize[myProc];       // ---- size of buffer on this proc
        hs_allprocsize[0] = (MPI_Offset)offsets[sortedGrids.size()];  // ---- size of buffer on all procs

        char dataname[1024];
        sprintf(dataname, "data:datatype=0");
        
        int dataspacedimmids[1];
        err = ncmpi_def_dim(vFile, "dataspace", hs_allprocsize[0], dataspacedimmids); ERR
        int dataset;
        err = ncmpi_def_var(vFile, (gL + "_" + dataname).c_str(), NC_DOUBLE, 1, dataspacedimmids, &dataset); ERR
        err = ncmpi_enddef(vFile); ERR  // Switch to data mode for writing the main data

        // ASim@lbl.gov CollectiveMetaData
	// commented out 10/3/2017
        //ret = H5Pset_all_coll_metadata_ops(dataset, true);
        //ret = H5Pset_coll_metadata_write(dataset, true);
        
        //select where in the file it will be written
        //H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, ch_offset, NULL,
        //                    hs_procsize, NULL);
        
        Vector<Real> a_buffer(procBufferSize[myProc], -1.0);
        long dataCount(0);
        for(MFIter mfi(*mf[level]); mfi.isValid(); ++mfi) {
            const Box &vbox    = mfi.validbox();
            const Real *dataPtr = (*mf[level])[mfi].dataPtr();
            for(int i(0); i < vbox.numPts() * nComp; ++i) {
                a_buffer[dataCount++] = dataPtr[i];
            }
        }
        if(ParallelDescriptor::IOProcessor()) {
            std::cout << "::---- calling NCPUT for the grid data on level " << level << std::endl;
        }
        
        /*
        // ASim@lbl.gov 4/10/2017 metadata sizing test 
        if(ParallelDescriptor::IOProcessor()) {
        // H5Fget_mdc_size(hid_t file_id, size_t *max_size_ptr, size_t *min_clean_size_ptr, size_t *cur_size_ptr, int *cur_num_entries_ptr)
        size_t max_size, min_clean_size, cur_size;
        int cur_num_entries;
        H5Fget_mdc_size(vFile, &max_size, &min_clean_size, &cur_size, &cur_num_entries);
        std::cout << "GET_MDC_SIZE = " << cur_size << std::endl;
        std::cout << "GET_MDC_SIZE_2 = " << max_size << " : " << min_clean_size << " : " << cur_num_entries << std::endl;
        }
        */
        
        // ASim@lbl.gov 6/15/2017
        double dPlotFileTime771(ParallelDescriptor::second());
        double dPlotFileTime781(dPlotFileTime771 - dPlotFileTime77);
        ParallelDescriptor::ReduceRealMax(dPlotFileTime781);
        if(ParallelDescriptor::IOProcessor()) {
            std::cout << "Write_NCATT_time_5 = " << dPlotFileTime781 << "  seconds." << std::endl;
            std::cout << "#%$: def_data_var_time: " << dPlotFileTime781 << std::endl;
        }

#ifdef H5INDEP
        // Switch to independent mode if we are writing independnetly
        err = ncmpi_begin_indep_data(vFile); ERR
#else
#endif
        
        // ASim@lbl.gov 6/15/2017
        double dPlotFileTime00(ParallelDescriptor::second());
        
        BL_PROFILE_VAR("NCVarPutGrids", h5dwg);

        // Write main data, PnetCDF has different sset of API for independnent and collective operation
#ifdef H5INDEP
        err = ncmpi_put_vara_double(vFile, dataset, ch_offset, hs_procsize, a_buffer.dataPtr());
#else
        err = ncmpi_put_vara_double_all(vFile, dataset, ch_offset, hs_procsize, a_buffer.dataPtr());
#endif
        BL_PROFILE_VAR_STOP(h5dwg);
        ERR
        //if(ret < 0) {
        //    std::cout << ParallelDescriptor::MyProc() << "_here 6:  ret = " << ret << std::endl;
        //}

#ifdef H5INDEP
        err = ncmpi_end_indep_data(vFile); ERR
#endif
        
        // ASim@lbl.gov 6/15/2017
        double dPlotFileTime11(ParallelDescriptor::second());
        double dPlotFileTime22(dPlotFileTime11 - dPlotFileTime00);
        ParallelDescriptor::ReduceRealMax(dPlotFileTime22);
        if(ParallelDescriptor::IOProcessor()) {
            std::cout << "Write_NCVARPUT_time = " << dPlotFileTime22 << "  seconds." << std::endl;
            std::cout << "Write_NCVARPUT_time_since = " << ParallelDescriptor::second() << std::endl;
            std::cout << "#%$: write_data_time: " << dPlotFileTime22 << std::endl;
        }
        
    	// ASim@lbl.gov 6/15/2017 for closing collective io
        //H5Pclose(dxfer_template); 

        BL_PROFILE_VAR_STOP(h5dwd);
    }
    
    err = ncmpi_close(vFile);
    }

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime791(ParallelDescriptor::second());

    // ASim@lbl.gov 6/15/2017
    double dPlotFileTime81(ParallelDescriptor::second());
    double dPlotFileTime82(dPlotFileTime81 - dPlotFileTime80);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime82);
    double dPlotFileTime792(dPlotFileTime81 - dPlotFileTime791);
    ParallelDescriptor::ReduceRealMax(dPlotFileTime791);
    if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Write_NC_time_7_closing = " << dPlotFileTime792 << "  seconds." << std::endl;
        std::cout << "Write_PNETCDF_time = " << dPlotFileTime82 << "  seconds." << std::endl;
        std::cout << "#%$: write_file_total_time: " << dPlotFileTime82 << std::endl;
    }    

}

#endif
