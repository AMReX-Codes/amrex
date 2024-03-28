#include <AMReX_VisMF.H>
#include <AMReX_AsyncOut.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FPC.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_ParmParse.H>

#include <fstream>
#include <iomanip>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#include <AMReX_PlotFileUtilOPENPMD.H>
#include <openPMD/openPMD.hpp>
#include <regex>

namespace amrex::openpmd_api {


    ////////////////////////////////////////
    //
    // Struct AMReX_VarNameParser
    //    parser var names to field and comp names
    //
    ////////////////////////////////////////

    AMReX_VarNameParser::AMReX_VarNameParser(std::string const& varname)
        :m_CompName(openPMD::MeshRecordComponent::SCALAR)
    {
        //auto [varname_no_mode, mode_index] = GetFieldNameModeInt(varname);
        GetFieldNameModeInt(varname);
        //bool var_in_theta_mode = mode_index != -1; // thetaMode or reconstructed Cartesian 2D slice
        //std::string m_FieldName = varname_no_mode;

        // assume fields are scalar unless they match the following match of known vector fields
    }

    void AMReX_VarNameParser::GetFieldNameModeInt (const std::string& varname)
    {
      // mode_index = -1 if varname isn't of form fieldName_mode_realOrImag
      // mode_index = 2 * mode - 1 + (realOrImag == 'imag')
      // in either case, there is a -1 in mode_index
      //int mode_index = -1;

      std::regex e_real_imag("(.*)_([0-9]*)_(real|imag)");
      std::smatch sm;
      std::regex_match(varname, sm, e_real_imag, std::regex_constants::match_default);

      if (sm.size() != 4 )
      {
          m_ThetaMode = (m_ModeIndex != -1); // thetaMode or reconstructed Cartesian 2D slice
          m_FieldName = varname;
      }
      else
      {
          // sm = [varname, field_name, mode, real_imag]
          int mode = std::stoi(sm[2]);
          if (mode == 0) {
              m_ModeIndex = 0;
          } else {
              if (sm[3] == "imag") {
                  m_ModeIndex += 1;
              }
              m_ModeIndex += 2 * mode;
          }
          m_ThetaMode = (m_ModeIndex != -1); // thetaMode or reconstructed Cartesian 2D slice
          m_FieldName = std::string(sm[1]);
      }
    }


    void AMReX_VarNameParser::GetMeshCompNames (int meshLevel)
    {
      std::string varname = m_FieldName;
      if (varname.size() >= 2U )
      {
          std::string const varname_1st = varname.substr(0U, 1U); // 1st character
          std::string const varname_2nd = varname.substr(1U, 1U); // 2nd character

          // Check if this field is a vector. If so, then extract the field name

          std::vector< std::string > const vector_fields = {"E", "B", "j"};
          std::vector< std::string > const field_components = getFieldComponentLabels();

          for( std::string const& vector_field : vector_fields )
          {
              for( std::string const& component : field_components )
              {
                  if( vector_field == varname_1st && component == varname_2nd )
                  {
                      m_FieldName = varname_1st + varname.substr(2); // Strip component
                      m_CompName  = varname_2nd;
                  }
              }
          }
      }


      if ( 0 == meshLevel ) {
          return;
      }

      m_FieldName += std::string("_lvl").append(std::to_string(meshLevel));
    }



    // create json options to pass to openpmd
    inline std::string
    getSeriesOptions (std::string const & operator_type,
                      std::map< std::string, std::string > const & operator_parameters,
                      std::string const & engine_type,
                      std::map< std::string, std::string > const & engine_parameters)
    {
      if (operator_type.empty() && engine_type.empty()) {
        return "{}";
      }

      std::string options;
      std::string top_block;
      std::string end_block;
      std::string op_block;
      std::string en_block;

      std::string op_parameters;
      for (const auto& kv : operator_parameters) {
        if (!op_parameters.empty()) { op_parameters.append(",\n"); }
        op_parameters.append(std::string(12, ' '))         /* just pretty alignment */
          .append("\"").append(kv.first).append("\": ")    /* key */
          .append("\"").append(kv.second).append("\""); /* value (as string) */
      }

      std::string en_parameters;
      for (const auto& kv : engine_parameters) {
        if (!en_parameters.empty()) { en_parameters.append(",\n"); }
        en_parameters.append(std::string(12, ' '))         /* just pretty alignment */
          .append("\"").append(kv.first).append("\": ")    /* key */
          .append("\"").append(kv.second).append("\""); /* value (as string) */
      }

      top_block = R"END(
{
  "adios2": {)END";
      end_block = R"END(
  }
})END";

      if (!operator_type.empty()) {
        op_block = R"END(
    "dataset": {
      "operators": [
        {
          "type": ")END";
        op_block += operator_type + "\"";
        if (!op_parameters.empty()) {
          op_block += R"END(,
          "parameters": {
)END";
          op_block += op_parameters +
            "\n          }";
        }
        op_block += R"END(
        }
      ]
    })END";

        if (!engine_type.empty() || !en_parameters.empty()) {
          op_block += ",";
        }
      }

      if (!engine_type.empty() || !en_parameters.empty())
        {
          en_block = R"END(
    "engine": {)END";
          if (!engine_type.empty()) {
            en_block += R"END(
      "type": ")END";
            en_block += engine_type + "\"";
            if(!en_parameters.empty()) {
              en_block += ",";
            }
          }
          if (!en_parameters.empty()) {
            en_block += R"END(
      "parameters": {
)END";
            en_block += en_parameters +
              "\n      }";
          }
          en_block += R"END(
    })END";
        }

      options = top_block + op_block + en_block + end_block;
      return options;
    }


    ////////////////////////////////////////
    //
    // Class AMReX_openPMDHandler
    //    files are saved as prefix/openpmd.bp
    //
    ////////////////////////////////////////

    AMReX_openPMDHandler::AMReX_openPMDHandler(std::string const& prefix) // NOLINT(modernize-pass-by-value) // match to diag_name in warpx
      :m_Writer(nullptr)
    {
      BL_PROFILE("AMReX_openPMDHandler::()");
      CreateWriter(prefix);
    }

    void AMReX_openPMDHandler::CreateWriter(const std::string& prefix)
    {
      ParmParse pp_prefix(prefix);

      // choose backend (e.g. ADIOS, ADIOS2 or HDF5). Default depends on openPMD-api configuration
      std::string openpmd_backend {"default"};
      pp_prefix.query("openpmd_backend", openpmd_backend);

      std::string  openpmd_encoding {"f"};
      pp_prefix.query("openpmd_encoding", openpmd_encoding);
      openPMD::IterationEncoding encoding = openPMD::IterationEncoding::groupBased;

      if ( openpmd_encoding == "v" ) {
        encoding = openPMD::IterationEncoding::variableBased;
      } else if ( openpmd_encoding == "g" ) {
        encoding = openPMD::IterationEncoding::groupBased;
      } else if ( openpmd_encoding == "f" ) {
        encoding = openPMD::IterationEncoding::fileBased;
      }

      auto lf_collect = [&](const char* key,
                            const std::string& parameter_tag,
                            std::string&  key_type,
                            std::map< std::string, std::string>& result)->void
      {
        //std::string key_type;
        pp_prefix.query(key, key_type);
        std::string const key_prefix = prefix + parameter_tag;
        ParmParse pp;
        auto entr = ParmParse::getEntries(key_prefix);

        auto const prefix_len = key_prefix.size() + 1;
        for (std::string k : entr) {
          std::string v;
          pp.get(k.c_str(), v);
          k.erase(0, prefix_len);
          result.insert({k, v});
        }
      };

      std::string  operator_type;
      std::map< std::string, std::string > operator_parameters;
      lf_collect("adios2_operator.type", ".adios2_operator.parameters", operator_type, operator_parameters);

      std::string  engine_type;
      std::map< std::string, std::string > engine_parameters;
      lf_collect("adios2_engine.type", ".adios2_engine.parameters", engine_type,  engine_parameters);

      std::string options=getSeriesOptions(operator_type, operator_parameters,
                                           engine_type, engine_parameters);

      m_Writer = std::make_unique<AMReX_openPMDWriter>(prefix, encoding, openpmd_backend, options);

      pp_prefix.query("file_min_digits", m_Writer->m_openPMDMinDigits);

    } // CreateWriter()


    void AMReX_openPMDHandler::SetWriter(amrex::openpmd_api::AMReX_openPMDWriter* w)
    {
       BL_ASSERT ( w != nullptr );

       // assuer that input key/values are inherited
       // so the openpmd filepath assigned from input file is still in use
       w->m_openPMDPrefix = m_Writer->m_openPMDPrefix;
       w->m_openPMDEncoding = m_Writer->m_openPMDEncoding;
       w->m_openPMDFileType = m_Writer->m_openPMDFileType;
       w->m_openPMDSeriesOptions = m_Writer->m_openPMDSeriesOptions;

       m_Writer.reset(w);
    }

    ////////////////////////////////////////
    //
    // Class AMReX_openPMDWriter
    //
    ////////////////////////////////////////

    AMReX_openPMDWriter::AMReX_openPMDWriter (std::string prefix,
                                              openPMD::IterationEncoding ie,
                                              std::string filetype,
                                              std::string options)
      :m_openPMDPrefix(std::move(prefix)),
       m_openPMDEncoding(ie),
       m_openPMDFileType(std::move(filetype)),
       m_openPMDSeriesOptions(std::move(options))
                                              //std::vector<bool> fieldPMLdirections // warpx specific
    {
        if( m_openPMDFileType == "default" ) {
#if openPMD_HAVE_ADIOS2==1
            m_openPMDFileType = "bp";
#elif openPMD_HAVE_ADIOS1==1
            m_openPMDFileType = "bp";
#elif openPMD_HAVE_HDF5==1
            m_openPMDFileType = "h5";
#else
            m_openPMDFileType = "json";
#endif
        }
    }

    AMReX_openPMDWriter::~AMReX_openPMDWriter ()
    {
      if( m_Series )
      {
          m_Series->flush();
          m_Series.reset( nullptr );
      }
    }

    void AMReX_openPMDWriter::SetStep(int ts)
    {
      m_CurrentStep = ts;

      Init(openPMD::Access::CREATE);
    }

    void AMReX_openPMDWriter::CloseStep(int /*ts*/)
    {
      if (m_Series) {
        GetIteration(m_CurrentStep).close();
      }
    }

    void AMReX_openPMDWriter::Init(openPMD::Access access)
    {
      std::string filepath = m_openPMDPrefix;
      GetFileName(filepath);

      if ( m_openPMDEncoding == openPMD::IterationEncoding::fileBased ) {
          m_Series = nullptr;
      } else if ( m_Series != nullptr ) {
          return;
      }

      if (amrex::ParallelDescriptor::NProcs() > 1)
        {
#if defined(AMREX_USE_MPI)
        m_Series = std::make_unique<openPMD::Series>(
                                                     filepath, access,
                                                     amrex::ParallelDescriptor::Communicator(),
                                                     m_openPMDSeriesOptions
                                                     );
#else
        amrex::Abort(Utils::TextMsg::Err("AMReX did not build with MPI support!"));
#endif
        }
      else
        {
          m_Series = std::make_unique<openPMD::Series>(filepath, access, m_openPMDSeriesOptions);
        }

      m_Series->setIterationEncoding( m_openPMDEncoding );

      m_Series->setMeshesPath( "fields" );
      // conform to ED-PIC extension of openPMD

      uint32_t const openPMD_ED_PIC = 1U;
      m_Series->setOpenPMDextension( openPMD_ED_PIC );
      // meta info

      m_Series->setSoftware( "AMReX", amrex::Version() );
    }

    void AMReX_openPMDWriter::CompSetup(int lev,
                                        openPMD::Container< openPMD::Mesh >& meshes,
                                        amrex::Geometry& full_geom,
                                        const std::vector<std::string>& varnames,
                                        const amrex::MultiFab* curr_mf) const
    {
      int const ncomp = curr_mf->nComp();
      for ( int icomp=0; icomp<ncomp; icomp++ )
        {
          std::string const & varname = varnames[icomp];
          amrex::openpmd_api::AMReX_VarNameParser curr(varname);
          curr.GetMeshCompNames( lev );
          {
            if (curr.m_CompName == openPMD::MeshRecordComponent::SCALAR)
              {
                if ( ! meshes.contains(curr.m_FieldName) )
                  {
                    auto mesh = meshes[curr.m_FieldName];
                    SetupMeshComp(  mesh, full_geom, *curr_mf, curr );
                  }
              }
            else
              {
                auto mesh = meshes[curr.m_FieldName];
                if ( ! mesh.contains(curr.m_CompName) ) {
                  SetupMeshComp(  mesh, full_geom, *curr_mf, curr );
                }
              }
          }
        } // icomp setup loop
    }

    void AMReX_openPMDWriter::CompStorage(int lev,
                                          openPMD::Container< openPMD::Mesh >& meshes,
                                          amrex::Geometry& full_geom,
                                          const std::vector<std::string>& varnames,
                                          const amrex::MultiFab* curr_mf)
    {
      int const ncomp = curr_mf->nComp();
      amrex::Box const & global_box = full_geom.Domain();


      for ( int icomp=0; icomp<ncomp; icomp++ )
        {
          std::string const & varname = varnames[icomp];

          amrex::openpmd_api::AMReX_VarNameParser curr(varname);
          curr.GetMeshCompNames( lev );

          auto mesh = meshes[curr.m_FieldName];
          auto mesh_comp = mesh[curr.m_CompName];

          for( amrex::MFIter mfi(*curr_mf); mfi.isValid(); ++mfi )
            {
              amrex::FArrayBox const& fab = (*curr_mf)[mfi];
              // TODO: fab.box() shows ghost cells, while validbox does not
              //       what should I use? dataPtr() covers ghost cell or not?
              //     NOTE that getReversedVec() turns everything into uint first.
              //amrex::Box const& local_box = fab.box();
              amrex::Box const& local_box = mfi.validbox();

              // Determine the offset and size of this chunk
              amrex::IntVect const box_offset = local_box.smallEnd() - global_box.smallEnd();
              auto chunk_offset = helper::getReversedVec( box_offset );
              auto chunk_size = helper::getReversedVec( local_box.size() );

              if (curr.m_ThetaMode)
                {
                  chunk_offset.emplace(chunk_offset.begin(), curr.m_ModeIndex);
                  chunk_size.emplace(chunk_size.begin(), 1);
                }
#ifdef AMREX_USE_GPU
              if (fab.arena()->isManaged() || fab.arena()->isDevice())
                {
                  amrex::BaseFab<amrex::Real> foo(local_box, 1, amrex::The_Pinned_Arena());
                  std::shared_ptr<amrex::Real> data_pinned(foo.release());
                  amrex::Gpu::dtoh_memcpy_async(data_pinned.get(), fab.dataPtr(icomp), local_box.numPts()*sizeof(amrex::Real));
                  // intentionally delayed until before we .flush(): amrex::Gpu::streamSynchronize();
                  mesh_comp.storeChunk(data_pinned, chunk_offset, chunk_size);
                }
              else
#endif
                {
                  amrex::Real const *local_data = fab.dataPtr(icomp);
                  mesh_comp.storeChunkRaw(local_data, chunk_offset, chunk_size);
                }
            }
        } // icomp store loop

    }

    void AMReX_openPMDWriter::WriteMesh(const std::vector<std::string>& varnames,
                                        const amrex::Vector<const amrex::MultiFab*>& mf,
                                        const amrex::Vector<amrex::Geometry>& geom,
                                        //const int iteration,
                                        const double time ) const

    {
      BL_PROFILE("AMReX_openPMDWriter::WriteMesh()");
      openPMD::Iteration series_iteration = GetIteration(m_CurrentStep);
      series_iteration.open();

      auto meshes = series_iteration.meshes;
      series_iteration.setTime( time );

      if ( varnames.empty() ) { return; }

      auto output_levels = int(geom.size());
      for (int lev=0; lev < output_levels; lev++)
        {
          amrex::Geometry full_geom = geom[lev];

          if ( 0 == lev ) {
            SetupFields(meshes, full_geom);
          }

          CompSetup(lev, meshes, full_geom, varnames, mf[lev]);
          CompStorage(lev, meshes, full_geom, varnames, mf[lev]);
#ifdef AMREX_USE_GPU
          amrex::Gpu::streamSynchronize();
#endif
          m_Series->flush();
      } // for lev loop
    }

    void AMReX_openPMDWriter::GetFileName(std::string& filepath) const
    {
      if (filepath.empty()) {
        filepath.append(".");
      }

      filepath.append("/");
      // transform paths for Windows
#ifdef _WIN32
      filepath = detail::replace_all(filepath, "/", "\\");
#endif

      std::string filename = "openpmd";
      //
      // OpenPMD supports timestepped names
      //
      if (m_openPMDEncoding == openPMD::IterationEncoding::fileBased)
        {
          std::string fileSuffix = std::string("_%0") + std::to_string(m_openPMDMinDigits) + std::string("T");
          filename = filename.append(fileSuffix);
        }
      filename.append(".").append(m_openPMDFileType);
      filepath.append(filename);
    }

    void AMReX_openPMDWriter::SetupFields (openPMD::Container< openPMD::Mesh >& meshes,
                                           amrex::Geometry& full_geom) const
    {
      //}
      // meta data for ED-PIC extension
      auto const period = full_geom.periodicity();

      std::vector<std::string> fieldBoundary(6, "reflecting");
      std::vector<std::string> particleBoundary(6, "absorbing");
      fieldBoundary.resize(AMREX_SPACEDIM * 2);
      particleBoundary.resize(AMREX_SPACEDIM * 2);

#if AMREX_SPACEDIM != 3
      fieldBoundary.resize(4);
      particleBoundary.resize(4);
#endif

      for (int i = 0; i < int(fieldBoundary.size() / 2); ++i) {
        if (period.isPeriodic(i)) {
          fieldBoundary.at(2 * i) = "periodic";
          fieldBoundary.at(2 * i + 1) = "periodic";
          particleBoundary.at(2 * i) = "periodic";
          particleBoundary.at(2 * i + 1) = "periodic";
        }
      }

      meshes.setAttribute("fieldBoundary", fieldBoundary);
      meshes.setAttribute("particleBoundary", particleBoundary);
    }


    void AMReX_openPMDWriter::SetupMeshComp (openPMD::Mesh& mesh,
                                             const amrex::Geometry& full_geom,
                                             amrex::MultiFab const& mf,
                                             const AMReX_VarNameParser& varName
                                             ) const
    {
      BL_PROFILE("SetupMeshComp(default)");

      auto mesh_comp = mesh[varName.m_CompName];
      amrex::Box const & global_box = full_geom.Domain();
      auto global_size = helper::getReversedVec(global_box.size() );

      // - Grid spacing
      std::vector<double> const grid_spacing = helper::getReversedVec(full_geom.CellSize());
      mesh.setGridSpacing(grid_spacing);

      // - Global offset
      std::vector<double> const global_offset = helper::getReversedVec(full_geom.ProbLo());
      mesh.setGridGlobalOffset(global_offset);

      // - AxisLabels
      std::vector<std::string> axis_labels = varName.getFieldAxisLabels();
      mesh.setAxisLabels(axis_labels);

      // Prepare the type of dataset that will be written
      openPMD::Datatype const datatype = openPMD::determineDatatype<amrex::Real>();
      auto const dataset = openPMD::Dataset(datatype, global_size);
      mesh.setDataOrder(openPMD::Mesh::DataOrder::C);

      if (varName.m_ThetaMode) {
        mesh.setGeometry("thetaMode");
      }

      mesh.setAttribute("fieldSmoothing", "none");
      mesh_comp.resetDataset(dataset);

      auto relative_cell_pos = helper::getRelativeCellPosition(mf);     // AMReX Fortran index order
      std::reverse( relative_cell_pos.begin(), relative_cell_pos.end() ); // now in C order
      mesh_comp.setPosition( relative_cell_pos );
    }


} // namespace amrex::openpmd_api
