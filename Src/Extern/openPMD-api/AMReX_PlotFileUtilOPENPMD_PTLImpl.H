#ifndef AMREX_PLOTFILE_OPENPMD_PTLS_API_H
#define AMREX_PLOTFILE_OPENPMD_PTLS_API_H

//extern  std::unique_ptr< amrex::openpmd_api::AMReX_openPMDHandler > amrex::openpmd_api::m_OpenPMDHandler;
namespace amrex::openpmd_api
{
    extern  std::unique_ptr< amrex::openpmd_api::AMReX_openPMDHandler > m_OpenPMDHandler;
    ////////////////////////////////////////////////////////////////////////
    //
    // Global function to write particle container
    //   each WriteParticles() call has two versions:
    //   with FUNC(for user supploed saveposid) or not (using default from the writer)
    //
    ////////////////////////////////////////////////////////////////////////

    template<class PC, typename FUNC_pc, typename FUNC_pit>
    void WriteParticles(PC const&  pc,
                        const std::string& specieName,
                        //int iteration,
                        const Vector<std::string>& real_comp_names,
                        const Vector<std::string>& int_comp_names,
                        const Vector<int>& do_write_real_comp,
                        const Vector<int>& do_write_int_comp,
                        FUNC_pc&& pcFunc,
                        FUNC_pit&& posIdFunc)
    {
      BL_PROFILE("amrex::openpmd_api::WriteParticles()");

      amrex::Print()<<"   Real: "<<pc.NStructReal<<" "<<pc.NumRealComps()<<"\n";
      amrex::Print()<<"   Int:  "<<pc.NStructInt<<" "<<pc.NumIntComps()<<"\n";

      if ((m_OpenPMDHandler == nullptr) || (m_OpenPMDHandler->m_Writer == nullptr)) {
        return;
      }

      AMREX_ASSERT( do_write_real_comp.size() <= real_comp_names.size() );
      AMREX_ASSERT( do_write_int_comp.size()  <= int_comp_names.size()  );
      AMREX_ASSERT( real_comp_names.size()    <= pc.NStructReal + pc.NumRealComps() );

      m_OpenPMDHandler->m_Writer->DumpParticles(pc, specieName, //iteration,
                                                do_write_real_comp, do_write_int_comp, real_comp_names, int_comp_names,
                                                pcFunc, posIdFunc);
    }

    // use default pos_id func
    template<class PC>
    void WriteParticles(PC const&  pc,
                        const std::string& specieName,
                        //int iteration,
                        const Vector<std::string>& real_comp_names,
                        const Vector<std::string>& int_comp_names,
                        const Vector<int>& do_write_real_comp,
                        const Vector<int>& do_write_int_comp)
    {
      BL_PROFILE("amrex::openpmd_api::WriteParticles(..)");

      amrex::Print()<<"   Real: "<<pc.NStructReal<<" "<<pc.NumRealComps()<<"\n";
      amrex::Print()<<"   Int:  "<<pc.NStructInt<<" "<<pc.NumIntComps()<<"\n";

      if ((m_OpenPMDHandler == nullptr) || (m_OpenPMDHandler->m_Writer == nullptr)) {
        return;
      }

      AMREX_ASSERT( do_write_real_comp.size() <= real_comp_names.size() );
      AMREX_ASSERT( do_write_int_comp.size()  <= int_comp_names.size()  );
      AMREX_ASSERT( real_comp_names.size()    <= pc.NStructReal + pc.NumRealComps() );

      m_OpenPMDHandler->m_Writer->DumpParticles(pc, specieName, //iteration,
                                                do_write_real_comp, do_write_int_comp, real_comp_names, int_comp_names);
    }

    template<class PC>
    void WriteParticles(PC const&  pc,
                        const std::string& specieName,
                        //int iteration,
                        const Vector<std::string>& real_comp_names,
                        const Vector<std::string>& int_comp_names)
    {
      Vector<int> write_real_comp;

      for (int i = 0; i < real_comp_names.size(); ++i )
        {
          write_real_comp.push_back(1);
        }

      Vector<int> write_int_comp;
      for (int i = 0; i < int_comp_names.size(); ++i )
        {
          write_int_comp.push_back(1);
        }

      WriteParticles(pc, specieName, //iteration,
                     real_comp_names, int_comp_names,
                     write_real_comp, write_int_comp);
    }


    template<class PC>
    void WriteParticles(PC const&  pc,
                        const std::string& specieName)
    {
      Vector<int> write_real_comp;
      Vector<std::string> real_comp_names;

      for (int i = 0; i < pc.NStructReal + pc.NumRealComps(); ++i )
        {
          write_real_comp.push_back(1);

          std::stringstream ss;
          if (i < pc.NStructReal) {
            ss << "aos-r_" << i;
          } else {
            ss << "SoA-r" << i;
          }

          real_comp_names.push_back(ss.str());
        }

      Vector<int> write_int_comp;
      Vector<std::string> int_comp_names;
      for (int i = 0; i < pc.NStructInt + pc.NumIntComps(); ++i )
        {
          write_int_comp.push_back(1);
          std::stringstream ss;

          if (i < pc.NStructInt) {
            ss << "aos-i_" << i;
          } else {
            ss << "SoA-i" << i;
          }

          int_comp_names.push_back(ss.str());
        }

      WriteParticles(pc, specieName, //iteration,
                     real_comp_names, int_comp_names,
                     write_real_comp, write_int_comp);
    }


    template<class PC, typename FUNC_pc, typename FUNC_pit>
    void WriteParticles(PC const&  pc,
                        const std::string& specieName,
                        FUNC_pc&& pcFunc,
                        FUNC_pit&& pitFunc
                        )
    {
      Vector<int> write_real_comp;
      Vector<std::string> real_comp_names;

      for (int i = 0; i < pc.NStructReal + pc.NumRealComps(); ++i )
        {
          write_real_comp.push_back(1);

          std::stringstream ss;
          if (i < pc.NStructReal) {
            ss << "aos-r_" << i;
          } else {
            ss << "SoA-r" << i;
          }

          real_comp_names.push_back(ss.str());
        }

      Vector<int> write_int_comp;
      Vector<std::string> int_comp_names;
      for (int i = 0; i < pc.NStructInt + pc.NumIntComps(); ++i )
        {
          write_int_comp.push_back(1);
          std::stringstream ss;

          if (i < pc.NStructInt) {
            ss << "aos-i_" << i;
          } else {
            ss << "SoA-i" << i;
          }
          int_comp_names.push_back(ss.str());
        }

      WriteParticles(pc, specieName, //iteration,
                     real_comp_names, int_comp_names,
                     write_real_comp, write_int_comp, pcFunc, pitFunc);
    }

    namespace helper {
      ////////////////////////////////////
      //
      // all helper functions
      //
      ////////////////////////////////////
      inline std::map< openPMD::UnitDimension, double >
      getUnitDimension ( std::string const & record_name )
      {

        if( record_name == "position" ) {
            return {{openPMD::UnitDimension::L,  1.}};
        } else if( record_name == "positionOffset" ) {
            return {{openPMD::UnitDimension::L,  1.}};
        } else if( record_name == "momentum" ) {
            return {{openPMD::UnitDimension::L,  1.},
                    {openPMD::UnitDimension::M,  1.},
                    {openPMD::UnitDimension::T, -1.}};
        } else if( record_name == "charge" ) {
            return {{openPMD::UnitDimension::T,  1.},
                    {openPMD::UnitDimension::I,  1.}};
        } else if( record_name == "mass" ) {
            return {{openPMD::UnitDimension::M,  1.}};
        } else if( record_name == "E" ) {
            return {{openPMD::UnitDimension::L,  1.},
                    {openPMD::UnitDimension::M,  1.},
                    {openPMD::UnitDimension::T, -3.},
                    {openPMD::UnitDimension::I, -1.}};
        } else if( record_name == "B" ) {
            return {{openPMD::UnitDimension::M,  1.},
                    {openPMD::UnitDimension::I, -1.},
                    {openPMD::UnitDimension::T, -2.}};
        } else {
            return {};
        }
      }

      //
      //
      //
      constexpr uint64_t
      localIDtoGlobal(int const id, int const cpu)
      {
        // from WarpXUtilIO::localIDtoGlobal
        static_assert(sizeof(int) * 2U <= sizeof(uint64_t),
                      "int size might cause collisions in global IDs");

        return uint64_t(id) | uint64_t(cpu) << 32U;
      }


       inline std::pair< std::string, std::string >
       name2openPMD ( std::string const& fullName )
       {
         std::string record_name = fullName;
         std::string component_name = openPMD::RecordComponent::SCALAR;

         // we use "_" as separator in names to group vector records
         std::size_t startComp = fullName.find_last_of('_');
         if( startComp != std::string::npos ) {  // non-scalar
           record_name = fullName.substr(0, startComp);
           component_name = fullName.substr(startComp + 1U);
         }
         return make_pair(record_name, component_name);
       }

      //
      // used by both setup/saverealproperty so better be here than a lambda function
      //
      inline
      auto  getComponentRecord (openPMD::ParticleSpecies& currSpecies,
                                std::string const& comp_name)
      {
        // handle scalar and non-scalar records by name
        const auto [record_name, component_name] = name2openPMD(comp_name);
        return currSpecies[record_name][component_name];
      }

    } // namespace helper






    template<class PIt>
    void
    AMReX_openPMDWriter::StoreAoS_Real(PIt& pti,
                                       openPMD::ParticleSpecies& currSpecies,
                                       const amrex::Vector<int>& write_real_comp,
                                       const amrex::Vector<std::string>& real_comp_names,
                                       unsigned long long offset) const
    {
      auto const& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile
      auto const numParticleOnTile = pti.numParticles();
      auto const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

      // first we concatenate the AoS into contiguous arrays
      {
        for( auto idx=0; idx<PIt::ContainerType::NStructReal; idx++ )
          {
            if( write_real_comp[idx] ) {
              const auto [record_name, component_name] = helper::name2openPMD(real_comp_names[idx]);
              auto currRecord = currSpecies[record_name];
              auto currRecordComp = currRecord[component_name];

              std::shared_ptr< amrex::ParticleReal > d(
                                                       new amrex::ParticleReal[numParticleOnTile],
                                                       [](amrex::ParticleReal const *p){ delete[] p; }
                                                       );

              for( auto kk=0; kk<numParticleOnTile; kk++ ) {
                d.get()[kk] = aos[kk].rdata(idx);
              }

              currRecordComp.storeChunk(d, {offset}, {numParticleOnTile64});
            }
          }
      }
    }


    template<class PIt>
    void
    AMReX_openPMDWriter::StoreAoS_Int(PIt& pti,
                                      openPMD::ParticleSpecies& currSpecies,
                                      const amrex::Vector<int>& write_int_comp,
                                      const amrex::Vector<std::string>& int_comp_names,
                                      unsigned long long offset) const
    {
      auto const& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile
      auto const numParticleOnTile = pti.numParticles();
      auto const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

      // first we concatenate the AoS into contiguous arrays
      {
        for( auto idx=0; idx<PIt::ContainerType::NStructInt; idx++ )   // lgtm [cpp/constant-comparison]
          {
            if( write_int_comp[idx] ) {
              const auto [record_name, component_name] = helper::name2openPMD(int_comp_names[idx]);
              auto currRecord = currSpecies[record_name];
              auto currRecordComp = currRecord[component_name];

              std::shared_ptr< int > d(
                                       new int[numParticleOnTile],
                                       [](int const *p){ delete[] p; }
                                       );

              for( auto kk=0; kk<numParticleOnTile; kk++ ) {
                d.get()[kk] = aos[kk].idata(idx);
              }

              currRecordComp.storeChunk(d, {offset}, {numParticleOnTile64});
            }
          }
      }
    }


    template<class PIt>
    void
    AMReX_openPMDWriter::StoreSoAReal(PIt& pti,
                                      openPMD::ParticleSpecies& currSpecies,
                                      const amrex::Vector<int>& write_real_comp,
                                      const amrex::Vector<std::string>& real_comp_names,
                                      unsigned long long offset) const
    {
      auto const& soa = pti.GetStructOfArrays();
      auto const numParticleOnTile = pti.numParticles();
      auto const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

      auto const real_counter = std::min(write_real_comp.size(), real_comp_names.size()) - PIt::ContainerType::NStructReal;

      for (auto idx=0; idx<real_counter; idx++)
        {
          auto ii = PIt::ContainerType::NStructReal + idx;  // jump over extra AoS names
          if (write_real_comp[ii])
            {
              //auto curr =  soa.GetRealData(idx);
              helper::getComponentRecord(currSpecies,real_comp_names[ii]).storeChunkRaw(soa.GetRealData(idx).dataPtr(),
                                                                                        {offset},
                                                                                        {numParticleOnTile64});
          }
        }
    }



    template<class PIt>
    void AMReX_openPMDWriter::StoreSoAInt(PIt& pti,
                                           openPMD::ParticleSpecies& currSpecies,
                                           const amrex::Vector<int>& write_int_comp,
                                           const amrex::Vector<std::string>& int_comp_names,
                                           unsigned long long offset) const
    {
      auto const& soa = pti.GetStructOfArrays();
      auto const numParticleOnTile = pti.numParticles();
      auto const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

      auto const int_counter = std::min(write_int_comp.size(), int_comp_names.size()) - PIt::ContainerType::NStructInt;
      for (auto idx=0; idx<int_counter; idx++)
        {
          auto ii = PIt::ContainerType::NStructInt + idx;  // jump over extra AoS names
          if (write_int_comp[ii]) {
            //auto curr =  soa.GetIntData(idx);
            helper::getComponentRecord(currSpecies,  int_comp_names[ii]).storeChunkRaw(soa.GetIntData(idx).dataPtr(),
                                                                                       {offset},
                                                                                       {numParticleOnTile64});
          }
        }
    }


    template<class PIt>
    void
    AMReX_openPMDWriter::SavePosId(PIt& pti,
                                   openPMD::ParticleSpecies& currSpecies,
                                   unsigned long long offset) const
    {
      BL_PROFILE("amrex::openpmd_api::SavePosId(default..)");
      auto const numParticleOnTile = pti.numParticles();
      auto const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

      // get position and particle ID from aos
      // note: this implementation iterates the AoS 4x...
      // if we flush late as we do now, we can also copy out the data in one go
      const auto& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile
      {
        // Save positions
        auto const positionComponents = /*helper::*/getParticlePositionComponentLabels();

        for (auto currDim = 0; currDim < AMREX_SPACEDIM; currDim++)
          {
            //amrex::AllPrint()<<aos.size()<<"  num ptls "<<numParticleOnTile<<"\n";
            std::shared_ptr<amrex::ParticleReal> curr(new amrex::ParticleReal[numParticleOnTile],
                                                      [](amrex::ParticleReal const *p) { delete[] p; }
                                                      );
            for (auto i = 0; i < numParticleOnTile; i++) {
              //amrex::AllPrint()<<"  i =  "<<i<<" of "<<numParticleOnTile<<"\n";
              //amrex::AllPrint()<<"  i =  "<<i<<" "<<aos[i]<<"\n";
              curr.get()[i] = aos[i].pos(currDim);
            }
            std::string const& positionComponent = positionComponents[currDim];
            currSpecies["position"][positionComponent].storeChunk(curr, {offset + GetGrandOffset()}, {numParticleOnTile64});
          }
      }
      SaveSpecieId(pti, currSpecies, offset);
    }

    template<class PIt>
    void
    AMReX_openPMDWriter::SaveSpecieId(PIt& pti,
                                      openPMD::ParticleSpecies& currSpecies,
                                      unsigned long long offset) const
    {
      BL_PROFILE("amrex::openpmd_api::SaveSpecieId");
      auto const numParticleOnTile = pti.numParticles();
      auto const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

      // get position and particle ID from aos
      // note: this implementation iterates the AoS 4x...
      // if we flush late as we do now, we can also copy out the data in one go
      const auto& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile
      {
        // save particle ID after converting it to a globally unique ID
        std::shared_ptr< uint64_t > ids(new uint64_t[numParticleOnTile],
                                        [](uint64_t const *p){ delete[] p; }
                                        );
        for (auto i=0; i<numParticleOnTile; i++)
          {
            ids.get()[i] = helper::localIDtoGlobal(aos[i].id(), aos[i].cpu());
          }
        auto const * const scalar = openPMD::RecordComponent::SCALAR;
        currSpecies["id"][scalar].storeChunk(ids, {offset + GetGrandOffset()}, {numParticleOnTile64});
      }
    }

    template<class PIt>
    void
    AMReX_openPMDWriter::SaveRealProperty (PIt& pti,
                                           openPMD::ParticleSpecies& currSpecies,
                                           unsigned long long offset,
                                           const amrex::Vector<int>& write_real_comp,
                                           const amrex::Vector<std::string>& real_comp_names,
                                           const amrex::Vector<int>& write_int_comp,
                                           const amrex::Vector<std::string>& int_comp_names) const
    {
      // note: WarpX does not yet use extra AoS Real attributes
      StoreAoS_Real(pti, currSpecies, write_real_comp, real_comp_names, offset);
      StoreAoS_Int(pti, currSpecies, write_int_comp, int_comp_names, offset);

      //auto const& soa = pti.GetStructOfArrays();
      StoreSoAReal(pti, currSpecies, write_real_comp, real_comp_names, offset);
      StoreSoAInt(pti, currSpecies, write_int_comp, int_comp_names, offset);
    }




    ////////////////////////////////////
    //
    // AMReX_openPMDWriter::DumpParticles()
    //
    ////////////////////////////////////

    template<class PC>
    void AMReX_openPMDWriter::DumpParticles(PC const& pc,
                                            const std::string& name,
                                            //int iteration,
                                            const amrex::Vector<int>& write_real_comp,
                                            const amrex::Vector<int>& write_int_comp,
                                            const amrex::Vector<std::string>& real_comp_names,
                                            const amrex::Vector<std::string>& int_comp_names
                                            ) const
    {
      DumpParticles(pc, name, write_real_comp, write_int_comp, real_comp_names, int_comp_names,
                    [=] (auto& /*pc*/, openPMD::ParticleSpecies& /*currSpecies*/, unsigned long long /*size*/)
                    {
                    }, // no extra pc level info to save
                    [=] (auto& pti, openPMD::ParticleSpecies& currSpecies, unsigned long long offset)
                    {
                      SavePosId(pti, currSpecies, offset); // default way of save pos & id
                    });
    }

    template<class PC, typename FUNC_pc, typename FUNC_pit>
    void AMReX_openPMDWriter::DumpParticles(PC const& pc,
                                            const std::string& name,
                                            //int iteration,
                                            const amrex::Vector<int>& write_real_comp,
                                            const amrex::Vector<int>& write_int_comp,
                                            const amrex::Vector<std::string>& real_comp_names,
                                            const amrex::Vector<std::string>& int_comp_names,
                                            FUNC_pc&& pc_func,
                                            FUNC_pit&& posId_func
                                            ) const
    {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Series != nullptr, "openPMD series must be initialized");

      // TODO have to Count before this function due to const restriction
      //     Count is not a const function
      //pc.CountParticles();

      openPMD::Iteration currIteration = GetIteration( m_CurrentStep );
      openPMD::ParticleSpecies currSpecies = currIteration.particles[name];

      //
      // define positions & offsets
      //
      /*AllocatePtlProperties(currSpecies, write_real_comp, real_comp_names, write_int_comp, int_comp_names,
                                pc.m_PtlCounter.GetTotalNumParticles()); */

      if ( AllocatePtlProperties(currSpecies, write_real_comp, real_comp_names, write_int_comp, int_comp_names,
                                 pc.m_PtlCounter.GetTotalNumParticles()) )
        {
          SetupRealProperties(pc, currSpecies, write_real_comp, real_comp_names, write_int_comp, int_comp_names,
                              pc.m_PtlCounter.GetTotalNumParticles());
          SetupConstant(currSpecies, pc.m_PtlCounter.GetTotalNumParticles());
        }

      pc_func(pc, currSpecies, pc.m_PtlCounter.GetTotalNumParticles());
      // open files from all processors, in case some will not contribute below
      m_Series->flush();

      bool emptyPC = true;
      for (auto currentLevel = 0; currentLevel <= pc.finestLevel(); currentLevel++)
        {
          auto offset = static_cast<uint64_t>( pc.m_PtlCounter.m_ParticleOffsetAtRank[currentLevel] );
          using ParIter = typename PC::ParConstIterType;

          //for (WarpXParIter pti(*pc, currentLevel); pti.isValid(); ++pti) {
          for (ParIter pti(pc, currentLevel); pti.isValid(); ++pti)
            {
              auto const numParticleOnTile = pti.numParticles();
              auto const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

              // Do not call storeChunk() with zero-sized particle tiles:
              //   https://github.com/openPMD/openPMD-api/issues/1147
              //   https://github.com/ECP-WarpX/WarpX/pull/1898#discussion_r745008290
              if (numParticleOnTile == 0) { continue; }

              emptyPC = false;

              //SavePosId(pti, currSpecies, offset);
              posId_func(pti, currSpecies, offset);
              //  save "extra" particle properties in AoS and SoA
              SaveRealProperty(pti,
                               currSpecies,
                               offset + GetGrandOffset(),
                               write_real_comp,real_comp_names,
                               write_int_comp, int_comp_names);

              offset += numParticleOnTile64;
          }
        }

      if (emptyPC && m_Series->backend() == "ADIOS2") {
        // TODO adiosWorkaround with empty rank reading error() for __BTD__
      }
      m_Series->flush();
    }


    template<class PC>
    void AMReX_openPMDWriter::SetupRealProperties (PC const& /*pc*/,
                                                   openPMD::ParticleSpecies& currSpecies,
                                                   const amrex::Vector<int>& write_real_comp,
                                                   const amrex::Vector<std::string>& real_comp_names,
                                                   const amrex::Vector<int>& write_int_comp,
                                                   const amrex::Vector<std::string>& int_comp_names,
                                                   const unsigned long long /*np*/) const
    {

      //AllocateRealProperties(currSpecies, write_real_comp, real_comp_names, write_int_comp, int_comp_names, np);

      // set attributes for some properties
      std::set< std::string > addedRecords; // add meta-data per record only once

      auto const real_counter = std::min(write_real_comp.size(), real_comp_names.size()) - PC::NStructReal;
      //for (auto idx=0; idx<pc.NumRealComps(); idx++)
      for (auto idx=0; idx<real_counter; idx++)
        {
          auto ii = PC::NStructReal + idx; // jump over extra AoS names
          if (write_real_comp[ii])
            {
              // handle scalar and non-scalar records by name
              const auto [record_name, component_name] = helper::name2openPMD(real_comp_names[ii]);
              auto currRecord = currSpecies[record_name];

              // meta data for ED-PIC extension
              [[maybe_unused]] const auto [_, newRecord] = addedRecords.insert(record_name);
              if( newRecord )
                {
                  currRecord.setUnitDimension( helper::getUnitDimension(record_name) );
                  if( record_name == "weighting" ) {
                    currRecord.setAttribute( "macroWeighted", 1U );
                  } else {
                    currRecord.setAttribute( "macroWeighted", 0U );
                  }
                  if( record_name == "momentum" || record_name == "weighting" ) {
                    currRecord.setAttribute( "weightingPower", 1.0 );
                  } else {
                    currRecord.setAttribute( "weightingPower", 0.0 );
                  }
                }
            }
        }

      auto const int_counter = std::min(write_int_comp.size(), int_comp_names.size()) - PC::NStructInt;
      for (auto idx=0; idx<int_counter; idx++)
        {
          auto ii = PC::NStructInt + idx; // jump over extra AoS names
          if (write_int_comp[ii])
            {
              // handle scalar and non-scalar records by name
              const auto [record_name, component_name] = helper::name2openPMD(int_comp_names[ii]);
              auto currRecord = currSpecies[record_name];

              // meta data for ED-PIC extension
              [[maybe_unused]] const auto [_, newRecord] = addedRecords.insert(record_name);
              if( newRecord )
                {
                  currRecord.setUnitDimension( helper::getUnitDimension(record_name) );
                  currRecord.setAttribute( "macroWeighted", 0U );
                  if( record_name == "momentum" || record_name == "weighting" ) {
                    currRecord.setAttribute( "weightingPower", 1.0 );
                  } else {
                    currRecord.setAttribute( "weightingPower", 0.0 );
                  }
                }
            }
        }
    }

} // namespace amrex::openpmd_api

#endif // AMREX_PLOTFILE_OPENPMD_API_H
