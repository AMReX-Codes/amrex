#ifndef _WritePlotFile_H_
#define _WritePlotFile_H_

#include <string>

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_AmrData.H>

// Write a MultiFab to a plotfile named "name" using the HyperCLawV1.1 parallel
// plotfile data format
void
writePlotFile (const char*                 name,
	       const amrex::MultiFab&            mf,
	       const amrex::Geometry&            geom,
	       const amrex::IntVect&             refRatio,
	       amrex::Real                       bgVal,
         const amrex::Vector<std::string>& names);

// Write an array of multifabs pulling all the data from the AmrData
void WritePlotFile(const amrex::Vector<amrex::MultiFab*>& mfa,
                   const amrex::Vector<amrex::Box>&       probDomain,
		   amrex::AmrData&                        amrdToMimic,
		   const std::string&                     oFile,
		   bool                                   verbose = false,
                   const amrex::Vector<std::string>&      varNames = amrex::Vector<std::string>());

// Write an array of multifabs pulling all the data from the AmrData
void WritePlotFile(const amrex::Vector<amrex::MultiFab*>& mfa,
		   amrex::AmrData&                        amrdToMimic,
		   const std::string&                     oFile,
		   bool                                   verbose = false,
                   const amrex::Vector<std::string>&      varNames = amrex::Vector<std::string>());

#endif

