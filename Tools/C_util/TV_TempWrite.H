#ifndef _TV_TempWrite_H_
#define _TV_TempWrite_H_

#include "iostream.h"
#include "MultiFab.H"
#include "VisMF.H"

extern "C" {
    
    // dump multifab to file
    bool writeMF (const MultiFab* mf,
		  const char*     file)
	{
	    
	    std::string FullPath = file;
	    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
		FullPath += '/';
	    
	    if (ParallelDescriptor::IOProcessor())
		if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
		    BoxLib::CreateDirectoryFailed(FullPath);
	    
	    static const std::string MultiFabBaseName("MultiFab");
	    FullPath += MultiFabBaseName;
	    
	    VisMF::Write(*mf,FullPath,VisMF::OneFilePerCPU);
	    
	    return true;
	}
    
    // read multifab from file
    bool readMF (MultiFab&   mf,
		 const char* file)
	{
	    
	    std::string FullPath = file;
	    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
		FullPath += '/';
	    
	    static const std::string MultiFabBaseName("MultiFab");
	    FullPath += MultiFabBaseName;

	    VisMF::Read(mf,FullPath);
	    
	    return true;
	}
}

void NeverCalledEXT()
{
    writeMF( NULL, "DORK" );
}





#endif
