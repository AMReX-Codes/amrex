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
	    
	    aString FullPath = file;
	    if (!FullPath.isNull() && FullPath[FullPath.length()-1] != '/')
		FullPath += '/';
	    
	    if (ParallelDescriptor::IOProcessor())
		if (!Utility::CreateDirectory(FullPath, 0755))
		    Utility::CreateDirectoryFailed(FullPath);
	    
	    static const aString MultiFabBaseName("MultiFab");
	    FullPath += MultiFabBaseName;
	    
	    VisMF::Write(*mf,FullPath,VisMF::OneFilePerCPU);
	    
	    return true;
	}
    
    // read multifab from file
    bool readMF (MultiFab&   mf,
		 const char* file)
	{
	    
	    aString FullPath = file;
	    if (!FullPath.isNull() && FullPath[FullPath.length()-1] != '/')
		FullPath += '/';
	    
	    static const aString MultiFabBaseName("MultiFab");
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
