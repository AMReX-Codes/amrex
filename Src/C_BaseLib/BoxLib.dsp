# Microsoft Developer Studio Project File - Name="BoxLib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=BoxLib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "BoxLib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "BoxLib.mak" CFG="BoxLib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "BoxLib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "BoxLib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "BoxLib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /iface:cref /libs:dll /nologo /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_LIB /D" /YX /FD _MBCS" /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "." /I "..\BoxLib" /D "NDEBUG" /D "_LIB" /D "BL_LANG_CC" /D "WIN32" /D "_MBCS" /D for="if(0);else for" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D "BL_USE_FLOAT" /YX /FD /c
# ADD BASE RSC /l 0x409
# ADD RSC /l 0x409
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "BoxLib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /debug:full /nologo /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /dbglibs /debug:full /iface:cref /libs:dll /nologo /threads /warn:argument_checking /warn:declarations /warn:nofileopt
# SUBTRACT F90 /traceback
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_LIB" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "." /I "..\BoxLib" /D "_LIB" /D "BL_LANG_CC" /D "_DEBUG" /D "BL_PROFILING" /D "BL_USE_FLOAT" /D "WIN32" /D "_MBCS" /D for="if(0);else for" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /FR /YX /FD /c
# ADD BASE RSC /l 0x409
# ADD RSC /l 0x409
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "BoxLib - Win32 Release"
# Name "BoxLib - Win32 Debug"
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\Arena.cpp
# End Source File
# Begin Source File

SOURCE=.\BArena.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BLMpi.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\BLProfiler.cpp
# End Source File
# Begin Source File

SOURCE=.\BLThread.cpp
# End Source File
# Begin Source File

SOURCE=.\BLWorkQueue.cpp
# End Source File
# Begin Source File

SOURCE=.\Box.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BoxArray.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BoxDomain.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BoxLib.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BoxList.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\CArena.cpp
# End Source File
# Begin Source File

SOURCE=.\DistributionMapping.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\FabArray.cpp
# End Source File
# Begin Source File

SOURCE=.\FabConv.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\FArrayBox.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\FPC.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\IndexType.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\IntVect.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\MultiFab.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\Orientation.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\ParallelDescriptor.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\ParmParse.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\UseCount.cpp
# End Source File
# Begin Source File

SOURCE=.\Utility.cpp
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\VisMF.cpp
# SUBTRACT CPP /I "."
# End Source File
# End Group
# Begin Group "C++ Headers"

# PROP Default_Filter "h"
# Begin Source File

SOURCE=.\Arena.H
# End Source File
# Begin Source File

SOURCE=.\Array.H
# End Source File
# Begin Source File

SOURCE=.\ArrayLim.H
# End Source File
# Begin Source File

SOURCE=.\BArena.H
# End Source File
# Begin Source File

SOURCE=.\BaseFab.H
# End Source File
# Begin Source File

SOURCE=.\BLassert.H
# End Source File
# Begin Source File

SOURCE=.\Blversion.h
# End Source File
# Begin Source File

SOURCE=.\Box.H
# End Source File
# Begin Source File

SOURCE=.\BoxArray.H
# End Source File
# Begin Source File

SOURCE=.\BoxDomain.H
# End Source File
# Begin Source File

SOURCE=.\BoxLib.H
# End Source File
# Begin Source File

SOURCE=.\BoxList.H
# End Source File
# Begin Source File

SOURCE=.\CArena.H
# End Source File
# Begin Source File

SOURCE=".\ccse-mpi.H"
# End Source File
# Begin Source File

SOURCE=.\Constants.h
# End Source File
# Begin Source File

SOURCE=.\DistributionMapping.H
# End Source File
# Begin Source File

SOURCE=.\FabArray.H
# End Source File
# Begin Source File

SOURCE=.\FabConv.H
# End Source File
# Begin Source File

SOURCE=.\FArrayBox.H
# End Source File
# Begin Source File

SOURCE=.\Fpc.h
# End Source File
# Begin Source File

SOURCE=.\IndexType.H
# End Source File
# Begin Source File

SOURCE=.\IntVect.H
# End Source File
# Begin Source File

SOURCE=.\Looping.H
# End Source File
# Begin Source File

SOURCE=.\MultiFab.H
# End Source File
# Begin Source File

SOURCE=.\Orientation.H
# End Source File
# Begin Source File

SOURCE=.\ParallelDescriptor.H
# End Source File
# Begin Source File

SOURCE=.\ParmParse.H
# End Source File
# Begin Source File

SOURCE=.\PArray.H
# End Source File
# Begin Source File

SOURCE=.\Pointers.H
# End Source File
# Begin Source File

SOURCE=.\Profiler.H
# End Source File
# Begin Source File

SOURCE=.\Real.h
# End Source File
# Begin Source File

SOURCE=.\Space.h
# End Source File
# Begin Source File

SOURCE=.\Space_f.h
# End Source File
# Begin Source File

SOURCE=.\Thread.H
# End Source File
# Begin Source File

SOURCE=.\Tuple.H
# End Source File
# Begin Source File

SOURCE=.\UseCount.H
# End Source File
# Begin Source File

SOURCE=.\Utility.H
# End Source File
# Begin Source File

SOURCE=.\VisMF.H
# End Source File
# Begin Source File

SOURCE=.\winstd.H
# End Source File
# Begin Source File

SOURCE=.\WorkQueue.H
# End Source File
# End Group
# End Target
# End Project
