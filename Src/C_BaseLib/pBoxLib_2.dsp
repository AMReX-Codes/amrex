# Microsoft Developer Studio Project File - Name="pBoxLib_2" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=pBoxLib_2 - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "pBoxLib_2.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "pBoxLib_2.mak" CFG="pBoxLib_2 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "pBoxLib_2 - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "pBoxLib_2 - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe

!IF  "$(CFG)" == "pBoxLib_2 - Win32 Release"

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
# ADD BASE F90 /include:"Release/" /compile_only /nologo /warn:nofileopt
# ADD F90 /include:"Release/" /compile_only /nologo /stand:f90 /iface:cref /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "." /D "NDEBUG" /D "_WINDOWS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D for="if(0);else for" /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "pBoxLib_2 - Win32 Debug"

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
# ADD BASE F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD F90 /browser /include:"Debug/" /compile_only /nologo /warn:declarations /debug:full /optimize:0 /check:bounds /warn:argument_checking /stand:f90 /iface:cref /threads /dbglibs /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /MTd /W3 /Gm /GX /Zi /Od /I "." /D "_WINDOWS" /D "_DEBUG" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D for="if(0);else for" /FR /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "pBoxLib_2 - Win32 Release"
# Name "pBoxLib_2 - Win32 Debug"
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\aString.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BArena.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\Box.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BoxArray.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BoxAssoc.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BoxDomain.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BoxLib.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\BoxList.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\CArena.cpp
# End Source File
# Begin Source File

SOURCE=.\DistributionMapping.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\FabConv.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\FArrayBox.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\FPC.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\IndexType.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\IntVect.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\MultiFab.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\Orientation.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\ParallelDescriptor.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\ParmParse.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\RunStats.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\tDir.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\Tracer.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\tVisMF.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\Utility.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# Begin Source File

SOURCE=.\VisMF.cpp
# ADD CPP /I "..\pBoxLib_2"
# SUBTRACT CPP /I "."
# End Source File
# End Group
# Begin Group "C++ Headers"

# PROP Default_Filter "h"
# Begin Source File

SOURCE=.\Arena.H
# End Source File
# Begin Source File

SOURCE=.\ArithFab.H
# End Source File
# Begin Source File

SOURCE=.\Array.H
# End Source File
# Begin Source File

SOURCE=.\ArrayLim.H
# End Source File
# Begin Source File

SOURCE=.\Assert.H
# End Source File
# Begin Source File

SOURCE=.\aString.H
# End Source File
# Begin Source File

SOURCE=.\BArena.H
# End Source File
# Begin Source File

SOURCE=.\BaseFab.H
# End Source File
# Begin Source File

SOURCE=.\Blversion.h
# End Source File
# Begin Source File

SOURCE=.\Boolean.H
# End Source File
# Begin Source File

SOURCE=.\Box.H
# End Source File
# Begin Source File

SOURCE=.\BoxArray.H
# End Source File
# Begin Source File

SOURCE=.\BoxAssoc.H
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

SOURCE=.\List.H
# End Source File
# Begin Source File

SOURCE=.\Looping.H
# End Source File
# Begin Source File

SOURCE=.\Misc.H
# End Source File
# Begin Source File

SOURCE=.\MultiFab.H
# End Source File
# Begin Source File

SOURCE=.\NormedFab.H
# End Source File
# Begin Source File

SOURCE=.\OrderedFab.H
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

SOURCE=.\Real.h
# End Source File
# Begin Source File

SOURCE=.\RunStats.H
# End Source File
# Begin Source File

SOURCE=.\Space.h
# End Source File
# Begin Source File

SOURCE=.\Space_f.h
# End Source File
# Begin Source File

SOURCE=.\Tracer.H
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
# End Group
# End Target
# End Project
