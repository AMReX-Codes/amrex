# Microsoft Developer Studio Project File - Name="amrlib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=amrlib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "amrlib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "amrlib.mak" CFG="amrlib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "amrlib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "amrlib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=xicl.exe
F90=df.exe

!IF  "$(CFG)" == "amrlib - Win32 Release"

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
# ADD F90 /include:"Release/" /compile_only /nologo /warn:declarations /iface:cref /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "." /I "..\bndrylib" /I "..\pBoxLib_2" /D "NDEBUG" /D "_WINDOWS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "amrlib - Win32 Debug"

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
# ADD F90 /include:"Debug/" /compile_only /nologo /warn:declarations /debug:full /optimize:0 /iface:cref /warn:nofileopt
# SUBTRACT F90 /dbglibs
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /I "." /I "..\bndrylib" /I "..\pBoxLib_2" /D "_DEBUG" /D "_WINDOWS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /FR /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "amrlib - Win32 Release"
# Name "amrlib - Win32 Debug"
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\Amr.cpp
# End Source File
# Begin Source File

SOURCE=.\AmrLevel.cpp
# End Source File
# Begin Source File

SOURCE=.\BCRec.cpp
# End Source File
# Begin Source File

SOURCE=.\Cluster.cpp
# End Source File
# Begin Source File

SOURCE=.\Derive.cpp
# End Source File
# Begin Source File

SOURCE=.\ErrorList.cpp
# End Source File
# Begin Source File

SOURCE=.\FluxRegister.cpp
# End Source File
# Begin Source File

SOURCE=.\Interpolater.cpp
# End Source File
# Begin Source File

SOURCE=.\StateData.cpp
# End Source File
# Begin Source File

SOURCE=.\StateDescriptor.cpp
# End Source File
# Begin Source File

SOURCE=.\TagBox.cpp
# End Source File
# End Group
# Begin Group "C++ Headers"

# PROP Default_Filter "h"
# Begin Source File

SOURCE=.\Amr.H
# End Source File
# Begin Source File

SOURCE=.\Amr_auxil.H
# End Source File
# Begin Source File

SOURCE=.\AmrLevel.H
# End Source File
# Begin Source File

SOURCE=.\Bc_types.h
# End Source File
# Begin Source File

SOURCE=.\BCRec.H
# End Source File
# Begin Source File

SOURCE=.\Cluster.H
# End Source File
# Begin Source File

SOURCE=.\Derive.H
# End Source File
# Begin Source File

SOURCE=.\Dims.h
# End Source File
# Begin Source File

SOURCE=.\ErrorList.H
# End Source File
# Begin Source File

SOURCE=.\Fluxreg_f.h
# End Source File
# Begin Source File

SOURCE=.\FluxRegister.H
# End Source File
# Begin Source File

SOURCE=.\Interp_f.h
# End Source File
# Begin Source File

SOURCE=.\Interpolater.H
# End Source File
# Begin Source File

SOURCE=.\LevelBld.H
# End Source File
# Begin Source File

SOURCE=.\Makeslice_f.h
# End Source File
# Begin Source File

SOURCE=.\Prob_amr_f.h
# End Source File
# Begin Source File

SOURCE=.\StateData.H
# End Source File
# Begin Source File

SOURCE=.\StateDescriptor.H
# End Source File
# Begin Source File

SOURCE=.\TagBox.H
# End Source File
# End Group
# Begin Group "Fortran"

# PROP Default_Filter ""
# Begin Group "2D"

# PROP Default_Filter ""
# End Group
# Begin Group "3D"

# PROP Default_Filter ""
# Begin Group "Temps"

# PROP Default_Filter "for"
# Begin Source File

SOURCE=.\Filcc_3d.for
# End Source File
# Begin Source File

SOURCE=.\FLUXREG_3D.for
# End Source File
# Begin Source File

SOURCE=.\INTERP_3D.for
# End Source File
# Begin Source File

SOURCE=.\MAKESLICE_3D.for
# End Source File
# End Group
# Begin Source File

SOURCE=.\FILCC_3D.F
# End Source File
# Begin Source File

SOURCE=.\FLUXREG_3D.F
# End Source File
# Begin Source File

SOURCE=.\INTERP_3D.F
# End Source File
# Begin Source File

SOURCE=.\MAKESLICE_3D.F
# End Source File
# End Group
# End Group
# End Target
# End Project
