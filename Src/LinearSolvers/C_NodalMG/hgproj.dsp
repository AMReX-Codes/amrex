# Microsoft Developer Studio Project File - Name="hgproj" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=hgproj - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "hgproj.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "hgproj.mak" CFG="hgproj - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "hgproj - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "hgproj - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "hgproj - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /include:"Release/" /compile_only /nologo /warn:nofileopt
# ADD F90 /include:"Release/" /compile_only /nologo /libs:dll /stand:f90 /iface:cref /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "..\pBoxLib_2" /I "." /D "_CONSOLE" /D "_MBCS" /D "NDEBUG" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D for="if(0);else for" /D "BL_USE_MPI" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib cvwmpi.lib /nologo /subsystem:console /profile /debug /machine:I386 /include:"__matherr"
# SUBTRACT LINK32 /verbose

!ELSEIF  "$(CFG)" == "hgproj - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD F90 /browser /include:"Debug/" /compile_only /nologo /warn:declarations /libs:dll /debug:full /optimize:0 /check:bounds /warn:argument_checking /fpe:0 /stand:f90 /iface:cref /threads /dbglibs /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /Zi /Od /I "..\pBoxLib_2" /I "." /D "_CONSOLE" /D "_MBCS" /D "_DEBUG" /D "HG_DEBUG" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D for="if(0);else for" /D "BL_USE_MPI" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib cvwmpi.lib /nologo /subsystem:console /debug /machine:I386 /include:"__matherr"
# SUBTRACT LINK32 /verbose /profile

!ENDIF 

# Begin Target

# Name "hgproj - Win32 Release"
# Name "hgproj - Win32 Debug"
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\amr_multi.cpp
# End Source File
# Begin Source File

SOURCE=.\boundary.cpp
# End Source File
# Begin Source File

SOURCE=.\fill_patch.cpp
# End Source File
# Begin Source File

SOURCE=.\hg_multi1.cpp
# End Source File
# Begin Source File

SOURCE=.\hg_multi2.cpp
# End Source File
# Begin Source File

SOURCE=.\hg_multi3.cpp
# End Source File
# Begin Source File

SOURCE=.\hg_projector.cpp
# End Source File
# Begin Source File

SOURCE=.\hgparallel.cpp
# End Source File
# Begin Source File

SOURCE=.\interface.cpp
# End Source File
# Begin Source File

SOURCE=.\interpolator.cpp
# End Source File
# Begin Source File

SOURCE=.\restrictor.cpp
# End Source File
# End Group
# Begin Group "C++ Headers"

# PROP Default_Filter "h"
# Begin Source File

SOURCE=.\amr_defs.H
# End Source File
# Begin Source File

SOURCE=.\amr_multi.H
# End Source File
# Begin Source File

SOURCE=.\boundary.H
# End Source File
# Begin Source File

SOURCE=.\fill_patch.H
# End Source File
# Begin Source File

SOURCE=.\hg_multi.H
# End Source File
# Begin Source File

SOURCE=.\hg_projector.H
# End Source File
# Begin Source File

SOURCE=.\hgparallel.h
# End Source File
# Begin Source File

SOURCE=.\interface.H
# End Source File
# Begin Source File

SOURCE=.\interpolator.H
# End Source File
# Begin Source File

SOURCE=.\RegType.H
# End Source File
# Begin Source File

SOURCE=.\restrictor.H
# End Source File
# End Group
# Begin Group "Fortran"

# PROP Default_Filter "F"
# Begin Source File

SOURCE=.\amr_real2d.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\amr_real3d.F

!IF  "$(CFG)" == "hgproj - Win32 Release"

!ELSEIF  "$(CFG)" == "hgproj - Win32 Debug"

# ADD F90 /warn:declarations /check:bounds

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_avg2d.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\hg_avg3d.F

!IF  "$(CFG)" == "hgproj - Win32 Release"

!ELSEIF  "$(CFG)" == "hgproj - Win32 Debug"

# ADD F90 /warn:declarations /check:bounds

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_multi2d.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\hg_multi3d.F

!IF  "$(CFG)" == "hgproj - Win32 Release"

!ELSEIF  "$(CFG)" == "hgproj - Win32 Debug"

# ADD F90 /warn:declarations /check:bounds

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_multi3d_terrain.F

!IF  "$(CFG)" == "hgproj - Win32 Release"

!ELSEIF  "$(CFG)" == "hgproj - Win32 Debug"

# ADD F90 /check:bounds

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_proj2d.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\hg_proj3d.F

!IF  "$(CFG)" == "hgproj - Win32 Release"

!ELSEIF  "$(CFG)" == "hgproj - Win32 Debug"

# ADD F90 /warn:declarations /check:bounds

!ENDIF 

# End Source File
# End Group
# Begin Source File

SOURCE=.\hgproj.pg
# End Source File
# Begin Source File

SOURCE=.\proj.cpp
# End Source File
# End Target
# End Project
