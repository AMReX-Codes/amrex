# Microsoft Developer Studio Project File - Name="hglib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=hglib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "hglib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "hglib.mak" CFG="hglib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "hglib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "hglib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=xicl.exe
F90=df.exe

!IF  "$(CFG)" == "hglib - Win32 Release"

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
# ADD F90 /include:"Release/" /compile_only /nologo /iface:cref /warn:nofileopt
# SUBTRACT F90 /browser
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "." /I "include\3d.v7" /I "..\pBoxLib_2" /D "NDEBUG" /D "_WINDOWS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /D "HG_CROSS_STENCIL" /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "hglib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "hglib___"
# PROP BASE Intermediate_Dir "hglib___"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "hglib___"
# PROP Intermediate_Dir "hglib___"
# PROP Target_Dir ""
# ADD BASE F90 /include:"hglib___/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD F90 /browser /include:"hglib___/" /compile_only /nologo /debug:full /optimize:0 /iface:cref /dbglibs /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /I "." /I "include\3d.v7" /I "..\pBoxLib_2" /D "_DEBUG" /D "_WINDOWS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /D "HG_CROSS_STENCIL" /FR /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "hglib - Win32 Release"
# Name "hglib - Win32 Debug"
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

SOURCE=.\cache.H
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

SOURCE=.\interface.H
# End Source File
# Begin Source File

SOURCE=.\interpolator.H
# End Source File
# Begin Source File

SOURCE=.\restrictor.H
# End Source File
# End Group
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\amr_multi.cpp
# End Source File
# Begin Source File

SOURCE=.\boundary.cpp
# End Source File
# Begin Source File

SOURCE=.\cache.cpp
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

SOURCE=.\interface.cpp
# End Source File
# Begin Source File

SOURCE=.\interpolator.cpp
# End Source File
# Begin Source File

SOURCE=.\restrictor.cpp
# End Source File
# End Group
# Begin Group "Fortran"

# PROP Default_Filter "F"
# Begin Source File

SOURCE=.\amr_real3d.F

!IF  "$(CFG)" == "hglib - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\amr_real3d.F
InputName=amr_real3d

"$(InputName).For" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /ansi /nologo /S. /Sinclude\3d.v7 /S..\pBoxLib_2 /DBL_LANG_FORT\
                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH $(InputName).F | perl\
                 ..\scripts\strip72 -c > $(InputName).for

# End Custom Build

!ELSEIF  "$(CFG)" == "hglib - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_avg3d.F

!IF  "$(CFG)" == "hglib - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\hg_avg3d.F
InputName=hg_avg3d

"$(InputName).For" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /ansi /nologo /S. /Sinclude\3d.v7 /S..\pBoxLib_2 /DBL_LANG_FORT\
                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH $(InputName).F | perl\
                 ..\scripts\strip72 -c > $(InputName).for

# End Custom Build

!ELSEIF  "$(CFG)" == "hglib - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_multi3d.F

!IF  "$(CFG)" == "hglib - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\hg_multi3d.F
InputName=hg_multi3d

"$(InputName).For" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /ansi /nologo /S. /Sinclude\3d.v7 /S..\pBoxLib_2 /DBL_LANG_FORT\
                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH $(InputName).F | perl\
                 ..\scripts\strip72 -c > $(InputName).for

# End Custom Build

!ELSEIF  "$(CFG)" == "hglib - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_multi3d_terrain.F
# End Source File
# Begin Source File

SOURCE=.\hg_proj3d.F

!IF  "$(CFG)" == "hglib - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\hg_proj3d.F
InputName=hg_proj3d

"$(InputName).For" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /ansi /nologo /S. /Sinclude\3d.v7 /S..\pBoxLib_2 /DBL_LANG_FORT\
                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH $(InputName).F | perl\
                 ..\scripts\strip72 -c > $(InputName).for

# End Custom Build

!ELSEIF  "$(CFG)" == "hglib - Win32 Debug"

!ENDIF 

# End Source File
# End Group
# End Target
# End Project
