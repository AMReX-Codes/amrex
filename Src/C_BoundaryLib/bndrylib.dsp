# Microsoft Developer Studio Project File - Name="bndrylib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=bndrylib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "bndrylib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "bndrylib.mak" CFG="bndrylib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "bndrylib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "bndrylib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=xicl.exe
F90=df.exe

!IF  "$(CFG)" == "bndrylib - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "." /I "..\pBoxLib_2" /D "NDEBUG" /D "_WINDOWS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "bndrylib - Win32 Debug"

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
# ADD F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /iface:cref /warn:nofileopt
# SUBTRACT F90 /dbglibs
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /I "." /I "..\pBoxLib_2" /D "_DEBUG" /D "_WINDOWS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /FR /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "bndrylib - Win32 Release"
# Name "bndrylib - Win32 Debug"
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\BndryRegister.cpp
# End Source File
# Begin Source File

SOURCE=.\CoordSys.cpp
# End Source File
# Begin Source File

SOURCE=.\FabSet.cpp
# End Source File
# Begin Source File

SOURCE=.\Geometry.cpp
# End Source File
# Begin Source File

SOURCE=.\RealBox.cpp
# End Source File
# End Group
# Begin Group "C++ Headers"

# PROP Default_Filter "h"
# Begin Source File

SOURCE=.\BndryRegister.H
# End Source File
# Begin Source File

SOURCE=.\CoordSys.H
# End Source File
# Begin Source File

SOURCE=.\Coordsys_f.h
# End Source File
# Begin Source File

SOURCE=.\FabSet.H
# End Source File
# Begin Source File

SOURCE=.\Geometry.H
# End Source File
# Begin Source File

SOURCE=.\RealBox.H
# End Source File
# End Group
# Begin Group "Fortran"

# PROP Default_Filter ""
# Begin Group "3D"

# PROP Default_Filter ""
# Begin Group "Temps"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\COORDSYS_3D.FOR
# End Source File
# End Group
# Begin Source File

SOURCE=.\COORDSYS_3D.F

!IF  "$(CFG)" == "bndrylib - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\COORDSYS_3D.F
InputName=COORDSYS_3D

"$(InputName).for" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /ansi /nologo /S. /S..\amrlib /S..\bndrylib /S..\pBoxLib_2 /DBL_LANG_FORT\
          /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH $(InputName).F | perl\
          ..\scripts\strip72 -c > $(InputName).for

# End Custom Build

!ELSEIF  "$(CFG)" == "bndrylib - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\COORDSYS_3D.F
InputName=COORDSYS_3D

"$(InputName).for" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /ansi /nologo /S. /S..\amrlib /S..\bndrylib /S..\pBoxLib_2 /DBL_LANG_FORT\
          /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH $(InputName).F | perl\
          ..\scripts\strip72 -c > $(InputName).for

# End Custom Build

!ENDIF 

# End Source File
# End Group
# End Group
# End Target
# End Project
