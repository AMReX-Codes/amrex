# Microsoft Developer Studio Project File - Name="mglib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=mglib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "mglib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "mglib.mak" CFG="mglib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "mglib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "mglib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=xicl.exe
F90=df.exe

!IF  "$(CFG)" == "mglib - Win32 Release"

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
# ADD CPP /nologo /W3 /GX /O2 /I "." /I "..\pBoxLib_2" /I "..\bndrylib" /I "..\amrlib" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "mglib - Win32 Debug"

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
# ADD F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /iface:cref /dbglibs /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /I "." /I "..\pBoxLib_2" /I "..\bndrylib" /I "..\amrlib" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "mglib - Win32 Release"
# Name "mglib - Win32 Debug"
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\ABecLaplacian.cpp
# End Source File
# Begin Source File

SOURCE=.\BndryData.cpp
# End Source File
# Begin Source File

SOURCE=.\CGSolver.cpp
# End Source File
# Begin Source File

SOURCE=.\InterpBndryData.cpp
# End Source File
# Begin Source File

SOURCE=.\Laplacian.cpp
# End Source File
# Begin Source File

SOURCE=.\LinOp.cpp
# End Source File
# Begin Source File

SOURCE=.\Mask.cpp
# End Source File
# Begin Source File

SOURCE=.\MultiGrid.cpp
# End Source File
# Begin Source File

SOURCE=.\WriteMultiFab.cpp
# End Source File
# End Group
# Begin Group "C++ Headers"

# PROP Default_Filter "H"
# Begin Source File

SOURCE=.\ABec_F.H
# End Source File
# Begin Source File

SOURCE=.\ABecLaplacian.H
# End Source File
# Begin Source File

SOURCE=.\BndryData.H
# End Source File
# Begin Source File

SOURCE=.\BoundCond.H
# End Source File
# Begin Source File

SOURCE=.\CG_F.H
# End Source File
# Begin Source File

SOURCE=.\CGSolver.H
# End Source File
# Begin Source File

SOURCE=.\InterpBndryData.H
# End Source File
# Begin Source File

SOURCE=.\INTERPBNDRYDATA_F.H
# End Source File
# Begin Source File

SOURCE=.\Laplacian.H
# End Source File
# Begin Source File

SOURCE=.\LinOp.H
# End Source File
# Begin Source File

SOURCE=.\LO_BCTYPES.H
# End Source File
# Begin Source File

SOURCE=.\LO_F.H
# End Source File
# Begin Source File

SOURCE=.\LP_F.H
# End Source File
# Begin Source File

SOURCE=.\Mask.H
# End Source File
# Begin Source File

SOURCE=.\MG_F.H
# End Source File
# Begin Source File

SOURCE=.\MultiGrid.H
# End Source File
# Begin Source File

SOURCE=.\WriteMultiFab.H
# End Source File
# End Group
# Begin Group "Fortran"

# PROP Default_Filter ""
# Begin Group "3D"

# PROP Default_Filter ""
# Begin Group "Temps"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\ABec_3D.FOR
# End Source File
# Begin Source File

SOURCE=.\Cg_3d.for
# End Source File
# Begin Source File

SOURCE=.\INTERPBNDRYDATA_3D.FOR
# End Source File
# Begin Source File

SOURCE=.\Lo_3d.for
# End Source File
# Begin Source File

SOURCE=.\Lp_3d.for
# End Source File
# Begin Source File

SOURCE=.\Mg_3d.for
# End Source File
# End Group
# Begin Source File

SOURCE=.\ABec_3D.F

!IF  "$(CFG)" == "mglib - Win32 Release"

# PROP Ignore_Default_Tool 1

!ELSEIF  "$(CFG)" == "mglib - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\ABec_3D.F
InputName=ABec_3D

"$(InputName).for" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo /S. /S..\amrlib /S..\bndrylib /S..\pBoxLib_2\
    /DBL_LANG_FORT                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE\
    /DBL_NO_FORT_FLUSH    $(InputName).F | perl          ..\scripts\strip72 -c >\
    $(InputName).for

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\CG_3D.F

!IF  "$(CFG)" == "mglib - Win32 Release"

# PROP Ignore_Default_Tool 1

!ELSEIF  "$(CFG)" == "mglib - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\CG_3D.F
InputName=CG_3D

"$(InputName).for" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo /S. /S..\amrlib /S..\bndrylib /S..\pBoxLib_2\
    /DBL_LANG_FORT                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE\
    /DBL_NO_FORT_FLUSH    $(InputName).F | perl          ..\scripts\strip72 -c >\
    $(InputName).for

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\INTERPBNDRYDATA_3D.F

!IF  "$(CFG)" == "mglib - Win32 Release"

# PROP Ignore_Default_Tool 1

!ELSEIF  "$(CFG)" == "mglib - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\INTERPBNDRYDATA_3D.F
InputName=INTERPBNDRYDATA_3D

"$(InputName).for" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo /S. /S..\amrlib /S..\bndrylib /S..\pBoxLib_2\
    /DBL_LANG_FORT                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE\
    /DBL_NO_FORT_FLUSH    $(InputName).F | perl          ..\scripts\strip72 -c >\
    $(InputName).for

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\LO_3D.F

!IF  "$(CFG)" == "mglib - Win32 Release"

# PROP Ignore_Default_Tool 1

!ELSEIF  "$(CFG)" == "mglib - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\LO_3D.F
InputName=LO_3D

"$(InputName).for" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo /S. /S..\amrlib /S..\bndrylib /S..\pBoxLib_2\
    /DBL_LANG_FORT                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE\
    /DBL_NO_FORT_FLUSH    $(InputName).F | perl          ..\scripts\strip72 -c >\
    $(InputName).for

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\LP_3D.F

!IF  "$(CFG)" == "mglib - Win32 Release"

# PROP Ignore_Default_Tool 1

!ELSEIF  "$(CFG)" == "mglib - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\LP_3D.F
InputName=LP_3D

"$(InputName).for" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo /S. /S..\amrlib /S..\bndrylib /S..\pBoxLib_2\
    /DBL_LANG_FORT                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE\
    /DBL_NO_FORT_FLUSH    $(InputName).F | perl          ..\scripts\strip72 -c >\
    $(InputName).for

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MG_3D.F

!IF  "$(CFG)" == "mglib - Win32 Release"

# PROP Ignore_Default_Tool 1

!ELSEIF  "$(CFG)" == "mglib - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\MG_3D.F
InputName=MG_3D

"$(InputName).for" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo /S. /S..\amrlib /S..\bndrylib /S..\pBoxLib_2\
    /DBL_LANG_FORT                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE\
    /DBL_NO_FORT_FLUSH    $(InputName).F | perl          ..\scripts\strip72 -c >\
    $(InputName).for

# End Custom Build

!ENDIF 

# End Source File
# End Group
# Begin Group "Temps No. 1"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Lo_util.for
# End Source File
# End Group
# Begin Source File

SOURCE=.\LO_UTIL.F

!IF  "$(CFG)" == "mglib - Win32 Release"

!ELSEIF  "$(CFG)" == "mglib - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\LO_UTIL.F
InputName=LO_UTIL

"$(InputName).for" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo /S. /S..\amrlib /S..\bndrylib /S..\pBoxLib_2\
    /DBL_LANG_FORT                 /DBL_SPACEDIM=3 /DBL_USE_DOUBLE\
    /DBL_NO_FORT_FLUSH    $(InputName).F | perl          ..\scripts\strip72 -c >\
    $(InputName).for

# End Custom Build

!ENDIF 

# End Source File
# End Group
# End Target
# End Project
