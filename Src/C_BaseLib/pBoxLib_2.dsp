# Microsoft Developer Studio Project File - Name="pboxlib_2" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=pboxlib_2 - Win32 Debug3D
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "pboxlib_2.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "pboxlib_2.mak" CFG="pboxlib_2 - Win32 Debug3D"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "pboxlib_2 - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "pboxlib_2 - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "pboxlib_2 - Win32 Debug3D" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=xicl.exe
F90=df.exe

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

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
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

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
# ADD F90 /browser /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /iface:cref /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /FR /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug3D"
# PROP BASE Intermediate_Dir "Debug3D"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug3D"
# PROP Intermediate_Dir "Debug3D"
# PROP Target_Dir ""
# ADD BASE F90 /include:"Debug3D/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD F90 /browser /include:"Debug3D/" /compile_only /nologo /debug:full /optimize:0 /iface:cref /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /Gi- /GX /Z7 /Od /D "_WINDOWS" /D "WIN32" /D "_DEBUG" /D BL_SPACEDIM=3 /FR /YX /FD /c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "pboxlib_2 - Win32 Release"
# Name "pboxlib_2 - Win32 Debug"
# Name "pboxlib_2 - Win32 Debug3D"
# Begin Group "C Sources"

# PROP Default_Filter "*.cpp"
# Begin Source File

SOURCE=.\aString.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\BArena.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Box.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\BoxArray.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\BoxAssoc.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\BoxDomain.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\BoxLib.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\BoxList.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\DistributionMapping.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FabConv.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FArrayBox.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FPC.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\IndexType.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\IntVect.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MultiFab.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Orientation.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\ParallelDescriptor.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\ParmParse.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Tracer.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Utility.cpp

!IF  "$(CFG)" == "pboxlib_2 - Win32 Release"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug"

# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=2

!ELSEIF  "$(CFG)" == "pboxlib_2 - Win32 Debug3D"

# ADD BASE CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I ".\\" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_LANG_CC" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for" /D BL_SPACEDIM=3

!ENDIF 

# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "*.H"
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

SOURCE=.\BLVERSION.H
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

SOURCE=.\CONSTANTS.H
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

SOURCE=.\FPC.H
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

SOURCE=.\REAL.H
# End Source File
# Begin Source File

SOURCE=.\SPACE.H
# End Source File
# Begin Source File

SOURCE=.\SPACE_F.H
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
# End Group
# End Target
# End Project
