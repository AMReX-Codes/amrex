# Microsoft Developer Studio Project File - Name="t_boxlib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=t_boxlib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "t_boxlib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "t_boxlib.mak" CFG="t_boxlib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "t_boxlib - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "t_boxlib - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "t_boxlib - Win32 Release"

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
# ADD F90 /compile_only /iface:cref /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib  kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib  kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "t_boxlib - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /iface:cref /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib  kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib  kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "t_boxlib - Win32 Release"
# Name "t_boxlib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\bl_constants.f90
NODEP_F90_BL_CO=\
	".\Debug\bl_types.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\bl_error.f90
NODEP_F90_BL_ER=\
	".\Debug\bl_types.mod"\
	".\Debug\parallel.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\bl_IO.f90
NODEP_F90_BL_IO=\
	".\Debug\bl_error_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\bl_kiss.f90
# End Source File
# Begin Source File

SOURCE=..\bl_mem_stat.f90
NODEP_F90_BL_ME=\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_string_module.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\parallel.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\bl_parmparse.f90
DEP_F90_BL_PA=\
	".\Debug\bl_error_module.mod"\
	".\Debug\f2kcli.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\bl_space.f90
# End Source File
# Begin Source File

SOURCE=..\bl_stream.f90
NODEP_F90_BL_ST=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_string_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\bl_string.f90
NODEP_F90_BL_STR=\
	".\Debug\bl_error_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\bl_timer.f90
NODEP_F90_BL_TI=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_string_module.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\parallel.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\bl_types.f90
# End Source File
# Begin Source File

SOURCE=..\bndry_reg.f90
NODEP_F90_BNDRY=\
	".\Debug\layout_module.mod"\
	".\Debug\multifab_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\box.f90
NODEP_F90_BOX_F=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_space.mod"\
	".\Debug\bl_stream_module.mod"\
	".\Debug\bl_string_module.mod"\
	".\Debug\bl_types.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\box_util.f90
NODEP_F90_BOX_U=\
	".\Debug\bl_IO_module.mod"\
	".\Debug\box_module.mod"\
	".\Debug\boxarray_module.mod"\
	".\Debug\mboxarray_module.mod"\
	".\Debug\mt19937_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\boxarray.f90
NODEP_F90_BOXAR=\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_mem_stat_module.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\box_module.mod"\
	".\Debug\list_box_module.mod"\
	".\Debug\sort_box_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\BoxLib.f90
NODEP_F90_BOXLI=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_space.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\parallel.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\cluster.f90
NODEP_F90_CLUST=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\box_module.mod"\
	".\Debug\list_box_module.mod"\
	".\Debug\multifab_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\f2kcli_win32.f90
# End Source File
# Begin Source File

SOURCE=..\f2kgetcl.c
# ADD CPP /D "UPPER"
# End Source File
# Begin Source File

SOURCE=..\fab.f90
NODEP_F90_FAB_F=\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_mem_stat_module.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\box_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\fabio.f90
NODEP_F90_FABIO=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_string_module.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\fab_module.mod"\
	".\Debug\multifab_module.mod"\
	".\Debug\parallel.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\fabio_c.c
# End Source File
# Begin Source File

SOURCE=..\filler.f90
NODEP_F90_FILLE=\
	".\Debug\plotfile_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\interp.f90
NODEP_F90_INTER=\
	".\Debug\bl_constants_module.mod"\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_types.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\knapsack.f90
NODEP_F90_KNAPS=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\sort_d_module.mod"\
	".\Debug\vector_i_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\layout.f90
NODEP_F90_LAYOU=\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_mem_stat_module.mod"\
	".\Debug\boxarray_module.mod"\
	".\Debug\knapsack_module.mod"\
	".\Debug\parallel.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\list_box.f90
NODEP_F90_LIST_=\
	".\Debug\bl_error_module.mod"\
	".\Debug\box_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\main.f90
NODEP_F90_MAIN_=\
	".\Debug\BoxLib.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\mboxarray.f90
NODEP_F90_MBOXA=\
	".\Debug\bl_IO_module.mod"\
	".\Debug\boxarray_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\mt19937ar.f90
NODEP_F90_MT199=\
	".\Debug\bl_error_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\multifab.f90
NODEP_F90_MULTI=\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_mem_stat_module.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\fab_module.mod"\
	".\Debug\layout_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\omp.f90
# End Source File
# Begin Source File

SOURCE=..\omp_stubs.f90
# End Source File
# Begin Source File

SOURCE=..\parallel_stubs.f90
NODEP_F90_PARAL=\
	".\Debug\bl_types.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\plotfile.f90
NODEP_F90_PLOTF=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_space.mod"\
	".\Debug\bl_stream_module.mod"\
	".\Debug\box_module.mod"\
	".\Debug\fabio_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\ppm_util.f90
NODEP_F90_PPM_U=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_string_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\ppm_util_c.c
# End Source File
# Begin Source File

SOURCE=..\sort_box.f90
NODEP_F90_SORT_=\
	".\Debug\box_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\sort_d.f90
NODEP_F90_SORT_D=\
	".\Debug\bl_types.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\sort_i.f90
# End Source File
# Begin Source File

SOURCE=.\t_main.f90
NODEP_F90_T_MAI=\
	".\Debug\bl_error_module.mod"\
	".\Debug\bl_IO_module.mod"\
	".\Debug\bl_kiss_module.mod"\
	".\Debug\bl_types.mod"\
	".\Debug\box_module.mod"\
	".\Debug\box_util_module.mod"\
	".\Debug\boxarray_module.mod"\
	".\Debug\cluster_module.mod"\
	".\Debug\f2kcli.mod"\
	".\Debug\fab_module.mod"\
	".\Debug\fabio_module.mod"\
	".\Debug\layout_module.mod"\
	".\Debug\list_box_module.mod"\
	".\Debug\mboxarray_module.mod"\
	".\Debug\mt19937_module.mod"\
	".\Debug\multifab_module.mod"\
	".\Debug\plotfile_module.mod"\
	".\Debug\sort_box_module.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\timer_c.c
# End Source File
# Begin Source File

SOURCE=..\vector_i.f90
NODEP_F90_VECTO=\
	".\Debug\bl_IO_module.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
