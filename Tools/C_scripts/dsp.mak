#!/bin/sh
USAGE="Usage: $myName [-p project] [-t paralleltopdir] [-v prversion] [-d spacedim] [-o outfile] [-O <0, 1>] [-f namefile] [-S] filenames"

myName=$0
PROJ=IAMRAll
TOP=..\..
DIM=2
OFILE=$PROJ.dsp
OLEVEL=0           # 0 - DEBUG, 1 - Release
VERSION=0
SERIAL=0
NFILE=

#
#  Parse command arguments, exit if incorrect usage
#
if [ $# = 0 ]
then
        echo $USAGE
        exit
fi
#set -- `getopt p:t:d:o:O:f:v:S $*`
if [ $? != 0 ]
  then
    echo $USAGE 
    exit
fi
for argument in $*
do
  case $argument in
  -p) shift;PROJ=$1;shift;;
  -t) shift;TOP=`echo $1 | sed 's/\\//\\\/g'`;shift;;
  -d) shift;DIM=$1;shift;;
  -o) shift;OFILE=$1;shift;;
  -O) shift;OLEVEL=$1;shift;;
  -f) shift;NFILE=$1;shift;;
  -v) shift;VERSION=$1;shift;;
  -S) shift;SERIAL=1;;
  --) shift;break;;
  esac
done

files=$*

vN=5
if [ ${VERSION} = v9 ]
  then
    vN=9
fi
CPROJVERS='/D "BL_PRVERSION='${vN}\"
FPROJVERS='/DBL_PRVERSION='${vN}
CPROJDEF=
FPROJDEF=
if [ ${SERIAL} = 1 ]
 then
   CPROJDEF='/D "BL_USE_HGPROJ_SERIAL"'
   FPROJDEF='/DBL_USE_HGPROJ_SERIAL'
fi

if [ -n $NFILE ]
  then
    files=`cat $NFILE` $files
fi

OptString=Debug
if [ $OLEVEL = 1 ] 
then
  OptString=Release
fi
#
# Split up the infiles
#
FEXE_sources=
CEXE_sources=
CEXE_headers=
fEXE_sources=
IncludeLocations=

for f in `cat $NFILE`
do
    if [ `echo $f | fgrep .F` ]
    then
        FEXE_sources="${FEXE_sources} $f"
    else if [ `echo $f | fgrep .cpp` ]
       then
           CEXE_sources="${CEXE_sources} $f"
       else if [ `echo $f | grep '\.[hH]'` ]
          then
              CEXE_headers="${CEXE_headers} $f"
          else if [ `echo $f | fgrep .f` ]
             then
                 fEXE_sources="${fEXE_sources} $f"
             else
                 dir=`echo $f | tr '\\/' '\\\'`
                 IncludeLocations="$IncludeLocations $dir"
             fi
          fi
       fi
    fi
done

cdirlist=
for d in $IncludeLocations
do
  cdirlist="$cdirlist /I \"$d\""
done

fdirlist=
for d in $IncludeLocations
do
  fdirlist="$fdirlist /S$d"
done

cat > $OFILE << EOF
# Microsoft Developer Studio Project File - Name="$PROJ" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=$PROJ - Win32 $OptString
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "$PROJ.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "$PROJ.mak" CFG="$PROJ - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "$PROJ - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "$PROJ - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "\$(CFG)" == "$PROJ - Win32 Release"

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
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /assume:noaccuracy_sensitive /compile_only /debug:none /iface:cref /include:"Release/" /math_library:fast /nologo /threads /tune:k7 /warn:nofileopt /unroll:4
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GR /GX /O2 $cdirlist /I "C:\WMPI\include" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D "BL_SPACEDIM=${DIM}" ${CPROJDEF} ${CPROJVERS} /D "BL_PARALLEL_IO" /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_USE_CHEM" /D for="if(0);else for" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib dformt.lib /nologo /subsystem:console /machine:I386 /include:"__matherr"

!ELSEIF  "\$(CFG)" == "$PROJ - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /dbglibs /debug:full /iface:cref /include:"Debug/" /libs:static /nologo /traceback /warn:argument_checking /optimize:0 /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MTd /W3 /Gm /GR /GX /ZI /Od $cdirlist /I "C:\WMPI\include" /D "_CONSOLE" /D "_MBCS" /D "_DEBUG" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D "BL_SPACEDIM=${DIM}" ${CPROJDEF} ${CPROJVERS} /D "BL_PARALLEL_IO" /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_USE_CHEM" /D for="if(0);else for" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib dformt.lib /nologo /subsystem:console /debug /machine:I386 /include:"__matherr" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "$PROJ - Win32 Release"
# Name "$PROJ - Win32 Debug"
EOF

#
# Do C++ sources
#
cat >> $OFILE <<EOF
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
EOF

for file in $CEXE_sources
 do
   cppfile=`echo $file | sed 's/\\//\\\/g'`
   cat >> $OFILE << EOF
# Begin Source File

SOURCE=$cppfile
# End Source File
EOF
 done
cat >> $OFILE << EOF
# End Group
EOF


#
# Do C++ headers
#
cat >> $OFILE <<EOF
# Begin Group "C++ Headers"

# PROP Default_Filter "H"
EOF

for file in $CEXE_headers
 do
   cpphdr=`echo $file | sed 's/\\//\\\/g'`
   cat >> $OFILE << EOF
# Begin Source File

SOURCE=$cpphdr
# End Source File
EOF
 done
cat >> $OFILE << EOF
# End Group
EOF

#
# Do processed fortran sources
#
cat >> $OFILE << EOF
# Begin Group "Fortran"

# PROP Default_Filter ""
EOF

for file in $fEXE_sources
do
   forfile=`echo $file | sed 's/\\//\\\/g'`
   cat >> $OFILE << EOF
# Begin Source File

SOURCE=$forfile
# End Source File
EOF
 done

cat >> $OFILE << EOF
# Begin Group "TEMPS"

# PROP Default_Filter "FOR"
EOF

for file in $FEXE_sources
 do
   filer=`basename $file .F`
   ffile=$filer.FOR
   forfile=`echo $ffile | sed 's/\\//\\\/g'`
   cat >> $OFILE << EOF
# Begin Source File

SOURCE=${OptString}\\$forfile
# End Source File
EOF
 done

cat >> $OFILE << EOF
# End Group
EOF

#
# Do preprocessing for Fortran sources
#
for file in $FEXE_sources
 do
   ifile=`basename $file .F`
   ffile=`echo $file | sed 's/\\//\\\/g'`
   cat >> $OFILE << EOF
# Begin Source File

SOURCE=$ffile

!IF  "\$(CFG)" == "$PROJ - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=$ffile
InputName=$ifile

"\$(IntDir)\\\$(InputName).FOR" : \$(SOURCE) "\$(INTDIR)" "\$(OUTDIR)"
	fpp /m /ansi /nologo $fdirlist /DBL_LANG_FORT  /DBL_SPACEDIM=${DIM} /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM ${FPROJDEF} ${FPROJVERS} \$(InputPath) | perl ${TOP}\scripts\strip72 -c > \$(IntDir)\\\$(InputName).FOR

# End Custom Build

!ELSEIF  "\$(CFG)" == "$PROJ - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=$ffile
InputName=$ifile

"\$(IntDir)\\\$(InputName).FOR" : \$(SOURCE) "\$(INTDIR)" "\$(OUTDIR)"
	fpp /m /ansi /nologo $fdirlist /DBL_LANG_FORT  /DBL_SPACEDIM=${DIM} /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM ${FPROJDEF} ${FPROJVERS} \$(InputPath) | perl ${TOP}\scripts\strip72 -c > \$(IntDir)\\\$(InputName).FOR

# End Custom Build

!ENDIF

# End Source File
EOF
 done

cat >> $OFILE << EOF
# End Group
EOF

cat >> $OFILE << EOF
# End Target
# End Project
EOF

exit 0

