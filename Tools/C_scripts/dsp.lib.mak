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

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=$PROJ - Win32 $OptString
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "$PROJ.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "$PROJ.mak" CFG="$PROJ - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "$PROJ - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "$PROJ - Win32 Debug" (based on "Win32 (x86) Static Library")
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
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /assume:noaccuracy_sensitive /compile_only /debug:none /iface:cref /include:"Release/" /math_library:fast /nologo /threads /tune:k7 /warn:nofileopt /unroll:4
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_LIB" /YX /FD /D "_MBCS" /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 $cdirlist /I "C:\WMPI\include" /D "NDEBUG" /D "_LIB" /D "BL_LANG_CC" /D "WIN32" /D "_MBCS" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_SPACEDIM=${DIM}" ${CPROJDEF} ${CPROJVERS} /D "BL_PARALLEL_IO" /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D for="if(0);else for" /YX /FD /c
# ADD BASE RSC /l 0x409
# ADD RSC /l 0x409
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

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
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /dbglibs /debug:full /iface:cref /libs:dll /nologo /threads /warn:argument_checking /warn:nofileopt
# SUBTRACT F90 /traceback
# ADD BASE CPP /nologo /W3 /GX /Od /D "WIN32" /D "_DEBUG" /D "_LIB" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MDd /W3 /Gm /Gi /GR /GX /ZI /Od $cdirlist /I "C:\WMPI\include" /D "_CONSOLE" /D "_MBCS" /D "_DEBUG" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_THREADS" /D "BL_SPACEDIM=${DIM}" ${CPROJDEF} ${CPROJVERS} /D "BL_PARALLEL_IO" /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D for="if(0);else for" /FR /YX /FD /c
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
# Begin Group "TEMPS-Debug"

# PROP Default_Filter "FOR"
EOF

for file in $FEXE_sources
 do
   filer=`basename $file .F`
   ffile=$filer.FOR
   forfile=`echo $ffile | sed 's/\\//\\\/g'`
   cat >> $OFILE << EOF
# Begin Source File

SOURCE=Debug\\Fort\\$forfile

!IF  "\$(CFG)" == "$PROJ - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "\$(CFG)" == "$PROJ - Win32 Debug"

!ENDIF 
# End Source File
EOF
 done

cat >> $OFILE << EOF
# End Group
EOF

cat >> $OFILE << EOF
# Begin Group "TEMPS-Release"

# PROP Default_Filter "FOR"
EOF

for file in $FEXE_sources
 do
   filer=`basename $file .F`
   ffile=$filer.FOR
   forfile=`echo $ffile | sed 's/\\//\\\/g'`
   cat >> $OFILE << EOF
# Begin Source File

SOURCE=Release\\Fort\\$forfile

!IF  "\$(CFG)" == "$PROJ - Win32 Release"

!ELSEIF  "\$(CFG)" == "$PROJ - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 
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

"Release\\Fort\\\$(InputName).FOR" : \$(SOURCE) "\$(INTDIR)" "\$(OUTDIR)"
	fpp /m /ansi /nologo $fdirlist /DBL_LANG_FORT  /DBL_SPACEDIM=${DIM} /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM ${FPROJDEF} ${FPROJVERS} \$(InputPath) | perl ${TOP}\scripts\strip72 -c > Release\\Fort\\\$(InputName).FOR

# End Custom Build

!ELSEIF  "\$(CFG)" == "$PROJ - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=$ffile
InputName=$ifile

"Debug\\Fort\\\$(InputName).FOR" : \$(SOURCE) "\$(INTDIR)" "\$(OUTDIR)"
	fpp /m /ansi /nologo $fdirlist /DBL_LANG_FORT  /DBL_SPACEDIM=${DIM} /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM ${FPROJDEF} ${FPROJVERS} \$(InputPath) | perl ${TOP}\scripts\strip72 -c > Debug\\Fort\\\$(InputName).FOR

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

# Add lf's needed by Windows junk
perl -e 'while(<>) {s/\n/\r\n/;print;}' < $OFILE > junk.$$
mv junk.$$ $OFILE
exit 0

