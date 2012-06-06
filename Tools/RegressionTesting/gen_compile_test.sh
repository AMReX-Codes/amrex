#!/bin/sh

# script to automatically generate the test parameters for compile tests
# for MAESTRO problems.  Run this in the MAESTRO/ directory.

for i in `find . -name GNUmakefile -print`
do 
   buildDir=`dirname $i`
   topDir=`dirname $buildDir`

   category=`basename $topDir`
   problem=`basename $buildDir`

   if [ "$category" != "docs" ]; then
       if [ "$problem" != "docs" ]; then
	   echo [$category-$problem]
	   echo buildDir = fParallel/MAESTRO/$buildDir
	   echo compileTest = 1
	   echo " "
       fi
   fi
done

