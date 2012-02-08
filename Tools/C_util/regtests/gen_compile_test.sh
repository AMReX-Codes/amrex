#!/bin/sh

for i in `find . -name GNUmakefile -print`
do 
   temp=`dirname $i`
   echo [`basename $temp`]
   echo buildDir = fParallel/MAESTRO/`dirname $i`
   echo compileTest = 1
   echo " "
done

