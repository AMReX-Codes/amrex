# this function borrowed from PlPlot, Copyright (C) 2006  Alan W. Irwin
function(TRANSFORM_VERSION numerical_result version)
# internal_version ignores everything in version after any character that
# is not 0-9 or ".".  This should take care of the case when there is
# some non-numerical data in the patch version.
message(STATUS "DEBUG: version = ${version}")
string(REGEX REPLACE "^([0-9.]+).*$" "\\1" internal_version ${version})

   # internal_version is normally a period-delimited triplet string of the form
   # "major.minor.patch", but patch and/or minor could be missing.
   # Transform internal_version into a numerical result that can be compared.
   string(REGEX REPLACE "^([0-9]*).+$" "\\1" major ${internal_version})
   string(REGEX REPLACE "^[0-9]*\\.([0-9]*).*$" "\\1" minor ${internal_version})
   string(REGEX REPLACE "^[0-9]*\\.[0-9]*\\.([0-9]*)$" "\\1" patch ${internal_version})

   if(NOT patch MATCHES "^[0-9]+$")
     set(patch 0)
   endif(NOT patch MATCHES "^[0-9]+$")

   if(NOT minor MATCHES "[0-9]+")
     set(minor 0)
   endif(NOT minor MATCHES "[0-9]+")

   if(NOT major MATCHES "[0-9]+")
     set(major 0)
   endif(NOT major MATCHES "[0-9]+")
   #message(STATUS "DEBUG: internal_version = ${internal_version}")
   #message(STATUS "DEBUG: major = ${major}")
   #message(STATUS "DEBUG: minor= ${minor}")
   #message(STATUS "DEBUG: patch = ${patch}")
   math(EXPR internal_numerical_result
     "${major}*1000000 + ${minor}*1000 + ${patch}"
     #"${major}*1000000 + ${minor}*1000"
     )
   #message(STATUS "DEBUG: ${numerical_result} = ${internal_numerical_result}")
   set(${numerical_result} ${internal_numerical_result} PARENT_SCOPE)
endfunction(TRANSFORM_VERSION)
