#!/bin/bash

function cerr() {
  echo $@ 1>&2
}

function extract() {
  local sargs='C:f:'
  local largs='chdir:,file:'
  parsed=$(getopt --options ${sargs} --longoptions ${largs} --name "extract" -- "$@")

  if [[ $? != 0 ]]; then
    return 1
  fi

  echo "$parsed"
  eval set -- "$parsed"
  extract_dir=$PWD
  file=
  while true; do
    case $1 in
      --chdir|-C)
          extract_dir=$2
          if [ ! -d $extract_dir ]; then
            mkdir -p $extract_dir
          fi
          extract_dir=$(realpath $extract_dir)
          shift 2
        ;;
      --file|-f)
          file=$2
          shift 2
        ;;
      --)
          shift
          break
        ;;
    esac
  done
  if [ -z $file ]; then
    if [ -f $1 ]; then
      file=$1
    else
      cerr "ERROR: No file to extract"
    fi
  fi
  file=$(realpath $file)

  local extract_command=
  case $file in
    *\.tar\.gz*)
      extract_command="tar -xzf"
      ;;
    *\.tar\.bz*)
      extract_command="tar -xJf"
      ;;
    *\.zip)
      extract_command="unzip"
      ;;
    *)
      cerr "ERROR: Cannot detect extract format"
      return 1
  esac
  echo "Entering ${extract_dir}"
  pushd $extract_dir &> /dev/null
  ${extract_command} $file
  echo "Leaving ${extract_dir}"
  popd &> /dev/null
}

function parse_version() {
  if [ -z ${1} ]; then
    exit 1
  fi
  version=$1

  if [ -z ${2} ]; then
    prefix=VERSION
  else
    prefix=$2
  fi

  eval set -- MAJOR MINOR PATCH PATCH_EXTRA --
  for vcomp in ${version//\./ }; do
    export ${prefix}_${1}=${vcomp}
    shift
  done
}

function parse_version_env() {
  parse_version ${!1} ${1}
}
