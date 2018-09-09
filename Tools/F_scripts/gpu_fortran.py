#!/usr/bin/env python

"""
Search the Fortran source code for procedures marked with !$gpu in their specification part:

  subroutine sub(a)

    integer :: a

    !$gpu

    ...

  end subroutine sub

For each one, generate a device copy by prepending attributes(device) prior to the procedure statement:

  attributes(device) subroutine sub_device(a)

    integer :: a

    !$gpu

    ...

  end subroutine sub_device

The host copy is retained with the same name.
"""

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for gpu_fortran.py")

if sys.version[0] == "2":
    reload(sys)
    sys.setdefaultencoding('utf8')

import os
import argparse



def get_procedure_type(procedure):

    procedure_type = "subroutine"

    if "function" in procedure and "subroutine" not in procedure:
        procedure_type = "function"

    return procedure_type



def get_procedure_name(procedure):

    procedure_type = get_procedure_type(procedure)

    sprocedure = procedure.split()

    procedure_name = sprocedure[sprocedure.index(procedure_type) + 1]

    # Now split out everything before the opening parenthesis, and truncate
    # out any spaces.

    procedure_name = procedure_name.split('(')[0].strip()

    return procedure_name



def append_device_to_line(line, subs):

    new_line = line

    # Look for statements of the form use module, only: In this case
    # we want to replace everything after the : but not before it.

    if 'use ' in line and 'only' in line and ':' in line:

        sline = line.split(':')

        new_line = sline[0] + ':'

        for var in sline[1].split(','):

            found_sub = False
            for sub_name in subs:
                if sub_name == var.strip():
                    new_line += var + '_device'
                    found_sub = True
                    break

            if not found_sub:
                new_line += var

            if var != sline[1].split(',')[-1]:
                new_line += ','

    else:

        for sub_name in subs:

            if 'call ' + sub_name + '(' in line and 'call ' + sub_name + '_device(' not in line:
                new_line = line.split('(')[0] + '_device('
                for word in line.split('(')[1:]:
                    new_line += word
                    if word != line.split('(')[-1]:
                        new_line += '('
                break

            elif sub_name in line.split() and sub_name + '_device' not in line.split():
                new_line = new_line.replace(sub_name, sub_name + '_device')
                break

            elif sub_name + '(' in line and sub_name + '_device' not in line:
                # Catch function calls here.
                new_line = new_line.replace(sub_name + '(', sub_name + '_device' + '(')
                break

    return new_line


def append_device(procedure):

    new_procedure = procedure
    
    called_subs = []
    for line in procedure.split('\n'):
        if "call " in line:
            sub_name = line.split('call ')[1].split('(')[0]
            if sub_name != "syncthreads":
                if "_device" not in sub_name:
                    called_subs.append(sub_name)

    in_specification_part = True

    for line in procedure.split('\n'):
        new_procedure = new_procedure.replace(line, append_device_to_line(line, called_subs))

    return new_procedure



def create_device_version(procedure):

    procedure_name = get_procedure_name(procedure)

    procedure_type = get_procedure_type(procedure)

    device_procedure = procedure

    # Append _device to the name before writing it to the file.
    # We want to cover this both in the subroutine / end subroutine
    # statement and in the bind(c) statement.

    for chunk in procedure.split(procedure_name):
        if procedure_type in chunk or "bind" in chunk:
            chunk += procedure_name
            new_chunk = chunk.replace(procedure_name, procedure_name + '_device')
            device_procedure = device_procedure.replace(chunk, new_chunk)

    # Make sure any procedure calls append _device to the procedure we are
    # calling. Note that this assumes that everything we are calling either
    # already has _device appended to it, or is being generated through this
    # script. We'll guard against known CUDA intrinsics.

    device_procedure = append_device(device_procedure)

    return device_procedure



def create_host_version(device_procedure):

    # Create a host version of the procedure as well.

    host_procedure = device_procedure.replace("_device","").replace("attributes(device)", "attributes(host)")

    # Some functions will have device code that is wrapped in #ifdef AMREX_USE_CUDA or
    # #if defined(AMREX_USE_CUDA), and host code otherwise. To deal with this, we'll undefine
    # AMREX_USE_CUDA inside the host version of the function.

    host_procedure = "#undef AMREX_USE_CUDA\n" + host_procedure + "\n#define AMREX_USE_CUDA\n"

    return host_procedure



def update_fortran_procedures(ffile):
    """For procedures marked up with !$gpu, generate a copy marked
    with attributes(device), and create a host copy as well, to
    avoid the CUF restriction on attributes(host, device) In order
    to avoid name collisions, we'll name one with a _device suffix
    and the other with a _host suffix, and then create an interface
    in the module so that they can both be addressed with the original
    procedure name."""

    gpu_targets = []

    # open the input Fortran file
    try:
        fin = open(ffile, "r")
    except IOError:
        sys.exit("Cannot open Fortran file {} for reading".format(ffile))

    # read in the file, close it, and then reopen it for writing
    lines = fin.readlines()

    fin.close()

    try:
        fout = open(ffile, "w")
    except IOError:
        sys.exit("Cannot open Fortran file {} for writing".format(ffile))

    # loop through the file and look for all procedures

    found_procedure = False
    device_procedure = ""
    host_procedure = ""
    procedure = ""

    for line in lines:

        if ("subroutine" in line or "function" in line) and not found_procedure and not "module" in line:

            # Check to make sure we are not in a comment.

            if '!' not in line:
                found_procedure = True
            else:
                if 'subroutine' in line and line.index('!') > line.index('subroutine'):
                    found_procedure = True
                elif 'function' in line and line.index('!') > line.index('function'):
                    found_procedure = True

        if ("end subroutine" in line or "end function" in line):

            procedure += line

            # Now if the procedure contains !$gpu, prepend
            # attributes(device) to the procedure.

            if "!$gpu" in procedure:

                # Some code may contain a legacy AMREX_DEVICE statement.

                if "AMREX_DEVICE" in procedure:
                    procedure = procedure.replace("AMREX_DEVICE", "attributes(device)", 1)
                else:
                    procedure = procedure.replace("subroutine", "attributes(device) subroutine", 1)
                    procedure = procedure.replace("function", "attributes(device) function", 1)

                # Capture the name of the procedure and add it to our targets list.

                if "subroutine" in line or "function" in line:

                    procedure_type = get_procedure_type(procedure)

                    procedure_name = get_procedure_name(procedure)

                    gpu_targets.append(procedure_name)

                    # Create device and host versions of the procedure.

                    device_procedure = create_device_version(procedure)

                    host_procedure = create_host_version(device_procedure)

                    # Blank out the original procedure; there is nothing left to do with it.

                    procedure = ""



            elif "AMREX_LAUNCH" in procedure:

                # Some code may contain a legacy AMREX_LAUNCH statement.
                # In this case all we want to do is update the calls
                # inside to append_device, but no other processing
                # needs to be done.

                procedure = append_device(procedure)



            # Now write out whatever we have generated.

            if host_procedure != "":
                fout.write(host_procedure)

            if device_procedure != "":
                fout.write(device_procedure)

            if procedure != "":
                fout.write(procedure)

            # Clear everything out for the next procedure.

            device_procedure = ""
            host_procedure = ""
            procedure = ""
            line = ""
            found_procedure = False

        if found_procedure:
            procedure += line
        else:
            fout.write(line)

    fout.close()



    # Now handle any remaining logistics that required a full pass through.

    fin = open(ffile, 'r')
    lines = fin.readlines()
    fin.close()

    fout = open(ffile, 'w')

    in_device_subroutine = False
    
    for line in lines:

        # Explicitly mark all newly created device targets as public.

        if "contains" in line:

            for target in gpu_targets:
                fout.write('public :: ' + target + '_device\n')

        # If there are any calls inside a device procedure to
        # a known device procedure but it is missing the _device
        # suffix, add that now. The most likely case will be
        # calls to functions, where it wasn't possible in the initial
        # processing to know what is a function (versus a variable).

        else:

            for target in gpu_targets:

                if "end function " + target + "_device" in line or "end subroutine " + target + "_device" in line:
                    in_device_subroutine = False
                    break

                elif "function " + target + "_device" in line or "subroutine " + target + "_device" in line:
                    in_device_subroutine = True
                    break

            if in_device_subroutine:

                line = append_device_to_line(line, gpu_targets)

        fout.write(line)

    fout.close()



if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--fortran",
                        help="the name of the Fortran file to update")

    args = parser.parse_args()

    update_fortran_procedures(args.fortran)
