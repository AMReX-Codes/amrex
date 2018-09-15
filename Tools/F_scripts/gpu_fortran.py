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
import re


re_detect_fun = re.compile('! *function\s')
re_detect_endcomment = re.compile('!')


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

        # Strip off any comments at the end of the line
        mcomment = re_detect_endcomment.search(line)
        comment = ''
        if mcomment:
            comment = line[mcomment.start():]
            line = line[:mcomment.start()]

        sline = line.split(':')

        new_line = sline[0] + ':'

        for var in sline[1].split(','):

            found_sub = False
            for sub_name in subs:
                if sub_name == var.strip():
                    new_line += var.strip() + '_device'
                    found_sub = True
                    break

            if not found_sub:
                new_line += var

            if var != sline[1].split(',')[-1]:
                new_line += ','

        new_line = new_line + " " + comment

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
                re_argument = re.compile("\S+ +" + sub_name + " +=")
                if not re_argument.search(line):
                    new_line = new_line.replace(sub_name, sub_name + '_device')
                    break
                else:
                    continue

            elif sub_name + '(' in line and sub_name + '_device' not in line:
                # Catch function calls here.
                # Make sure "bar(" is not any of "foobar(" or "foo % bar(" or "foo_bar("
                old_fun = sub_name + '('
                new_fun = sub_name + '_device' + '('
                re_old_fun = re.compile(sub_name + "\(")
                re_not_fun = re.compile("((\w)|(% *))" + sub_name + "\(")
                actual_new_line = ''
                while(new_line):
                    m_old_fun = re_old_fun.search(new_line)
                    if m_old_fun:
                        m_not_fun = re_not_fun.search(new_line)
                        if m_not_fun and m_old_fun.end() == m_not_fun.end():
                            # reject this replacement, it's not a function call
                            actual_new_line = actual_new_line + new_line[:m_not_fun.end()]
                            new_line = new_line[m_not_fun.end():]
                        else:
                            # do the replacement, it's a function call
                            substring = new_line[:m_old_fun.end()].replace(old_fun, new_fun)
                            actual_new_line = actual_new_line + substring
                            new_line = new_line[m_old_fun.end():]
                    else:
                        actual_new_line = actual_new_line + new_line
                        break
                new_line = actual_new_line
                # continue because there could be multiple functions calls on the same line
                # or functions that satisfy the elif but not necessarily the regex criteria.
                continue

    return new_line


def get_function_uses(line):
    # Device functions in another module will
    # be detected expecting that the 'use/only' line contains
    # nothing but functions and ends with "! function"

    imported_functions = []

    mfun = re_detect_fun.search(line)
    if 'use ' in line and 'only' in line and ':' in line and mfun:
        line = line[:mfun.start()]
        sline = line.split(':')
        for var in sline[1].split(','):
            imported_functions.append(var.strip())

    print(line)
    print(imported_functions)

    return imported_functions


def append_device(procedure):

    new_procedure = procedure
    
    called_subs = []
    for line in procedure.split('\n'):
        called_subs = called_subs + get_function_uses(line)
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
    and the other with a _host suffix."""

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
                    if get_procedure_type(procedure) == "subroutine":
                        procedure = procedure.replace("subroutine", "attributes(device) subroutine", 1)
                    elif get_procedure_type(procedure) == "function":
                        procedure = procedure.replace("function", "attributes(device) function", 1)
                    else:
                        sys.exit("Error - Procedure is not a subroutine or function!")

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

    module_imported_functions = []
    subroutine_imported_functions = []
    in_module_header = True
    
    for line in lines:

        if in_module_header:
            # Detect any used functions from other modules
            module_imported_functions = module_imported_functions + get_function_uses(line)

        if in_device_subroutine:
            # Detect any used functions from other modules
            subroutine_imported_functions = subroutine_imported_functions + get_function_uses(line)

        # Explicitly mark all newly created device targets as public.

        if "contains" in line:

            in_module_header = False

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
                    subroutine_imported_functions = []
                    break

                elif "function " + target + "_device" in line or "subroutine " + target + "_device" in line:
                    print(line)
                    in_device_subroutine = True
                    break

            if in_device_subroutine:
                print(gpu_targets + module_imported_functions + subroutine_imported_functions)
                line = append_device_to_line(line, gpu_targets + module_imported_functions + subroutine_imported_functions)

        fout.write(line)

    fout.close()



if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--fortran",
                        help="the name of the Fortran file to update")

    args = parser.parse_args()

    update_fortran_procedures(args.fortran)
