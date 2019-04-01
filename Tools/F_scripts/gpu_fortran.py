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

Limitations:

- At present, this detects function imports in use statements if there is a use/only
  statement that imports only functions and ends in "! function". Cannot span multiple lines.

- In general use/only statements that continue to multiple lines via "&" are not handled properly

- Unrestricted private statements are disallowed, since we do not generate public statements for
  the newly created device functions. Doing this generation would require us to keep track of which
  ifdefs should be applied, since the IBM compiler won't allow you to generate a public statement
  if the corresponding function is ifdef'd out.
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

def case_insensitive_replace(s, old, new, count=-1, on_condition=None):
    # Does a case-insensitive replacement of string old
    # by string new in the string s.
    # Does count number of replacements.
    if not old or old==new:
        # We were passed the empty string for old.
        # or we're replacing old with itself.
        return s
    n_replaced = 0
    ridx_start = 0
    old = old.lower()
    new = new.lower()
    while(n_replaced < count or count == -1):
        if on_condition:
            istart, iend, found = on_condition(s[ridx_start:])
            if found:
                istart += ridx_start
                iend += ridx_start
            else:
                break
        else:
            try:
                istart = s.lower().index(old, ridx_start)
            except ValueError:
                break
            else:
                iend = istart + len(old)
        if new.lower() == old.lower() + '_device':
            actual_new = s[istart:iend] + '_device'
        else:
            actual_new = new
        s = s[:istart] + actual_new + s[iend:]
        ridx_start = istart + len(actual_new)
        n_replaced = n_replaced + 1
    return s

def get_procedure_type(procedure):

    procedure_type = "subroutine"

    if "function" in procedure.lower() and "subroutine" not in procedure.lower():
        procedure_type = "function"

    return procedure_type



def get_procedure_name(procedure):

    procedure_type = get_procedure_type(procedure)

    sprocedure_lower = procedure.lower().split()
    sprocedure = procedure.split()

    procedure_name = sprocedure[sprocedure_lower.index(procedure_type) + 1]

    # Now split out everything before the opening parenthesis, and truncate
    # out any spaces.

    procedure_name = procedure_name.split('(')[0].strip()

    return procedure_name



def append_device_to_line(line, subs):

    new_line = line

    # Look for statements of the form use module, only: In this case
    # we want to replace everything after the : but not before it.

    if 'use ' in line.lower() and 'only' in line.lower() and ':' in line:

        # Strip off any comments at the end of the line
        mcomment = re_detect_endcomment.search(line)
        comment = ''
        if mcomment:
            comment = line[mcomment.start():]
            line = line[:mcomment.start()]

        sline = line.split(':')

        new_line = sline[0] + ': '

        for var in sline[1].split(','):

            found_sub = False
            for sub_name in subs:
                if sub_name.lower() == var.strip().lower():
                    new_line += var.strip().lower() + '_device'
                    found_sub = True
                    break

            if not found_sub:
                new_line += var

            if var != sline[1].split(',')[-1]:
                new_line += ', '

        new_line = " ".join([new_line, comment])

    else:

        for sub_name in subs:
            re_subroutine_call = re.compile('(call ' + sub_name.lower() + '\s|\()')
            re_subroutine_device_call = re.compile('call ' + sub_name.lower() + '_device\s|\(')
            m_subroutine_call = re_subroutine_call.search(new_line.lower())
            m_subroutine_device_call = re_subroutine_device_call.search(new_line.lower())
            if m_subroutine_call and not m_subroutine_device_call:
                actual_new_line = new_line[:m_subroutine_call.start()]
                actual_new_line = actual_new_line + m_subroutine_call.group(1).replace(sub_name.lower(), sub_name + '_device')
                actual_new_line = actual_new_line + new_line[m_subroutine_call.end():]
                new_line = actual_new_line
                break

            if sub_name.lower() in line.lower().split() and sub_name.lower() + '_device' not in line.lower().split():
                new_line = case_insensitive_replace(new_line, sub_name, sub_name + '_device', 
                                                    on_condition=lambda x: get_replace_procedure(x, sub_name.lower()))

            elif sub_name in line.lower() and sub_name + '_device' not in line.lower():
                # Catch function calls here.
                # Make sure "bar(" is not any of "foobar(" or "foo % bar(" or "foo_bar("
                old_fun = sub_name.lower() + '('
                new_fun = sub_name.lower() + '_device' + '('
                re_old_fun = re.compile(sub_name.lower() + " *\(")
                re_not_fun = re.compile("((\w)|(% *))" + sub_name.lower() + " *\(")
                actual_new_line = ''
                while(new_line):
                    m_old_fun = re_old_fun.search(new_line.lower())
                    if m_old_fun:
                        m_not_fun = re_not_fun.search(new_line.lower())
                        if m_not_fun and m_old_fun.end() == m_not_fun.end():
                            # reject this replacement, it's not a function call
                            actual_new_line = actual_new_line + new_line[:m_not_fun.end()]
                            new_line = new_line[m_not_fun.end():]
                        else:
                            # do the replacement, it's a function call
                            substring = case_insensitive_replace(new_line[:m_old_fun.end()], old_fun, new_fun)
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

    mfun = re_detect_fun.search(line.lower())
    if 'use ' in line.lower() and 'only' in line.lower() and ':' in line and mfun:
        line = line[:mfun.start()]
        sline = line.split(':')
        for var in sline[1].split(','):
            if not '_device' in var.strip().lower() and not var.strip().lower() in imported_functions:
                imported_functions.append(var.strip())

    return imported_functions


def append_device(procedure):

    new_procedure = procedure
    
    called_subs = []
    for line in procedure.split('\n'):
        called_subs = called_subs + get_function_uses(line)
        if "call " in line.lower():
            try:
                ex_index = line.index('!')
            except:
                pass
            else:
                if ex_index > line.lower().index('call'):
                    pass
                else:
                    continue
            sub_name = line.lower().split('call ')[1].split('(')[0].strip()
            if sub_name != "syncthreads":
                if not "_device" in sub_name and not sub_name.lower() in lowerlist(called_subs):
                    called_subs.append(sub_name)

    new_procedure = [append_device_to_line(line, called_subs) for line in procedure.split('\n')]
    new_procedure = '\n'.join(new_procedure)

    return new_procedure


def lowerlist(l):
    x = [li.lower() for li in l]
    return x


def get_replace_procedure(line, pname):
    m = re.search('(?P<head>(\w|(% *))*)(?P<proc>' + pname.lower() + ')(?P<tail>((\w| *%| *=)*))', line.lower())
    istart = -1
    iend = -1
    found = False
    if m and m.group('proc') and not m.group('head') and not m.group('tail'):
        found = True
        istart = m.start('proc')
        iend = m.end('proc')
    return istart, iend, found


def create_device_version(procedure):

    procedure_name = get_procedure_name(procedure).lower()

    procedure_type = get_procedure_type(procedure)

    # Append _device to the name before writing it to the file.
    # We want to cover this both in the subroutine / end subroutine
    # statement and in the bind(c) statement.

    device_procedure = []

    re_subroutine_dec = re.compile(procedure_type + ' +' + procedure_name)

    for line in procedure.split('\n'):
        m_subroutine_dec = re_subroutine_dec.search(line.lower())
        if (m_subroutine_dec or "bind" in line.lower()) and not line.strip().startswith('!'):
            new_line = case_insensitive_replace(line, procedure_name, procedure_name + '_device', on_condition=lambda x: get_replace_procedure(x, procedure_name))
        else:
            new_line = line
        device_procedure.append(new_line)

    device_procedure = '\n'.join(device_procedure)

    # Make sure any procedure calls append _device to the procedure we are
    # calling. Note that this assumes that everything we are calling either
    # already has _device appended to it, or is being generated through this
    # script. We'll guard against known CUDA intrinsics.

    device_procedure = append_device(device_procedure)

    return device_procedure



def create_host_version(device_procedure):

    # Create a host version of the procedure as well.

    host_procedure = device_procedure.replace("_device","").replace("attributes(device)", "attributes(host)")
    assert(host_procedure.count('_device') == 0)

    # Some functions will have device code that is wrapped in #ifdef AMREX_USE_GPU_PRAGMA or
    # #if defined(AMREX_USE_GPU_PRAGMA), and host code otherwise. To deal with this, we'll undefine
    # AMREX_USE_GPU_PRAGMA inside the host version of the function.

    host_procedure = "#undef AMREX_USE_GPU_PRAGMA\n" + host_procedure + "\n#define AMREX_USE_GPU_PRAGMA\n"

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

        if ("subroutine" in line.lower() or "function" in line.lower()) and not found_procedure and not "module" in line.lower():

            # Check to make sure we are not in a comment.

            if '!' not in line:
                found_procedure = True
            else:
                if 'subroutine' in line.lower() and line.index('!') > line.lower().index('subroutine'):
                    found_procedure = True
                elif 'function' in line.lower() and line.index('!') > line.lower().index('function'):
                    found_procedure = True

        if ("end subroutine" in line.lower() or "end function" in line.lower()):

            procedure += line

            # Now if the procedure contains !$gpu, prepend
            # attributes(device) to the procedure.

            if "!$gpu" in procedure:

                # Some code may contain a legacy AMREX_DEVICE statement.

                if "AMREX_DEVICE" in procedure:
                    procedure = procedure.replace("AMREX_DEVICE", "attributes(device)", 1)
                else:
                    if get_procedure_type(procedure) == "subroutine":
                        procedure = case_insensitive_replace(procedure, "subroutine", "attributes(device) subroutine", 1)
                    elif get_procedure_type(procedure) == "function":
                        procedure = case_insensitive_replace(procedure, "function", "attributes(device) function", 1)
                    else:
                        sys.exit("Error - Procedure is not a subroutine or function!")

                # Capture the name of the procedure and add it to our targets list.

                if "subroutine" in line.lower() or "function" in line.lower():

                    procedure_type = get_procedure_type(procedure)

                    procedure_name = get_procedure_name(procedure)

                    if not procedure_name.lower() in gpu_targets:
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

            if host_procedure:
                fout.write(host_procedure)

            if device_procedure:
                fout.write(device_procedure)

            if procedure:
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

    subroutine_imported_functions = []
    
    for line in lines:

        if in_device_subroutine:
            # Detect any used functions from other modules
            subroutine_imported_functions = subroutine_imported_functions + get_function_uses(line)

        # If there are any calls inside a device procedure to
        # a known device procedure but it is missing the _device
        # suffix, add that now. The most likely case will be
        # calls to functions, where it wasn't possible in the initial
        # processing to know what is a function (versus a variable).

        for target in gpu_targets:
            if "end function " + target.lower() + "_device" in line.lower() or "end subroutine " + target.lower() + "_device" in line.lower():
                in_device_subroutine = False
                subroutine_imported_functions = []
                break

            elif "function " + target.lower() + "_device" in line.lower() or "subroutine " + target.lower() + "_device" in line.lower():
                in_device_subroutine = True
                break

        if in_device_subroutine:
            line = append_device_to_line(line, set(gpu_targets + subroutine_imported_functions))

        fout.write(line)

    fout.close()



if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--fortran",
                        help="the name of the Fortran file to update")

    args = parser.parse_args()

    update_fortran_procedures(args.fortran)
