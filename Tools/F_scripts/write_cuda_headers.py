#!/usr/bin/env python

"""
Search the Fortran source code for subroutines marked as:

  AMREX_DEVICE subroutine sub(a)

    ...

  end subroutine

and maintain a list of these.

Then copy the C++ headers for Fortran files (typically *_F.H) into a
temp directory, modifying any subroutines marked with #pragma gpu to
have both a device and host signature.
"""

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for dep.py")

if sys.version[0] == "2":
    reload(sys)
    sys.setdefaultencoding('utf8')

import os
import re
import argparse

import find_files_vpath as ffv
import preprocess


TEMPLATE = """
__global__ static void cuda_{}
{{
   int blo[3];
   int bhi[3];
   for (int k = lo_3 + blockIdx.z * blockDim.z + threadIdx.z; k <= hi_3; k += blockDim.z * gridDim.z) {{
     blo[2] = k;
     bhi[2] = k;
     for (int j = lo_2 + blockIdx.y * blockDim.y + threadIdx.y; j <= hi_2; j += blockDim.y * gridDim.y) {{
       blo[1] = j;
       bhi[1] = j;
       for (int i = lo_1 + blockIdx.x * blockDim.x + threadIdx.x; i <= hi_1; i += blockDim.x * gridDim.x) {{
         blo[0] = i;
         bhi[0] = i;
         {};
       }}
     }}
   }}
}}
"""

# for finding just the variable definitions in the function signature (between the ())
decls_re = re.compile("(.*?)(\\()(.*)(\\))", re.IGNORECASE|re.DOTALL)

# for finding a fortran function subroutine marked with AMREX_DEVICE or attributes(device)
fortran_re = re.compile("(AMREX_DEVICE)(\\s+)(subroutine)(\\s+)((?:[a-z][a-z_]+))",
                        re.IGNORECASE|re.DOTALL)
fortran_attributes_re = re.compile("(attributes\\(device\\))(\\s+)(subroutine)(\\s+)((?:[a-z][a-z_]+))",
                                   re.IGNORECASE|re.DOTALL)

# for finding a header entry for a function binding to a fortran subroutine
fortran_binding_re = re.compile("(void)(\\s+)((?:[a-z][a-z_]+))",
                                re.IGNORECASE|re.DOTALL)

class HeaderFile(object):
    """ hold information about one of the headers """

    def __init__(self, filename):

        self.name = filename

        # when we preprocess, the output has a different name
        self.cpp_name = None


def find_fortran_targets(fortran_names):
    """read through the Fortran files and look for those marked up with
    attributes(device) or AMREX_DEVICE"""

    targets = []

    for f in fortran_names.split():
        # open the Fortran file
        try:
            fin = open(f, "r")
        except IOError:
            sys.exit("Cannot open Fortran file {}".format(f))

        # loop through the file and look for the target subroutines
        line = fin.readline()
        while line:
            m = fortran_re.search(line)
            if m is not None:
                targets.append(m.group(5).lower().replace("_device",""))
            else:
                m = fortran_attributes_re.search(line)
                if m is not None:
                    targets.append(m.group(5).lower().replace("_device",""))

            line = fin.readline()

    return targets


def convert_headers(outdir, fortran_targets, header_files, cpp):
    """rewrite the C++ headers that contain the Fortran routines"""

    print('looking for targets: {}'.format(fortran_targets))
    print('looking in header files: {}'.format(header_files))

    # first preprocess all the headers and store them in a temporary
    # location.  The preprocessed headers will only be used for the
    # search for the signature, not as the basis for writing the final
    # CUDA header
    pheaders = []

    for h in header_files:
        hdr = "/".join([h[1], h[0]])
        hf = HeaderFile(hdr)

        # preprocess -- this will create a new file in our temp_dir that
        # was run through cpp and has the name CPP-filename
        cpp.preprocess(hf, add_name="CPP")

        pheaders.append(hf)

    # now scan the preprocessed headers and find any of our function
    # signatures and output to a new unpreprocessed header
    for h in pheaders:

        # open the preprocessed header file -- this is what we'll scan
        try:
            hin = open(h.cpp_name, "r")
        except IOError:
            sys.exit("Cannot open header {}".format(h.cpp_name))

        # we'll keep track of the signatures that we need to mangle
        signatures = {}

        line = hin.readline()
        while line:

            # if the line does not start a function signature that
            # matches one of our targets, then we ignore it.
            # Otherwise, we need to capture the function signature
            found = None

            # strip comments
            idx = line.find("//")
            tline = line[:idx]

            for target in fortran_targets:
                fort_target_match = fortran_binding_re.search(tline.lower())
                if fort_target_match:
                    if target == fort_target_match.group(3):
                        found = target
                        print('found target {} in header {}'.format(target, h.cpp_name))
                        break

            # we found a target function, so capture the entire
            # signature, which may span multiple lines
            if found is not None:
                launch_sig = ""
                sig_end = False
                while not line.strip().endswith(";"):
                    launch_sig += line
                    line = hin.readline()
                launch_sig += line

                signatures[found] = launch_sig

            line = hin.readline()

        hin.close()

        # we've now finished going through the header. Note: there may
        # be more signatures here than we really need, because some may
        # have come in via '#includes' in the preprocessing.


        # Now we'll go back to the original file, parse it, making note
        # of any of the signatures we find, but using the preprocessed
        # version in the final output.

        # open the CUDA header for output
        _, tail = os.path.split(h.name)
        ofile = os.path.join(outdir, tail)
        try:
            hout = open(ofile, "w")
        except IOError:
            sys.exit("Cannot open output file {}".format(ofile))

        # and back to the original file (not preprocessed) for the input
        try:
            hin = open(h.name, "r")
        except IOError:
            sys.exit("Cannot open output file {}".format(ofile))

        signatures_needed = []

        line = hin.readline()
        while line:

            # if the line does not start a function signature that we
            # need, then we ignore it
            found = None

            # strip comments
            idx = line.find("//")
            tline = line[:idx]

            # if the line is not a function signature that we already
            # captured then we just write it out
            for target in signatures:

                fort_target_match = fortran_binding_re.search(tline.lower())
                if fort_target_match:
                    if target == fort_target_match.group(3):
                        found = target
                        signatures_needed.append(found)

                        print('found target {} in unprocessed header {}'.format(target, h.name))
                        break

            if found is not None:
                launch_sig = "" + line
                sig_end = False
                while not sig_end:
                    line = hin.readline()
                    launch_sig += line
                    if line.strip().endswith(";"):
                        sig_end = True

            else:
                # this was not one of our device headers
                hout.write(line)

            line = hin.readline()

        # we are done with the pass through the header and we know all
        # of the signatures that need to be CUDAed

        # remove any dupes in the signatures needed
        signatures_needed = list(set(signatures_needed))

        # now do the CUDA signatures
        print('signatures needed: {}'.format(signatures_needed))

        hout.write("\n")
        hout.write("#include <AMReX_ArrayLim.H>\n")
        hout.write("#include <AMReX_BLFort.H>\n")
        hout.write("#include <AMReX_Device.H>\n")
        hout.write("\n")

        hdrmh = os.path.basename(h.name).strip(".H")

        # Add an include guard -- do we still need this?
        hout.write("#ifndef _cuda_" + hdrmh + "_\n")
        hout.write("#define _cuda_" + hdrmh + "_\n\n")

        # Wrap the device declarations in extern "C"
        hout.write("#ifdef AMREX_USE_CUDA\n")
        hout.write("extern \"C\" {\n\n")

        for name in signatures_needed:

            func_sig = signatures[name]

            # First write out the device signature
            device_sig = "__device__ {};\n\n".format(func_sig)

            idx = func_sig.lower().find(name)

            # here's the case-sensitive name
            case_name = func_sig[idx:idx+len(name)]

            # Add _device to the function name.

            device_sig = device_sig.replace(case_name, case_name + "_device")

            device_sig = device_sig.replace("AMREX_ARLIM_VAL", "AMREX_ARLIM_REP")
            hout.write(device_sig)

            # Now write out the global signature. This involves
            # getting rid of the data type definitions and also
            # replacing the lo and hi (which must be in the function
            # definition) with blo and bhi.
            dd = decls_re.search(func_sig)
            vars = []

            has_lo = False
            has_hi = False

            for n, v in enumerate(dd.group(3).split(",")):

                # we will assume that our function signatures _always_ include
                # the name of the variable
                _tmp = v.split()
                var = _tmp[-1].replace("*", "").replace("&", "").strip()

                # Replace AMReX Fortran macros
                var = var.replace("BL_FORT_FAB_ARG_3D", "BL_FORT_FAB_VAL_3D")
                var = var.replace("BL_FORT_IFAB_ARG_3D", "BL_FORT_FAB_VAL_3D")

                if var == "AMREX_ARLIM_VAL(lo)":
                    var = "blo"
                    has_lo = True
                elif var == "lo":
                    var = "blo"
                    has_lo = True
                    func_sig = func_sig.replace("const int* lo", "AMREX_ARLIM_VAL(lo)")

                if var == "AMREX_ARLIM_VAL(hi)":
                    var = "bhi"
                    has_hi = True
                elif var == "hi":
                    var = "bhi"
                    has_hi = True
                    func_sig = func_sig.replace("const int* hi", "AMREX_ARLIM_VAL(hi)")

                vars.append(var)

            if not has_lo or not has_hi:
                sys.exit("ERROR: function signature must have variables lo and hi defined:\n--- function name:\n {} \n--- function signature:\n {}\n---".format(name,func_sig))

            # reassemble the function sig
            all_vars = ", ".join(vars)
            new_call = "{}({})".format(case_name + "_device", all_vars)


            hout.write(TEMPLATE.format(func_sig[idx:].replace(';',''), new_call))
            hout.write("\n")


        # Close out the extern "C" region
        hout.write("\n}\n")
        hout.write("#endif\n")

        # Close out the include guard
        hout.write("\n")
        hout.write("#endif\n")

        hin.close()
        hout.close()


def convert_cxx(outdir, cxx_files):
    """look through the C++ files for "#pragma gpu" and switch it
    to the appropriate CUDA launch macro"""

    print('looking in C++ files: {}'.format(cxx_files))

    for c in cxx_files:
        cxx = "/".join([c[1], c[0]])

        # open the original C++ file
        try:
            hin = open(cxx, "r")
        except IOError:
            sys.exit("Cannot open header {}".format(cxx))

        # open the C++ file for output
        _, tail = os.path.split(cxx)
        ofile = os.path.join(outdir, tail)
        try:
            hout = open(ofile, "w")
        except IOError:
            sys.exit("Cannot open output file {}".format(ofile))

        # look for the appropriate pragma, and once found, capture the
        # function call following it
        line = hin.readline()
        while line:

            # if the line starts with "#pragma gpu", then we need
            # to take action
            if line.startswith("#pragma gpu"):
                # we don't need to reproduce the pragma line in the
                # output, but we need to capture the whole function
                # call that follows
                func_call = ""
                line = hin.readline()
                while not line.strip().endswith(";"):
                    func_call += line
                    line = hin.readline()
                # last line -- remove the semi-colon
                func_call += line.rstrip()[:-1]

                # now split it into the function name and the
                # arguments
                print(func_call)
                dd = decls_re.search(func_call)
                func_name = dd.group(1).strip().replace(" ", "")
                args = dd.group(3)

                # finally output the code in the form we want, with
                # the device launch
                hout.write("dim3 {}numBlocks, {}numThreads;\n" \
                            "Device::grid_stride_threads_and_blocks({}numBlocks, {}numThreads);\n" \
                            "#if ((__CUDACC_VER_MAJOR__ > 9) || (__CUDACC_VER_MAJOR__ == 9 && __CUDACC_VER_MINOR__ >= 1))\n" \
                            "CudaAPICheck(cudaFuncSetAttribute(&cuda_{}, cudaFuncAttributePreferredSharedMemoryCarveout, 0));\n" \
                            "#endif\n" \
                            "cuda_{}<<<{}numBlocks, {}numThreads, 0, Device::cudaStream()>>>\n    ({});\n".format(
                                func_name, func_name, func_name, func_name, func_name, func_name, func_name, func_name, args))

            else:
                # we didn't find a pragma
                hout.write(line)

            line = hin.readline()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--vpath",
                        help="the VPATH to search for files")
    parser.add_argument("--fortran",
                        help="the names of the fortran files to search")
    parser.add_argument("--headers",
                        help="the names of the header files to convert")
    parser.add_argument("--cxx",
                        help="the names of the C++ files to process pragmas")
    parser.add_argument("--output_dir",
                        help="where to write the new header files",
                        default="")
    parser.add_argument("--cpp",
                        help="command to run C preprocessor on .F90 files first.  If omitted, then no preprocessing is done",
                        default="")
    parser.add_argument("--defines",
                        help="defines to send to preprocess the files",
                        default="")
    parser.add_argument("--exclude_defines",
                        help="space separated string of directives to remove from defines",
                        default="")


    args = parser.parse_args()

    defines = args.defines

    if args.exclude_defines != "":
        excludes = args.exclude_defines.split()
        for ex in excludes:
            defines = defines.replace(ex, "")

    print("defines: ", defines)

    if args.cpp != "":
        cpp_pass = preprocess.Preprocessor(temp_dir=args.output_dir, cpp_cmd=args.cpp,
                                           defines=defines)
    else:
        cpp_pass = None

    # part I: for each Fortran routine marked with
    # !$gpu, we need to append a new header in the
    # corresponding *_F.H file

    # find the names of the Fortran subroutines that are marked as
    # device
    targets = find_fortran_targets(args.fortran)

    # find the location of the headers
    headers, _ = ffv.find_files(args.vpath, args.headers)

    # copy the headers to the output directory, replacing the
    # signatures of the target Fortran routines with the CUDA pair
    convert_headers(args.output_dir, targets, headers, cpp_pass)


    # part II: for each C++ file, we need to expand the `#pragma gpu`
    cxx, _ = ffv.find_files(args.vpath, args.cxx)

    convert_cxx(args.output_dir, cxx)
