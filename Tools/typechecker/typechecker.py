#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse
from pycparser import parse_file, c_ast

def typechecker(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--output", required=True)
    script_args = parser.parse_args()

    debug = True

    funcs_to_check = []

    # (1) let's collect all c functions we are going to check into a list
    class FuncDeclVisitor(c_ast.NodeVisitor):
        def visit_FuncDecl(self, node):
            funcs_to_check.append(node.type.declname.rstrip('_'))

    for f in os.listdir(script_args.workdir):
        if f.endswith('-cppd.h'):
            fname = os.path.join(script_args.workdir,f)
            if debug:
                print("Checking {0}...".format(f))
            ast = parse_file(fname)
            v = FuncDeclVisitor()
            v.visit(ast)

    # (2) find out which Fortran files have the functions
    func_src = {}
    findFortranSources(funcs_to_check, func_src, script_args.workdir)

    # (3) we visit each function declaration again
    fout = open(script_args.output,'w')

    aux_info = {'numerrors':0, 'numfuncs':0, 'current_c_header':''}

    class FuncDeclVisitor(c_ast.NodeVisitor):
        def visit_FuncDecl(self, node):
            check_doit(node, script_args.workdir, func_src, fout, aux_info)

    for f in os.listdir(script_args.workdir):
        if f.endswith('-cppd.h'):
            fname = os.path.join(script_args.workdir,f)
            fout.write("\nChecking "+fname+"...\n")
            aux_info['current_c_header'] = os.path.basename(fname.replace('-cppd.h','.H'))
            ast = parse_file(fname)
            v = FuncDeclVisitor()
            v.visit(ast)

    fout.close()

    print("{0} functions checked, {1} error(s) found.  More details can be found in {2}."
          .format(aux_info['numfuncs'],aux_info['numerrors'], script_args.output))


def check_doit(node, workdir, func_src, fout, aux_info):
    c_funcname = node.type.declname.rstrip('_')
    if not c_funcname in func_src.keys():
        return
    f_filename = func_src[c_funcname]
    f_ret_type, f_arg_type = getFortranArg(c_funcname, os.path.join(workdir,f_filename))

    if f_ret_type:
        # found a fortran function with that name
        c_ret_type = c_ast_get_type(node.type)
        c_arg_type = []
        if node.args:
            for param in node.args.params:
                typ = c_ast_get_type(param.type).split(' ')
                if len(typ) == 1:
                    c_arg_type.append([typ[0], 'value'])
                elif len(typ) == 2:
                    c_arg_type.append(typ)
                else:
                    print("check_doit: how did this happen? ", typ)
                    sys.exit(1)

        error_msg = []

        error_msg.append("\nFunction "+c_funcname+" in "+aux_info['current_c_header']+
                         " vs. Fortran procedure in "+f_filename.replace(".orig",'')+"\n")
        fout.write(error_msg[-1])

        if c_to_f_type(c_ret_type) == f_ret_type:
            fout.write("    C return type {0} matches Fortran {1}.\n"
                       .format(c_ret_type,f_ret_type))
        else:
            error_msg.append("    C return type {0} does NOT match Fortran {1}.\n"
                             .format(c_ret_type,f_ret_type))
            fout.write(error_msg[-1])
                
        if len(c_arg_type) == len(f_arg_type):
            fout.write("    number of arguments {0} matches {1}.\n"
                       .format(len(c_arg_type),len(f_arg_type)))
        else:
            error_msg.append("    number of arguments {0} does NOT match {1}.\n"
                             .format(len(c_arg_type),len(f_arg_type)))
            fout.write(error_msg[-1])

        cnt = 0
        for ctyp,ftyp in zip(c_arg_type, f_arg_type):
            cnt = cnt+1
            # exact match, or * == type(c_ptr), value, or void* == any Fortran 'reference'
            if c_to_f_type(ctyp[0]) == ftyp[0] and ctyp[1] == ftyp[1] or \
               ctyp[1] == "pointer" and ftyp[0] == "DERIVED c_ptr" and ftyp[1] == "value" or \
               ctyp[0] == "void" and ctyp[1] == "pointer" and ftyp[1] == "pointer":
                fout.write("    arg #{0}: C type {1} matches Fortran type {2}.\n"
                           .format(cnt, ctyp, ftyp))
            else:
                error_msg.append("    arg #{0}: C type {1} does NOT match Fortran type {2}.\n"
                                 .format(cnt, ctyp, ftyp))
                fout.write(error_msg[-1])

        if len(error_msg) > 1:
            aux_info['numerrors'] = aux_info['numerrors']+1
            for msg in error_msg:
                print(msg,end='')
        # note some C declarations may not have Fortran function defined
        # those don't count in numfuncs
        aux_info['numfuncs'] = aux_info['numfuncs']+1


def c_to_f_type(ctyp):
    # supported C types: char, short, int, long, float, double, void, _Bool
    if ctyp == 'char':
        return 'CHARACTER 1 1'
    elif ctyp == 'short':
        return 'INTEGER 2'
    elif ctyp == 'int':
        return 'INTEGER 4'
    elif ctyp == 'long':
        return 'INTEGER 8'  # yes, this is our assumption
    elif ctyp == 'float':
        return 'REAL 4'
    elif ctyp == 'double':
        return 'REAL 8'
    elif ctyp == 'amrex_real':
        return 'REAL 8'
    elif ctyp == 'void':
        return 'void'
    elif ctyp == '_Bool':
        return 'LOGICAL 1'


def c_ast_get_type(decl):
    typ = type(decl)
    if typ == c_ast.IdentifierType:
        #return ' '.join(decl.names)
        return decl.names[0]
    elif typ == c_ast.PtrDecl or typ == c_ast.ArrayDecl:
        return c_ast_get_type(decl.type) + ' pointer'
    else:
        return c_ast_get_type(decl.type)


def findFortranSources(funcs, func_src, workdir):
    for fname in os.listdir(workdir):
        if fname.endswith('.orig'):
            f = open(os.path.join(workdir,fname), 'r')
            for line in f.readlines():
                if "procedure name = " in line:
                    pname = line.split('=')[1].strip()
                    if pname in funcs:
                        if pname in func_src.keys():
                            print("two sources for", pname, "in", func_src[pname],
                                  "and", fname)
                            sys.exit(1)
                        else:
                            func_src[pname] = fname

    
def getFortranArg(funcname, fortranfile):
    """given a function name and a gortran tree original file, return a
        tuple of function return type and list of arguments.  For
        example, (INTEGER 4, [('INTERGER 4', 'value'), ('REAL 8',
        'pointer')]).  Note that pointer here is not Fortran pointer.
    """

    debug = False

    module_tok = "symtree: '" + funcname + "'"
    proc_tok = "procedure name = " + funcname + "\n"
    this_func = []
    f = open(fortranfile, 'r')
    # let's collect the text of this function into a list of strings
    in_module_block = False
    ws_module_block = 0
    in_proced_block = False
    ws_proced_block = 0
    numblocks = 0
    for line in f.readlines():
        if line.isspace(): continue
        num_white_spaces = len(line) - len(line.lstrip())
        if module_tok in line and not in_proced_block:
            this_func.append(line)
            in_module_block = True
            ws_module_block = num_white_spaces
            in_proced_block = False
            ws_proced_block = 0
            numblocks = numblocks + 1
        elif proc_tok in line and not in_module_block:
            this_func.append(line)
            in_proced_block = True
            ws_proced_block = num_white_spaces
            in_module_block = False
            ws_module_block = 0
            numblocks = numblocks + 1
        elif in_module_block:
            if ws_module_block < num_white_spaces:
                this_func.append(line)
            else:
                in_module_block = False
                ws_module_block = 0
                if numblocks == 2:
                    break
        elif in_proced_block:
            if line.strip() != "code:" and ws_proced_block < num_white_spaces:
                this_func.append(line)
            else:
                in_proced_block = False
                ws_proced_block = 0
                if numblocks == 2:
                    break
    f.close()            
    if not this_func:
        print(fortranfile, "doesn't contain", function)

    if debug:
        print(funcname, "this_func...")
        for line in this_func:
            print(line,)
            
    return_type = ''
    arguments_type = []
    found = False
    ws = 0
    arglist = []

    def parse_type_spec(line):
        # The line may look like, type spec : (REAL 8 C_INTEROP)
        #                     or  type spec : (UNKNOWN 0 C_INTEROP)
        #                     or  type spec : (INTEGER 4)
        #                     or  type spec : (CHARACTER 1 1 C_INTEROP ISO_C)
        if "UNKNOWN" in line:
            type_spec = 'void'
        else:
            type_spec = line.split('(')[1].split(')')[0].split('C_INTEROP')[0].strip()
        return type_spec

    for line in this_func:
        num_white_spaces = len(line) - len(line.lstrip())
        # searching for blocks like
        #   symtree: 'fort_fab_copy'|| symbol: 'fort_fab_copy' 
        #     type spec : (UNKNOWN 0)
        #     attributes: (PROCEDURE MODULE-PROC  BIND(C) SUBROUTINE)
        #     Formal arglist: lo hi dst dlo dhi src slo shi sblo ncomp
        if "symtree: '"+funcname+"'" in line: # e.g. 
            found = True
            ws = num_white_spaces
        elif found:
            if ws < num_white_spaces:
                if "type spec :" in line:
                    return_type = parse_type_spec(line)
                elif "Formal arglist:" in line:
                    arglist = line.replace("Formal arglist:", "").split()
            else:
                break

    if arglist:
        func_args = {}
        found = False
        argname = ''
        argtype = ''
        for line in this_func:
            # search for blocks like
            #   symtree: 'dhi'         || symbol: 'dhi'          
            #     type spec : (INTEGER 4)
            #     attributes: (VARIABLE  DIMENSION DUMMY(IN))
            #     Array spec:(1 [0] AS_EXPLICIT 1 3 )
            if "symtree: '" in line:
                argname = line.split("'")[1]
                if argname in arglist:
                    found = True
            elif found:
                if "type spec :" in line:
                    argtype = parse_type_spec(line)
                elif "attributes:" in line: # attributes: (VARIABLE  VALUE DUMMY(IN))
                    if "VALUE" in line:
                        argattrib = 'value'
                    else:
                        argattrib = 'pointer'
                    found = False
                    func_args[argname] = (argtype, argattrib, argname)
                    if len(func_args) == len(arglist):
                        break
        if len(func_args) != len(arglist):
            print("getFortranArg: len != len", funcname, func_args, arglist)
            sys.exit(1)

        if debug:
            print ("-------")
            print(arglist)
            print(func_args)
            print ("-------")
            
        for a in arglist:
            arguments_type.append(func_args[a])
        
    return return_type, arguments_type

if __name__ == "__main__":
    typechecker(sys.argv)
