#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse
from pycparser import parse_file, c_ast

error_found = False

def typechecker(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--output", required=True)
    script_args = parser.parse_args()

    fout = open(script_args.output,'w')

    class FuncDeclVisitor(c_ast.NodeVisitor):
        def visit_FuncDecl(self, node):
            check_doit(node, script_args.workdir, fout)

    for f in os.listdir(script_args.workdir):
        if f.endswith('-cppd.h'):
            fname = os.path.join(script_args.workdir,f)
            fout.write("\nChecking "+fname+"...\n")
            ast = parse_file(fname)
            v = FuncDeclVisitor()
            v.visit(ast)

    fout.close()

    global error_found
    if error_found:
        print("\nError found")
    else:
        print("\nNo error found")

def check_doit(node, workdir, fout):
    c_funcname = node.type.declname.rstrip('_')
    f_ret_type, f_arg_type = getFortranArg(c_funcname, workdir)
    if f_ret_type:
        # found a fortran function with that name
        c_ret_type = c_ast_get_type(node.type)
        c_arg_type = []
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

        error_msg.append("\nFunction "+c_funcname+": C vs. Fortran\n")
        fout.write(error_msg[-1])

        if c_to_f_type(c_ret_type) == f_ret_type:
            fout.write("    C return type {0} matches F {1}.\n"
                       .format(c_ret_type,f_ret_type))
        else:
            error_msg.append("    C return type {0} does NOT match F {1}.\n"
                             .format(c_ret_type,f_ret_type))
            fout.write(error_msg[-1])
                
        if len(c_arg_type) == len(f_arg_type):
            fout.write("    number of arguments {0} matches {1}.\n"
                       .format(len(c_arg_type),len(f_arg_type)))
        else:
            error_msg.append("    number of arguments {0} does NOT matche {1}.\n"
                             .format(len(c_arg_type),len(f_arg_type)))
            fout.write(error_msg[-1])

        for ct,ft in zip(c_arg_type, f_arg_type):
            if c_to_f_type(ct[0]) == ft[0] and ct[1] == ft[1]:
                fout.write("    C arg type {0} matches F arg type {1}.\n"
                           .format(ct, ft))
            else:
                error_msg.append("    C arg type {0} does NOT match F arg type {1}.\n"
                                 .format(ct, ft))
                fout.write(error_msg[-1])

        global error_found
        if len(error_msg) > 1:
            error_found = True
            for msg in error_msg:
                print(msg,end='')


def c_to_f_type(ctyp):
    # supported C types: char, short, int, long, float, double, void
    if ctyp == 'char':
        return 'INTEGER 1'
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


def c_ast_get_type(decl):
    typ = type(decl)
    if typ == c_ast.IdentifierType:
        return ' '.join(decl.names)
    elif typ == c_ast.PtrDecl:
        return c_ast_get_type(decl.type) + ' pointer'
    else:
        return c_ast_get_type(decl.type)
    
def getFortranArg(funcname, workdir):
    """ given a function name and a directory for searching, return a tuple
        of function return type and list of arguments.  For example,
        (INTEGER 4, [('INTERGER 4', 'value'), ('REAL 8', 'pointer')])."""

    debug = False

    module_tok = "symtree: '" + funcname + "'"
    proc_tok = "procedure name = " + funcname + "\n"
    this_func = []
    for fname in os.listdir(workdir):
        if fname.endswith('.orig'):
            f = open(os.path.join(workdir,fname), 'r')
            # let's collect the text of this function into a list of strings
            in_module_block = False
            ws_module_block = 0
            in_proced_block = False
            ws_proced_block = 0
            for line in f.readlines():
                num_white_spaces = len(line) - len(line.lstrip())
                if module_tok in line and not in_proced_block:
                    this_func.append(line)
                    in_module_block = True
                    ws_module_block = num_white_spaces
                    in_proced_block = False
                    ws_proced_block = 0
                    found = True
                elif proc_tok in line and not in_module_block:
                    this_func.append(line)
                    in_proced_block = True
                    ws_proced_block = num_white_spaces
                    in_module_block = False
                    ws_module_block = 0
                    found = True
                elif in_module_block:
                    if ws_module_block < num_white_spaces:
                        this_func.append(line)
                    else:
                        in_module_block = False
                        ws_module_block = 0
                elif in_proced_block:
                    if ws_proced_block < num_white_spaces:
                        this_func.append(line)
                    else:
                        in_proced_block = False
                        ws_proced_block = 0
            f.close()            
            if this_func:
                break

    if debug:
        print(funcname, "this_func...")
        print(this_func)
            
    return_type = ''
    arguments_type = []
    found = False
    ws = 0
    arglist = []
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
                    if "UNKNOWN" in line:
                        return_type = 'void'
                    else:
                        return_type = line.split('(')[1].split(')')[0]
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
                    argtype = line.split('(')[1].split(')')[0]
                elif "attributes:" in line: # attributes: (VARIABLE  VALUE DUMMY(IN))
                    if "VALUE" in line:
                        argattrib = 'value'
                    else:
                        argattrib = 'pointer'
                    found = False
                    func_args[argname] = (argtype, argattrib)
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
