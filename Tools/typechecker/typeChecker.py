import os
import sys
from pycparser import parse_file, c_parser, c_ast, c_generator


#we will add more types into the dictionaries later
c_primitive_types= ["char", "short", "int", "long", "float", "double"]
fortran_primitive_types= ["CHARACTER", "INTEGER", "LOGICAL", "REAL", "COMPLEX", "DOUBLE PRECISION"]

def may_type_safe(ctype, cdim, ftype, fdim):
  if cdim==0: return False #we can't pass by value 
  if not ctype in c_primitive_types: return True
  if not ftype in fortran_primitive_types: return true
  #now cdim always >=1
  if cdim!=fdim: 
    if ctype=="float" or ctype == "double":return True #this is an exception in Boxlib when floating point data in C is represented as 1D array
    elif cdim==1 and fdim==0: return True #passing by reference
    else: return False
  if ctype== "char": return ftype=="CHARACTER"
  elif ctype=="short": return False
  elif ctype=="int": return ftype=="INTEGER" or ftype=="LOGICAL"
  elif ctype=="long": return ftype=="INTEGER" or ftype=="LOGICAL"
  elif ctype=="float": return ftype=="REAL" or ftype=="COMPLEX"
  elif ctype=="double": return ftype=="DOUBLE PRECISION" 

def checkFortranArguments(funcName, astNode):
  str= "symtree: '"+funcName+"'"
  for file in os.listdir("."):
    if file.endswith("f90.orig"):
      f = open(file, 'r')
      a_line = f.readline()
      arglist = ""
      state=0
      while a_line:
        if str in a_line:
          if(state==0): #first occurence
            f.readline()
            f.readline()
            a_line = f.readline()
            if not "arglist" in a_line:
              a_line = f.readline()
            arglist= a_line.replace("Formal arglist:", "")
            state = state+1
            str= "procedure name = "+funcName
          else: #state==1 --- we can get the type or arguments
            argarray= arglist.split()
            f_argcount= len(argarray)
            if f_argcount != len(astNode.args.params):
              print "Argument mismatch: %s has %d parameters but received %d arguments" %(funcName, len(astNode.args.params), f_argcount)
              return False
            a_line = f.readline()
            while a_line:
              if "procedure name" in a_line: break
              if("symtree" in a_line):
                argname = a_line.split()[1] 
                argname = argname[1:]
                argname = argname[:-1]
                if(argname in argarray): #this local variable is an argument
                  #find the type of this argument
                  a_line = f.readline()
                  word_array = a_line.split()
                  f_type = word_array[3]
                  f_type = f_type[1:]
                  f.readline() #skip 1 line
                  f_dimension =0
		  a_line = f.readline()
                  if "Array spec" in a_line: 
		    word_array= a_line.split()[1]
		    f_dimension= word_array[-1:] 
                  c_dimension=0
                  c_type =astNode.args.params[argarray.index(argname)].type
                  while isinstance(c_type,c_ast.PtrDecl)==True:
                    c_dimension = c_dimension+1
                    c_type= c_type.type
                  if may_type_safe(c_type.type.names[0], c_dimension, f_type, int(f_dimension)) == False:
	            print "Possible type mismatch when calling routine %s with argument number %d" %(funcName, argarray.index(argname))
                  continue
              a_line = f.readline()
            state= state+1
        if state==2: break
        a_line = f.readline()
  return True

class FuncDeclVisitor(c_ast.NodeVisitor):
  def visit_FuncDecl(self, node):
    checkFortranArguments(node.type.declname, node)

def preprocessCHeaders():
  for file in os.listdir("."):
    if file.endswith(".h"):
      f = open(file, 'rw')
      if os.path.exists(file[:-1]+"i"):
        print("Error")
        break
      else:
        f1 = open(file[:-1]+"i", 'w')
      line = f.readline()
      while line:
        line= line.replace("extern \"C\"","\nextern \"C\"\n")
        line= line.replace("{","\n{\n")
        line= line.replace("}","\n}\n")
        f1.write(line)
        line = f.readline()
  
#create a new file without namespace, extern "C" {} so pycparser can process
os.chdir(sys.argv[1])
preprocessCHeaders()
for file in os.listdir("."):
  if file.endswith(".i"):
    f = open(file, 'r')
    if os.path.exists(file[:-1]+"j"):
      print("Error")
      break
    else:
      f1 = open(file[:-1]+"j", 'w')
    a_line = f.readline()
    state_stack = []
    while a_line:
      if "extern \"C\"" in a_line or "namespace" in a_line:
        state_stack.append(1) #remove extern "C" and namespace
      elif "{" in a_line:
        if(len(state_stack) >0): state= state_stack.pop()
        else: state=0
        if state==1: #remove this bracket
          state_stack.append(1) 
          state_stack.append(2) 
        else: 
          state_stack.append(0) 
          f1.write(a_line)
      elif "}" in a_line:
        state= state_stack.pop()
        if state==2: 
          state_stack.pop()
        else: 
          f1.write(a_line)
      else: 
        if(a_line[0] != '\n'): 
          if(len(state_stack) >0): state= state_stack.pop()
          else: state=-1
          if state==1: pass  #extern C isn't followed by {
          elif state!=-1: state_stack.append(state)
          replaced_line= a_line.replace('&','')  #remove all & occurences
          f1.write(replaced_line)
      a_line = f.readline()
    f.close()
    f1.close()
    os.remove(file)

for file in os.listdir("."):
  if file.endswith(".j"):
    f = open(file, 'r')
    #print ("Processing file "+file)
    ast = parse_file(file, use_cpp=True)
    v = FuncDeclVisitor()
    v.visit(ast)
    f.close()
    os.remove(file)
  


