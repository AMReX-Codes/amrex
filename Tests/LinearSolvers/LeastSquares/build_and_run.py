import argparse
import sys
import os
import glob

parser = argparse.ArgumentParser(description='Build and run the Least Squares tests.')
parser.add_argument('-d', '--dim', type=int, default=3, help='Dimensionality')
parser.add_argument('-c', '--clean', action='store_true', default=False, help='Clean build')
parser.add_argument('-f', '--fcompare', action='store_true', default=False, help='Run fcompare on the plot files')
args = parser.parse_args()

ALL_TESTS = []

if args.dim == 2:
   ALL_TESTS = glob.glob('inputs.2d.*')
   EXE = './main2d.gnu.DEBUG.MPI.ex'
   MAKE_CMD = 'make DIM=2 DEBUG=TRUE -j10'
elif args.dim == 3:
   ALL_TESTS = glob.glob('inputs.3d.*')
   EXE = './main3d.gnu.DEBUG.MPI.ex'
   MAKE_CMD = 'make DIM=3 DEBUG=TRUE -j10'

if args.clean:
   print('Clearing plot files and building clean')
   os.system('rm -r plot*')
   os.system('make clean')

print("Building the test executable")
err = os.system(MAKE_CMD);
if not err == 0:
   sys.exit(1)

for t in ALL_TESTS:
   print('################################################################')
   print('Running ' + t)
   err = os.system(EXE + ' ' + t);
   if not err == 0:
      sys.exit(1)
   print('Done ' + t)

if args.fcompare:
   all_pfs = glob.glob('plot*')
   numerical_pfs = []
   analytic_pfs = []

   for f in all_pfs:
      if('analytic') in f:
         s = f.replace('-analytic', '')
         analytic_pfs.append(s)
      else:
         numerical_pfs.append(f)

   for f in numerical_pfs:
      if f in analytic_pfs:
         print('################################################################')
         print('Running fcompare for ' + f)
         err = os.system('fcompare ' + f + ' ' + f + '-analytic')
         print('Done fcompare for ' + f)
