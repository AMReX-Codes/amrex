# Copyright 2018-2020 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import os, shutil, datetime
import argparse, re, time
from functions_perftest import *

# This script runs automated performance tests for WarpX.
# It runs tests in list test_list defined below, and write
# results in file performance_log.txt in warpx/performance_tests/

# ---- User's manual ----
# Before running performance tests, make sure you have the latest version
# of performance_log.txt

# ---- Running a custom set of performance tests ----
# > python run_alltests_1node.py --no-recompile --compiler=intel
# > --architecture=knl --mode=run --input_file=uniform_plasma
# > --n_node=1 --log_file='my_performance_log.txt'

# ---- Running the pre-drefined automated tests ----
# Compile and run:
# > python run_alltests_1node.py --automated --recompile
# Just run:
# > python run_alltests_1node.py --automated

# To add a new test item, extent the test_list with a line like
# test_list.extend([['my_input_file', n_node, n_mpi, n_omp]]*n_repeat)
# - my_input_file must be in warpx/performance_tests

# ---- Developer's manual ----
# This script can run in two modes:
# - 'run' mode: for each test item, a batch job is executed.
#     create folder '$SCRATCH/performance_warpx/'
#     recompile the code if option --recompile is used
#     loop over test_list and submit one batch script per item
#     Submit a batch job that executes the script in read mode
#     This last job runs once all others are completed
# - 'read' mode: Get performance data from all test items
#     create performance log file if does not exist
#     loop over test_file
#         read initialization time and step time
#         write data into the performance log file
#         push file performance_log.txt on the repo

# Read command-line arguments
# ---------------------------
# Create parser and read arguments
parser = argparse.ArgumentParser(
    description='Run performance tests and write results in files')
parser.add_argument('--recompile', dest='recompile', action='store_true', default=False)
parser.add_argument('--no-recompile', dest='recompile', action='store_false', default=False)
parser.add_argument('--commit', dest='commit', action='store_true', default=False)
parser.add_argument( '--compiler', choices=['gnu', 'intel'], default='intel',
    help='which compiler to use')
parser.add_argument( '--architecture', choices=['cpu', 'knl'], default='knl',
    help='which architecture to cross-compile for NERSC machines')
parser.add_argument( '--mode', choices=['run', 'read'], default='run',
    help='whether to run perftests or read their perf output. run calls read')
parser.add_argument( '--log_file', dest = 'log_file', default='my_performance_log.txt',
    help='name of log file where data will be written. ignored if option --commit is used')
parser.add_argument('--n_node', dest='n_node', default=1, help='nomber of nodes for the runs')
parser.add_argument('--input_file', dest='input_file', default='input_file.pixr',
    type=str, help='input file to run')
parser.add_argument('--automated', dest='automated', action='store_true', default=False,
                    help='Use to run the automated test list')

args = parser.parse_args()
log_file = args.log_file
do_commit = args.commit
run_name = args.input_file

# list of tests to run and analyse.
# Note: This is overwritten if option --automated is used
# each element of test_list contains
# [str input_file, int n_node, int n_mpi PER NODE, int n_omp]
test_list = []
n_repeat = 2
filename1 = args.input_file
test_list.extend([[filename1, 1, 128, 1]]*n_repeat)
test_list.extend([[filename1, 1, 64, 2]]*n_repeat)

# Nothing should be changed after this line
# if flag --automated is used, test_list and do_commit are
# overwritten

if args.automated == True:
    test_list = []
    n_repeat = 2
    test_list.extend([['automated_test_1_uniform_rest_32ppc', 1, 16, 8]]*n_repeat)
    test_list.extend([['automated_test_2_uniform_rest_1ppc',  1, 16, 8]]*n_repeat)
    test_list.extend([['automated_test_3_uniform_drift_4ppc', 1, 16, 8]]*n_repeat)
    test_list.extend([['automated_test_4_labdiags_2ppc',      1, 16, 8]]*n_repeat)
    test_list.extend([['automated_test_5_loadimbalance',      1, 16, 8]]*n_repeat)
    test_list.extend([['automated_test_6_output_2ppc',        1, 16, 8]]*n_repeat)
    do_commit = False
    run_name = 'automated_tests'

n_tests   = len(test_list)
if do_commit == True:
    log_file = 'performance_log.txt'

# Dictionaries
# compiler names. Used for WarpX executable name
compiler_name = {'intel': 'intel', 'gnu': 'gcc'}
# architecture. Used for WarpX executable name
module_name = {'cpu': 'haswell', 'knl': 'mic-knl'}
# architecture. Used in batch scripts
module_Cname = {'cpu': 'haswell', 'knl': 'knl,quad,cache'}
# Define environment variables
cwd = os.getcwd() + '/'
res_dir_base = os.environ['SCRATCH'] + '/performance_warpx/'
bin_dir = cwd + 'Bin/'
bin_name = 'perf_tests3d.' + args.compiler + '.' + module_name[args.architecture] + '.TPROF.MPI.OMP.ex'
log_dir  = cwd

day = time.strftime('%d')
month = time.strftime('%m')
year = time.strftime('%Y')
n_node   = int(args.n_node)

perf_database_file = cwd + 'perf_database_warpx.h5'

# Initialize tests
# ----------------
if args.mode == 'run':
# Set default options for compilation and execution
    config_command = ''
    config_command += 'module unload darshan;'
    config_command += 'module load craype-hugepages4M;'
    if args.architecture == 'knl':
        if args.compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif args.compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-haswell;'
        config_command += 'module load craype-mic-knl;'
    elif args.architecture == 'cpu':
        if args.compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif args.compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-mic-knl;'
        config_command += 'module load craype-haswell;'
    # Create main result directory if does not exist
    if not os.path.exists(res_dir_base):
        os.mkdir(res_dir_base)

# Recompile if requested
if args.recompile == True:
    with open(cwd + 'GNUmakefile_perftest') as makefile_handler:
        makefile_text = makefile_handler.read()
    makefile_text = re.sub('\nCOMP.*', '\nCOMP=%s' %compiler_name[args.compiler], makefile_text)
    with open(cwd + 'GNUmakefile_perftest', 'w') as makefile_handler:
        makefile_handler.write( makefile_text )
    os.system(config_command + " make -f GNUmakefile_perftest realclean ; " + " rm -r tmp_build_dir *.mod; make -j 8 -f GNUmakefile_perftest")

# This function runs a batch script with dependencies to perform the analysis
# when performance runs are done.
def process_analysis():
    dependencies = ''
    f_log = open(cwd + 'log_jobids_tmp.txt','r')
    line = f_log.readline()
    print(line)
    dependencies += line.split()[3] + ':'
    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#SBATCH --job-name=warpx_1node_read\n'
    batch_string += '#SBATCH --time=00:05:00\n'
    batch_string += '#SBATCH -C haswell\n'
    batch_string += '#SBATCH -N 1\n'
    batch_string += '#SBATCH -S 4\n'
    batch_string += '#SBATCH -q regular\n'
    batch_string += '#SBATCH -e read_error.txt\n'
    batch_string += '#SBATCH -o read_output.txt\n'
    batch_string += '#SBATCH --mail-type=end\n'
    batch_string += '#SBATCH --account=m2852\n'
    batch_string += 'python ' + __file__ + ' --no-recompile --compiler=' + \
                    args.compiler + ' --architecture=' + args.architecture + \
                    ' --mode=read' + ' --log_file=' + log_file + \
                ' --input_file=' + args.input_file
    if do_commit == True:
        batch_string += ' --commit'
    if args.automated == True:
        batch_string += ' --automated'
    batch_string += '\n'
    batch_file = 'slurm_perfread'
    f_exe = open(batch_file,'w')
    f_exe.write(batch_string)
    f_exe.close()
    os.system('chmod 700 ' + batch_file)
    os.system('sbatch  --dependency afterok:' + dependencies[0:-1] + ' ' + batch_file)
    return 0

# Loop over the tests and return run time + details
# -------------------------------------------------
if args.mode == 'run':
    # Remove file log_jobids_tmp.txt if exists.
    # This file contains the jobid of every perf test
    # It is used to manage the analysis script dependencies
    if os.path.isfile(cwd + 'log_jobids_tmp.txt'):
        os.remove(cwd + 'log_jobids_tmp.txt')
    res_dir = res_dir_base
    res_dir += '_'.join([run_name, args.compiler,\
                         args.architecture, str(n_node)]) + '/'
    # Run the simulation.
    run_batch_nnode(test_list, res_dir, bin_name, config_command,\
                    architecture=args.architecture, Cname=module_Cname[args.architecture], \
                    n_node=n_node)
    os.chdir(cwd)
    process_analysis()

if args.mode == 'read':
    # Create log_file for performance tests if does not exist
    if not os.path.isfile(log_dir + log_file):
        log_line = '## year month day input_file compiler architecture n_node n_mpi ' +\
                   'n_omp time_initialization time_one_iteration Redistribute '+\
                   'FillBoundary ParallelCopy CurrentDeposition FieldGather '+\
                   'ParthiclePush Copy EvolveEM Checkpoint '+\
                   'WriteParticles Write_FabArray '+\
                   'WriteMultiLevelPlotfile(unit: second) '+\
                   'RedistributeMPI\n'
        f_log = open(log_dir + log_file, 'a')
        f_log.write(log_line)
        f_log.close()
    for count, current_run in enumerate(test_list):
        # Results folder
        print('read ' + str(current_run))
        input_file = current_run[0]
        # Do not read n_node = current_run[1], it is an external parameter
        n_mpi    = current_run[2]
        n_omp    = current_run[3]
        n_steps  = get_nsteps(cwd  + input_file)
        print('n_steps = ' + str(n_steps))
        res_dir = res_dir_base
        res_dir += '_'.join([run_name, args.compiler,\
                             args.architecture, str(n_node)]) + '/'
        # Read performance data from the output file
        output_filename = 'out_' + '_'.join([input_file, str(n_node), str(n_mpi), str(n_omp), str(count)]) + '.txt'
        timing_list = read_run_perf(res_dir + output_filename, n_steps)
        # Write performance data to the performance log file
        log_line = ' '.join([year, month, day, input_file, args.compiler,\
                             args.architecture, str(n_node), str(n_mpi),\
                             str(n_omp)] +  timing_list + ['\n'])
        write_perf_logfile(log_dir + log_file, log_line)
        # Read data for all test to put in hdf5 a database
        # ------------------------------------------------
        # This is an hdf5 file containing ALL the simulation parameters and results. Might be too large for a repo
        df_newline = extract_dataframe(res_dir + output_filename, n_steps)
        # Add all simulation parameters to the dataframe
        df_newline['run_name'] = run_name
        df_newline['n_node'] = n_node
        df_newline['n_mpi'] = n_mpi
        df_newline['n_omp'] = n_omp
        df_newline['n_steps'] = n_steps
        df_newline['rep'] = count
        df_newline['date'] = datetime.datetime.now()
        input_file_open = open(cwd + input_file, 'r')
        input_file_content = input_file_open.read()
        input_file_open.close()
        df_newline['inputs_content'] = input_file_content
        if os.path.exists(perf_database_file):
            df_base = pd.read_hdf(perf_database_file, 'all_data')
            updated_df = df_base.append(df_newline, ignore_index=True)
        else:
            updated_df = df_newline
        updated_df.to_hdf(perf_database_file, key='all_data', mode='w')

    # Store test parameters fot record
    dir_record_base = './perf_warpx_record/'
    if not os.path.exists(dir_record_base):
        os.mkdir(dir_record_base)
    count = 0
    dir_record = dir_record_base + '_'.join([year, month, day]) + '_0'
    while os.path.exists(dir_record):
        count += 1
        dir_record = dir_record[:-1] + str(count)
    os.mkdir(dir_record)
    shutil.copy(__file__, dir_record)
    shutil.copy(log_dir + log_file, dir_record)
    for count, current_run in enumerate(test_list):
        shutil.copy(current_run[0], dir_record)

    # Rename directory with precise date for archive purpose
    res_dir_arch = res_dir_base
    res_dir_arch += '_'.join([year, month, day, run_name, args.compiler,\
                              args.architecture, str(n_node)]) + '/'
    os.rename(res_dir, res_dir_arch)

    # Commit results to the Repo
    if do_commit == True:
        os.system('git add ' + log_dir + log_file + ';'\
                  'git commit -m "performance tests";'\
                  'git push -u origin master')

    # Plot file
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    filename0 = 'my_performance_log'
    filename = filename0 + '.txt'
    fontsize = 14
    matplotlib.rcParams.update({'font.size': fontsize})
    nsteps = 100.
    nrepeat = 4
    legends = [ 'n_node', 'n_mpi', 'n_omp', 'time_initialization', 'time_one_iteration', \
                'Redistribute', 'FillBoundary', 'ParallelCopy', 'CurrentDeposition', \
                'FieldGather', 'ParthiclePush', 'Copy', 'EvolveEM', 'Checkpoint', \
                'WriteParticles', 'Write_FabArray', 'WriteMultiLevelPlotfile', \
                'RedistributeMPI']
    date = np.loadtxt( filename, usecols = np.arange(0, 3 ))
    data = np.loadtxt( filename, usecols = np.arange(6, 6+len(legends)) )
    # Read run name
    with open(filename) as f:
        namelist_tmp = zip(*[line.split() for line in f])[3]
        # Remove first line = comments
        namelist = list(namelist_tmp[1:])
    selector_list = ['automated_test_1_uniform_rest_32ppc',\
                     'automated_test_2_uniform_rest_1ppc',\
                     'automated_test_3_uniform_drift_4ppc',\
                     'automated_test_4_labdiags_2ppc',\
                     'automated_test_5_loadimbalance',\
                     'automated_test_6_output_2ppc']
    selector_string = selector_list[0]
    selector = [idx for idx in range(len(namelist)) if selector_string in namelist[idx]]
    lin_date = date[:,0]+date[:,1]/12.+date[:,2]/366.
    unique_lin_date = np.unique(lin_date)
    my_xticks = unique_lin_date
#     cmap = plt.get_cmap("tab20")
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for selector_string in selector_list:
        selector = [idx for idx in range(len(namelist)) if selector_string in namelist[idx]]
        plt.figure(num=0, figsize=(8,4))
        plt.clf()
        plt.title('warpx ' + selector_string)
        for i in np.arange(data.shape[1]):
            icolors = i-3
            if i>3 and (data[selector,i] > 5./100*data[selector,4]).any():
                plt.plot(lin_date[selector], data[selector,i],'+', ms=6, \
                         mew=2, label=legends[i] )
                # plt.plot(lin_date[selector], data[selector,i],'+', ms=6, \
                #          mew=2, label=legends[i], color=cmap(i) )
            plt.xlabel('date')
            plt.ylabel('time/step (s)')
            plt.grid()
            plt.legend(loc='best')
            plt.legend(bbox_to_anchor=(1.1, 1.05))
            plt.savefig( selector_string + '.pdf', bbox_inches='tight')
            plt.savefig( selector_string + '.png', bbox_inches='tight')
