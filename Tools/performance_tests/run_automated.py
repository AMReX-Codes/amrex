#!/usr/common/software/python/2.7-anaconda-4.4/bin/python

import os, sys, shutil, datetime
import argparse, re, time, copy
import pandas as pd
from functions_perftest import *

# typical use: python run_automated.py --n_node_list='1,8,16,32' --automated
# Assume warpx, picsar, amrex and perf_logs repos ar in the same directory and
# environment variable AUTOMATED_PERF_TESTS contains the path to this directory

# Handle parser
###############
parser = argparse.ArgumentParser( description='Run performance tests and write results in files' )
parser.add_argument('--recompile',
                    dest='recompile',
                    action='store_true',
                    default=False)
parser.add_argument('--commit',
                    dest='commit',
                    action='store_true',
                    default=False)
parser.add_argument('--automated',
                    dest='automated',
                    action='store_true',
                    default=False, 
                    help='Use to run the automated test list')
parser.add_argument('--n_node_list',
                    dest='n_node_list',
                    default=[],
                    help='list ofnumber of nodes for the runs', type=str)
parser.add_argument('--start_date',
                    dest='start_date' )
parser.add_argument('--compiler',
                    choices=['gnu', 'intel'],
                    default='intel', 
                    help='which compiler to use')
parser.add_argument('--architecture',
                    choices=['cpu', 'knl'],
                    default='knl',
                    help='which architecture to cross-compile for NERSC machines')
parser.add_argument('--mode',
                    choices=['run', 'read', 'browse_output_files', 'write_csv'],
                    default='run', 
                    help='whether to run perftests or read their perf output. run calls read')
args = parser.parse_args()
n_node_list_string   = args.n_node_list.split(',')
n_node_list = [int(i) for i in n_node_list_string]
start_date = args.start_date

# Set behavior variables 
########################
write_csv = False
browse_output_files = False
if args.mode == 'write_csv':
    write_csv = True
if args.mode == 'browse_output_files':
    browse_output_file = True
if args.mode == 'read':
    write_csv = True
    browse_output_files = True
recompile = args.recompile
perf_database_file = 'my_tests_database.h5'
if args.automated == True:
    run_name = 'automated_tests'
    perf_database_file = 'automated_tests_database.h5'
    rename_archive = True
    store_full_input = False
    update_perf_log_repo = True
    push_on_perf_log_repo = False
    pull_3_repos = True
    recompile = True

# Each instance of this class contains information for a single test.
class test_element():
    def __init__(self, input_file=None, n_node=None, n_mpi_per_node=None, 
                 n_omp=None, n_cell=None, n_step=None):
        self.input_file = input_file
        self.n_node = n_node
        self.n_mpi_per_node = n_mpi_per_node
        self.n_omp = n_omp
        self.n_cell = n_cell
        self.n_step = n_step

    def scale_n_cell(self, n_node=0):
        n_cell_scaled = copy.deepcopy(self.n_cell)
        index_dim = 0
        while n_node > 1:
            n_cell_scaled[index_dim] *= 2
            n_node /= 2
            index_dim = (index_dim+1) % 3
        self.n_cell = n_cell_scaled

# List of tests to perform
# ------------------------
test_list_unq = []
# Each test runs n_repeat times
n_repeat = 2
# n_node is kept to None and passed in functions as an external argument
# That way, several test_element_instance run with the same n_node on the same batch job
test_list_unq.append( test_element(input_file='automated_test_1_uniform_rest_32ppc', 
                                   n_mpi_per_node=8, 
                                   n_omp=8, 
                                   n_cell=[128, 128, 128], 
                                   n_step=10) )
test_list_unq.append( test_element(input_file='automated_test_2_uniform_rest_1ppc', 
                                   n_mpi_per_node=8, 
                                   n_omp=8, 
                                   n_cell=[256, 256, 512], 
                                   n_step=10) )
test_list_unq.append( test_element(input_file='automated_test_3_uniform_drift_4ppc', 
                                   n_mpi_per_node=8, 
                                   n_omp=8, 
                                   n_cell=[128, 128, 128], 
                                   n_step=10) )
test_list_unq.append( test_element(input_file='automated_test_4_labdiags_2ppc', 
                                   n_mpi_per_node=8, 
                                   n_omp=8, 
                                   n_cell=[64, 64, 128], 
                                   n_step=50) )
test_list_unq.append( test_element(input_file='automated_test_5_loadimbalance', 
                                   n_mpi_per_node=8, 
                                   n_omp=8, 
                                   n_cell=[128, 128, 128], 
                                   n_step=10) )
test_list_unq.append( test_element(input_file='automated_test_6_output_2ppc', 
                                   n_mpi_per_node=8, 
                                   n_omp=8, 
                                   n_cell=[128, 256, 256], 
                                   n_step=0) )
test_list = [copy.deepcopy(item) for item in test_list_unq for _ in range(n_repeat) ]

# Define directories
# ------------------
source_dir_base = os.environ['AUTOMATED_PERF_TESTS']
warpx_dir = source_dir_base + '/WarpX/'
picsar_dir = source_dir_base + '/picsar/'
amrex_dir = source_dir_base + '/amrex/'
res_dir_base = os.environ['SCRATCH'] + '/performance_warpx/'
perf_logs_repo = source_dir_base + 'perf_logs/'

# Define dictionaries
# -------------------
compiler_name = {'intel': 'intel', 'gnu': 'gcc'}
module_name = {'cpu': 'haswell', 'knl': 'mic-knl'}
module_Cname = {'cpu': 'haswell', 'knl': 'knl,quad,cache'}
cwd = os.getcwd() + '/'
bin_dir = cwd + 'Bin/'
bin_name = 'perf_tests3d.' + args.compiler + '.' + module_name[args.architecture] + '.TPROF.MPI.OMP.ex'
log_dir  = cwd
perf_database_file = cwd + perf_database_file
day = time.strftime('%d')
month = time.strftime('%m')
year = time.strftime('%Y')

# Initialize tests
# ----------------
if args.mode == 'run':
    start_date = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
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
    # ----------------------
    if recompile == True:
        if pull_3_repos == True:
            git_repo = git.cmd.Git( picsar_dir )
            git_repo.pull()
            git_repo = git.cmd.Git( amrex_dir  )
            git_repo.pull()
            git_repo = git.cmd.Git( warpx_dir  )
            git_repo.pull()
        with open(cwd + 'GNUmakefile_perftest') as makefile_handler:
            makefile_text = makefile_handler.read()
        makefile_text = re.sub('\nCOMP.*', '\nCOMP=%s' %compiler_name[args.compiler], makefile_text)
        with open(cwd + 'GNUmakefile_perftest', 'w') as makefile_handler:
            makefile_handler.write( makefile_text )
        os.system(config_command + " make -f GNUmakefile_perftest realclean ; " + " rm -r tmp_build_dir *.mod; make -j 8 -f GNUmakefile_perftest")
        if os.path.exists( cwd + 'store_git_hashes.txt' ):
            os.remove( cwd + 'store_git_hashes.txt' )
        store_git_hash(repo_path=picsar_dir, filename=cwd + 'store_git_hashes.txt', name='picsar')
        store_git_hash(repo_path=amrex_dir , filename=cwd + 'store_git_hashes.txt', name='amrex' )
        store_git_hash(repo_path=warpx_dir , filename=cwd + 'store_git_hashes.txt', name='warpx' )

# This function runs a batch script with 
# dependencies to perform the analysis 
# after all performance tests are done.
def process_analysis():
    dependencies = ''
    f_log = open(cwd + 'log_jobids_tmp.txt' ,'r')
    for line in f_log.readlines():
        dependencies += line.split()[3] + ':'
    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#SBATCH --job-name=warpx_1node_read\n'
    batch_string += '#SBATCH --time=00:07:00\n'
    batch_string += '#SBATCH -C knl\n'
    batch_string += '#SBATCH -N 1\n'
    batch_string += '#SBATCH -S 4\n'
    batch_string += '#SBATCH -q regular\n'
    batch_string += '#SBATCH -e read_error.txt\n'
    batch_string += '#SBATCH -o read_output.txt\n'
    batch_string += '#SBATCH --mail-type=end\n'
    batch_string += '#SBATCH --account=m2852\n'
    batch_string += 'python ' + __file__ + ' --compiler=' + \
                    args.compiler + ' --architecture=' + args.architecture + \
                    ' --mode=read' + \
                ' --n_node_list=' + '"' + args.n_node_list + '"' + \
                ' --start_date=' + start_date
    if args.automated == True:
        batch_string += ' --automated'
    batch_string += '\n'
    batch_file = 'slurm_perfread'
    f_exe = open(batch_file,'w')
    f_exe.write(batch_string)
    f_exe.close()
    os.system('chmod 700 ' + batch_file)
    print( 'process_analysis line:  ' + 'sbatch  --dependency afterok:' + dependencies[0:-1] + ' ' + batch_file)
    os.system('sbatch  --dependency afterok:' + dependencies[0:-1] + ' ' + batch_file)
    return 0

# Loop over the tests and run all simulations:
# One batch job submitted per n_node. Several
# tests run within the same batch job. 
# --------------------------------------------
if args.mode == 'run':
    if os.path.exists( 'log_jobids_tmp.txt' ):
        os.remove( 'log_jobids_tmp.txt' )
    # loop on n_node. One batch script per n_node
    for n_node in n_node_list:
        res_dir = res_dir_base
        res_dir += '_'.join([run_name, args.compiler, args.architecture, str(n_node)]) + '/'
        runtime_param_list = []
        # Deep copy as we change the attribute n_cell of
        # each instance of class test_element
        test_list_n_node = copy.deepcopy(test_list)
        # Loop on tests
        for current_run in test_list_n_node:
            current_run.scale_n_cell(n_node)
            runtime_param_string  = ' amr.n_cell=' + ' '.join(str(i) for i in current_run.n_cell)
            runtime_param_string += ' max_step=' + str( current_run.n_step )
            runtime_param_list.append( runtime_param_string )
        # Run the simulations.
        run_batch_nnode(test_list_n_node, res_dir, bin_name, config_command,\
                        architecture=args.architecture, Cname=module_Cname[args.architecture], \
                        n_node=n_node, runtime_param_list=runtime_param_list)
    os.chdir(cwd)
    # submit batch for analysis
    process_analysis()

# read the output file from each test and store timers in
# hdf5 file with pandas format
# -------------------------------------------------------
for n_node in n_node_list:
    print(n_node)
    if browse_output_files:
        for count, current_run in enumerate(test_list):
            res_dir = res_dir_base
            res_dir += '_'.join([run_name, args.compiler,\
                                 args.architecture, str(n_node)]) + '/'
            # Read performance data from the output file
            output_filename = 'out_' + '_'.join([current_run.input_file, str(n_node), str(current_run.n_mpi_per_node), str(current_run.n_omp), str(count)]) + '.txt'
            # Read data for all test to put in hdf5 a database
            # This is an hdf5 file containing ALL the simulation
            # parameters and results. Might be too large for a repo
            df_newline = extract_dataframe(res_dir + output_filename, current_run.n_step)
            # Add all simulation parameters to the dataframe
            df_newline['git_hashes'] = get_file_content(filename=cwd+'store_git_hashes.txt')
            df_newline['start_date'] = start_date
            df_newline['run_name'] = run_name
            df_newline['input_file'] = current_run.input_file
            df_newline['n_node'] = n_node
            df_newline['n_mpi_per_node'] = current_run.n_mpi_per_node
            df_newline['n_omp'] = current_run.n_omp
            df_newline['n_steps'] = current_run.n_step
            df_newline['rep'] = count%n_repeat
            df_newline['date'] = datetime.datetime.now()
            if store_full_input:
                df_newline['inputs_content'] = get_file_content( filename=cwd+current_run.input_file )
            # Load file perf_database_file if exists, and
            # append with results from this scan
            if os.path.exists(perf_database_file):
                df_base = pd.read_hdf(perf_database_file, 'all_data')
                updated_df = df_base.append(df_newline, ignore_index=True)
            else:
                updated_df = df_newline
            # Write dataframe to file perf_database_file 
            # (overwrite if file exists)
            updated_df.to_hdf(perf_database_file, key='all_data', mode='w')
 
        # Rename directory with precise date+hour for archive purpose
        if rename_archive == True:
            loc_counter = 0
            res_dir_arch = res_dir_base
            res_dir_arch += '_'.join([year, month, day, run_name, args.compiler,\
                                      args.architecture, str(n_node), str(loc_counter)]) + '/'
            while os.path.exists( res_dir_arch ):
                loc_counter += 1
                res_dir_arch = res_dir_base
                res_dir_arch += '_'.join([year, month, day, run_name, args.compiler,\
                                          args.architecture, str(n_node), str(loc_counter)]) + '/'
            os.rename( res_dir, res_dir_arch )

# Extract sub-set of pandas data frame, write it to
# csv file and copy this file to perf_logs repo
# -------------------------------------------------
if write_csv:
    # Extract small data from data frame and write them to 
    # First, generate csv files
    df = pd.read_hdf( perf_database_file )
    # One large file
    df.loc[:,'step_time'] = pd.Series(df['time_running']/df['n_steps'], index=df.index)
    # Make smaller dataframe with only data to be written to csv file
    df_small = df.copy()
    df_small.loc[ df_small['input_file']=='automated_test_6_output_2ppc', 'step_time'] = \
        df_small[ df_small['input_file']=='automated_test_6_output_2ppc' ]['time_WritePlotFile']
    df_small = df_small.loc[:, ['date', 'input_file', 'git_hashes', 'n_node', 'n_mpi_per_node', 'n_omp', 'rep', 'start_date', 'time_initialization', 'step_time'] ]
    # Write to csv
    df_small.to_csv( 'cori_knl.csv' )
    # Errors may occur depending on the version of pandas. I had errors with v0.21.0 solved with 0.23.0
    # Second, move files to perf_logs repo
    if update_perf_log_repo:
        git_repo = git.Repo( perf_logs_repo )
        if push_on_perf_log_repo:
            git_repo.git.stash('save')
            git_repo.git.pull()
        shutil.move( 'cori_knl.csv', perf_logs_repo + '/logs_csv/cori_knl.csv' )
        os.chdir( perf_logs_repo )
        sys.path.append('./')
        import generate_index_html
        git_repo.git.add('./index.html')
        git_repo.git.add('./logs_csv/cori_knl.csv')
        index = git_repo.index
        index.commit("automated tests")
