# Copyright 2018-2019 Axel Huebl, Luca Fedeli, Maxence Thevenet
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import os, sys, shutil, datetime, git
import argparse, time, copy
import pandas as pd
from functions_perftest import store_git_hash, get_file_content, \
    run_batch_nnode, extract_dataframe

# Get name of supercomputer and import configuration functions from
# machine-specific file
if os.getenv("LMOD_SYSTEM_NAME") == 'summit':
    machine = 'summit'
    from summit import executable_name, process_analysis, \
        get_config_command, time_min, get_submit_job_command, \
        get_batch_string, get_run_string, get_test_list
if os.getenv("NERSC_HOST") == 'cori':
    machine = 'cori'
    from cori import executable_name, process_analysis, \
        get_config_command, time_min, get_submit_job_command, \
        get_batch_string, get_run_string, get_test_list

# typical use: python run_automated.py --n_node_list='1,8,16,32' --automated
# Assume warpx, picsar, amrex and perf_logs repos ar in the same directory and
# environment variable AUTOMATED_PERF_TESTS contains the path to this directory

# requirements:
# - python packages: gitpython and pandas
# - AUTOMATED_PERF_TESTS: environment variables where warpx,
#   amrex and picsar are installed ($AUTOMATED_PERF_TESTS/warpx etc.)
# - SCRATCH: environment variable where performance results are written.
#   This script will create folder $SCRATCH/performance_warpx/

if "AUTOMATED_PERF_TESTS" not in os.environ:
    raise ValueError("environment variable AUTOMATED_PERF_TESTS is not defined.\n"
                     "It should contain the path to the directory where WarpX, "
                     "AMReX and PICSAR repos are.")
if "SCRATCH" not in os.environ:
    raise ValueError("environment variable SCRATCH is not defined.\n"
                     "This script will create $SCRATCH/performance_warpx/ "
                     "to store performance results.")
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
                    choices=['gnu', 'intel', 'pgi'],
                    default='intel',
                    help='which compiler to use')
parser.add_argument('--architecture',
                    choices=['cpu', 'knl', 'gpu'],
                    default='knl',
                    help='which architecture to cross-compile for NERSC machines')
parser.add_argument('--mode',
                    choices=['run', 'read', 'browse_output_files'],
                    default='run',
                    help='whether to run perftests or read their perf output. run calls read')
parser.add_argument('--path_source',
                    default=None,
                    help='path to parent folder containing amrex, picsar and warpx folders')
parser.add_argument('--path_results',
                    default=None,
                    help='path to result directory, where simulations run')

args = parser.parse_args()
n_node_list_string   = args.n_node_list.split(',')
n_node_list = [int(i) for i in n_node_list_string]
start_date = args.start_date

# Set behavior variables
########################
run_name = 'custom_perftest'
perf_database_file = 'my_tests_database.h5'
rename_archive = False
store_full_input = False
update_perf_log_repo = False
push_on_perf_log_repo = False
recompile = args.recompile
pull_3_repos = False
recompile = True
compiler = args.compiler
architecture = args.architecture
source_dir_base = args.path_source
res_dir_base = args.path_results

browse_output_files = False
if args.mode == 'browse_output_files':
    browse_output_file = True
if args.mode == 'read':
    browse_output_files = True

if args.automated == True:
    run_name = 'automated_tests'
    perf_database_file = machine + '_results.h5'
    rename_archive = True
    store_full_input = False
    update_perf_log_repo = True
    push_on_perf_log_repo = False
    pull_3_repos = True
    recompile = True
    source_dir_base = os.environ['AUTOMATED_PERF_TESTS']
    res_dir_base = os.environ['SCRATCH'] + '/performance_warpx/'
    if machine == 'summit':
        compiler = 'gnu'
        architecture = 'gpu'

# List of tests to perform
# ------------------------
# Each test runs n_repeat times
n_repeat = 2
# test_list is machine-specific
test_list = get_test_list(n_repeat)

# Define directories
# ------------------
warpx_dir = source_dir_base + '/warpx/'
picsar_dir = source_dir_base + '/picsar/'
amrex_dir = source_dir_base + '/amrex/'
perf_logs_repo = source_dir_base + 'perf_logs/'

# Define dictionaries
# -------------------
compiler_name = {'intel': 'intel', 'gnu': 'gcc', 'pgi':'pgi'}
module_Cname = {'cpu': 'haswell', 'knl': 'knl,quad,cache', 'gpu':''}
csv_file = {'cori':'cori_knl.csv', 'summit':'summit.csv'}
# cwd = os.getcwd() + '/'
cwd = warpx_dir + 'Tools/performance_tests/'

path_hdf5 = cwd
if args.automated:
    path_hdf5 = perf_logs_repo + '/logs_hdf5/'

bin_dir = cwd + 'Bin/'
bin_name = executable_name(compiler, architecture)

log_dir  = cwd
day = time.strftime('%d')
month = time.strftime('%m')
year = time.strftime('%Y')

# Initialize tests
# ----------------
if args.mode == 'run':
    start_date = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    # Set default options for compilation and execution
    config_command = get_config_command(compiler, architecture)
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

        # Copy WarpX/GNUmakefile to current directory and recompile
        # with specific options for automated performance tests.
        # This way, performance test compilation does not mess with user's
        # compilation
        shutil.copyfile("../../GNUmakefile","./GNUmakefile")
        make_realclean_command = " make realclean WARPX_HOME=../.. " \
            "AMREX_HOME=../../../amrex/ PICSAR_HOME=../../../picsar/ " \
            "EBASE=perf_tests COMP=%s" %compiler_name[compiler] + ";"
        make_command = "make -j 16 WARPX_HOME=../.. " \
            "AMREX_HOME=../../../amrex/ PICSAR_HOME=../../../picsar/ " \
            "EBASE=perf_tests COMP=%s" %compiler_name[compiler]
        if machine == 'summit':
            make_command += ' USE_GPU=TRUE '
        os.system(config_command + make_realclean_command + \
                  "rm -r tmp_build_dir *.mod; " + make_command )

        # Store git hashes for WarpX, AMReX and PICSAR into file, so that
        # they can be read when running the analysis.
        if os.path.exists( cwd + 'store_git_hashes.txt' ):
            os.remove( cwd + 'store_git_hashes.txt' )
        store_git_hash(repo_path=picsar_dir, filename=cwd + 'store_git_hashes.txt', name='picsar')
        store_git_hash(repo_path=amrex_dir , filename=cwd + 'store_git_hashes.txt', name='amrex' )
        store_git_hash(repo_path=warpx_dir , filename=cwd + 'store_git_hashes.txt', name='warpx' )

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
        res_dir += '_'.join([run_name, compiler, architecture, str(n_node)]) + '/'
        runtime_param_list = []
        # Deep copy as we change the attribute n_cell of
        # each instance of class test_element
        test_list_n_node = copy.deepcopy(test_list)
        job_time_min = time_min(len(test_list))
        batch_string = get_batch_string(test_list_n_node, job_time_min, module_Cname[architecture], n_node)
        # Loop on tests
        for count, current_run in enumerate(test_list_n_node):
            current_run.scale_n_cell(n_node)
            runtime_param_string  = ' amr.n_cell=' + ' '.join(str(i) for i in current_run.n_cell)
            runtime_param_string += ' amr.max_grid_size=' + str(current_run.max_grid_size)
            runtime_param_string += ' amr.blocking_factor=' + str(current_run.blocking_factor)
            runtime_param_string += ' max_step=' + str( current_run.n_step )
            # runtime_param_list.append( runtime_param_string )
            run_string = get_run_string(current_run, architecture, n_node, count, bin_name, runtime_param_string)
            batch_string += run_string
            batch_string += 'rm -rf plotfiles lab_frame_data diags\n'

        submit_job_command = get_submit_job_command()
        # Run the simulations.
        run_batch_nnode(test_list_n_node, res_dir, cwd, bin_name, config_command, batch_string, submit_job_command)
    os.chdir(cwd)
    # submit batch for analysis
    if os.path.exists( 'read_error.txt' ):
        os.remove( 'read_error.txt' )
    if os.path.exists( 'read_output.txt' ):
        os.remove( 'read_output.txt' )
    process_analysis(args.automated, cwd, compiler, architecture,
                     args.n_node_list, start_date, source_dir_base, res_dir_base)

# read the output file from each test and store timers in
# hdf5 file with pandas format
# -------------------------------------------------------
for n_node in n_node_list:
    print(n_node)
    if browse_output_files:
        res_dir = res_dir_base
        res_dir += '_'.join([run_name, compiler,\
                             architecture, str(n_node)]) + '/'
        for count, current_run in enumerate(test_list):
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
            if os.path.exists(path_hdf5 + perf_database_file):
                df_base = pd.read_hdf(path_hdf5 + perf_database_file, 'all_data')
                updated_df = df_base.append(df_newline, ignore_index=True)
            else:
                updated_df = df_newline
            # Write dataframe to file perf_database_file
            # (overwrite if file exists)
            updated_df.to_hdf(path_hdf5 + perf_database_file, key='all_data', mode='w', format='table')

# Extract sub-set of pandas data frame, write it to
# csv file and copy this file to perf_logs repo
# -------------------------------------------------
if args.mode=='read' and update_perf_log_repo:
    # get perf_logs repo
    git_repo = git.Repo( perf_logs_repo )
    if push_on_perf_log_repo:
        git_repo.git.stash('save')
        git_repo.git.pull()
    os.chdir( perf_logs_repo )
    sys.path.append('./')
    import write_csv
    git_repo.git.add('./logs_csv/' + csv_file[machine])
    git_repo.git.add('./logs_hdf5/' + perf_database_file)
    index = git_repo.index
    index.commit("automated tests")

# Rename all result directories for archiving purposes:
# include date in the name, and a counter to avoid over-writing
for n_node in n_node_list:
    if browse_output_files:
        res_dir = res_dir_base
        res_dir += '_'.join([run_name, compiler,\
                             architecture, str(n_node)]) + '/'
        # Rename directory with precise date+hour for archive purpose
        if rename_archive == True:
            loc_counter = 0
            res_dir_arch = res_dir_base
            res_dir_arch += '_'.join([year, month, day, run_name, compiler,\
                                      architecture, str(n_node), str(loc_counter)]) + '/'
            while os.path.exists( res_dir_arch ):
                loc_counter += 1
                res_dir_arch = res_dir_base
                res_dir_arch += '_'.join([year, month, day, run_name, compiler,\
                                          architecture, str(n_node), str(loc_counter)]) + '/'
            print("renaming " + res_dir + " -> " + res_dir_arch)
            os.rename( res_dir, res_dir_arch )
