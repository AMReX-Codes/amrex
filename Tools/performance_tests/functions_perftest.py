import os, shutil, re
import pandas as pd
import numpy as np
import git

def scale_n_cell(ncell, n_node):
     ncell_scaled = ncell[:]
     index_dim = 0
     while n_node > 1:
         ncell_scaled[index_dim] *= 2
         n_node /= 2
         index_dim = (index_dim+1) % 3
     return ncell_scaled

def store_git_hash(repo_path=None, filename=None, name=None):
    repo = git.Repo(path=repo_path)
    sha = repo.head.object.hexsha
    file_handler = open( filename, 'a+' )
    file_handler.write( name + ':' + sha + ' ')
    file_handler.close()

def get_file_content(filename=None):
    file_handler = open( filename, 'r' )
    file_content = file_handler.read()
    file_handler.close()
    return file_content

def run_batch_nnode(test_list, res_dir, bin_name, config_command, architecture='knl', Cname='knl', n_node=1, runtime_param_list=[]):
    # Clean res_dir
    if os.path.exists(res_dir):
        shutil.rmtree(res_dir, ignore_errors=True)
    os.makedirs(res_dir)
    # Copy files to res_dir
    cwd = os.environ['AUTOMATED_PERF_TESTS'] + '/WarpX/Tools/performance_tests/'
    bin_dir = cwd + 'Bin/'
    shutil.copy(bin_dir + bin_name, res_dir)
    os.chdir(res_dir)
    # Calculate simulation time. Take 5 min + 2 min / simulation
    job_time_min = 5. + len(test_list)*5.
    job_time_str = str(int(job_time_min/60)) + ':' + str(int(job_time_min%60)) + ':00'
    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#SBATCH --job-name=' + test_list[0].input_file + '\n'
    batch_string += '#SBATCH --time=' + job_time_str + '\n'
    batch_string += '#SBATCH -C ' + Cname + '\n'
    batch_string += '#SBATCH -N ' + str(n_node) + '\n'
    batch_string += '#SBATCH -q regular\n'
    batch_string += '#SBATCH -e error.txt\n'
    batch_string += '#SBATCH --account=m2852\n'
    
    for count, current_test in enumerate(test_list):
        shutil.copy(cwd + current_test.input_file, res_dir)
        srun_string = ''
        srun_string += 'export OMP_NUM_THREADS=' + str(current_test.n_omp) + '\n'
        # number of logical cores per MPI process
        if architecture == 'cpu':
            cflag_value = max(1, int(32/current_test.n_mpi_per_node) * 2) # Follow NERSC directives
        elif architecture == 'knl':
            cflag_value = max(1, int(64/current_test.n_mpi_per_node) * 4) # Follow NERSC directives
        output_filename = 'out_' + '_'.join([current_test.input_file, str(n_node), str(current_test.n_mpi_per_node), str(current_test.n_omp), str(count)]) + '.txt'
        srun_string += 'srun --cpu_bind=cores '+ \
                       ' -n ' + str(n_node*current_test.n_mpi_per_node) + \
                       ' -c ' + str(cflag_value)   + \
                       ' ./'  + bin_name + \
                       ' ' + current_test.input_file + \
                       runtime_param_list[ count ] + \
                       ' > ' + output_filename + '\n'
        batch_string += srun_string
        batch_string += 'rm -rf plotfiles ; rm -rf lab_frame_data\n'
    batch_file = 'slurm'
    f_exe = open(batch_file,'w')
    f_exe.write(batch_string)
    f_exe.close()
    os.system('chmod 700 ' + bin_name)
    os.system(config_command + 'sbatch ' + batch_file + ' >> ' + cwd + 'log_jobids_tmp.txt')
    return 0

def run_batch(run_name, res_dir, bin_name, config_command, architecture='knl',\
              Cname='knl', n_node=1, n_mpi=1, n_omp=1):
    # Clean res_dir
    if os.path.exists(res_dir):
        shutil.rmtree(res_dir)
    os.makedirs(res_dir)
    # Copy files to res_dir
    # Copy files to res_dir
    cwd = os.environ['WARPX'] + '/Tools/performance_tests/'
    bin_dir = cwd + 'Bin/'
    shutil.copy(bin_dir + bin_name, res_dir)
    shutil.copyfile(cwd + run_name, res_dir + 'inputs')
    os.chdir(res_dir)
    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#SBATCH --job-name=' + run_name + str(n_node) + str(n_mpi) + str(n_omp) + '\n'
    batch_string += '#SBATCH --time=00:23:00\n'
    batch_string += '#SBATCH -C ' + Cname + '\n'
    batch_string += '#SBATCH -N ' + str(n_node) + '\n'
    batch_string += '#SBATCH -q regular\n'
    batch_string += '#SBATCH -e error.txt\n'
    batch_string += '#SBATCH --account=m2852\n'
    batch_string += 'export OMP_NUM_THREADS=' + str(n_omp) + '\n'
    if architecture == 'cpu':
        cflag_value = max(1, int(32/n_mpi) * 2) # Follow NERSC directives
        batch_string += 'srun --cpu_bind=cores '+ \
                    ' -n ' + str(n_node*n_mpi) + \
                    ' -c ' + str(cflag_value)   + \
                    ' ./'  + bin_name + ' inputs > perf_output.txt'
    elif architecture == 'knl':
        # number of logical cores per MPI process
        cflag_value = max(1, int(64/n_mpi) * 4) # Follow NERSC directives
        batch_string += 'srun --cpu_bind=cores '     + \
                        ' -n ' + str(n_node*n_mpi) + \
                        ' -c ' + str(cflag_value)   + \
                        ' ./'  + bin_name + ' inputs > perf_output.txt\n'
    batch_file = 'slurm'
    f_exe = open(batch_file,'w')
    f_exe.write(batch_string)
    f_exe.close()
    os.system('chmod 700 ' + bin_name)
    os.system(config_command + 'sbatch ' + batch_file + ' >> ' + cwd + 'log_jobids_tmp.txt')
    return 0

# Read output file and return init time and 1-step time
def read_run_perf(filename, n_steps):
    timing_list = []
    # Search inclusive time to get simulation step time
    partition_limit = 'NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %'
    with open(filename) as file_handler:
        output_text = file_handler.read()
    # Get total simulation time
    line_match_totaltime = re.search('TinyProfiler total time across processes.*', output_text)
    total_time = float(line_match_totaltime.group(0).split()[8])
    search_area = output_text.partition(partition_limit)[2]
    line_match_looptime = re.search('\nWarpX::Evolve().*', search_area)
    time_wo_initialization = float(line_match_looptime.group(0).split()[3])
    timing_list += [str(total_time - time_wo_initialization)]
    timing_list += [str(time_wo_initialization/n_steps)]
    partition_limit1 = 'NCalls  Excl. Min  Excl. Avg  Excl. Max   Max %'
    partition_limit2 = 'NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %'
    file_handler.close()
    with open(filename) as file_handler:
        output_text = file_handler.read()
    # Search EXCLISUSIVE routine timings
    search_area = output_text.partition(partition_limit1)[2].partition(partition_limit2)[0]
    pattern_list = ['\nParticleContainer::Redistribute().*',\
                    '\nFabArray::FillBoundary().*',\
                    '\nFabArray::ParallelCopy().*',\
                    '\nPPC::CurrentDeposition.*',\
                    '\nPPC::FieldGather.*',\
                    '\nPPC::ParticlePush.*',\
                    '\nPPC::Evolve::Copy.*',\
                    '\nWarpX::EvolveEM().*',\
                    'Checkpoint().*',\
                    'WriteParticles().*',\
                    '\nVisMF::Write(FabArray).*',\
                    '\nWriteMultiLevelPlotfile().*',\
                    '\nParticleContainer::RedistributeMPI().*']
    for pattern in pattern_list:
        timing = '0'
        line_match = re.search(pattern, search_area)
        if line_match is not None:
            timing = [str(float(line_match.group(0).split()[3])/n_steps)]
        timing_list += timing
    return timing_list

# Write time into logfile
def write_perf_logfile(log_file, log_line):
    f_log = open(log_file, 'a')
    f_log.write(log_line)
    f_log.close()
    return 0

def get_nsteps(run_name):
    with open(run_name) as file_handler:
        run_name_text = file_handler.read()
    line_match_nsteps = re.search('\nmax_step.*', run_name_text)
    nsteps = float(line_match_nsteps.group(0).split()[2])
    return nsteps

def extract_dataframe(filename, n_steps):
    # Get init time and total time through Inclusive time                                                                                                                            
    partition_limit_start = 'NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %'
    with open(filename) as file_handler:
        output_text = file_handler.read()
    # get total simulation time                                                                                                                                                      
    line_match_totaltime = re.search('TinyProfiler total time across processes.*', output_text)
    total_time = float(line_match_totaltime.group(0).split()[8])
    # get time performing steps as Inclusive WarpX::Evolve() time                                                                                                                    
    search_area = output_text.partition(partition_limit_start)[2]
    line_match_looptime = re.search('\nWarpX::Evolve().*', search_area)
    time_wo_initialization = float(line_match_looptime.group(0).split()[3])
    # New, might break something
    line_match_WritePlotFile = re.search('\nWarpX::WritePlotFile().*', search_area)
    if line_match_WritePlotFile is not None:
         time_WritePlotFile = float(line_match_WritePlotFile.group(0).split()[3])
    else:
         time_WritePlotFile = 0.
    # Get timers for all routines                                                                                                                                                    
    # Where to start and stop in the output_file                                                                                                                                     
    partition_limit_start = 'NCalls  Excl. Min  Excl. Avg  Excl. Max   Max %'
    partition_limit_end   = 'NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %'
    # Put file content in a string                                                                                                                                                   
    with open(filename) as file_handler:
        output_text = file_handler.read()
    # Keep only profiling data                                                                                                                                                       
    search_area = output_text.partition(partition_limit_start)[2]\
                             .partition(partition_limit_end)[0]
    list_string = search_area.split('\n')[2:-4]
    time_array = np.zeros(len(list_string))
    column_list= []
    for i in np.arange(len(list_string)):
        column_list.append(list_string[i].split()[0])
        time_array[i] = float(list_string[i].split()[3])
    df = pd.DataFrame(columns=column_list)
    df.loc[0] = time_array
    df['time_initialization'] = total_time - time_wo_initialization
    df['time_running'] = time_wo_initialization
    df['time_WritePlotFile'] = time_WritePlotFile
    # df['string_output'] = partition_limit_start + '\n' + search_area
    return df

# Run a performance test in an interactive allocation                                                                                                                                                                                                                            
# def run_interactive(run_name, res_dir, n_node=1, n_mpi=1, n_omp=1):
#     # Clean res_dir                                                                                                                                                                                                                                                           #  
#     if os.path.exists(res_dir):
#         shutil.rmtree(res_dir)
#     os.makedirs(res_dir)
#     # Copy files to res_dir                                                                                                                                                                                                                                                   #  
#     shutil.copyfile(bin_dir + bin_name, res_dir + bin_name)
#     shutil.copyfile(cwd  + run_name, res_dir + 'inputs')
#     os.chdir(res_dir)
#     if args.architecture == 'cpu':
#         cflag_value = max(1, int(32/n_mpi) * 2) # Follow NERSC directives                                                                                                                                                                                                     #  
#         exec_command = 'export OMP_NUM_THREADS=' + str(n_omp) + ';' +\
#                        'srun --cpu_bind=cores '     + \
#                        ' -n ' + str(n_node*n_mpi) + \
#                        ' -c ' + str(cflag_value)   + \
#                        ' ./'  + bin_name + ' inputs > perf_output.txt'
#     elif args.architecture == 'knl':
#         # number of logical cores per MPI process                                                                                                                                                                                                                             #  
#         cflag_value = max(1,int(68/n_mpi) * 4) # Follow NERSC directives                                                                                                                                                                                                      #  
#         exec_command = 'export OMP_NUM_THREADS=' + str(n_omp) + ';' +\
#                        'srun --cpu_bind=cores '     + \
#                        ' -n ' + str(n_node*n_mpi) + \
#                        ' -c ' + str(cflag_value)   + \
#                        ' ./'  + bin_name + ' inputs > perf_output.txt'
#     os.system('chmod 700 ' + bin_name)
#     os.system(config_command + exec_command)
#     return 0
