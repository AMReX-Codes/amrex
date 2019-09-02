
module_name = {'cpu': 'haswell.', 'knl': 'mic-knl.', 'gpu':'.'}

def executable_name(compiler, architecture):
    return 'perf_tests3d.' + compiler + \
        '.' + module_name[architecture] + 'TPROF.MPI.OMP.ex'

def get_config_command(compiler, architecture):
    config_command = ''
    config_command += 'module unload darshan;' 
    config_command += 'module load craype-hugepages4M;'
    if architecture == 'knl':
        if compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-haswell;'
        config_command += 'module load craype-mic-knl;'
    elif architecture == 'cpu':
        if compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-mic-knl;'
        config_command += 'module load craype-haswell;'
    return config_command

# This function runs a batch script with 
# dependencies to perform the analysis 
# after all performance tests are done.
def process_analysis(cwd, compiler, architecture, n_node_list, start_date):
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
    batch_string += 'module load h5py-parallel\n'
    batch_string += 'python run_automated.py --compiler=' + \
                    compiler + ' --architecture=' + architecture + \
                    ' --mode=read' + \
                ' --n_node_list=' + '"' + n_node_list + '"' + \
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

# Calculate simulation time. Take 5 min + 5 min / simulation
def time_min(nb_simulations):
    return 5. + len(test_list)*5.

def get_batch_string(test_list, job_time_str, Cname, n_node):
    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#SBATCH --job-name=' + test_list[0].input_file + '\n'
    batch_string += '#SBATCH --time=' + job_time_str + '\n'
    batch_string += '#SBATCH -C ' + Cname + '\n'
    batch_string += '#SBATCH -N ' + str(n_node) + '\n'
    batch_string += '#SBATCH -q regular\n'
    batch_string += '#SBATCH -e error.txt\n'
    batch_string += '#SBATCH --account=m2852\n'
    return batch_string

def get_run_string(current_test, architecture, n_node, count, bin_name, runtime_param_list):
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
    return srun_string
