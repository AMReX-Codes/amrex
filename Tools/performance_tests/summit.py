

def executable_name(compiler,architecture):
    return 'perf_tests3d.' + compiler + 'TPROF.MPI.ACC.CUDA.ex'

def get_config_command(compiler, architecture):
    config_command = ''
    config_command += 'module load pgi;'
    config_command += 'module load cuda;'
    return config_command

# This function runs a batch script with 
# dependencies to perform the analysis 
# after all performance tests are done.
def process_analysis(cwd, compiler, architecture, n_node_list, start_date):
    print("process analysis...")

# Calculate simulation time. Take 2 min + 2 min / simulation
def time_min(nb_simulations):
    return 2. + len(test_list)*2.
