from subprocess import call

nthreads = 2
nprocs   = 2
steps    = 10000

cmd  = "OMP_NUM_THREADS=2 mpirun -np 2 "
cmd += "./main2d.gnu.TPROF.MPI.OMP.ex inputs max_step=%d do_nl=0 &> run_gold.out" % steps
print(cmd)
call(cmd, shell=True)

cmd = "mv particles%.5d particles%.5d" % (steps, steps) + "gold" 
print(cmd)
call(cmd, shell=True)

cmd  = "OMP_NUM_THREADS=2 mpirun -np 2 "
cmd += "./main2d.gnu.TPROF.MPI.OMP.ex inputs max_step=%d do_nl=1 &> run_nl.out" % steps
print(cmd)
call(cmd, shell=True)

cmd = "mv particles%.5d particles%.5d" % (steps, steps) + "nl" 
print(cmd)
call(cmd, shell=True)

with open("particles%.5d" % steps + "gold", "r") as f:
    lines_gold = f.readlines()
    lines_gold.sort()

with open("particles%.5d" % steps + "nl", "r") as f:
    lines_nl = f.readlines()
    lines_nl.sort()

for l1, l2 in zip(lines_gold, lines_nl):
    assert(l1 == l2)
