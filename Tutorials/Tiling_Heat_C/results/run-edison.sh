#PBS -q debug
#PBS -l walltime=00:30:00
#PBS -l mppwidth=24
#PBS -j eo

cd $PBS_O_WORKDIR

# aprun -n 1 -N 1 -S 1 -ss ./Heat3d.Linux.Intel.Intel.ex inputs_3d fabarray.mfiter_tile_size="1024 4 4" &> run-4.out
# aprun -n 1 -N 1 -S 1 -ss ./Heat3d.Linux.Intel.Intel.ex inputs_3d fabarray.mfiter_tile_size="1024 8 8" &> run-8.out
# aprun -n 1 -N 1 -S 1 -ss ./Heat3d.Linux.Intel.Intel.ex inputs_3d fabarray.mfiter_tile_size="1024 16 16" &> run-16.out
# aprun -n 1 -N 1 -S 1 -ss ./Heat3d.Linux.Intel.Intel.ex inputs_3d fabarray.mfiter_tile_size="1024 32 32" &> run-32.out
# aprun -n 1 -N 1 -S 1 -ss ./Heat3d.Linux.Intel.Intel.ex inputs_3d fabarray.mfiter_tile_size="1024 64 64" &> run-64.out
# aprun -n 1 -N 1 -S 1 -ss ./Heat3d.Linux.Intel.Intel.ex inputs_3d do_tiling=0 &> run-notilig.out

export KMP_AFFINITY=verbose,granularity=core,explicit,proclist=[0]
export OMP_NUM_THREADS=1
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=0 &> run-nt${OMP_NUM_THREADS}-dt0.out
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=1 &> run-nt${OMP_NUM_THREADS}-dt1.out

export KMP_AFFINITY=verbose,granularity=core,explicit,proclist=[0,1]
export OMP_NUM_THREADS=2
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=0 &> run-nt${OMP_NUM_THREADS}-dt0.out
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=1 &> run-nt${OMP_NUM_THREADS}-dt1.out

export KMP_AFFINITY=verbose,granularity=core,explicit,proclist=[0,1,2,3]
export OMP_NUM_THREADS=4
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=0 &> run-nt${OMP_NUM_THREADS}-dt0.out
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=1 &> run-nt${OMP_NUM_THREADS}-dt1.out

export KMP_AFFINITY=verbose,granularity=core,explicit,proclist=[0,1,2,3,4,5,6,7]
export OMP_NUM_THREADS=8
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=0 &> run-nt${OMP_NUM_THREADS}-dt0.out
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=1 &> run-nt${OMP_NUM_THREADS}-dt1.out

export KMP_AFFINITY=granularity=core,explicit,proclist=[0,1,2,3,4,5,6,7,8,9,10,11]
export OMP_NUM_THREADS=12
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=0 &> run-nt${OMP_NUM_THREADS}-dt0.out
aprun -n 1 -N 1 -S 1 -d ${OMP_NUM_THREADS} -cc none ./Heat3d.Linux.Intel.Intel.OMP.ex inputs_3d do_tiling=1 &> run-nt${OMP_NUM_THREADS}-dt1.out

