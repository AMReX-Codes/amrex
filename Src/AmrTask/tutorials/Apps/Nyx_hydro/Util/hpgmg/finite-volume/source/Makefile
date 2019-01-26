#CC = cc
#CFLAGS = -g -O2 -I. -openmp -DUSE_BICGSTAB=1 -DUSE_SUBCOMM=1 -DUSE_FCYCLES=1 -DUSE_CHEBY=1 -DUSE_MPI=1 -DUNLIMIT_FMG_VCYCLES=1
#LFLAGS = -g -O2 -openmp

#CC = cc
#CFLAGS = -g -O2 -I. -fopenmp -DUSE_BICGSTAB=1 -DUSE_SUBCOMM=1  -DUSE_CHEBY=1 -DUSE_MPI=1 -DUSE_FCYCLES=1 -DUNLIMIT_FMG_VCYCLES=1
#LFLAGS = -g -O2 -fopenmp

CC = cc
CFLAGS = -g -O0 -I. -fopenmp -DUSE_BICGSTAB=1 -DUSE_SUBCOMM=1 -DUSE_FCYCLES=1 -DUSE_CHEBY=1 -DUSE_MPI=1 -DUNLIMIT_FMG_VCYCLES=1
LFLAGS = -g -O0 -fopenmp

#CC = cc
#CFLAGS = -g -O0 -I. -DUSE_BICGSTAB=1 -DUSE_SUBCOMM=1 -DUSE_FCYCLES=1 -DUSE_CHEBY=1 -DUSE_MPI=1 -DUNLIMIT_FMG_VCYCLES=1
#LFLAGS = -g -O0

#CC = cc
#CFLAGS = -g -O2 -I. -DUSE_BICGSTAB=1 -DUSE_SUBCOMM=1 -DUSE_FCYCLES=1 -DUSE_CHEBY=1 -DUSE_MPI=1
#LFLAGS = -g -O2

#CC = icc
#CFLAGS = -g -O2 -I. -DUSE_BICGSTAB=1 -DUSE_FCYCLES=1 -DUSE_CHEBY=1
#LFLAGS = -g -O2

#CC = mpiicc
#CFLAGS = -g -O2 -ip -xHost -I. -DUSE_BICGSTAB=1 -DUSE_SUBCOMM=1 -DUSE_FCYCLES=1 -DUSE_CHEBY=1 -DUSE_MPI=1
#LFLAGS = -g -O2 -ip -xHost

#CC = mpiicc
#CFLAGS = -g -O0 -I. -DUSE_BICGSTAB=1 -DUSE_SUBCOMM=1 -DUSE_FCYCLES=1 -DUSE_CHEBY=1 -DUSE_MPI=1 -DUSE_PERIODIC_BC=1 -DUNLIMIT_FMG_VCYCLES=1
#LFLAGS = -g -O0

OBJ = timers.o level.o operators.7pt.o mg.o solvers.o hpgmg_setup.o
MAIN = call_hpgmg_setup.o

%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) -o $@ $<

libhpgmg_test.a: $(OBJ)
	ar rcs $@ $^

clean:
	$(RM) *.o *.a

lib: libhpgmg_test.a

all: $(OBJ) $(MAIN)
	$(CC) -o hpgmg $^ $(LFLAGS)
