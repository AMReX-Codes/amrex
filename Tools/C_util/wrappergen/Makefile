# Generated automatically from Makefile.in by configure.
#  (C) 1993 by Argonne National Laboratory and Mississipi State University.
#

MPI_PROTO = '-DPROTO_FILE="./mpi_proto"'
CPPFLAGS = $(INCLUDE_DIR)
CFLAGS = -g

default: wrappergen

all: $(EXECS)

readproto.o: readproto.c
	$(CC) $(CFLAGS) $(MPI_PROTO) -c $*.c

wrappergen: readproto.o wrappergen.o args.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o wrappergen $(CFLAGS) readproto.o wrappergen.o args.o

write_proto: prof_wrapper.o driver.o doc.o write_proto.o system.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o write_proto driver.o prof_wrapper.o doc.o write_proto.o system.o -lm

args.o: args.h expandingList.h
wrappergen.o: expandingList.h wrappergen.h
readproto.o: args.h expandingList.h wrappergen.h

clean:
	/bin/rm -rf $(EXECS) *.o *~
