# Generated automatically from Makefile.in by configure.
#  (C) 1993 by Argonne National Laboratory and Mississipi State University.
#  (C) 2001 by Lawrence Berkely National Laboratory
#

CC=gcc -Wall
CPPFLAGS = $(INCLUDE_DIR)
CFLAGS = -g

default: wrappergen

all: $(EXECS)

wrappergen: readproto.o wrappergen.o args.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o wrappergen $(CFLAGS) readproto.o wrappergen.o args.o

write_proto: prof_wrapper.o driver.o doc.o write_proto.o system.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o write_proto driver.o prof_wrapper.o doc.o write_proto.o system.o -lm

args.o: args.h expandingList.h
wrappergen.o: expandingList.h wrappergen.h
readproto.o: args.h expandingList.h wrappergen.h

prof_wrapper.o: expandingList.h wrappergen.h expandingList.h
driver.o: tools.h doc.h expandingList.h
doc.o: doc.h
write_proto.o: wrappergen.h expandingList.h
sytem.o: system.h tools.h


clean:
	/bin/rm -rf $(EXECS) *.o *~
