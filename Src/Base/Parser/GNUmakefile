default:
	bison --no-lines --defines=amrex_parser.tab.h --output=amrex_parser.tab.cpp amrex_parser.y
	flex --noline --header-file=amrex_parser.lex.h --outfile=amrex_parser.lex.cpp amrex_parser.l
	sed -i '1s/^/#pragma GCC diagnostic ignored "-Wnull-dereference"\n#pragma GCC diagnostic ignored "-Wunreachable-code"\n#pragma GCC diagnostic ignored "-Wfloat-conversion"\n/' amrex_parser.lex.cpp
