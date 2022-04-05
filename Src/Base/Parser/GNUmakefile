default:
	bison --no-lines --defines=amrex_parser.tab.h --output=amrex_parser.tab.cpp amrex_parser.y
	flex --noline --header-file=amrex_parser.lex.h --outfile=amrex_parser.lex.cpp amrex_parser.l
	sed -i '1s/^/#ifdef __GNUC__\n#pragma GCC diagnostic ignored "-Wnull-dereference"\n#pragma GCC diagnostic ignored "-Wunreachable-code"\n#pragma GCC diagnostic ignored "-Wfloat-conversion"\n#elif defined(__clang__)\n#pragma clang diagnostic ignored "-Wnull-dereference"\n#pragma clang diagnostic ignored "-Wunreachable-code"\n#pragma clang diagnostic ignored "-Wfloat-conversion"\n#endif\n/' amrex_parser.lex.cpp
	bison --no-lines --defines=amrex_iparser.tab.h --output=amrex_iparser.tab.cpp amrex_iparser.y
	flex --noline --header-file=amrex_iparser.lex.h --outfile=amrex_iparser.lex.cpp amrex_iparser.l
	sed -i '1s/^/#ifdef __GNUC__\n#pragma GCC diagnostic ignored "-Wnull-dereference"\n#pragma GCC diagnostic ignored "-Wunreachable-code"\n#elif defined(__clang__)\n#pragma clang diagnostic ignored "-Wnull-dereference"\n#pragma clang diagnostic ignored "-Wunreachable-code"\n#endif\n/' amrex_iparser.lex.cpp
