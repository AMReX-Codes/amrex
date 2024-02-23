default:
	bison --no-lines --defines=amrex_parser.tab.h --output=amrex_parser.tab.nolint.H amrex_parser.y
	flex --noline --header-file=amrex_parser.lex.h --outfile=amrex_parser.lex.nolint.H amrex_parser.l
	bison --no-lines --defines=amrex_iparser.tab.h --output=amrex_iparser.tab.nolint.H amrex_iparser.y
	flex --noline --header-file=amrex_iparser.lex.h --outfile=amrex_iparser.lex.nolint.H amrex_iparser.l
