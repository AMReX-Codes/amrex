#!/usr/bin/perl
# ------------------------------------------------------------------
# Perl script to generate a "make"-style dependency list for a
# C, C++, source file with CPP include directives.
#
# Usage:  mkdep --debug [--I <dir>]*  filename ...
# Notes:  *  --I <path> defines a search path for include files
#         *  --debug turn on debug flag
#         *   searches current directory only if -I. is in search
#             path or #include directive uses double quotes rather
#             than angle brackets.
#         *   dependency list is sent to standard output
#         *   follows all #include directives, even those protected
#             by #if and #ifdef directives.
#         *   ignores include files not found in search path
#         *   No duplications in dependency list
#
# Author:  Michael Welcome
#          4/26/95
#          Lawrence Livermore National Laboratory
# Revised: CCSE
# 	   5/29/04
# ------------------------------------------------------------------

use 5.6.0;
use strict;
use warnings;
use Getopt::Long;

my $debug = 0;
my @incdirs;

my $odir = ".";
my $obj_ext = 'o';
my $pathsep = "/";
my $opt_help;
my $dirsep = ":";

if ( "$^O" eq "MSWin32" )
  {
    $obj_ext = 'obj';
    $pathsep = "\\";
    $dirsep = ";";
  }

# search command line for -I and -X options

my $usage = "usage: $0 [--odir dir] [--I ldir]\n";

GetOptions(
    "help"	=> \$opt_help,
    "debug"	=> \$debug,
    "objext"	=> \$obj_ext,
    "odir=s"	=> \$odir,
    "I=s@"	=> \@incdirs) or die $usage;

@incdirs = (".", split(/$dirsep/, join($dirsep,@incdirs)));

$odir = $odir . $pathsep;

my @suffixes = qw( .c .cpp );
use File::Basename;

foreach my $ifile (@ARGV)
  {
    print "# PARSING FILE: $ifile\n" if $debug;
    die "cannot read $ifile\n" if (!(-e $ifile && -r _));

    my $base = basename($ifile, @suffixes);
    my $ofile = $odir . $base . "." . $obj_ext;
    my @searchlist = ("$ifile\"");
    my %usedfiles = ();
    my %deptable = ();
    $usedfiles{ $ifile } = 1;

    while (@searchlist)
      {
        # get next file off search list
	my $file = shift(@searchlist);

        # NOTE: the last char in $file is either a double quote (") or
        #       a right angle bracket (>).  Strip off this character and
	#       save it.  If it is a quote, search current directory first.
	#       If its a right angle bracket, only search directories
	#       specified with -I options.
	my $incltype = chop $file;

	# NOTE: if the first char in $file is a "/" indicating an absolute
	#       path, do not search directory list
	my $abspath = ($file =~ /^\// ? 1 : 0);

	foreach my $d (@incdirs)
	  {
	    if ($d ne "" && $abspath)
	      {
		# this is an absolute path, dont search current directory
		next;
	      }
	    if ($d eq "" && $incltype eq ">")
	      {
		# dont search current directory
		next;
	      }
	    my $dep = "$d" . $pathsep . "$file";
	    print "# Will search $d for $file : $dep\n" if $debug;
	    if (-e $dep)
	      {

		# file found, build dependency, enter in table
	        # print "$ofile: $dep\n";
		$deptable{ $dep } = 1;

		# grep file for includes and, if not seen before, add
		# to end of search list
		open(FL,"$dep") || die "cant open $dep\n";
		while (<FL>)
		  {
		    if (/^\s*#\s*include\s+["<]\s*([^">]+[">])/)
		      {
			if ($usedfiles{ $1 }++ == 0)
			  {
			    print "# including  $1\n" if $debug;
			    # this file not searched yet, add to list
			    # NOTE: last char is either double quote
			    #       or right angle bracket.
			    push(@searchlist,($1));
			}
		    }
		}

		# since file was found in search path, jump out of loop
		last;
	    }
	}
	# print "# @searchlist\n" id $debug;
    }

    # now generate dependency list
    for my $dep (keys %deptable)
      {
	print "$ofile: $dep\n";
    }
}
