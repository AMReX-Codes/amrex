use 5.6.0; 
use strict;
use warnings;
use Getopt::Long;

my $noglob;
BEGIN {
  eval "use File::Glob";
  if ($@) {
    $noglob =1;
  }
}

my $dirsep = ":";
my $is_win32 = 0;
if ( "$^O" eq "MSWin32" ) {
  $dirsep = ";";
  $is_win32 = 1;
}

my $usage = "usage: $0 [--debug] [--output file] files ... \n";

my $debug;
my $opt_help;
my $the_outfile = "-";

GetOptions(
	   "help"     => \$opt_help,
	   "debug"    => \$debug,
	   "output=s" => \$the_outfile,
	  ) or die $usage;

open( OUT, ">$the_outfile") || die "Couldn't write to \"$the_outfile\": $!.\n";

unless ( $noglob && ! $is_win32 ) {
  @ARGV = map {
    my @g = File::Glob::glob($_) if /[*?]/;
    @g ? @g : $_;
  } @ARGV;
}

if ( @ARGV == 0 ) {
  push @ARGV, "-";
}

foreach my $file (@ARGV) {
  open(FILE, $file) || die "Can't open $file: $!\n";
  my $line;
  while ( defined($line=<FILE>) ) {
    chomp $line;
    if ( $line =~ /^\?\?/ ) {
      $line =~ s/^\?\? *//;
      $line =~ s/!.*$/ /;
      $line =~ s/ *$//;
      if ( $line =~ s/&$// ) {
	my $nline = <FILE>;
	$nline =~ s/^\?\? *//;
	$line .= $nline;
	redo unless eof(FILE);
      }
      print OUT "line = \"$line\"\n";
    }
  }
  close(FILE);
}
