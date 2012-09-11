use 5.6.0;
use strict;
use warnings;
use Getopt::Long;

my $noglob;
BEGIN {
    eval "use File::Glob";
    if ($@) { $noglob =1 }
}

my $tsort;
my $debug;
my @incdirs;
my $opt_help;
my $def_fixed_source = 0;
my @forexts;
my @vpath = qw( . );
my $odir = ".";

my $obj_ext = 'o';
my $pathsep = "/";
my $dirsep = ":";
my @srcpaths;
my $is_win32;
if ( "$^O" eq "MSWin32" ) {
    $obj_ext = 'obj';
    $pathsep = "\\";
    $dirsep  = ";";
    $is_win32 = 1;
}

my $usage = "usage: $0 [--debug] [--odir dir] [--tsort] [--objext ext] [--fixed] [--I ldir]\n";

GetOptions(
    "help" 	=> \$opt_help,
    "debug"	=> \$debug,
    "fixed"	=> \$def_fixed_source,
    "tsort"	=> \$tsort,
    "objext=s"	=> \$obj_ext,
    "odir=s"    => \$odir,
    "forexts=s@" => \@forexts,
    "srcpath=s@" => \@srcpaths,
    "I=s@"	=> \@incdirs) or die $usage;

@incdirs = (".", split(/$dirsep/, join($dirsep,@incdirs)));

#print "incdirs = @incdirs\n";
$odir = $odir . $pathsep;

my @suffixes = qw( .f90 .f .f95 .for );
push @suffixes, @forexts;
push @vpath, @srcpaths;


# needed for Win32 expansion of file names.
unless ( $noglob && ! $is_win32 ) {
    @ARGV = map {
	    my @g = File::Glob::glob($_) if /[*?]/;
	    @g ? @g : $_;
	} @ARGV;
}

my %sources;

foreach my $file (@ARGV) {
    if ( $file =~ m/--fixed/ ) {
	$def_fixed_source = 1;
	next;
    }
    if ( $file =~ m/--free/ ) {
        $def_fixed_source = 0;
    	next;
    }
    my $myfile = $file;

    if ( defined($myfile) ) {
      $sources{$myfile} = new Source_File($myfile);
      $sources{$myfile}->find_uses();
    }
    else
    {
        foreach my $dir (@vpath) {
            my $t = $dir . $pathsep . $file;
            if ( -f $t ) { $myfile = $t; last; }
        }
        if ( defined($myfile) ) {
            $sources{$myfile} = new Source_File($myfile);
            $sources{$myfile}->find_uses();
        }
    }
}

foreach my $target (sort keys(%sources)) {
    $sources{$target}->print();
}

package Source_File;

use File::Basename;

my %mod_files = ();

sub new {
    my ($type, $filename) = @_;
    my $class = ref($type) || $type;
    my $self = {};
    $self->{'Source_File'} = $filename;
    $self->{'uses'}	   = {};
    $self->{'modules'}     = {};
    $self->{'includes'}    = {};
    bless $self, $class;
}

sub find_uses {
    my $self = shift;
    my $file = $self->{'Source_File'};
    local(*FILE);
    local($_);

    if (-f $file) {
	open(FILE, $file) || warn "Can't open $file: $!\n";
    } else {
	return;
    }
    # assume that .f, .for files are f77 files.
    my $fixed_source = 0;
    $fixed_source = 1 if $file =~ m/\.f(or)?$/i || $def_fixed_source;
    print "#\tchecking $file with  fixed_source = $fixed_source\n" if $debug;
    my $line;
    my $nline;
    while ( defined($line=<FILE>) ) {
	$nline = "";
	chomp $line;
	# eliminate comments
	$line =~ s/!.*$/ /;
	unless ( $fixed_source ) {
	    # eat up free-format continuation lines
	    if ( $line =~ s/&$/ / ) {
	        $line .= <FILE>;
	        redo unless eof(FILE);
	    }
	} else {
	    # eat up fixed-format continuation lines
	    next if $line =~ m/^[c*]/i;
	    $nline = <FILE> unless eof(FILE);
	    if ( $nline =~ s/^     [^ ]/ / ) {
		$line .= $nline;
		redo unless eof(FILE);
	    }
	}
	# allow the module statement after a semi-colon;
    	if ( $line =~ /(^|;)\s*module\s+(\w+)/i) {
	    my $module = $2;
	    next if (lc($module) eq "procedure");
	    $mod_files{$module} = $file;
	    $self->{'modules'}{$module} = 1;
	    print "#\t\tdefines module = $module\n" if $debug;
	}
	# allow more than one 'use on a line'
	while ( $line =~ /(^|;)\s*use\s+(\w+)/ig ) {
	    my $use = $2;
	    next if (lc($use) eq "iso_c_binding");
	    $self->{'uses'}{$use} = 1;
	    print "#\t\tuses = $use\n" if $debug;
	}
	# fortran style includes must be the only 'thing' on a line.
	# they are not 'statements'
	if ( $line =~ /^\s*include\s+['"]([^'"]*)['"]/i) {
	    my $inc = $1;
	    $self->{'includes'}{$inc} = 1;
	    print "#\t\tincludes = $inc\n" if $debug;
	}
	if ( $nline ) {
	    $line = $nline;
	    redo;
	}
    }
    close(FILE);
}

sub print {
    my $self = shift;
    my $source = $self->{'Source_File'};
    my $base = basename($source, @suffixes);
    my $target = $odir . $base . "." . $obj_ext;
    if ( $tsort ) {
	print "$source $source\n";
    } else {
	print "$target: $source\n";
    }
    print "#\tchecking source = $source\n" if $debug;
    my $mod;
    foreach $mod ( keys %{$self->{'uses'}}) {
	print "#\t\tchecking for self-usage $mod\n" if $debug;
	if ( ${$self->{'modules'}}{$mod} ) {
	    print "#\t\t\tdeleting $mod\n" if $debug;
	    delete ${$self->{'uses'}}{$mod};
	}
    }
    print "#\tprinting $target dependencies\n" if $debug;
    my %dep_done;
    foreach $mod ( keys %{$self->{'uses'}}) {
	print "#\t\tmod($mod)\n" if $debug;
	if ( !exists($mod_files{$mod}) ) {
	    warn "#\t\tfile not found for module $mod";
	} elsif ( !exists($dep_done{$mod}) ) {
	    my $file = $mod_files{$mod};
	    my $base = basename($file, @suffixes);
	    print "#\t\tfile($file) --> base($base)\n" if $debug;
	    my $dep = $odir . $base . "." . $obj_ext;
	    if ( $tsort ) {
		print "$file $source\n";
	    } else {
		print "$target: $dep\n";
	    }
	    $dep_done{$mod} = 1;
	}
    }
    foreach my $file (keys %{$self->{'includes'}}) {
	print "#\tinc file = $file\n" if $debug;
	foreach my $dir (@incdirs) {
	    my $path = $dir . $pathsep . $file;
	    print "#\t\ttesting path = $path\n" if $debug;
	    if ( -f "$path" ) {
		print "#\t\t\t$path works\n" if $debug;
		print "$target: $path\n";
		last;
	    }
	}
    }
}
__END__

=pod

=head1 NAME

moddep.pl - Determine MODULE/INCLUDE dependencies for Fortran.

=head1 AUTHOR

LBNL/CCSE.

=cut
