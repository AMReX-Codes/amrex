require 5.6.0;
use strict;
use warnings;

use vars qw($opt_b $opt_d $opt_r);
use Getopt::Std;
my $usage = "usage: $0 [-b|-d] [-r] [filename]\n";
getopts("bdr") or die $usage;
die $usage if ($opt_b && $opt_d);

my %pairs;	# all pairs ($l, $r)
my %npred;	# number of predecessors
my %succ;	# list of successors
my @selfs;

while (<>) {
    my ($l, $r) = my @l = split;
    next unless @l == 2;
    if ( $l eq $r ) {
        push @selfs, $l;
	next;
    }
    next if defined $pairs{$l}{$r};
    $pairs{$l}{$r}++;
    $npred {$l} += 0;
    ++$npred{$r};
    push @{$succ{$l}}, $r;
}

# create a list of nodes without predecessors

my @list = grep {!$npred{$_}} keys %npred;

my @olist;
while (@list) {
    $_ = pop @list;
    push @olist, $_;
    foreach my $child (@{$succ{$_}}) {
	if ($opt_b) {	# breadth-first
	    unshift @list, $child unless --$npred{$child};
	} else {	# depth-first (default)
	    push @list, $child unless --$npred{$child};
	}

    }
}

warn "cycle detected\n" if grep {$npred{$_}} keys %npred;

my @oself = grep {!$succ{$_}} @selfs;

@oself = (@oself, @olist);

@oself = reverse @oself if $opt_r;

foreach my $oli (@oself) {
    print "$oli\n";
}

__END__

=pod

=head1 NAME

TCSORT - Toplogical sort

=head1 DESCRIPTION

tcsort: Reads lists of pairs of node names of directed arcs in a graph.
Prints nodes in topological order.
Cycles indicate an error.

=head1 AUTHER

Jeffrey S. Haemer, unknown affiliation

Modifications, CCSE/LBNL.

=cut
