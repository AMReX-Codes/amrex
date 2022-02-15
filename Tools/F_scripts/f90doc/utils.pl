package utils;

use strict;

sub copy_list {
   my ($listref) = @_;
   my @list;
   @list = @$listref;
   \@list;
}

sub copy_hash {
   my ($hashref) = @_;
   my %hash;
   %hash = %$hashref;
   \%hash;
}

sub hash2str {
   my ($hash) = @_;
   my ($key, $s);
   $s = "{\n";
   foreach $key (keys %$hash) {
      $s .= "   $key => $hash->{$key}\n";
   }
   $s .= "}";
}

sub trim {
   my ($s) = @_;
   $s =~ s/^\s*//;
   $s =~ s/\s*$//;
   $s;
}

# balsplit (sep, string) splits string into pieces divided by sep when
# sep is "outside" ()s.  Returns a list just like split.
sub balsplit {
   my ($sep, $str) = @_;
   my ($i, $c);
   my ($len, $level, $left) = (length ($str), 0, 0);
   my (@list) = ();

   for ($i = 0; $i < $len; $i++) {
      $c = substr ($str, $i, 1);
      if ($c eq "(") {
         $level++;
      } elsif ($c eq ")") {
         $level--;
         die "balsplit: Unbalanced parens (too many )'s)" if $level < 0;
      } elsif ($c eq $sep && $level == 0) {
         push (@list, substr ($str, $left, $i-$left));
         $left = $i + 1;
      }
   }

   push (@list, substr ($str, $left));
   return @list;
}

# Takes the first word of each element of the list.
sub leftword {
   my ($listref) = @_;
   my @out = ();
   my ($x);
   foreach $x (@$listref) {
      $x =~ s/^\s*//;
      $x =~ /^\w*/;
      push (@out, $&);
   }
   @out;
}

sub remove_blanks {
   my ($listref) = @_;
   my @out = ();
   my ($x);
   foreach $x (@$listref) {
      push (@out, $x) unless $x =~ /^\s*$/;
   }
   @out;
}

sub do_nothing {
}

1;
