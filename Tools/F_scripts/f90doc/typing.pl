package typing;

use strict;

# Stores the type of each variable.
$typing::typeof = "";
# Stack: one typeof per scope.
@typing::typeofs = ();

# Stores the definition of each type.
$typing::typedef = "";
# Stack: one typedef per scope.
@typing::typedefs = ();

# Stores the definition of each function/operator.
$typing::code = "";
# Stack: one code per scope.
@typing::codes = ();


# DOUBLE PRECISION type.
$typing::double_precision = typing::make_type ('real', 8, "double precision");

# Default character kind.
$typing::default_character_kind = 1;

# Default types.
%typing::default_type = (
  'complex' => typing::make_type ('complex', 8, "complex"),
  'integer' => typing::make_type ('integer', 4, "integer"),
  'logical' => typing::make_type ('logical', 1, "logical"),
  'real'    => typing::make_type ('real', 4, "real"),
);
$typing::default_type{'character'} = typing::make_character_type ();

# Types with wild sub and any other info (just a base defined).
$typing::wild_type = {
   'complex'   => typing::make_type ('complex'),
   'real'      => typing::make_type ('real'),
   'integer'   => typing::make_type ('integer'),
   'logical'   => typing::make_type ('logical'),
   'character' => typing::make_type ('character')
};


# Precedence of operations; based on that which is in expr_parse.y.
# Higher precedence indicated by larger number.
$typing::precedence = {
  '.eqv.'  => 1,
  '.neqv.' => 1,
  '.or.'   => 2,
  '.and.'  => 3,
  '.not.'  => 4,
  '<'      => 5,
  '>'      => 5,
  '<='     => 5,
  '>='     => 5,
  '=='     => 5,
  '/='     => 5,
  '//'     => 6,
  '+'      => 7,
  '-'      => 7,
  'u+'     => 8,
  'u-'     => 8,
  '*'      => 9,
  '/'      => 9,
  '**'     => 10,
  '%'      => 11,
  '%call'  => 11,
  '%colon' => 30, # this is a guess
  '%namedarg' => 30, # this is a guess
  '%array' => 40,    # as in "forty days and forty nights," which means
  '%const' => 40,    #    "a long time," here we use 40 as an approx. to infty.
  '%var'   => 40,
  '%do'    => 40,
};

#####
# Starts a new scope.  If this is a top-level scope, initializes the codes
# to intrinsics and the like.
#####
sub new_scope {
   my ($newtypeof, $newtypedef, $newcode);

   if (@typing::typeofs) {
      $typing::typeof = utils::copy_hash ($typing::typeof);
      $typing::typedef = utils::copy_hash ($typing::typedef);
      $typing::code = utils::copy_hash ($typing::code);
   } else {
      $typing::typeof = {};
      $typing::typedef = {};
      $typing::code = {};
      $typing::code{"//"} = [ {
         'parms' => [ $typing::wild_type{'character'},
                      $typing::wild_type{'character'} ],
         'return' => $typing::wild_type{'character'}
      } ];
      my ($int, $real, $logical, $char) = ( $typing::wild_type{'integer'},
         $typing::wild_type{'real'}, $typing::wild_type{'logical'},
         $typing::wild_type{'character'} );
      my ($op);
      foreach $op ("+", "-", "*", "/") {
         $typing::code->{$op} = [
            { 'parms' => [ $int, $int ], 'return' => $int },
            { 'parms' => [ $real, $int ], 'return' => $real },
            { 'parms' => [ $int, $real ], 'return' => $real },
            { 'parms' => [ $real, $real ], 'return' => $real }
         ];
      }
      $typing::code->{"**"} = [
         { 'parms' => [ $int, $int ], 'return' => $int },
         { 'parms' => [ $real, $int ], 'return' => $real },
         { 'parms' => [ $int, $real ], 'return' => $real },
         { 'parms' => [ $real, $real ], 'return' => $real },
      ];
      foreach $op ("u+", "u-") {
         $typing::code->{$op} = [
            { 'parms' => [ $int ], 'return' => $int },
            { 'parms' => [ $real ], 'return' => $real }
         ];
      }
      foreach $op ("<", "<=", "==", "/=", ">", ">=") {
         $typing::code->{$op} = [
            { 'parms' => [ $int, $int ], 'return' => $logical },
            { 'parms' => [ $real, $int ], 'return' => $logical },
            { 'parms' => [ $int, $real ], 'return' => $logical },
            { 'parms' => [ $real, $real ], 'return' => $logical },
            { 'parms' => [ $char, $char ], 'return' => $logical }
         ];
      }
      foreach $op (".or.", ".and.", ".eqv.", ".neqv.") {
         $typing::code->{$op} = [
            { 'parms' => [ $logical, $logical ], 'return' => $logical }
         ];
      }
      $typing::code->{".not."} = [
         { 'parms' => [ $logical ], 'return' => $logical }
      ];
      $typing::code->{"//"} = [
         { 'parms' => [ $char, $char ], 'return' => $char }
      ];
   }

   push @typing::typeofs, $typing::typeof;
   push @typing::typedefs, $typing::typedef;
   push @typing::codes, $typing::code;
}

#####
# Ends an old scope.
#####
sub end_scope {
   pop @typing::typeofs;
   pop @typing::typedefs;
   pop @typing::codes;

   if ($typing::typeofs) {
      $typing::typeof = $typing::typeofs[$#typing::typeofs];
      $typing::typedef = $typing::typedefs[$#typing::typedefs];
      $typing::code = $typing::codes[$#typing::codes];
   }
}

#####
# Creates a new type with specified base and sub.
# Note that sub corresponds to kind for built-in types.
# sub can be left out for a wild type.
# A third argument, print, can specify how the type should print.  Used for
# default types, double precision, etc.
#####
sub make_type {
  my ($base, $sub, $print) = @_;
  my $type = { 'base' => $base };
  $type->{'sub'} = $sub if $sub;
  $type->{'print'} = $print;
  return $type;
}

#####
# Creates a new complex type with specified types of "sides."
#####
sub make_complex_type {
  my ($type1, $type2) = @_;
  my ($base1, $base2) = ($type1->{'base'}, $type2->{'base'});
  die "Complex constant must have real and/or integer parts, but I found types $base1 and $base2"
    unless ($base1 eq 'integer' || $base1 eq 'real') &&
           ($base2 eq 'integer' || $base2 eq 'real');
  my $which;
  # From Metcalf and Reed's Fortran 90 Explained, if one of the types is an
  # integer then the kind of the complex is the kind of the other type.
  if ($base1 eq 'integer') {
    $which = $type2;
  } elsif ($base2 eq 'integer') {
    $which = $type1;
  } else {
    if ($type1->{'sub'} > $type2->{'sub'}) {
      $which = $type1;
    } else {
      $which = $type2;
    }
  }
  return {
    'base'    => 'complex',
    'sub'     => $which
  };
}

#####
# Creates a new character type with specified sub (kind) and len.
#####
sub make_character_type {
  my ($sub, $len) = @_;
  $sub = $typing::default_character_kind unless defined $sub;
  $sub = [ "%const", $typing::default_type{'integer'}, $sub ] unless ref $sub;
  $len = "1" unless defined $len;
  $len = [ "%const", $typing::default_type{'integer'}, $len ]
    unless ref $len || $len eq "*";
  return {
    'base' => 'character',
    'sub'  => $sub,
    'len'  => $len
  };
}

#####
# Returns true iff the given type was created to be the default of its kind.
# This has no meaning for compound types (hence it returns false).  For
# characters, there's a slight bug in that it will say that the type was
# created default even if you specify the default explicitly.  No biggie.
# Note that the defaultness is only for the KIND, not the LENGTH.
# 
# I could fix the above-mentioned problem by storing a 'default' entry just for
# the default types.  Then is_default_kind just translates to an exists test.
# This is much simpler and avoids the wierd checks for double precision numbers
# (0.0d0 ==> don't show a kind.  This is really "default").  This would be
# kinda nice but 'default' is probably the wrong word.
#####
sub is_default_kind {
   my ($type) = @_;

   if ($type->{'base'} eq "character") {
     my ($top, @rest) = @{$type->{'sub'}};
     return ($top eq "%const" && $rest[0] eq $typing::default_type{'integer'}
          && $rest[1] == $typing::default_character_kind);
   } else {
      return (exists $typing::default_type{$type->{'base'}} && $typing::default_type{$type->{'base'}} eq $type);
   }
}

#####
# Converts the given type to a string, written in Fortran 90 code.
# Only displays the kind if it was specified explicitly.  Slight bug:
# if you say character (kind=1) :: c, then it will print character :: c.
# (This is only for characters with default kind.  For other types with
# default kind explicitly specified, it is printed.)
#####
sub type_to_f90 {
  my ($type) = @_;

  # This covers the case where the kind is the default, except for characters.
  return $type->{'print'} if defined $type->{'print'};

  my $mods = "";
  if ($type->{'base'} eq "character") {
    if ($type->{'len'} eq "*") {
      $mods = "len=*";
    } elsif ($type->{'len'}->[0] ne "%const" ||
             $type->{'len'}->[1] != $typing::default_type{'integer'} ||
             $type->{'len'}->[2] ne "1") {
      $mods = "len=" . expr_to_f90 ($type->{'len'});
    }
    unless (is_default_kind ($type)) {
      $mods .= ", " unless $mods eq '';
      $mods .= "kind=" . expr_to_f90 ($type->{'sub'});
    }
  } elsif ($type->{'base'} eq "type") {
    $mods = "$type->{'sub'}";
  } else {
    $mods = "kind=" . expr_to_f90 ($type->{'sub'});
  }
  $mods = " ($mods)" unless $mods eq '';
  return $type->{'base'} . $mods;
}

#####
# Converts an expression right back to a string, doing "no" conversion (i.e.,
# output is in Fortran 90).  Optionally returns the precedence of the outmost
# operation in the expression (see $typing::precedence).
#####
sub expr_to_f90 {
  my ($exprptr) = @_;
  my ($op, @children) = @$exprptr;

  die "Unrecognized operation $op",%$op," (has no precedence?)"
    unless exists $typing::precedence->{$op};
  my $prec = $typing::precedence->{$op};

  my $answer;
  if ($op eq "%") {
    my ($struct, $elem) = @children;
    my ($s, $sprec) = expr_to_f90 ($struct);
    $s = "($s)" if $prec > $sprec;
    $answer = "$s%$elem";
  } elsif ($op eq "%var") {
    $answer = $children[0];
  } elsif ($op eq "%const") {
    my ($type, $val) = @children;
    if ($type->{'base'} eq 'complex') {
      if (!is_default_kind ($type->{'sub'})) {
        my ($k1, $k2) = ("", "");
        $k1 = "_$type->{'sub'}->{'sub'}" unless $val->[0] =~ /D[+-]?\d+$/i;
        $k2 = "_$type->{'sub'}->{'sub'}" unless $val->[1] =~ /D[+-]?\d+$/i;
        $answer = "($val->[0]$k1, $val->[1]$k2)";
      } else {
        $answer = "($val->[0], $val->[1])";
      }
    } elsif (is_default_kind ($type) || $val =~ /D[+-]?\d+$/i) {
      $answer = $val;
    } else {
      $answer = "${val}_$type->{'sub'}";
    }
  } elsif ($op eq "%array") {
    $answer = "(/ " . join (", ", map { (expr_to_f90 ($_))[0] } @children)
            . " /)";
  } elsif ($op eq "%colon") {
    my ($left, $right) = @children;
    $left = (expr_to_f90 ($left))[0] if $left ne '';
    $right = (expr_to_f90 ($right))[0] if $right ne '';
    $answer = $left . ":" . $right;  # : has ultimately low precedence
  } elsif ($op eq "%namedarg") {
    my ($left, $right) = @children;
    $answer = $left . " = " .
              (expr_to_f90 ($right))[0];  # = has ultimately low precedence
  } elsif ($op eq "%do") {
    my ($child, $var, @args) = @children;
    $answer = "(" . expr_to_f90 ($child) . ", " . $var . " = " .
              join (", ", map { (expr_to_f90 ($_))[0] } @args) . ")";
  } elsif ($op eq "%call") {
    ($op, @children) = @children;
    my ($s, $sprec) = expr_to_f90 ($op);
    $s = "($s)" if $prec > $sprec;
    $answer = "$s (" . join (", ", map ((expr_to_f90 ($_))[0], @children))
      . ")";
  } elsif (scalar @children == 1) {
    $op = substr ($op, 1) if substr ($op, 0, 1) eq 'u';
    my ($s, $sprec) = expr_to_f90 ($children[0]);
    $s = "($s)" if $prec > $sprec;
    $answer = "$op$s";
  } elsif (scalar @children == 2) {
    my ($s1, $sprec1) = expr_to_f90 ($children[0]);
    $s1 = "($s1)" if $prec > $sprec1;
    my ($s2, $sprec2) = expr_to_f90 ($children[1]);
    $s2 = "($s2)" if $prec > $sprec2;
    $answer = "$s1 $op $s2";
  } else {
    die "expr_to_f90: Unrecognized operation $op with " . (scalar @children) .
      " children";
  }

  if (wantarray) {
    return ($answer, $prec);
  } else {
    return $answer;
  }
}

#####
# Computes the type of the given expression (which is passed by reference).
# Returns a reference to the actual type.
#####
sub expr_type {
   my ($exprptr) = @_;
   my ($op, @children) = @$exprptr;

   if ($op eq "%") {
      my ($struct, $elem) = @children;
      my ($type) = expr_type ($struct);
      die "expr_type: \%$elem failed: left part is not a compound type" unless $type->{'base'} eq "type";
      my ($typedef) = $typing::typedef->{$type->{'sub'}};
      my ($elemtype) = $typedef->{$elem};
      die "expr_type: \%$elem failed: left part does not include $elem" unless $elemtype;
      return $elemtype;
   } elsif ($op eq "%var") {
      my ($var) = @children;
      my ($vartype) = $typing::typeof->{$var};
      die "expr_type: Variable $var undefined" unless $vartype;
      return $vartype;
   } elsif ($op eq "%const") {
      my ($type, $val) = @children;
      return $type;
   } elsif ($op eq "%array") {
      # HERE
   } elsif ($op eq "%colon") {
      my ($string, $left, $right) = @children;
      my ($stringtype) = expr_type ($string);
      die "expr_type: colon notation for non-character string" if $stringtype->{'base'} ne "character";
      die "expr_type: colon notation for character array" if $stringtype->{'dimension'};
      return typing::make_character_type ($stringtype->{'sub'}, "*");
   } elsif ($op eq "%call") {
      ($op, @children) = @children;
      my ($subop, @subchildren) = @$op;
      if ($subop eq "%var") {
         ($op) = @subchildren;
         # Fall through: we allow overloaded function name in this special case.
      } else {
         # Function call without overloading or an array reference.
         my ($optype) = expr_type ($op);

         if ($optype->{'dimension'}) {  # array reference
            return make_type ($optype->{'base'}, $optype->{'sub'});
         } else {
            die "expr_type: Array/function call for something that is neither" unless $optype->{'base'} eq "interface";
            # HERE function call without overloading.
         }
      }
   }

   my ($opcodes) = $typing::code->{$op};
   die "Operation/function $op undefined" unless $opcodes;
   my (@childtypes) = ();
   my ($child);
   foreach $child (@children) {
      print "childtypes was: @childtypes\n";
      print "type of $child is ", expr_type ($child), "\n";
      push @childtypes, expr_type ($child);
      print "childtypes is now: @childtypes\n";
   }
   my ($opcode);
   foreach $opcode (@$opcodes) {
      print "children: @children\n";
      print "childtypes: @childtypes\n";
      if (typing::subtypes_list (\@childtypes, $opcode->{'parms'})) {
         my ($parm);
         my ($ret) = $opcode->{'return'};
         if ($ret->{'base'} eq "character" && ! $ret->{'len'}) {
            $ret->{'len'} = 0;
find_len:
            foreach $parm (@$opcode->{'parms'}) {
               if ($parm->{'base'} eq $ret->{'base'}) {
                  if ($parm->{'len'} eq "*") {
                     $ret->{'len'} = "*";
                     last find_len;
                  } else {
                     $ret->{'len'} += $parm->{'len'};
                  }
               }
            }
         }
         if ($ret->{'sub'}) {
            return $ret;
         } else {
            # Make intrinsic type's kind: look for all parameters with the same
            # base type, and use the maximum kind out of those.
            my ($maxkind) = -1;
            foreach $parm (@$opcode->{'parms'}) {
               if ($parm->{'base'} eq $ret->{'base'}) {
                  $maxkind = $parm->{'sub'} if $maxkind < $parm->{'sub'};
               }
            }
            die "expr_type: Internal error caused by new_scope" if $maxkind < 0;
            return { %$ret, 'sub' => $maxkind };
         }
      }
   }
   die "Operation/function $op defined but not for this (these) type(s)";
}

#####
# Returns if first type is a subtype of the second type.
# This currently only supports intrinsic types (integer*4 subtypes integer*?).
#####
sub subtypes {
   my ($t1, $t2) = @_;
   return 0 if $t1->{'base'} ne $t2->{'base'};
   if ($t1->{'base'} eq "type") {
      return 0 if $t1->{'sub'} eq $t2->{'sub'};
   } else {
      if ($t1->{'base'} eq "character") {
         if ($t1->{'len'}) {
            return 0 unless $t1->{'len'};
            return 0 if $t2->{'len'} != $t1->{'len'};
         }
      }
      if ($t1->{'base'} eq "interface") {
         # HERE fill this in when I do function types ("interface").
      }
      if ($t1->{'sub'}) {
         return 0 unless $t1->{'sub'};
         return 0 if $t2->{'sub'} ne $t1->{'sub'};
      }
   }
   return 1;
}

#####
# Returns if first type is a subtype of the second type, where the first
# and second type are (conceptually) tuples.  That is, the lengths must be
# equal, and each element must subtype the corresponding element.
# The lists are passed as references.
#####
sub subtypes_list {
   my ($l1ptr, $l2ptr) = @_;
   my (@l1) = @$l1ptr;
   my (@l2) = @$l2ptr;
   return 0 if $#l1 != $#l2;

   print "l1 is: @l1\n";
   print "l2 is: @l2\n";

   my ($i);
   for ($i = 0; $i <= $#l1; $i++) {
      print "calling subtypes with $l1[$i] and $l2[$i]\n";
      return 0 unless typing::subtypes ($l1[$i], $l2[$i]);
   }
   return 1;
}
