package stmts;

use strict;

require "expr_parse.pl";
require "typing.pl";
require "utils.pl";

#########################################################################
# PUBLIC GLOBALS

# Set to a reference to a routine to take !! comments if !! comments are
# to be caught.
$stmts::bangbang = "";

# Set to a reference to a routine to return accumulated comments if !! comments
# are caught.  You should reset them after each time you call read_line or
# read_stmt.
$stmts::comments = "";

# Set this to disable warnings.  Don't use this for a compiler!  Suitable for
# something like f90doc though.  This shouldn't be used once stmts supports
# all Fortran 90 statements and attributes; until then, it's pretty much
# needed; after then, it should be removed.
$stmts::diable_warns = 0;

# Set this to use fixed-form Fortran, like good old Fortran 77.
$stmts::fixed_form = 0;

#########################################################################
# PRIVATE GLOBALS

# A "left-over" piece of a statement is stored here when semi-colons are
# encountered.
$stmts::leftover = "";

# Number of opened files.
$stmts::nfile = 0;

# List of string's values.
@stmts::strings = ();

# List of structure pointers that we're currently nested in.
# topnest stores the top of the stack.
@stmts::nesting = ();
$stmts::topnest = undef;

# List of structure pointers that we're currently nested in, but for a
# specified type.
%stmts::nesting_by = ();

#########################################################################
# ROUTINES

#####
# Reads an entire file, and returns all the top-level structures found.
# If specified, a given function will be called after every statement
# (usually this is for resetting !! comments and such).
#####
sub read_file {
  my ($filename, $every_stmt) = @_;
  stmts::open_file ($filename);

  my ($stmt, $struct, @rval);
  my @toplevel = ();
  while ((@rval = stmts::read_stmt ()) [0]) {
    push @toplevel, $rval[1] if !defined $stmts::topnest && ref $rval[1];
    &$every_stmt () if defined $every_stmt;
  }

  return @toplevel;
}

#####
# Starts reading the specified filename.
#####
sub open_file {
   my ($filename) = @_;
   $stmts::FILE = "";

   open IN, $filename
     or die "Couldn't open $filename";
   $stmts::{'FILE' . $stmts::nfile} = $stmts::{'IN'};
}

#####
# Cleans up from reading the current file.
# This is automatically called by read_line, so most don't have to worry
# about it.
# Returns false if there are no files left.
#####
sub close_file {
   close IN;
   $stmts::nfile--;
   if ($stmts::nfile > 0) {
      # CHECK--does this still do the desired thing, in light of open_file?
      $stmts::{'IN'} = $stmts::{'FILE' . $stmts::nfile};
      return 1;
   } else {
      # Clean up strings.
      @stmts::strings = ();
      return 0;
   }
}

#####
# Reads a line of Fortran 90 doing whatever it takes.  This may involve
# reading multiple lines from the current file, walking into files, etc.
# INCLUDE is parsed at this level.
# Note that the returned string may have various cases (lc isn't called).
#####
sub read_line {

ALLOVERAGAIN:
  my $line;
  if ($stmts::leftover ne '') {
    $line = $stmts::leftover;
    $stmts::leftover = '';
  } else {
    $line = <IN>;
    until (defined $line) {
      return "" unless close_file ();
      $line = <IN>;
    }
    chomp $line;

    substr ($line, 0, 1) = '!' if $stmts::fixed_form && $line =~ /^\S/;
  }

  # This is used for fixed-form continuations.
  my $lastlen = length $line;

  my $continue = 0;

  while (1) {
    # Grab doubled comments (!!) if requested.
    if ($stmts::bangbang && $line =~ /^([^"'!]|('[^']*')|("[^"]*"))*(!!.*)$/) {
      $line = substr ($line, 0, length ($line) - length ($4));
      &$stmts::bangbang ($4);
    }

    # Delete comments.
    elsif ($line =~ /^([^"'!]|(\'[^']*')|("[^"]*"))*(!.*)$/) {
      $line = substr ($line, 0, length ($line) - length ($4));
    }

    # Fixed-form continuations.
    if ($stmts::fixed_form) {

      # Check next line for continuation mark.
      $stmts::leftover = <IN>;
      $stmts::leftover = '' unless defined $stmts::leftover;
      chomp $stmts::leftover;
      substr ($stmts::leftover, 0, 1) = '!' if $stmts::leftover =~ /^\S/;
      if ($stmts::leftover =~ /^\s....\S/) {

        # Pad previous line with spaces if it had less than 72 characters.
        $line .= ' ' x (72-$lastlen) if $lastlen < 72;

        # Add next (continuation) line to the line.
        $line .= substr ($stmts::leftover, 6);
        $lastlen = length $stmts::leftover;
        
        # Continue on to check the next line.
        $stmts::leftover = '';
        next;
      }
      
    # Free-form continuations.
    } elsif ($continue || $line =~ /&\s*$/) {
      $line = $` if $line =~ /&\s*$/;
      my $rest = <IN>;
      chomp $rest;
      $rest = $' if $rest =~ /^\s*&/;
      $line = "$line$rest";
      # Blank lines don't stop the continuation.
      $continue = ($rest =~ /^\s*(?:!.*)?$/);
      next;
    }

    last;
  }

  # Semicolons.
  if ($line =~ /^([^;]*);(.*)$/) {
    $line = $1;
    if ($stmts::leftover eq '') {
      $stmts::leftover = $2;
    } else {
      $stmts::leftover .= ";$2";
    }
  }

  # Replace strings to avoid confusion.
  my @quotes;
  while ($line =~ / " ([^"]|"")* " | ' ([^']|'')* ' /xg) {
    push @quotes, [length $`, length $&, $&];
  }
  for my $quote (reverse @quotes) {
    ## Process in reverse order so that $start is preserved despite replacement
    my ($start, $length, $string) = @$quote;
    push @stmts::strings, $string;
    substr ($line, $start, $length) = "\'" . $#stmts::strings . "\'";
  }

  # Get rid of spaces on either end.
  $line = utils::trim ($line);

  goto ALLOVERAGAIN if $line eq '';

  #print "read line `$line'\n";

  return $line;
}

#####
# Returns the physical value for the given string number.
#####
sub get_string {
   my ($n) = @_;
   return $stmts::strings[$n];
}

#####
# Reads a Fortran 90 statement from the current input.
# Checks for proper nesting, etc., and keeps tracks of what's in what.
# Possible results:
#    ('?', $the_line)
#    ('program', \%structure)
#    ('endprogram', \%structure)
#    ('module', \%structure)
#    ('endmodule', \%structure)
#    ('subroutine', \%structure)
#    ('endsubroutine', \%structure)
#    ('function', \%structure)
#    ('endfunction', \%structure)
#    ('program', \%structure)
#    ('endprogram', \%structure)
#    ('type', \%structure)
#    ('endtype', \%structure)
#    ('interface', \%structure)
#    ('endinterface', \%structure)
#    ('var', \%struct1, \%struct2, ...)
#    ('contains', \%parent)
#    ('public', $name1, $name2, ...)          empty means global default
#    ('private', $name1, $name2, ...)         empty means global default
#    ('optional', $name1, $name2, ...)
#    ('call', $arg1, $arg2, ...)              currently args are unparsed
#####
sub read_stmt {
   my ($line) = read_line ();
   if (! $line) {
      die "File ended while still nested" if @stmts::nesting;
      return ("", "");
   }

   # MODULE PROCEDURE (must be before module)
   if ($line =~ /^module\s+procedure\s+(\w.*)$/i) {
      die "module procedure outside of interface block" unless defined $stmts::topnest && $stmts::topnest->{'type'} eq "interface" && $stmts::topnest->{'name'} ne "";
      my (@list) = split (/\s*,\s*/, utils::trim ($1));
      my ($p);
      foreach $p (@list) {
         die "Invalid module procedure `$p'" unless $p =~ /^\w+$/;
         new_struct ({
            'type'   => "mprocedure",
            'name'   => $p,
            hashed_comments ()
         });
      }
      return ("mprocedure", @list);
   }

   # MODULE/PROGRAM
   elsif ($line =~ /^(module|program)(?:\s+(\w+))?$/i) {
      die "$1 begun not at top level" if defined $stmts::topnest;
      return new_nest ({
         'type' => lc $1,
         'name' => (defined $2 ? $2 : ''),
         hashed_comments ()
      });
   }

   # END MODULE/SUBROUTINE/FUNCTION/PROGRAM/TYPE/INTERFACE, or general END
   elsif ($line =~ /^end\s*(?:(module|subroutine|function|program|type|interface)(?:\s+(\w+))?)?$/i) {
      die "END statement outside of any nesting" unless defined $stmts::topnest;
      my $top = $stmts::topnest;

      # We do some special "fixing up" for modules, which resolves named
      # references (module procedures) and computes publicity.
      #
      # Note that end_nest will ensure that the type of thing ended matches
      # the thing the user says it is ending, so we don't have to worry about
      # that.
      if ($top->{'type'} eq "module") {

        # Set publicity (visibility) of objects within the module.

        # First, the explicitly set ones.
        my $name;
        foreach $name (@{$top->{'publiclist'}}) {
          do_attrib ($name, "vis", 'public', "visibility");
        }
        foreach $name (@{$top->{'privatelist'}}) {
          do_attrib ($name, "vis", 'private', "visibility");
        }

        # Second, the globally set ones (those obeying the default).
        my $obj;
        $top->{'defaultvis'} = "public" unless exists $top->{'defaultvis'};
        foreach $obj (@{$top->{'ocontains'}}) {
          $obj->{'vis'} = $top->{'defaultvis'} unless exists $obj->{'vis'};
        }

        # Traverse (arbitrarily deeply) nested structures.
        sub traverse {
          my ($node) = @_;
          my $top = $stmts::topnest;   # HAVE NO IDEA WHY THIS IS NEEDED
          
          # Graduate nested MODULE PROCEDURE (mprocedure) to point to the
          # appropriate thing (either a function or a subroutine with that
          # name).
          if ($node->{'type'} eq "mprocedure") {
            die "Couldn't find module procedure $node->{'name'} (nothing with that name in module $top->{'name'})"
              unless exists $top->{'contains'}->{lc $node->{'name'}};
            
            my ($possibles) =
              $top->{'contains'}->{lc $node->{'name'}};
            die "Couldn't find module procedure $node->{'name'} in module $top->{'name'} (wrong type)"
              if !exists $possibles->{'subroutine'}
              && !exists $possibles->{'function'};
            die "Found both a subroutine and function to match module procedure $node->{'name'} in module $top->{'name'}"
              if exists $possibles->{'subroutine'}
              && exists $possibles->{'function'};
            
            if (exists $possibles->{'subroutine'}) {
              $node->{'bind'} = $possibles->{'subroutine'};
            } else {
              $node->{'bind'} = $possibles->{'function'};
            }
          }

          # Recurse.
          map { traverse ($_) } @{$node->{'ocontains'}}
          if exists $node->{'ocontains'};
        }
        map { traverse ($_) } @{$top->{'ocontains'}};
      }

      my @return_val = end_nest ($1, $2);

      # Subroutines and functions in interface blocks must be noted at the
      # top level.  We do this with "interface" structures with the names
      # of the actual contained routines (unless this is already the
      # case).  Make sense?
      if ($top->{'type'} eq "interface" && $top->{'name'} eq "") {
          my $sub;
          foreach $sub (@{$top->{'ocontains'}}) {
              next if $sub->{'name'} eq $top->{'name'} ||
                      $sub->{'type'} eq "mprocedure";

              my %copy = %$top;
              $copy{'name'} = $sub->{'name'};
              new_nest (\%copy);
              my $old_within = $sub->{'within'};
              new_struct ($sub);
              $sub->{'within'} = $old_within;
              end_nest ('interface', $sub->{'name'});
          }
      }

      return @return_val;
   }

   # SUBROUTINE/FUNCTION
   elsif ($line =~ /^(?:(.+?)\s+)?(subroutine|function)\s+(\w+)\s*(\([^()]*\))?(?:\s*result\s*\(\s*(\w+)\s*\))?$/i) {
      my ($type, $name, $parmstr, $rtype, $result) =
         (lc $2, $3,    $4,       $1,     $5);

      die "Start of $type $name before `contains' section of $stmts::topnest->{'type'} $stmts::topnest->{'name'}"
          if defined $stmts::topnest && ! $stmts::topnest->{'incontains'} &&
             $stmts::topnest->{'type'} ne "interface";
      if (exists $stmts::nesting_by{'subroutine'} ||
          exists $stmts::nesting_by{'function'}) {
         my $n = 0;
         $n += scalar @{$stmts::nesting_by{'subroutine'}}
            if exists $stmts::nesting_by{'subroutine'};
         $n += scalar @{$stmts::nesting_by{'function'}}
            if exists $stmts::nesting_by{'function'};
#FIXME  #die "Routine nested in routine nested in routine" if $n > 1;
      }

      $parmstr = "()" unless defined $parmstr;
      $parmstr = utils::trim (substr ($parmstr, 1, length ($parmstr) - 2));
      my (@parms);
      if ($parmstr) {
         @parms = split (/\s*,\s*/, $parmstr);
         my ($parm);
         foreach $parm (@parms) {
            die "Parameter `$parm' is not just a word or *"
              unless $parm =~ /^\w+|\*$/;
            ## * as a final argument allows the calling to specify a statement
            ## to jump as an alternative return address.  (Legacy Fortran!)
            ## Thanks to Art Olin for this info.
         }
      } else {
         @parms = ();
      }

      my $struct = {
         'type'      => $type,
         'name'      => $name,
         'parms'     => \@parms,
         hashed_comments ()
      };
      new_nest ($struct);

      $struct->{'result'} = $result if defined $result;

      $rtype = "" unless defined $rtype;
      while ($rtype =~ /(?:^|\s+)(recursive|pure|elemental)$/i ||
             $rtype =~ /^(recursive|pure|elemental)(?:\s+|$)/i) {
        $rtype = $` . $'; # actually whichever is not blank
        $struct->{lc $1} = 1;
      }
      if ($rtype ne '') {
        $struct->{'rtype'} = parse_type ($rtype);
        new_struct ({
          'type'        => 'var',
          'name'        => (defined $result ? $result : $name),
          'vartype'     => $struct->{'rtype'},
          'comments'    => ''
        });
      }

      return ($type, $struct);
   }

   # TYPE definition (must go before variable declarations)
   elsif ($line =~ /^type(?:\s+|\s*(,.*)?::\s*)(\w+)$/i) {
     my $struct = new_nest ({
       'type' => 'type',
       'name' => $2,
       hashed_comments ()
     });
     if (defined $1) {
       my $attrib = utils::trim (substr ($1, 1));
       if ($attrib =~ /^(public|private)$/i) {
         $struct->{'vis'} = lc $attrib;
       } elsif ($attrib) {
         warn "Invalid attribute `$attrib' for derived-type declaration--should be just public or private";
       }
     }
     return $struct;
   }

   # INTERFACE block (for overloading) or statement (for definition of external)
   elsif ($line =~ /^interface(?:\s+(\S.+))?$/i) {
       return new_nest ({
           'type' => 'interface',
           'name' => (defined $1 ? $1 : ""),
           hashed_comments ()
       });
   }

   # CONTAINS
   elsif ($line =~ /^contains$/i) {
      die "`contains' found at top level" unless defined $stmts::topnest;
      die "`contains' found in $stmts::topnest->{'type'} $stmts::topnest->{'name'}" unless exists $stmts::topnest->{'incontains'};
      die "Multiple `contains' found in same scope"
         if $stmts::topnest->{'incontains'};
      die "`contains' found in interface definition"
         if $stmts::topnest->{'interface'};
      $stmts::topnest->{'incontains'} = 1;
      return ("contains", $stmts::topnest);
   }

   # PUBLIC/PRIVATE/SEQUENCE
   elsif ($line =~ /^(public|private|sequence)(?=\s+[^=(]|::|$)(\s*::\s*)?/i) {
     my ($what, $rest) = (lc $1, $');

     if (defined $stmts::topnest && $stmts::topnest->{'type'} eq "type") {
       die "public statement not allowed in a type declaration"
         if $what eq 'public';
       die "$1 cannot be qualified inside type declaration" if $rest;
       $stmts::topnest->{$what . 'type'} = 1;
       return ($what);
     } else {
       die "sequence statement only allowed immediately inside type declaration"
         if $1 eq 'sequence';

       die "$1 statement not immediately inside a module or type declaration"
         unless defined $stmts::topnest && $stmts::topnest->{'type'} eq "module";
       if ($rest eq "") {  # Unqualified
         die "Unqualified $what in addition to unqualified " .
           $stmts::topnest->{'defaultvis'}
         if exists $stmts::topnest->{'defaultvis'};
         $stmts::topnest->{'defaultvis'} = $what;
         return ($what);
         
       } else {  # Qualified
         my @namelist = map {
           die "Invalid name `$_' specified in $what statement"
             unless /^\s*(\w+)(?:\s*(\([^()]+\)))?\s*$/i;
           $1 . (defined $2 ? $2 : "");
         } (split ',', $rest);
         push @{$stmts::topnest->{"${what}list"}}, @namelist;
         return ($what, @namelist);
       }
     }
   }

    # OPTIONAL
    elsif ($line =~ /^optional(\s+|\s*::\s*)((\w|\s|,)+)$/i) {
        my $name;
        my @namelist = split (/\s*,\s*/, utils::trim ($2));
        foreach $name (@namelist) {
            do_attrib ($name, "optional", 1, "optional attribute");
        }
        return ('optional', @namelist);
    }

   # Variable declarations
   elsif ($line =~ /^(integer|real|double\s*precision|character|complex|logical|type)\s*(\(|\s\w|[:,*])/i) {
      my ($vartype, $rest) = parse_part_as_type ($line);
      my (@attribs, @right);
      if ($rest =~ /^(.*)\:\:(.*)/) {
         my ($a, $b) = ($1, $2);
         @attribs = map (( utils::trim ($_) ), utils::balsplit (",", $a));
         @right = map (( utils::trim ($_) ), utils::balsplit (",", $b));
      } else {
         @attribs = ();
         @right = map (( &utils::trim ($_) ), utils::balsplit (",", $rest));
      }
      my ($r, @structs);
      foreach $r (@right) {
          my ($rl, $rassign) = &utils::balsplit ("=", $r);
          my ($rll, $starpart) = &utils::balsplit ("*", $rl);
          if (defined $starpart) {
            die "Sorry, I don't support 'character var*kind' yet; use 'character*kind var' instead";
          }
          $rll =~ /^ (\w+) (\s* \(.*\))? \s* $/x
              or die "Invalid variable declaration `$rll'";
          my ($name, $dimension) = ($1, $2);
          my ($initop, $initial);
          if (defined $rassign) {
            # implicit lead =
            $rassign =~ /^ (>?) \s* (.*) $/x
              or die "Invalid variable initialization `= $rassign'";
            ($initop, $initial) = ("=" . $1, $2);
          }

          my $struct;
          $struct = {
              'type'        => 'var',
              'name'        => $name,
              'vartype'     => $vartype,
              hashed_comments ()
          };
          if (defined $initial) {
            $struct->{'initop'} = $initop;
            $struct->{'initial'} = expr_parse::parse_expr ($initial);
          }
          new_struct ($struct);
          push @structs, $struct;

          my @attribs_copy = @attribs;
          push @attribs_copy, "dimension $dimension" if defined $dimension;

          my ($attrib, @tempattribs);
          foreach $attrib (@attribs_copy) {
              if ($attrib =~ /^(public|private)$/i) {
                  $attrib = lc $attrib;
                  $struct->{'vis'} = $attrib;
              } elsif ($attrib =~ /^optional$/i) {
                  $attrib = lc $attrib;
                  $struct->{$attrib} = 1;
              } elsif ($attrib) {
                  warn "Unrecognized attribute `$attrib'"
                      unless $stmts::disable_warns;
                  push @tempattribs, $attrib;
              }
          }

          $struct->{'tempattribs'} = \@tempattribs;
      }

      return ('var', @structs);
   }

   # USE
   elsif ($line =~ /^use\s+(\w+)($|,\s*)/i) {
      die "`use' found at top level" unless defined $stmts::topnest;
      die "`use' found in $stmts::topnest->{'type'} $stmts::topnest->{'name'}" unless exists $stmts::topnest->{'uses'};
      my $extra = length $' ? $' : undef;
      push @{$stmts::topnest->{'uses'}}, [$1, $extra];

      return ('use', $1, $extra);
   }
   
   # CALL or IF (...) CALL [hack--xxx]
   elsif ($line =~ /^(?:if\s*\(.*\)\s*)?call\s+(\w+)\s*(?:\(\s*(.*?)\s*\))?$/i) {
      die "`call' found at top level" unless defined $stmts::topnest;
      die "`call' found in $stmts::topnest->{'type'} $stmts::topnest->{'name'}" unless exists $stmts::topnest->{'calls'};
      $stmts::topnest->{'calls'}->{$1} = 1;
      my @args = ();
      @args = split /\s*,\s*/, $2 if defined $2;
      return ('call', @args);
   }
   
   # Unrecognized statement
   else {
      if ($line =~ /^\w+/) {
         warn "Unrecognized statement beginning with word $&" unless $stmts::disable_warns;
      } else {
         warn "Unrecognized statement" unless $stmts::disable_warns;
      }
      return ('?', $line);
   }
}

#####
# Returns a list that would fit right into a hash table you're making.  If
# there are no comments, returns the empty list.  The entry is called
# 'comments'.
#####
sub hashed_comments {
   if ($stmts::comments) {
      return ( 'comments', &$stmts::comments () );
   } else {
      return ();
   }
}

#####
# Makes note of a new structure.  Called by new_nest, for example.
#####
sub new_struct {
   my ($struct) = @_;
   my $type = $struct->{'type'};

   die "Basic structure must be found at a nesting level"
     unless defined $stmts::topnest;

   if (exists ($stmts::topnest->{'contains'}->{lc $struct->{'name'}})) {
      die "Redefinition of $type $struct->{'name'} in $stmts::topnest->{'type'} $stmts::topnest->{'name'}"
         if exists ($stmts::topnest->{'contains'}->{lc $struct->{'name'}}->{$type});
      $stmts::topnest->{'contains'}->{lc $struct->{'name'}}->{$type} = $struct;
   } else {
      $stmts::topnest->{'contains'}->{lc $struct->{'name'}} =
         { $type => $struct };
   }
   push @{$stmts::topnest->{'ocontains'}}, $struct;
   $struct->{'within'} = $stmts::topnest;
}

#####
# Starts a new nesting level represented by the given structure.  The
# structure must define the 'type' and 'name' entries.  You should not
# define the 'contains' or 'defaultvis' entry.
#####
sub new_nest {
   my ($struct) = @_;
   my ($type) = $struct->{'type'};

   $struct->{'contains'} = { };
   $struct->{'ocontains'} = [ ];

   # Program unit
   if ($type eq "subroutine" || $type eq "function" || $type eq "module" || $type eq "program") {
     $struct->{'incontains'} = 0;
     $struct->{'uses'} = [ ];
     $struct->{'interface'} = 0 if $type eq "subroutine" || $type eq "function";
   }

   # Program unit with code
   if ($type eq "subroutine" || $type eq "function" || $type eq "program") {
     $struct->{'calls'} = { };
   }

   if (defined $stmts::topnest) {
      my ($toptype) = $stmts::topnest->{'type'};
      if ($toptype eq "interface" && ($struct->{'type'} eq "subroutine" || $struct->{'type'} eq "function")) {
         $struct->{'interface'} = 1;
      } else {
         die "Nesting in $toptype not allowed" unless $toptype eq "subroutine" || $toptype eq "function" || $toptype eq "module" || $toptype eq "program";
      }
      new_struct ($struct) unless $struct->{'name'} eq "";
   }
   push @stmts::nesting, $struct;
   if (exists ($stmts::nesting_by{$type})) {
      push @{$stmts::nesting_by{$type}}, $struct;
   } else {
      $stmts::nesting_by{$type} = [ $struct ];
   }
   $stmts::topnest = $struct;
   return ( $type, $struct );
}

#####
# Ends the current nesting level.  Optionally, you can pass the 'type' that
# it's supposed to be as the first argument.  Optionally, you can pass the
# 'name' it should have after that (as the second argument).
#####
sub end_nest {
  my ($type, $name) = @_;
  $type = lc $type if defined $type;
  unless (defined $stmts::topnest) {
    if (defined $name && defined $type) {
      die "Ended $type $name at top level";
    } elsif (defined $type) {
      die "Ended unnamed $type at top level";
    } else {
      die "END statement at top level";
    }
  }
  my ($struct) = pop @stmts::nesting;
  die "Ended $type while in $struct->{'type'} $struct->{'name'}"
    if defined $type && $type ne $struct->{'type'};
  die "Ended $name while in $struct->{'type'} $struct->{'name'}"
    if defined $name && $name !~ /^\Q$struct->{'name'}\E$/i;
  if (@stmts::nesting) {
    $stmts::topnest = $stmts::nesting[$#stmts::nesting];
  } else {
    $stmts::topnest = undef;
  }
  pop @{$stmts::nesting_by{$struct->{'type'}}};
  return ( "end" . (defined $type ? $type : ''), $struct );
}

#####
# Parses the basic type that prefixes the given string.
# Returns (parsed type, string portion remaining).
#####
sub parse_part_as_type {
  my ($str) = @_;

  $str =~ /^integer|real|double\s*precision|character|complex|logical|type/i
    or die "parse_part_as_type: Invalid input `$str'";
  my ($base, $rest) = ($&, $');

  my $level = 0;
  ## Wait till we are outside of all parens and see a letter, colon, or comma.
  while ($rest =~ /[()a-zA-Z_:,]/g) {
    if ($& eq '(') {
      $level++;
    } elsif ($& eq ')') {
      $level--;
      die "Unbalanced parens (too many )'s)" if $level < 0;
    } elsif ($level == 0) {
      return (parse_type ($base . $`), $& . $');
    }
  }
  
  die "Couldn't split into type and rest for `$str'";

# Some old, presumably less-efficient code:
#  my ($level, $len) = (0, length ($str));
#  my ($i, $c);
#  for ($i = length ($&); $i < $len; $i++) {
#    $c = substr ($str, $i, 1);
#    if ($c eq "(") {
#      $level++;
#    } elsif ($c eq ")") {
#      $level--;
#      die "Unbalanced parens (too many )'s)" if $level < 0;
#    } elsif ($level == 0 && $c =~ /^\w|:|,$/) {
#      last;
#    }
#  }
#  return (parse_type (substr ($str, 0, $i)), substr ($str, $i));
}

#####
# Parses a basic type, creating a type structure for it:
#     integer [( [kind=] kind_val )]
#     real [( [kind=] kind_val )]
#     double precision                  (no kind is allowed)
#     complex [( [kind=] kind_val )]
#     character [( char_stuff )]
#     logical [( [kind=] kind_val )]
#     type (type_name)
#
# integer*number, real*number, complex*number, and logical*number are also
# supported as nonstandard Fortran extensions for kind specification.
# "number" can either be a direct integer or an expression in parentheses.
# 
# char_stuff is empty or (stuff), where stuff is one of:
#     len_val [, [kind=] kind_val]
#     kind=kind_val [, [len=] len_val]
#     len=len_val [, kind=kind_val]
# kind_val and len_val are expressions; len_val can also be just `*'.
# 
# The length can also be specified using the nonstandard Fortran extension
# character*number.  If number is `*', it must be in parentheses (indeed,
# any expression other than a number must be in parentheses).
#####
sub parse_type {
  my ($str) = @_;

  # print "Parsing type: $str\n";

  $str = utils::trim ($str);
  $str =~ /^(integer|real|double\s*precision|complex|character|logical|type)
    \s* (?: \( (.*) \) | \* \s* (\d+ | \(.*\)) )?$/ix
    or die "Invalid type `$str'";
  my $base = lc $1;

  if ($base =~ /^double\s*precision$/) {
    die "double precision cannot have kind specification"
      if defined $2 || defined $3;
    return $typing::double_precision;
  }

  if (defined $2 || defined $3) {
    my $star = defined $3;
    my $args = utils::trim ($star ? $3 : $2);

    if ($base eq 'type') {
      die "type$args invalid--use type($args)" if $star;
      die "type(w) for non-word w" unless $args =~ /^\w+$/;
      return typing::make_type ($base, $args);
    } elsif ($base eq 'character') {
      my ($kind, $len, $rest);
      if ($star) {
        if ($args =~ /^\(\s*\*\s*\)$/) {
          $len = '*';
        } else {
          $len = expr_parse::parse_expr ($args);
        }
      } elsif ($args =~ /^kind\s*=\s*/i) {
        $args = substr ($args, length ($&));
        ($kind, $rest) = expr_parse::parse_part_as_expr ($args);
        if (defined $rest) {
          $rest = utils::trim ($rest);
          $rest =~ s/^len\s*=\s*//i;
          $len = ($rest eq '*' ? '*' : expr_parse::parse_expr ($rest));
        }
      } elsif ($args =~ /^len\s*=\s*/i) {
        $args = substr ($args, length ($&));
        if (substr ($args, 0, 1) eq '*') {
          $len = '*';
          $rest = $args;
          $rest =~ s/^\*\s*,// or $rest = undef;
        } else {
          ($len, $rest) = expr_parse::parse_part_as_expr ($args);
        }
        if (defined $rest) {
          $rest = utils::trim ($rest);
          $rest =~ /^kind\s*=\s*/
            or die "kind= specifier needed when len= specifier is given";
          $rest = substr ($rest, length ($&));
          $kind = expr_parse::parse_expr ($rest);
        }
      } else {  # len
        if (substr ($args, 0, 1) eq '*') {
          $len = "*";
          $rest = $args;
          $rest =~ s/^\*\s*,// or $rest = undef;
        } else {
          ($len, $rest) = expr_parse::parse_part_as_expr ($args);
        }
        if (defined $rest) {
          $rest = utils::trim ($rest);
          $rest = substr ($rest, length ($&)) if $rest =~ /^kind\s*=\s*/i;
          $kind = expr_parse::parse_expr ($rest);
        }
      }
      return typing::make_character_type ($kind, $len);
    } else {
      $args =~ s/^kind\s*=\s*//i unless $star;
      return typing::make_type ($base, expr_parse::parse_expr ($args));
    }
  } else {
    die "type without (type-name) after it" if $base eq 'type';
    die "No default type for `$base'"
      unless exists $typing::default_type{$base};
    return $typing::default_type{$base};
  }
}

sub do_attrib {
    my ($name, $attrib, $val, $attribname) = @_;
    my ($struct);
    foreach $struct (values %{$stmts::topnest->{'contains'}->{lc $name}}) {
        die "Redefining $attribname of $struct->{'type'} $name from " .
            "$struct->{$attrib} to $val" if exists $struct->{$attrib};
        $struct->{$attrib} = $val;
    }
}

1;
