#!/usr/bin/perl -w

$CKFILE    = "";
$TRANFILE  = "";
$THERMFILE = "";

@species   = ();
@elements  = ();
@reactions = ();
@thermo    = ();

$LowestT  = "";
$CommonT  = "";
$HighestT = "";

%Weight = ('H'  => 1.00797,
           'HE' => 4.0026,
           'LI' => 6.939,
           'BE' => 9.01220,
           'B'  => 10.811,
           'C'  => 12.01115,
           'N'  => 14.0067,
           'O'  => 15.9994,
           'F'  => 18.9984,
           'NE' => 20.183,
           'NA' => 22.9898,
           'MG' => 24.312,
           'AL' => 26.9815,
           'SI' => 28.086,
           'P'  => 30.9738,
           'S'  => 32.064,
           'CL' => 35.453,
           'AR' => 39.948,
           'K'  => 39.102,
           'CA' => 40.08,
           'SC' => 44.956,
           'TI' => 47.9,
           'V'  => 50.942,
           'CR' => 51.996,
           'MN' => 54.938,
           'FE' => 55.847,
           'CO' => 58.9332,
           'NI' => 58.71,
           'CU' => 63.54,
           'ZN' => 65.37,
           'GA' => 69.72,
           'GE' => 72.59,
           'AS' => 74.9216,
           'SE' => 78.96,
           'BR' => 79.9009,
           'KR' => 83.8,
           'RB' => 85.47,
           'SR' => 87.62,
           'Y'  => 88.905,
           'ZR' => 91.22,
           'NB' => 92.9064,
           'MO' => 95.94,
           'TC' => 98,
           'RU' => 101.07,
           'RH' => 102.906,
           'PD' => 106.42,
           'AG' => 107.868,
           'CD' => 112.41,
           'IN' => 114.82,
           'SN' => 118.71,
           'SB' => 121.75,
           'TE' => 127.6,
           'I'  => 126.905,
           'XE' => 131.29,
           'CS' => 132.905,
           'BA' => 137.33,
           'LA' => 138.906,
           'CE' => 140.12,
           'PR' => 140.908,
           'ND' => 144.24,
           'PM' => 145,
           'SM' => 150.36,
           'EU' => 151.96,
           'GD' => 157.25,
           'TB' => 158.925,
           'DY' => 162.5,
           'HO' => 164.93,
           'ER' => 167.26,
           'TM' => 168.934,
           'YB' => 173.04,
           'LU' => 174.967,
           'HF' => 178.49,
           'TA' => 180.948,
           'W'  => 183.85,
           'RE' => 186.207,
           'OS' => 190.2,
           'IR' => 192.22,
           'PT' => 195.08,
           'AU' => 196.967,
           'HG' => 200.59,
           'TL' => 204.383,
           'PB' => 207.2,
           'BI' => 208.98,
           'PO' => 209,
           'AT' => 210,
           'RN' => 222,
           'FR' => 223,
           'RA' => 226.025,
           'AC' => 227.028,
           'TH' => 232.038,
           'PA' => 231.036,
           'U'  => 238.029,
           'NP' => 237.048,
           'PU' => 244,
           'AM' => 243,
           'CM' => 247,
           'BK' => 247,
           'CF' => 251,
           'EI' => 252,
           'FM' => 257,
           'MD' => 258,
           'NO' => 259,
           'LW' => 269,
           'D'  => 2.0141,
           'E'  => 5.45E-4);
    
$energy   = "CAL/MOLE";
$quantity = "MOLES";

sub get_options
{
    use Getopt::Long;

    GetOptions("chemkin=s" => \$CKFILE,
               "therm=s"   => \$THERMFILE,
               "tran=s"    => \$TRANFILE);

    if (!$CKFILE || !$THERMFILE || !$TRANFILE)
    {
        die "Usage: ./ckread.pl --chemkin=file1 --therm=file2 --tran=file3\n";
    }
}

sub is_integer
{
    my $str = shift;

    if ($str =~ /^[+-]?\d+$/o)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

sub is_numeric
{
    use POSIX qw(strtod);

    my $str = shift;

    $str =~ s/^\s+//o;
    $str =~ s/\s+$//o;

    $! = 0;

    my($num,$unparsed) = strtod($str);

    if (($str eq '') || ($unparsed != 0) || $!)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

sub is_species
{
    my $str = shift;

    for $f (@species)
    {
        return 1 if ($f =~ /$str/i);
    }

    return 0;
}

sub get_elements
{
    open(CKFILE) || die "Couldn't open $CKFILE\n";

    my $got_one = 0;

    print "Reading elements from $CKFILE ... ";

  TOP:
    while (<CKFILE>)
    {
        chop; s/^\s+//o; s/\s+$//o; tr/a-z/A-Z/;

        if (/[!]/o)
        {
            my ($l) = split('!',$_); $_ = $l;
        }

        last if (/^SPEC/o);

        if (/^END/o)
        {
            $got_one = 0; next;
        }

        if (/^ELEM/o || $got_one == 1)
        {
            my @tokens = split(/[ \/]+/,$_);

            if (/^ELEM/o)
            {
                if (!($tokens[0] =~ /(ELEM)|(ELEMENTS)/o))
                {
                    die "malformed line: $_\n";
                }

                $got_one = 1;

                shift @tokens;
            }

            while ($#tokens >= 0)
            {
                if ($tokens[0] eq "END")
                {
                    $got_one = 0; next TOP;
                }

                die "Invalid element: $tokens[0]\n" if (!exists($Weight{$tokens[0]}));

                my $duplicate = 0;

                for $f (@elements)
                {
                    $duplicate = 1 if ($f eq $tokens[0]);
                }

                push @elements, $tokens[0] if !$duplicate;

                my $element = shift @tokens;

                if ($#tokens >= 0 && is_numeric($tokens[0]))
                {
                    $Weight{$element} = $tokens[0];

                    shift @tokens;
                }
            }
        }
    }

    close(CKFILE);

    die "No elements found!!!\n" if $#elements < 0;

    print "done\n";

    print "ELEMENTS(", ($#elements+1), "): ";
    for $f (@elements)
    {
        print $f, " ";
    }
    print "\n";
}

sub get_species
{
    open(CKFILE) || die "Couldn't open $CKFILE\n";

    print "Reading species from $CKFILE ... ";

    my $got_one = 0;

  TOP:
    while (<CKFILE>)
    {
        chop; s/^\s+//o; s/\s+$//o; tr/a-z/A-Z/;

        if (/[!]/o)
        {
            my ($l) = split('!',$_); $_ = $l;
        }

        last if (/^REAC/o);

        if (/^END/o)
        {
            $got_one = 0; next;
        }

        if (/^SPEC/o || $got_one == 1)
        {
            my @tokens = split(' ',$_);

            if (/^SPEC/o)
            {
                if (!($tokens[0] =~ /(SPEC)|(SPECIES)/o))
                {
                    die "malformed line: $_\n";
                }

                $got_one = 1;

                shift @tokens;
            }

            while ($#tokens >= 0)
            {
                if ($tokens[0] eq "END")
                {
                    $got_one = 0; next TOP;
                }

                die "Species names must be less that 17 symbols: $tokens[0]\n" if (length($tokens[0]) > 16);

                die "Species names must not begin with a +, =, or a number: $tokens[0]\n" if ($tokens[0] =~ /^[+=0-9]/o);

                die "Ionic species are not supported: $tokens[0]\n" if ($tokens[0] =~ /[-+]$/o);

                my $duplicate = 0;

                for $f (@species)
                {
                    $duplicate = 1 if ($f eq $tokens[0]);
                }

                push @species, $tokens[0] if !$duplicate;

                shift @tokens;
            }
        }
    }

    close(CKFILE);

    die "No species found!!!\n" if $#species < 0;

    print "done\n";

    print "SPECIES(", ($#species+1), "): ";
    for $f (@species)
    {
        print $f, " ";
    }
    print "\n";
}

sub get_reactions
{
    open(CKFILE) || die "Couldn't open $CKFILE\n";

    my $got_one = 0;

    print "Reading reactions ... ";

  TOP:
    while (<CKFILE>)
    {
        chop; s/^\s+//o; s/\s+$//o; tr/a-z/A-Z/;

        if (/[!]/o)
        {
            my ($l) = split('!',$_); $_ = $l;
        }

        if (/^END/o)
        {
            $got_one = 0; next;
        }

        if (/^REAC/o || $got_one == 1)
        {
            while (/[&]$/o)
            {
                s/[&]$/ /o;
                $_ .= <CKFILE>;
            }

            my @tokens = split(' ',$_);

            next if $#tokens < 0;

            if (/^REAC/o)
            {
                if (!($tokens[0] =~ /(REAC)|(REACTIONS)/o))
                {
                    die "malformed line: $_\n";
                }

                $got_one = 1;

                shift @tokens;

                while ($#tokens >= 0)
                {
                    if ($tokens[0] =~ /(CAL\/MOLE)|(KCAL\/MOLE)|(JOULES\/MOLE)|(KJOULES\/MOLE)|(KELVINS)/o)
                    {
                        $energy = $tokens[0];
                    }
                    elsif ($tokens[0] =~ /(MOLES)|(MOLECULES)/o)
                    {
                        $quantity = $tokens[0];
                    }

                    shift @tokens;
                }

                next TOP;
            }

            if (/[=]/o)
            {
                die "malformed reaction line: $_\n" if ($#tokens < 3);

                my $line = "";

                for ($i = 0; $i <= $#tokens - 3; $i++)
                {
                    $line .= $tokens[$i];
                }

                $line .= " "; $line .= $tokens[$#tokens-2];
                $line .= " "; $line .= $tokens[$#tokens-1];
                $line .= " "; $line .= $tokens[$#tokens-0];

                my $duplicate = 0;

                for $f (@reactions)
                {
                    $duplicate = 1 if ($f eq $line);
                }

                die "Duplicate reaction lines not allowed: $line\n" if $duplicate;

                push @reactions, $line;
            }
            else
            {
                die "No preceding reaction line: $_\n" if ($#reactions < 0);

                @tokens = split(/[ \/]+/,$_);

                my $line = $reactions[$#reactions];

                while ($#tokens >= 0)
                {
                    if (is_species($tokens[0]))
                    {
                        die "Missing 3rd body efficiency: $_\n" if ($#tokens < 1);

                        $line .= " "; $line .= $tokens[0]; shift @tokens;
                        $line .= " "; $line .= $tokens[0]; shift @tokens;
                    }
                    elsif ($tokens[0] =~ /(DUP)|(DUPLICATE)/o)
                    {
                        shift @tokens;
                    }
                    elsif ($tokens[0] =~ /(REV)|(REVERSE)/o)
                    {
                        die "Missing reverse rate parameter(s): $_\n" if ($#tokens < 3);

                        $line .= " "; $line .= $tokens[0]; shift @tokens;
                        $line .= " "; $line .= $tokens[0]; shift @tokens;
                        $line .= " "; $line .= $tokens[0]; shift @tokens;
                        $line .= " "; $line .= $tokens[0]; shift @tokens;
                    }
                    elsif ($tokens[0] =~ /(LOW)|(HIGH)|(TROE)/o)
                    {
                        my $keyword = $tokens[0];

                        die "Missing pressure dependency info: $_\n" if ($#tokens < 3);

                        $line .= " "; $line .= $tokens[0]; shift @tokens;
                        $line .= " "; $line .= $tokens[0]; shift @tokens;
                        $line .= " "; $line .= $tokens[0]; shift @tokens;
                        $line .= " "; $line .= $tokens[0]; shift @tokens;

                        if ($#tokens >= 0 && $keyword =~ /TROE/o && is_numeric($tokens[0]))
                        {
                            $line .= " "; $line .= $tokens[0]; shift @tokens;
                        }
                    }
                    else
                    {
                        die "Unsupported auxiliary data: $tokens[0]\n";
                    }
                }

                my $duplicate = 0;

                for ($i = 0; $i < $#reactions; $i++)
                {
                    $duplicate = 1 if ($i eq $line);
                }

                die "Duplicate reaction lines not allowed: $line\n" if $duplicate;

                $reactions[$#reactions] = $line;
            }
        }
    }

    close(CKFILE);

    die "No reactions found!!!" if $#reactions < 0;

    print "done\n";

    print "REACTIONS(", ($#reactions+1), "):\n";
    for $f (@reactions)
    {
        print $f, "\n";
    }
}

sub get_therm
{
    #
    # THERM data can be in either the CHEMKIN file or its own them.dat file.
    #
    open THERM, "$_[0]" || die "Couldn't open $_[0]\n";

    my $got_one           = 0;
    my $expect_temp_range = 0;

    print "Reading thermodynamic data out of $_[0] ... ";

  TOP:
    while (<THERM>)
    {
        chop; tr/a-z/A-Z/;

        if (/[!]/o)
        {
            my ($l) = split('!',$_); $_ = $l;
        }
        #
        # If data is in CHEMKIN file it must precede reaction data.
        # Note that thermo data in CHEMKIN files overrides that in therm file.
        #
        last if (/^REAC/o);

        if (/^END/o)
        {
            $got_one = 0; next;
        }

        if (/^THER/o || $got_one == 1)
        {
            my @tokens = split(' ',$_);

            if (/^THER/o)
            {
                die "malformed line: $_\n" if (!($tokens[0] =~ /(THER)|(THERMO)/o));

                $got_one = 1; shift @tokens;

                $expect_temp_range = 1 if (($#tokens >= 0 && $tokens[0] =~ /ALL/o) || $_[0] eq $THERMFILE);

                next;
            }

            if ($expect_temp_range)
            {
                if (/(.{10})(.{10})(.{1,10})/o)
                {
                    $LowestT  = $1;
                    $CommonT  = $2;
                    $HighestT = $3;

                    $LowestT  =~ s/[ ]+//og;
                    $CommonT  =~ s/[ ]+//og;
                    $HighestT =~ s/[ ]+//og;
                }
                else
                {
                    die "Expected temperature range coefficients (3F10.0): $_\n";
                }

                $expect_temp_range = 0; next;
            }

            if ($#tokens >= 0)
            {
                my $line = "";

                if ($tokens[0] eq "END")
                {
                    $got_one = 0; next TOP;
                }

                if (/1$/o)
                {
                    chop;

                    if (/(.{18}).{6}(.)(.)(...)(.)(.)(...)(.)(.)(...)(.)(.)(...)(.)(.{10})(.{10})(.{8})(.)(.)(...)/o)
                    {
                        my @results = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20);

                        for ($i = 0; $i <= $#results; $i++)
                        {
                            #
                            # Use '@' to represent a blank field.
                            #
                            $results[$i] = "@" if ($results[$i] =~ /^[ ]+$/o);
                            #
                            # Remove blanks. Fortran allow blanks in numbers.
                            #
                            $results[$i] =~ s/[ ]+//og;
                        }

                        die "Phase must be S, L or G: $_\n" if ($results[13] !~ /^[SLG]$/o);

                        die "Lo temperature must be numeric: $_\n" if (!is_numeric($results[14]));
                        die "Hi temperature must be numeric: $_\n" if (!is_numeric($results[15]));

                        die "Common temperature must be numeric: $_\n" if ($results[16] ne "@" && !is_numeric($results[16]));

                        for ($i = 0; $i <= $#results; $i++)
                        {
                            $line .= $results[$i]; $line .= " ";
                        }
                        #
                        # Remove old entry for this species in THERMO data.
                        #
                        for ($i = 0; $i <= $#thermo; $i++)
                        {
                            my ($spec) = split(' ',$thermo[$i]);

                            if ($results[0] eq $spec)
                            {
                                my @newthermo = ();

                                for ($j = 0; $j <= $#thermo; $j++)
                                {
                                    push @newthermo, $thermo[$j] if ($j != $i);
                                }

                                @thermo = @newthermo;

                                last;
                            }
                        }
                        #
                        # Always push entries in progress on the end of thermo.
                        #
                        push @thermo, $line;

                        next TOP;
                    }
                    else
                    {
                        die "malformed thermo line (type 1): $_\n";
                    }
                }
                elsif (/[23]$/o)
                {
                    chop;

                    if (/(.{15})(.{15})(.{15})(.{15})(.{15})/o)
                    {
                        my @results = ($1,$2,$3,$4,$5);

                        for ($i = 0; $i <= $#results; $i++)
                        {
                            #
                            # Fortran allows spaces in numbers on input :-(
                            #
                            $results[$i] =~ s/[ ]+//og;

                            die "Coefficients must be numeric: $_\n" if (!is_numeric($results[$i]));
                        }

                        for ($i = 0; $i <= $#results; $i++)
                        {
                             $line .= $results[$i]; $line .= " ";
                        }

                        $thermo[$#thermo] .= $line;

                        next TOP;
                    }
                    else
                    {
                        die "malformed thermo line (type 2 or 3): $_\n";
                    }
                }
                elsif (/4$/o)
                {
                    chop;

                    if (/(.{15})(.{15})(.{15})(.{15})/o)
                    {
                        my @results = ($1,$2,$3,$4);

                        for ($i = 0; $i <= $#results; $i++)
                        {
                            #
                            # Fortran allows spaces in numbers on input :-(
                            #
                            $results[$i] =~ s/[ ]+//og;

                            die "Coefficients must be numeric: $_\n" if (!is_numeric($results[$i]));
                        }

                        for ($i = 0; $i <= $#results-1; $i++)
                        {
                            $line .= $results[$i]; $line .= " ";
                        }
                        $line .= $results[$#results];

                        $thermo[$#thermo] .= $line;

                        next TOP;
                    }
                    else
                    {
                        die "malformed thermo line (type 4): $_\n";
                    }
                }
                else
                {
                    die "malformed line: $_\n";
                }
            }
        }
    }

    close(THERM);

    die "No THERMO data found!!!\n" if $#thermo < 0;

    print "done\n";

    if ($_[0] eq $CKFILE)
    {
        print "THERMO data:\n";
        for $f (@thermo)
        {
            print $f, "\n";
        }
    }
}

get_options;
get_elements;
get_species;
get_reactions;
get_therm($THERMFILE);
get_therm($CKFILE);

if (!$LowestT || !$CommonT || !$HighestT)
{
    die "Lowest T, Common T & Highest T MUST be set in THERMO data\n";
}
#
# Make sure each species has a thermodynamic entry.
#
for $f (@species)
{
    my $found = 0;

    for $t (@thermo)
    {
        my ($spec) = split(' ',$t);

        if ($spec eq $f)
        {
            $found = 1; last;
        }
    }

    die "Species $f does not have THERMO data.\n" if !$found;
}
