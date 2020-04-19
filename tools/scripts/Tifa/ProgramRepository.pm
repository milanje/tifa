package Tifa::ProgramRepository;

# Copyright (C) 2011 CNRS - Ecole Polytechnique - INRIA.
#
# This file is part of TIFA.
#
# TIFA is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# TIFA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
# more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.

#------------------------------------------------------------------------------
#                  A repository of TIFA programs.
#------------------------------------------------------------------------------
# File    : ProgramRepository.pm
# Author  : Jerome Milan
# Version : 0.1
#
# This is a collection of available TIFA::Program's with a very simple
# interface to access them.
#------------------------------------------------------------------------------

use 5.006002;
use strict;
use warnings;
use Class::Struct;
use File::Basename;
use Carp;

use Tifa::Program;

require Exporter;

our @ISA        = qw(Exporter);
our @EXPORT_OK  = qw();
our $VERSION    = '0.1';

my %algo_to_program = ();


my @param_names   = ();
my @param_descrs  = ();
my @param_types   = ();
my $cmdline       = '';
my $help          = '';

#------------------------------------------------------------------------------
#                           CFRAC algorithm
#------------------------------------------------------------------------------
my $cfrac_program = new Tifa::Program();

$algo_to_program{"cfrac"} = $cfrac_program;

$cfrac_program->set_algo("cfrac");
$cfrac_program->set_descr("The Continued FRACtion algorithm");
$cfrac_program->set_exe("cfrac");

@param_names = (
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_res",
    "filter_method",
    "nsteps_early_abort",
    "nrelations",
    "linalg_method",
    "use_large_primes",
    "nprimes_tdiv"
);
@param_descrs = (
    "Number of primes in factor base",
    "Number of primes to trial divide by the residues",
    "Method of smooth residue detection",
    "Number of steps in early abort strategy",
    "Number of relations",
    "Method of linear system resolution",
    "Whether to use the large prime variation or not",
    "Number of primes to trial divide the number to factor"
);
@param_types = (
    "int",
    "int",
    "int",
    "int",
    "int",
    "int",
    "switch",
    "int"
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_res",
    "filter_method",
    "nsteps_early_abort",
    "nrelations",
    "linalg_method",
    "use_large_primes",
    "nprimes_tdiv"
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$cfrac_program->add_mode("no_defaults",
                         "Default mode: specify all parameter values",
                         \@param_names, \@param_descrs, \@param_types,
                         $cmdline);

$cfrac_program->set_default_mode("no_defaults");
$cfrac_program->set_mode("no_defaults");

@param_names  = (
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes to trial divide the number to factor",
);
@param_types  = (
    "int",
);
$cmdline =  "exe nprimes_tdiv";
#
# use_defaults mode: Let CFRAC use the precomputed optimal parameter values
#                    depending on the size of the number to factor
#
$cfrac_program->add_mode("use_defaults",
                         "'Best' mode: let CFRAC use optimal values",
                         \@param_names, \@param_descrs, \@param_types,
                         $cmdline);
$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to cfrac _and_ the exe
      parameter is not defined or left empty.
      Default: cfrac_program

  --nprimes_in_factor_base=i
      (mandatory in 'no_defaults' mode)
      Number of primes of the factor base.

  --nprimes_tdiv_smooth_res=i
      (mandatory in 'no_defaults' mode)
      Number of primes to consider in the trial division of the
      smooth residues.

  --filter_method=i
      (mandatory in 'no_defaults' mode)
      Method to use for smooth residue detection.
      0: Trial division
      1: Trial division + early abort
      2: Bernstein's batch

  --nsteps_early_abort=i
      Number of steps in early abort strategy for residues factorization.
      No early abort is performed if nsteps_early_abort = 0.

  --nrelations=i
      (mandatory in 'no_defaults' mode)
      Number of relations to find during the linear resolution
      phase of the CFRAC algorithm.

  --linalg_method=i
      (mandatory in 'no_defaults' mode)
      Linear system resolution method to use.
      Accepted values:
          0 : Smart gaussian elimination

  --use_large_primes
      (mandatory in 'no_defaults' mode)
      Use the large prime variation.

  --mode=s
      (optional, default value used if none provided)
      Sets the program mode to use.

      use_defaults: In this mode, the CFRAC program will determine the values
                    of the aforementioned parameters. These values are choosen
                    according to the size of the number to factor and are more
                    or less optimal if this size is between 60 and 200 bits.

      no_defaults:  In this mode, all parameters should be specified.

      Default value: no_defaults

  --nprimes_tdiv=i
      (mandatory in all modes)
      Number of primes to use to trial divide the integer to factor.

EOS

$cfrac_program->set_help($help);

#------------------------------------------------------------------------------
#                           ECM algorithm
#------------------------------------------------------------------------------
my $ecm_program = new Tifa::Program();

$algo_to_program{"ecm"} = $ecm_program;

$ecm_program->set_algo("ecm");
$ecm_program->set_descr("The Elliptic Curve Method");
$ecm_program->set_exe("ecm");

@param_names = (
    "nprimes_tdiv",
    "b1",
    "b2",
    "ncurves"
);
@param_descrs = (
    "Number of primes to trial divide the number to factor",
    "Bound for phase 1",
    "Bound for phase 2",
    "Number of curves to try before giving up"
);
@param_types = (
    "int",
    "int",
    "int",
    "int"
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_tdiv",
    "b1",
    "b2",
    "ncurves"
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$ecm_program->add_mode("no_defaults",
                       "Default mode: specify all parameter values",
                        \@param_names, \@param_descrs, \@param_types,
                        $cmdline);

$ecm_program->set_default_mode("no_defaults");
$ecm_program->set_mode("no_defaults");

@param_names  = (
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes to trial divide the number to factor",
);
@param_types  = (
    "int",
);
$cmdline =  "exe nprimes_tdiv";
#
# use_defaults mode: Let ECM use the precomputed optimal parameter values
#                    depending on the size of the number to factor
#
$ecm_program->add_mode("use_defaults",
                       "'Best' mode: let ECM use optimal values",
                        \@param_names, \@param_descrs, \@param_types,
                        $cmdline);
$help = <<EOS;
Parameter(s):

  --b1=i
      (mandatory in 'no_defaults' mode)
      Bound used in phase 1 of the ECM.

  --b2=i
      (mandatory in 'no_defaults' mode)
      Bound used in phase 2 of the ECM.

  --ncurves=i
      (mandatory in 'no_defaults' mode)
      Number of curves to try before giving up.

  --mode=s
      (optional, default value used if none provided)
      Sets the program mode to use.

      use_defaults: In this mode, the ECM program will determine the values
                    of the aforementioned parameters. These values are choosen
                    according to the size of the number to factor and are more
                    or less optimal if this size is between 60 and 200 bits.

      no_defaults:  In this mode, all parameters should be specified.

      Default value: no_defaults

  --nprimes_tdiv=i
      (mandatory in all modes)
      Number of primes to use to trial divide the integer to factor.

EOS

$ecm_program->set_help($help);

#------------------------------------------------------------------------------
#                       Fermat's algorithm (McKee's speedup)
#------------------------------------------------------------------------------
my $fermat_program = new Tifa::Program();

$algo_to_program{"fermat"} = $fermat_program;

$fermat_program->set_algo("fermat");
$fermat_program->set_descr("Fermat's algorithm (McKee's speedup)");
$fermat_program->set_exe("fermat");

@param_names = (
    "nprimes_tdiv"
);
@param_descrs = (
    "Number of primes to trial divide the number to factor"
);
@param_types = (
    "int"
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_tdiv"
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$fermat_program->add_mode("no_defaults",
                          "Default mode: specify all parameter values",
                          \@param_names, \@param_descrs, \@param_types,
                          $cmdline);

$fermat_program->set_default_mode("no_defaults");
$fermat_program->set_mode("no_defaults");

$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to fermat _and_ the exe
      parameter is not defined or left empty.
      Default: fermat_program

  --nprimes_tdiv=i
      (mandatory)
      Number of primes to use to trial divide the integer to factor.

EOS

$fermat_program->set_help($help);
#------------------------------------------------------------------------------
#                           SIQS algorithm
#------------------------------------------------------------------------------
my $siqs_program = new Tifa::Program();

$algo_to_program{"siqs"} = $siqs_program;

$siqs_program->set_algo("siqs");
$siqs_program->set_descr("The Self-Initializing Quadratic Sieve algorithm");
$siqs_program->set_exe("siqs");

@param_names = (
    "sieve_half_width",
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_res",
    "nrelations",
    "linalg_method",
    "use_large_primes",
    "nprimes_tdiv",
);
@param_descrs = (
    "Half width of the sieving interval",
    "Number of primes in factor base",
    "Number of primes to trial divide by the residues",
    "Number of relations",
    "Method of linear system resolution",
    "Whether to use the large prime variation or not",
    "Number of primes to trial divide the number to factor",
);
@param_types = (
    "int",
    "int",
    "int",
    "int",
    "int",
    "switch",
    "int",
);
$cmdline = join(
    ' ',
    "exe",
    "sieve_half_width",
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_res",
    "nrelations",
    "linalg_method",
    "use_large_primes",
    "nprimes_tdiv",
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$siqs_program->add_mode("no_defaults",
                        "Default mode: specify all parameter values",
                        \@param_names, \@param_descrs, \@param_types,
                        $cmdline);

$siqs_program->set_default_mode("no_defaults");
$siqs_program->set_mode("no_defaults");

@param_names  = (
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes to trial divide the number to factor",
);
@param_types  = (
    "int",
);
$cmdline =  "exe nprimes_tdiv";
#
# use_defaults mode: Let SIQS use the precomputed optimal parameter values
#                    depending on the size of the number to factor
#
$siqs_program->add_mode("use_defaults",
                        "'Best' mode: let QS use optimal values",
                        \@param_names, \@param_descrs, \@param_types,
                        $cmdline);

$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to siqs _and_ the exe
      parameter is not defined or left empty.
      Default: siqs_program

  --sieve_half_width=i
      (mandatory)
      Half width of the sieving interval. The complete sieving interval is
      thus given by [-sieve_half_width, +sieve_half_width].

  --nprimes_in_factor_base=i
      (mandatory)
      Number of primes of the factor base.

  --nprimes_tdiv_smooth_res=i
      (mandatory)
      Number of primes to consider in the trial division of the
      smooth residues.

  --nrelations=i
      (mandatory)
      Number of relations to find during the linear resolution
      phase of the SIQS algorithms.

  --linalg_method=i
      (mandatory)
      Linear system resolution method to use.
      Accepted values:
          0 : Smart gaussian elimination

  --use_large_primes
      (optional)
      Use the large prime variation

  --nprimes_tdiv=i
      (mandatory)
      Number of primes to use to trial divide the integer to factor.
EOS
$siqs_program->set_help($help);
#------------------------------------------------------------------------------
#                          SQUFOF algorithm
#------------------------------------------------------------------------------
my $squfof_program = new Tifa::Program();

$algo_to_program{"squfof"} = $squfof_program;

$squfof_program->set_algo("squfof");
$squfof_program->set_descr("The SQUare FOrm Factorization algorithm");
$squfof_program->set_exe("squfof");

@param_names = (
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes to trial divide the number to factor",
);
@param_types = (
    "int",
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_tdiv",
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$squfof_program->add_mode("no_defaults",
                          "Default mode: specify all parameter values",
                          \@param_names, \@param_descrs, \@param_types,
                          $cmdline);

$squfof_program->set_default_mode("no_defaults");
$squfof_program->set_mode("no_defaults");

$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to squfof _and_ the exe
      parameter is not defined or left empty.
      Default: squfof_program

  --nprimes_tdiv=i
      (mandatory)
      Number of primes to use to trial divide the integer to factor.

EOS
$squfof_program->set_help($help);
#------------------------------------------------------------------------------
#                          Trial division algorithm
#------------------------------------------------------------------------------
my $tdiv_program = new Tifa::Program();

$algo_to_program{"tdiv"} = $tdiv_program;

$tdiv_program->set_algo("tdiv");
$tdiv_program->set_descr("The naive trial division algorithm");
$tdiv_program->set_exe("tdiv");

@param_names = (
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes used to trial divide the number to factor",
);
@param_types = (
    "int",
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_tdiv",
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$tdiv_program->add_mode("no_defaults",
                        "Default mode: specify all parameter values",
                        \@param_names, \@param_descrs, \@param_types,
                        $cmdline);

$tdiv_program->set_default_mode("no_defaults");
$tdiv_program->set_mode("no_defaults");

$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to tdiv _and_ the exe
      parameter is not defined or left empty.
      Default: tdiv_program

  --nprimes_tdiv=i
      (mandatory)
      Number of primes to use to trial divide the integer to factor.
EOS
$tdiv_program->set_help($help);
#------------------------------------------------------------------------------
#                          Pollard's rho algorithm
#------------------------------------------------------------------------------
my $rho_program = new Tifa::Program();

$algo_to_program{"rho"} = $rho_program;

$rho_program->set_algo("rho");
$rho_program->set_descr("Pollard's rho algorithm");
$rho_program->set_exe("rho");

@param_names = (
    "nprimes_tdiv",
    "a",
    "b",
    "max_niterations",
    "ndeltas_to_accumulate",
    "cycle_finding_method"
);
@param_descrs = (
    "Number of primes to trial divide the number to factor",
    "Iteration function given by (a * x^2 + b) mod n",
    "Iteration function given by (a * x^2 + b) mod n",
    "Maximum number of iteration to perform before giving up",
    "Number of (x_{2i} - x_{i}) to multiply before computing a gcd",
    "Cycle finding method to use"
);
@param_types = (
    "int",
    "int",
    "int",
    "int",
    "int",
    "int"
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_tdiv",
    "a",
    "b",
    "max_niterations",
    "ndeltas_to_accumulate",
    "cycle_finding_method"
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$rho_program->add_mode("no_defaults",
                       "Default mode: specify all parameter values",
                       \@param_names, \@param_descrs, \@param_types, $cmdline);

$rho_program->set_default_mode("no_defaults");
$rho_program->set_mode("no_defaults");

@param_names  = (
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes to trial divide the number to factor",
);
@param_types  = (
    "int",
);
$cmdline =  "exe nprimes_tdiv";
#
# use_defaults mode
#
$rho_program->add_mode("use_defaults",
                       "use default values",
                       \@param_names, \@param_descrs, \@param_types, $cmdline);
$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to rho _and_ the exe
      parameter is not defined or left empty.
      Default: rho_program

  --nprimes_tdiv=i
      (mandatory in all modes)
      Number of primes to use to trial divide the integer to factor.

  --a=i
      (mandatory in 'no_defaults' mode)
      The 'a' coefficient in the iteration function given by
      (a * x^2 + b) mod n, where n is the number to factor.

  --b=i
      (mandatory in 'no_defaults' mode)
      The 'b' coefficient in the iteration function given by
      (a * x^2 + b) mod n, where n is the number to factor.

  --max_niterations=i
      (mandatory in 'no_defaults' mode)
      Maximum number of iterations to perform before giving up the
      factorization.

  --ndeltas_to_accumulate=i
      (mandatory in 'no_defaults' mode)
      Number of (x_{j} - x_{i}) to multiply together before computing
      gcd(prod (x_{j} - x_{i})), n).

  --cycle_finding_method=i
      (mandatory in 'no_defaults' mode)
      Cycle finding method to use.
      Accepted values:
          0 : Floyd's cycle finding method.
          1 : Brent's cycle finding method.

EOS

$rho_program->set_help($help);
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
sub get_algo_to_program_hash {
    return %algo_to_program;
}
#------------------------------------------------------------------------------
sub get_algo_list {
    return keys %algo_to_program;
}
#------------------------------------------------------------------------------
sub get_program_list {
    return values %algo_to_program;
}
#------------------------------------------------------------------------------
sub get_all_param_names {
    my @all_params = ();
    foreach my $algo (keys %algo_to_program) {
        push(@all_params, $algo_to_program{$algo}->get_all_param_names());
    }
    my %hash = map {$_, 1} @all_params; # get rid of duplicate values

    return keys %hash;
}
#------------------------------------------------------------------------------
sub get_all_param_to_descrs_hash {
    my %hash = ();
    foreach my $algo (keys %algo_to_program) {
        my %subhash = $algo_to_program{$algo}->get_all_param_to_descrs_hash();
        %hash = (%hash, %subhash);
    }
    return %hash;
}
#------------------------------------------------------------------------------
sub get_all_param_to_types_hash {
    my %hash = ();
    foreach my $algo (keys %algo_to_program) {
        my %subhash = $algo_to_program{$algo}->get_all_param_to_types_hash();
        %hash = (%hash, %subhash);
    }
    return %hash;
}
#------------------------------------------------------------------------------
sub get_getopt_strings {
    my $self = shift;

    my %hash = get_all_param_to_types_hash();

    my @getopt_strings = ();
    my $tmp = "";

    foreach my $param (keys %hash) {
        $tmp = $param;
        if ($hash{$param} ne "switch") {
            if ($hash{$param} eq "string") {
                $tmp .= '=s';
                push(@getopt_strings, $tmp);
            }
            if ($hash{$param} eq "int") {
                $tmp .= '=i';
                push(@getopt_strings, $tmp);
            }
            if ($hash{$param} eq "extint") {
                $tmp .= '=o';
                push(@getopt_strings, $tmp);
            }
            if ($hash{$param} eq "float") {
                $tmp .= '=f';
                push(@getopt_strings, $tmp);
            }
        } else {
            push(@getopt_strings, $tmp);
        }
    }
    push(@getopt_strings, "algo=s");

    return @getopt_strings;
}
#------------------------------------------------------------------------------
sub get_program_from_name {
    my $name     = basename(shift);

    foreach my $prog (values %algo_to_program) {
        if ( $name eq basename($prog->get_exe()) ) {
            return $prog;
        }
    }
    carp("ERROR: ProgramRepository: $name is not a known program");
    return undef;
}
#------------------------------------------------------------------------------

1;

__END__

#------------------------------------------------------------------------------
# Stub documentation for this module.
#------------------------------------------------------------------------------

=head1 NAME

Tifa::ProgramRepository - A repository of available Tifa::Program's

=head1 SYNOPSIS

 use Tifa::ProgramRepository;

 %algo_to_program = Tifa::ProgramRepository::get_algo_to_program_hash();

=head1 REQUIRE

Perl 5.006002, Carp, Class::Struct, Exporter and Tifa::Program.

=head1 SUMMARY

The Tifa::ProgramRepository module acts as a repository of all available
Tifa::Program objects.

=head1 DESCRIPTION

The Tifa::ProgramRepository module acts as a repository of all available
(factoring) Tifa::Program's. The listed Tifa::Program objects can then be used
in other scripts or modules.

=head2 Available methods

    get_algo_to_program_hash()
    get_algo_list()
    get_program_list()
    get_program_from_name($progname)
    get_all_param_names()
    get_all_param_to_descrs_hash()
    get_all_param_to_types_hash()

=head2 Methods description

    get_algo_to_program_hash()
        Returns a hashtable mapping algorithm names to the Tifa::Program
        objects implementing it.

    get_algo_list()
        Returns an array listing the names of the implemented algorithm.

    get_program_list()
        Returns an array listing the Tifa::Program objects available.

    get_program_from_name($progname)
        Returns a reference to the Tifa::Program object whose default
        program name is given by $progname (or the base name of $progname
        if $progname is a path).

    get_all_param_names()
        Returns an array listing all the parameter names from all
        programs.

    get_all_param_to_descrs_hash()
        Returns a hashtable mapping all the parameter names (from all
        programs and regardless of their current modes), to their
        descriptions.

    get_all_param_to_types_hash()
        Returns a hashtable mapping all the parameter types (from all
        programs and regardless of their current modes), to their
        descriptions.

=head1 EXAMPLE

This following code snippet gives an example of how the Tifa::ProgramRepository
module can be used.

 %algo_to_program = Tifa::ProgramRepository::get_algo_to_program_hash();

 $program = $algo_to_program{"squfof"};

 %param_vals = {
     "exe"   => "./squfof",
     "ntdiv" => 128,
 };

 $program->execute(\%param_vals, 816379, "> trace.txt");

=head1 EXPORT

No functions are exported from this package by default.

=head1 SEE ALSO

The Tifa::Program module's man page.

=head1 AUTHOR

Jerome Milan, E<lt>milanje at gmail dot comE<gt>

=head1 VERSION

0.1.0

=head1 COPYRIGHT AND LICENSE

INCLUDE_LICENSE_AS_TEXT

=cut
