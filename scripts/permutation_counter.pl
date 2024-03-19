#/usr/bin/perl

#     permutation_counter
#     Copyright (C) 2024  Egon A. Ozer

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to the Free Software Foundation, Inc.,
#     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


use warnings;
use strict;
use bignum;
no warnings 'recursion';

my $usage = "
permutation_counter.pl <n> <r> <order important?> <repetition allowed?> <core_and_pangenome?>

Inputs:
  n:                    Types to choose from. Positive integer.
  r:                    Number chosen. Positive integer.
  order important?:     Y or N
  repetition allowed?:  Y or N
  core_and_pangenome?:  Y or N.  If Y, ignores the order and repetition questions

Outputs:
  Total number of possible permutations. Prints to STDOUT.

";

die $usage unless @ARGV >= 4;

my ($n, $r, $order, $rep, $cp) = @ARGV;
die $usage unless $n =~ m/^\d+$/; #is it an integer?
die $usage unless $r =~ m/^\d+$/; #is it an integer?
die $usage unless $n > 0;
die $usage unless $r > 0;
$order = uc($order);
$rep = uc($rep);
die $usage unless $order =~ m/[YN]/;
die $usage unless $rep =~ m/[YN]/;
if ($cp){
    die $usage unless $cp =~ m/[YN]/;
}

my $perms = 0;
if ($cp and $cp eq "N"){
    if ($order eq "N"){
        if ($rep eq "N"){
            $perms = fac($n) / (fac($n - $r) * fac($r));
        } else {
            $perms = (fac($n + $r -1) / (fac($r) * fac($n - 1)));
        }
    } else {
        if ($rep eq "N"){
            $perms = fac($n) / fac($n - $r);
        } else {
            $perms = $n ** $r;
        }
    }
} elsif ($cp and $cp eq "Y") {
    $perms = fac($n) / (fac($r - 1) * fac($n - $r));
}
print "$perms\n";

#-----------------
sub fac {
  $_[0]>1?$_[0]*fac($_[0]-1):1;
}
