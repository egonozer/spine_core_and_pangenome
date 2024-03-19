# Core and Pangenome Permutation Calculator

**spine_core_and_pangenome.pl** - Uses position counts output by Spine (when not using the --mini option) to calculate sizes of the pangenome, core genome, and novel genomic sequence for every combination of genomes from 1 to n where n is the total number of genomes given to Spine.  Also calculates the core genome size of every possible combination of the n genomes at every possible definition of core genome, i.e. from requiring that a base be present in all n genomes to be considered core (absent from 0), to that a base can be present in only 1 of the n genomes to be considered core (absent fron n - 1).

## Install
Make sure the **scripts** directory stays in the same directory as the **spine_core_and_pangenome.pl** script.

## Usage

`perl spine_core_and_pangenome.pl` [options] \<position_counts.txt>

## Options
*  `-s`    processing steps to run:
            1: Only run core and pangenome calculations,
            2: Only run core at each possible a-value calculations,
            3: Run both step 1 and step 2 [default: 3]
*  `-p`    minumum percentage of genomes in which a region can be present and
        still be considered \"core\". For example, if you wish to consider a
        region that is present in at least 9 of 10 genomes to be \"core\",
        enter \"90\" as the value for this option. A true core (i.e. regions
        present in 100% of genomes) will also be output regardless of the
        value entered here.
        (default: 90)
*  `-a`    maximum number of genomes from which a region can be absent and still
        be counted as core. For example, if the value \"2\" is given, any
        sequence present in at least 8 of 10 or at least 15 of 17 genomes will
        be considered \"core\". A true core (i.e. regions
        present in 100% of genomes) will also be output regardless of the
        value entered here.
        Any value entered here will take precedence over -p, i.e. if -p is
        given as 90 and -a is given as 3, the secondary core genome output will
        be based on a maximum of 3 missing genomes at each position.
*  `-m`    minimum number of permutations to perform at each number of genomes
        during step 1 (pan and core genome at increasing # of input genomes).
        This setting will only to apply to genome numbers with at least the set
        number of possible permutations.
        (default: 0)
        A setting of \"0\" will indicate no minimum permutations
*  `-n`    minimum number of permutations to perform at each number of genomes
        during step 2 (core genome at increasingly flexible definitions of core)
        Can usually go higher with this number than for step 1 as it is a much
        faster calculation.
        (default: same as what is given for -m)
        A setting of \"0\" will indicate no minimum permutations. The number of
        possible permutations for your dataset can be calculated by \"n!\"
        where n = total number of input genomes
*  `-d`    maximum delta around the average, in percent
        Moves on to next genome number when the average of the calculations
        changes by less than this setting as new calculations are added.
        Only applies if a setting > 0 is given to -m
        (default: 0.01)
*  `-o`    output prefix [default: \"out\"]
*  `-t`    number of parallel processes to run [default: 4]
*  `-1`    starting genome number
        (default: 1)
*  `-2`    ending genome number
        (default: total number of genomes in the set)

## License
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License or LICENSE file for more details.
