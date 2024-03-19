#!/usr/bin/perl

my $version = "0.5";
## Changes from v0.4
## Moved calculation of numbers of permutations to an external script. This is because the "use bignum" dependency for large factorial calculations slows down ALL calculations in the script.

## Changes from v0.3:
## Increased efficiency by instead of looping through 2 arrays for core and new, instead loops through array once and collects counts for both core and new
## Clean up STDERR output a bit
## Calculation and output of pan, core, and accessory now done on the fly rather than after complete calculation of all raw values. Should save time and memory.
## Increased efficiency of pan and core calculations by using regex for all values rather than one at a time. Overall this process now is running about twice as fast as before
## Massively simplified the second part: calculation of core genome size at every definition of core (a-value).  Now runs ~ 155x faster (seriously).
## Fixed problem where less than the minimum number of permutations were being output for both parts 1 and 2
## No longer waits for unneeded processes to finish after reaching the minimum permutation count and delta of the averages is less than required. Now just kills any remaining processes and moves on.
## Added option to set separate minimum number of permutations for step 2 (option -n)

## Changes from v0.2:
## For core genome calculations, permutations are calculated immediately upon generation rather than waiting to populate an array with all of the permutations first.  Saves time and memory.
## Added ability to run each calculation step individually (option -s)
## Removed dependency on Algorithm::Permute

## Changes from v0.1:
## Added processing limits for large numbers of genomes (i.e. minimum calculations, maximum delta around averages)
## Added option to produce core either using minimum percentage of genomes (-p) or by using a maximum number of genomes from which a region can be absent (-a)

use strict;
use warnings;
no warnings 'recursion';

my $perm_count = "scripts/permutation_counter.pl"; #since "bignum" that is needed to count permutations is sooooo sloooooow, will do these calculations externally.

my $usage = "
spine_core_and_pangenome.pl - Uses position counts output by Spine (when not
using the --mini option) to calculate sizes of the pangenome, core genome, 
and novel genomic sequence for every combination of genomes from 1 to n where
n is the total number of genomes given to Spine.  Also calculates the core
genome size of every possible combination of the n genomes at every possible
definition of core genome, i.e. from requiring that a base be present in all n
genomes to be considered core (absent from 0), to that a base can be present in
only 1 of the n genomes to be considered core (absent fron n - 1).

version: $version

usage:
perl $0 [options] <position_counts.txt>

options:
  -s    processing steps to run:
            1 - Only run core and pangenome calculations
            2 - Only run core at each possible a-value calculations
            3 - Run both step 1 and step 2 [default]
  -p    minumum percentage of genomes in which a region can be present and
        still be considered \"core\". For example, if you wish to consider a
        region that is present in at least 9 of 10 genomes to be \"core\",
        enter \"90\" as the value for this option. A true core (i.e. regions
        present in 100% of genomes) will also be output regardless of the
        value entered here.
        (default: 90)
  -a    maximum number of genomes from which a region can be absent and still
        be counted as core. For example, if the value \"2\" is given, any
        sequence present in at least 8 of 10 or at least 15 of 17 genomes will
        be considered \"core\". A true core (i.e. regions
        present in 100% of genomes) will also be output regardless of the
        value entered here.
        Any value entered here will take precedence over -p, i.e. if -p is
        given as 90 and -a is given as 3, the secondary core genome output will
        be based on a maximum of 3 missing genomes at each position.
  -m    minimum number of permutations to perform at each number of genomes
        during step 1 (pan and core genome at increasing # of input genomes).
        This setting will only to apply to genome numbers with at least the set
        number of possible permutations.
        (default: 0)
        A setting of \"0\" will indicate no minimum permutations
  -n    minimum number of permutations to perform at each number of genomes
        during step 2 (core genome at increasingly flexible definitions of core)
        Can usually go higher with this number than for step 1 as it is a much
        faster calculation.
        (default: same as what is given for -m)
        A setting of \"0\" will indicate no minimum permutations. The number of
        possible permutations for your dataset can be calculated by \"n!\"
        where n = total number of input genomes
  -d    maximum delta around the average, in percent
        Moves on to next genome number when the average of the calculations
        changes by less than this setting as new calculations are added.
        Only applies if a setting > 0 is given to -m
        (default: 0.01)
  -o    output prefix [default: \"out\"]
  -t    number of parallel processes to run [default: 4]

  -1    starting genome number
        (default: 1)
  -2    ending genome number
        (default: total number of genomes in the set)

";

use Getopt::Std;
our ($opt_s, $opt_p, $opt_a, $opt_m, $opt_d, $opt_o, $opt_t, $opt_n, $opt_1, $opt_2);
getopts('s:p:a:m:d:o:t:n:1:2:');
die $usage unless (@ARGV);

my $step    = $opt_s ? $opt_s : 3;
my $pctcore = $opt_p ? $opt_p : 90;
my $maxmiss = $opt_a ? $opt_a : 0;
my $minperm = $opt_m ? $opt_m : 0;
my $maxdelt = $opt_d ? $opt_d : 0.01;
my $pref    = $opt_o ? $opt_o : "out";
my $threads = $opt_t ? $opt_t : 4;
my $c_minperm = $opt_n ? $opt_n : $minperm;
my $startnum    = $opt_1 if $opt_1;
my $stopnum     = $opt_2 if $opt_2;

my $core_frac = (100 - $pctcore) / 100;

die "Option -s must be set to 1, 2, or 3.\n" unless ($step == 1 or $step == 2 or $step == 3);

print "Loading data ...\n";
my %perm_hash;
my %gen_size; #by collecting the genome sizes at this step, saves cycling through the hash each time the pangeonme calculation needs to start
open (my $in, "< $ARGV[0]") or die "Can't open input file $ARGV[0]: $!\n";
my $g_count = 0;
while (my $line = <$in>){
    chomp $line;
    next if $line =~ m/^ref_genome/;
    my ($ref, $keep_count, $missing, $val) = split('\t', $line);
    $perm_hash{$ref}{$keep_count}{$missing} = $val;
    $g_count = $keep_count if $keep_count > $g_count;
    $gen_size{$ref} += $val;
}
close ($in);

print STDERR "g_count: $g_count\n";

if ($startnum){
    die "ERROR: -1 must be set to an integer between 1 and $g_count\n" unless $startnum >= 1 and $startnum <= $g_count and $startnum =~ m/^\d+$/;
} else {
    $startnum = 1;
}
if ($stopnum){
    die "ERROR: -2 must be set to an integer between 1 and $g_count\n" unless $stopnum >= 1 and $stopnum <= $g_count and $stopnum =~ m/^\d+$/;
} else {
    $stopnum = $g_count;
}
die "ERROR: Value for -2 cannot be less than value for -1\n" if $stopnum < $startnum;

print STDERR "startnum: $startnum stopnum: $stopnum\n";

if ($step == 1 or $step == 3){
    if ($opt_a){
        die "Option -a cannot exceed the number of query genomes.\n" if $maxmiss >= $g_count;
    }

    print "Calculating pangenome, core genome, and new genome sequence at each permutation ...\n";
    print "\tPerforming pan/core/new calculations ...\n";
    ## Determine pangenome, core genome, and new genome from permutations

    open (my $pan_raw_out, "> $pref\_pangenome_raw.txt");
    if ($opt_a){
        print $pan_raw_out "n_genomes\tgenome_codes\tpangenome_size\tcore_genome_size\ta$maxmiss\_core_genome_size\tnew_genome_size\n";
    } else {
        print $pan_raw_out "n_genomes\tgenome_codes\tpangenome_size\tcore_genome_size\tcore_$pctcore\_pct_core_genome_size\tnew_genome_size\n";
    }
    open (my $pan_avg_out, "> $pref\_pangenome_avg.txt");
    if ($opt_a){
        print $pan_avg_out "n_genomes\tavg_pangenome_size\tavg_core_genome_size\ta$maxmiss\_core_genome_size\tavg_new_genome_size\n";
    } else {
        print $pan_avg_out "n_genomes\tavg_pangenome_size\tavg_core_genome_size\tcore_$pctcore\_pct_core_genome_size\tavg_new_genome_size\n";
    }
    my @raw_array;
    my $num_threads_running = 0;
    for my $array_size ($startnum .. $stopnum){
        my %raw_results;
        my $dp = `perl $perm_count $g_count $array_size N N Y`;
        #my $dp = fac($g_count) / (fac($array_size - 1) * fac($g_count - $array_size));
        print STDERR "Number of potential permutations of $array_size genome(s): $dp\n";
        my %seen_array;
        my $stop = $dp;
        $stop = $minperm if ($minperm != 0);
        my $stop_delta = 0;
        my ($last_pan_avg, $last_core_avg, $last_new_avg) = (0) x 3;
        my $stat_run;
        my $finished_recs = 0;
        my %running_pids;
        while (keys %seen_array < $dp){
            last if ($finished_recs >= $stop and $stop_delta);
            my @array;
            my $index = int rand ($g_count) + 1;
            if ($array_size == 1){
                redo if $seen_array{$index}++;
                push @array, $index;
                #print "** $index\n";
            } else {
                my %seen_cand;
                for (2 .. $array_size){
                    my $cand = int rand($g_count) + 1;
                    redo if $cand == $index;
                    redo if $seen_cand{$cand}++;
                    push @array, $cand;
                }
                @array = sort {$a <=> $b} @array;
                unshift @array, $index;
                my $string = join(',', @array);
                redo if $seen_array{$string}++;
                #print "** $string\n";
            }
            my $recs = keys %seen_array;
            #print STDERR join(',', @array), "\n";

            if ($num_threads_running < $threads){
                my $pid = fork;
                if (0 == $pid){
                    open (my $out, "> $$.tmp.txt");
                    my ($pan, $core, $core_pct, $new) = (0) x 4;
                    my $pan_miss_string;
                    my @test_array = @array;
                    shift @test_array;
                    my $corenew_miss_string = join("|", @test_array);
                    for my $j (0 .. $#array){
                        my $ref = $array[$j];
                        $pan += $gen_size{$ref} if $j == 0;
                        my $last = $g_count - $j;
                        for my $k (1 .. $last){
                            next if !$perm_hash{$ref}{$k};
                            my %test_hash = %{$perm_hash{$ref}{$k}};
                            if ($j == 0){
                                if ($array_size > 1){
                                    foreach my $missing (keys %test_hash){
                                        my $count = 0;
                                        $count++ while ($missing =~ /\b(?:$corenew_miss_string)\b/g ); #count how many of the array values are in $missing
                                        my $val;
                                        if ($count == 0 or $count <= $maxmiss or $count <= ($core_frac * $array_size) or $count == ($array_size - 1)){ #only looks up $val if it is needed. Maybe save some time.
                                            $val = $test_hash{$missing};
                                            $core += $val if $count == 0;
                                            if ($opt_a){
                                                $core_pct += $val if ($count <= ($maxmiss) or $count == 0);
                                            } else {
                                                $core_pct += $val if $count <= ($core_frac * $array_size);
                                            }
                                            $new += $val if $count == ($array_size - 1);
                                        }
                                    }
                                }
                            } else {

                                ##test
                                #foreach my $missing (grep /$pan_miss_string/, keys %test_hash){
                                #    $pan += $test_hash{$missing};
                                #}
                                ###end test

                                ##test
                                #foreach my $missing (keys %test_hash){
                                #    if ($missing =~ m/$pan_miss_string/){
                                #        $pan += $test_hash{$missing};
                                #    }
                                #}
                                ###end test

                                ##Benchmarked, this is almost twice as fast as the two methods above (grep with no array, or cycling through all keys). Can't exactly figure out why...
                                my @matches = grep /$pan_miss_string/, keys %test_hash;
                                if (@matches){
                                    for my $l (0 .. $#matches){
                                        $pan += $test_hash{$matches[$l]};
                                    }
                                }
                            }
                        }
                        $pan_miss_string .= "(?=.*,$ref,)";
                    }
                    my $p_string = join(",", @array);
                    print $out "$p_string\t$pan\t$core\t$new\t$core_pct\n";
                    #print STDERR "\t$p_string\t$pan\t$core\t$new\t$core_pct\n" if $array_size == 1; #debug
                    close $out;
                    exit(0);
                }
                $num_threads_running++;
                $running_pids{$pid}++;
            }
            if ($num_threads_running == $threads){ #stop to wait for a thread to finish before starting another
                my $pid = wait;
                open (my $in, "< $pid.tmp.txt") or die "Can't open $pid.tmp.txt: $!\n";
                $finished_recs++;
                print STDERR "\r\t$array_size: $finished_recs";
                while (my $line = <$in>){
                    chomp $line;
                    my ($string, $pan, $core, $new, $core_pct) = split('\t', $line);
                    print $pan_raw_out "$array_size\t$string\t$pan\t$core\t$core_pct\t$new\n";
                    push @{$raw_results{'pan'}}, $pan;
                    push @{$raw_results{'core'}}, $core;
                    push @{$raw_results{'new'}}, $new;
                    push @{$raw_results{'corepct'}}, $core_pct;
                    if ($finished_recs >= $stop){
                        my $pan_avg = average(\@{$raw_results{'pan'}});
                        my $core_avg = average(\@{$raw_results{'core'}});
                        my $new_avg = average(\@{$raw_results{'new'}});
                        my $d_count = 0;
                        $d_count++ if pct_diff($pan_avg, $last_pan_avg) < $maxdelt;
                        if ($array_size > 1){
                            $d_count++ if pct_diff($core_avg, $last_core_avg) < $maxdelt;
                            $d_count++ if pct_diff($new_avg, $last_new_avg) < $maxdelt;
                            print STDERR "\t".sprintf("%.1f", $pan_avg)." (".sprintf("%.4f", pct_diff($pan_avg, $last_pan_avg))."%)\t".sprintf("%.1f", $core_avg)." (".sprintf("%.4f", pct_diff($core_avg, $last_core_avg))."%)\t".sprintf("%.1f", $new_avg)." (".sprintf("%.4f", pct_diff($new_avg, $last_new_avg))."%)" if $stat_run;
                            $stat_run = 1;
                        }
                        $stop_delta = 1 if $d_count == 3;
                        $stop_delta = 1 if $d_count == 1 and $array_size == 1;
                        ($last_pan_avg, $last_core_avg, $last_new_avg) = ($pan_avg, $core_avg, $new_avg);
                    }
                }
                close ($in);
                unlink ("$pid.tmp.txt");
                delete $running_pids{$pid};
                $num_threads_running--;
            }
        }
        if ($stop_delta == 1){
            foreach my $pid (keys %running_pids){
                #print STDERR "\tKilling $pid\n";
                kill ('KILL', $pid) if kill ('ZERO', $pid);
                unlink "$pid.tmp.txt" if -e "$pid.tmp.txt";
                delete $running_pids{$pid};
            }
            while (my $pid = wait){ #round up stragglers
                last if $pid < 0;
                unlink "$pid.tmp.txt" if -e "$pid.tmp.txt";
            }
        } else {
            while (my $pid = wait){
                last if $pid < 0;
                if ($minperm == 0 or $dp < $minperm){
                    open (my $in, "< $pid.tmp.txt") or die "Can't open $pid.tmp.txt: $!\n";
                    $finished_recs++;
                    print STDERR "\r\t$array_size: $finished_recs";
                    while (my $line = <$in>){
                        chomp $line;
                        my ($string, $pan, $core, $new, $core_pct) = split('\t', $line);
                        print $pan_raw_out "$array_size\t$string\t$pan\t$core\t$core_pct\t$new\n";
                        push @{$raw_results{'pan'}}, $pan;
                        push @{$raw_results{'core'}}, $core;
                        push @{$raw_results{'new'}}, $new;
                        push @{$raw_results{'corepct'}}, $core_pct;
                    }
                    close ($in);
                }
                unlink ("$pid.tmp.txt");
            }
        }
        print STDERR "\n";
        $num_threads_running = 0;
        print "\tOutputting average values ... \n";
        #output average values
        print $pan_avg_out "$array_size\t";
        my @p_lengs = @{$raw_results{'pan'}};
        my ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq, $stdev, $ste) = stats(\@p_lengs);
        print $pan_avg_out "$rounded_mean\t";
        my @c_lengs = @{$raw_results{'core'}};
        ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq, $stdev, $ste) = stats(\@c_lengs);
        print $pan_avg_out "$rounded_mean\t";
        my @cp_lengs = @{$raw_results{'corepct'}};
        ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq, $stdev, $ste) = stats(\@cp_lengs);
        print $pan_avg_out "$rounded_mean\t";
        my @n_lengs = @{$raw_results{'new'}};
        ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq, $stdev, $ste) = stats(\@n_lengs);
        print $pan_avg_out "$rounded_mean\n";
    }
    close ($pan_raw_out);
    close ($pan_avg_out);
}

###Let's just start from scratch here...
if ($step == 2 or $step == 3){
    print "Calculating core genome sizes at every definition of core using $threads threads ...\n";
    my $num_threads_running = 0;
    my %all_values;
    open (my $perm_out, "> $pref\_cores_raw.txt");
    print $perm_out "a_val\tref_genomes\tcore_size\n";
    my $dp = fac($g_count);
    print STDERR "$dp possible permutations of $g_count genomes\n";
    my %seen_array;
    my $stop = $dp;
    $stop = $c_minperm if ($c_minperm != 0);
    my $stop_delta = 0;
    my $last_pan_avg = 0;
    my $finished_recs = 0;
    my %running_pids;
    while (keys %seen_array < $dp){
        last if ($finished_recs >= $stop and $stop_delta);
        my @perm;
        my %seen_cand;
        for my $j (1 .. $g_count){
            my $cand = int rand($g_count) + 1;
            redo if $seen_cand{$cand}++;
            push @perm, $cand;
        }
        my $string = join(',', @perm);
        redo if $seen_array{$string}++;
        my $recs = keys %seen_array;
        if ($num_threads_running < $threads){
            my $pid = fork;
            if (0 == $pid){
                my $outfile = "$$.tmp_core.txt";
                open (my $out, ">", $outfile) or die "Can't open $outfile: $!\n";
                my %core;
                my @test_array = @perm;
                shift @test_array;
                my $corenew_miss_string = join("|", @test_array);
                for my $j (0 .. 0){
                    my $ref = $perm[$j];
                    my $last = $g_count - $j;
                    for my $k (1 .. $last){
                        next if !$perm_hash{$ref}{$k};
                        my %test_hash = %{$perm_hash{$ref}{$k}};
                        foreach my $missing (keys %test_hash){
                            my $count = 0;
                            $count++ while ($missing =~ /\b(?:$corenew_miss_string)\b/g);
                            my $val = $test_hash{$missing};
                            for my $l ($count .. ($g_count - 1)){ #generate core size at every possible "a-value" for this permutation
                                $core{"$l"} += $val;
                            }
                        }
                    }
                }
                for my $i (0 .. ($g_count - 1)){
                    if ($core{"$i"}){
                        my $val = $core{"$i"};
                        my $gens = join(",", @perm[0..$i]);
                        print $out "$i:$gens:$val\n";
                    }
                }
                close $out;
                exit (0);
            }
            $num_threads_running++;
            $running_pids{$pid}++;
        }
        if ($num_threads_running == $threads){ #stop to wait for a thread to finish before starting another
            my $pid = wait;
            open (my $in, "< $pid.tmp_core.txt") or die "Can't open $pid.tmp_core.txt: $!\n";
            $finished_recs++;
            print STDERR "\r\tpermutation: $finished_recs";
            delete $running_pids{$pid};
            while (my $line = <$in>){
                chomp $line;
                my ($aval, $gens, $core) = split(':', $line);
                print $perm_out "$aval\t$gens\t$core\n";
                push @{$all_values{"$aval"}}, $core;
                if ($finished_recs >= $stop){
                    my $test = $g_count - 1; #we'll check the delta around the pangenome (a-value = $g_count - 1);
                    my $pan_avg = average(\@{$all_values{"$test"}});
                    my $d_count = 0;
                    $stop_delta = 1 if pct_diff($pan_avg, $last_pan_avg) < $maxdelt;
                    $last_pan_avg = $pan_avg;
                }
            }
            close ($in);
            unlink ("$pid.tmp_core.txt");
            $num_threads_running--;
        }
    }
    if ($stop_delta == 1){
        foreach my $pid (keys %running_pids){
            #print STDERR "\tkilling $pid\n";
            kill ('KILL', $pid) if kill ('ZERO', $pid);
            unlink "$pid.tmp_core.txt" if -e "$pid.tmp_core.txt";
            delete $running_pids{$pid};
        }
        while (my $pid = wait){ #round up stragglers
            last if $pid < 0;
            unlink "$pid.tmp.txt" if -e "$pid.tmp.txt";
        }
    } else {
        while (my $pid = wait){
            last if $pid < 0;
            if ($c_minperm == 0 or $dp < $c_minperm){
                open (my $in, "< $pid.tmp_core.txt") or die "Can't open $pid.tmp_core.txt: $!\n";
                $finished_recs++;
                print STDERR "\r\tpermutation: $finished_recs";
                while (my $line = <$in>){
                    chomp $line;
                    my ($aval, $gens, $core) = split(':', $line);
                    print $perm_out "$aval\t$gens\t$core\n";
                    push @{$all_values{"$aval"}}, $core;
                }
                close ($in);
            }
            unlink ("$pid.tmp_core.txt");
        }
    }
    print STDERR "\n";
    $num_threads_running = 0;
    print STDERR "\tOutputting average values ... \n";
    open (my $perm_a_out, "> $pref\_cores_avg.txt");
    print $perm_a_out "a_val\tnum_permutations\taverage_size\tstd_error_average_size\n";
    for my $i (0 .. ($g_count - 1)){
        my @lengs = @{$all_values{"$i"}};
        my ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq, $stdev, $ste) = stats(\@lengs);
        #print "Average: $i\t$num\t$rounded_mean\t$ste\n";
        print $perm_a_out "". $i ."\t$num\t$rounded_mean\t$ste\n";
    }
    close ($perm_a_out);
}



#if ($step == 2 or $step == 3){
#    print "Calculating all possible core genome permutations using $threads threads ...\n";
#    ##### GENERATE ALL POSSIBLE PERMUTATIONS OF @to_use
#    my @new_order = (1..$g_count);
#    my $num_threads_running = 0;
#    my %all_values;
#    open (my $perm_out, "> $pref\_cores_raw.txt");
#    print $perm_out "a_val\tref_genomes\tcore_size\n";
#    #use Algorithm::Permute;
#    #my @permutations;
#    for my $i (1 .. $g_count){
#        my $dp = fac($g_count) / fac($g_count - $i);
#        my %seen_array;
#        my ($total, $count) = (0) x 2;
#        my $stop = $dp;
#        $stop = $minperm if ($minperm != 0);
#        my $stop_delta = 0;
#        while (keys %seen_array < $dp){
#            last if (keys %seen_array >= $stop and $stop_delta);
#            my @perm;
#            my %seen_cand;
#            for my $j (1 .. $i){
#                my $cand = int rand($g_count) + 1;
#                redo if $seen_cand{$cand}++;
#                push @perm, $cand;
#            }
#            my $string = join(',', @perm);
#            redo if $seen_array{$string}++;
#            my $recs = keys %seen_array;
#            if ($num_threads_running < $threads){
#                my $pid = fork;
#                if (0 == $pid){
#                    my $outfile = "$$.tmp_perm.txt";
#                    my $status = permute(\@perm, $outfile);
#                    exit($status);
#                }
#                $num_threads_running++;
#            } else {
#                my $pid = wait;
#                open (my $in, "< $pid.tmp_perm.txt") or die "Can't open $pid.tmp_perm.txt: $!\n";
#                my ($num, $miss, $core_size);
#                while (my $line = <$in>){
#                    chomp $line;
#                    ($num, $miss, $core_size) = split('\t', $line);
#                    push @{$all_values{$num}}, $core_size;
#                    print $perm_out "". $num - 1 ."\t$miss\t$core_size\n";
#                    #print STDERR "\r$miss\t$core_size";
#                }
#                close ($in);
#                unlink ("$pid.tmp_perm.txt");
#                $num_threads_running--;
#
#                my $old_avg = 0;
#                if ($total > 0){
#                    $old_avg = $total / $count;
#                }
#                $total += $core_size;
#                $count ++;
#                print STDERR sprintf("\r$i %10s\t$miss", $count);
#                my $avg = $total / $count;
#                if ($count >= $stop){
#                    my $delt = pct_diff($old_avg, $avg);
#                    print STDERR "\n\t$avg ($delt%)\n";
#                    if ($delt <= $maxdelt){
#                        $stop_delta = 1;
#                        next;
#                    }
#                }
#                $pid = fork;
#                if (0 == $pid){
#                    my $outfile = "$$.tmp_perm.txt";
#                    my $status = permute(\@perm, $outfile);
#                    exit($status);
#                }
#                $num_threads_running++;
#            }
#        }
#        while (my $pid = wait){
#            last if $pid < 0;
#            if ($stop_delta == 0){
#                my ($num, $miss, $core_size);
#                open (my $in, "< $pid.tmp_perm.txt") or die "Can't open $pid.tmp_perm.txt: $!\n";
#                while (my $line = <$in>){
#                    chomp $line;
#                    ($num, $miss, $core_size) = split('\t', $line);
#                    push @{$all_values{$num}}, $core_size;
#                    print $perm_out "". $num - 1 ."\t$miss\t$core_size\n";
#                    #print STDERR "\r$miss\t$core_size";
#                }
#                close ($in);
#                my $old_avg = 0;
#                if ($total > 0){
#                    $old_avg = $total / $count;
#                }
#                $total += $core_size;
#                $count ++;
#                print STDERR sprintf("\r$i %10s\t$miss", $count);
#                my $avg = $total / $count;
#                if ($count >= $stop){
#                    my $delt = pct_diff($old_avg, $avg);
#                    print STDERR "\n\t$avg ($delt%)\n";
#                    if ($delt <= $maxdelt){
#                        $stop_delta = 1;
#                        next;
#                    }
#                }
#            }
#            unlink ("$pid.tmp_perm.txt");
#        }
#        $num_threads_running = 0;
#        print STDERR "\n";
#    }
#    print STDERR "\n";
#
#    ##calculate averages
#    open (my $perm_a_out, "> $pref\_cores_avg.txt");
#    print $perm_a_out "a_val\tnum_permutations\taverage_size\tstd_error_average_size\n";
#    for my $i (1 .. $g_count){
#        my @lengs = @{$all_values{$i}};
#        my ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq, $stdev, $ste) = stats(\@lengs);
#        #print "Average: $i\t$num\t$rounded_mean\t$ste\n";
#        print $perm_a_out "". $i - 1 ."\t$num\t$rounded_mean\t$ste\n";
#    }
#    close ($perm_a_out);
#}
print STDERR "Done\n";

#------------------------------------------------------------------------------
sub fac {
  $_[0]>1?$_[0]*fac($_[0]-1):1;
}

sub permute {
    my @pmuts = @{ $_[0] };
    my $ofile = $_[1];
    open (my $out, ">", $ofile);
    #for my $i (0 .. $#pmuts){
        my @to_use = @pmuts;
        my $use_count = scalar @to_use;
        my $depth = $g_count - $use_count + 1;

        my %keep_hash;
        my @p_build;
        my $core_size = 0;
        for my $m (0 .. $#to_use){
            my $ref = $to_use[$m];
            for my $j ($depth .. ($g_count - $m)){
                my $last_val = 1;
                $last_val = 2 if $ref == 1;
                push @p_build, $last_val;
                my $size = scalar @p_build;
                if ($size == ($g_count - $j)){
                    my $keep_count = 0;
                    for my $k (0 .. $#p_build){
                        $keep_count++ if $keep_hash{$p_build[$k]};
                    }
                    if ($keep_count == keys %keep_hash){
                        my $missing = ",".join(",", @p_build).",";
                        my $check_size = $g_count - scalar @p_build;
                        $core_size += $perm_hash{$ref}{$check_size}{$missing} if ($perm_hash{$ref}{$check_size}{$missing});
                        #debug
                        #print "\t$i\t$ref\t$check_size\t$missing\t$core_size\n";
                        #\debug
                    }
                }
                my $final = $g_count;
                $final = $g_count - 1 if $ref == $final;
                my ($stop, $next);
                while (!$stop){
                    while ($last_val < $final){
                        $next = $last_val + 1;
                        $next = $last_val + 2 if $next == $ref;
                        push @p_build, $next;
                        my $size = scalar @p_build;
                        if ($size == ($g_count - $j)){
                            my $keep_count = 0;
                            for my $k (0 .. $#p_build){
                                $keep_count++ if $keep_hash{$p_build[$k]};
                            }
                            if ($keep_count == keys %keep_hash){
                                my $missing = ",".join(",", @p_build).",";
                                my $check_size = $g_count - scalar @p_build;
                                $core_size += $perm_hash{$ref}{$check_size}{$missing} if ($perm_hash{$ref}{$check_size}{$missing});
                                #debug
                                #print "\t$i\t$ref\t$check_size\t$missing\t$core_size\n";
                                #\debug
                            }
                        }
                        $stop = 1 if (scalar @p_build == 1 and $p_build[0] == $final);
                        $last_val = $next;
                    }
                    pop (@p_build);
                    $last_val = pop(@p_build);
                }
            }
            if ($m == 0){
                $core_size += $perm_hash{$ref}{$g_count}{",0,"} if ($perm_hash{$ref}{$g_count}{",0,"});
                #debug
                #print "\t$i\t$ref\t$g_count\t,0,\t$core_size\n";
                #\debug
            }
            $keep_hash{$ref}++;
        }
        #print scalar @to_use, "\t", join(",", @to_use), "\t$core_size\n";
        print $out scalar @to_use, "\t", join(",", @to_use), "\t$core_size\n";
    #}
    close ($out);
    return (0);
}

sub average {
    my @lengths = @{$_[0]};
    my $num = scalar @lengths;
    my $sum = 0;
    for my $i (0 .. $#lengths){
        $sum += $lengths[$i];
    }
    my $avg = $sum / $num;
    return ($avg);
}

sub pct_diff{
    my ($one, $two) = ($_[0], $_[1]);
    my $diff = abs ($two - $one);
    my $avg = ($one + $two) / 2;
    my $pct_diff = 100 * ($diff / $avg);
    return ($pct_diff);
}


sub stats{
	my @lengths = @{$_[0]};
	my ($sum, $num, $min, $maxi, $mean, $median, $mode, $mode_freq, $stdev, $ste);
	my %seen;
	my @sorted_leng = sort {$a <=> $b} @lengths;
	for my $i (0 .. $#sorted_leng){
		$sum += $sorted_leng[$i];
		$seen{$sorted_leng[$i]}++;
	}
	$num = $#sorted_leng + 1;
	$min = $sorted_leng[0];
	$maxi = $sorted_leng[$#sorted_leng];
	$mean = $sum/$num;
	my $rounded_mean = sprintf("%.2f", $mean);
	my @modes;
	foreach my $leng (sort {$seen{$b} <=> $seen{$a}} keys %seen) {
		push @modes, ([$leng, $seen{$leng}]);
	}
	$mode = $modes[0][0];
	$mode_freq = $modes[0][1];
	my $mid = int @sorted_leng/2;
	if (@sorted_leng % 2){
		$median = $sorted_leng[$mid];
	} else {
		$median = ($sorted_leng[$mid-1] + $sorted_leng[$mid])/2;
	}
        my $sqtotal = 0;
        foreach (@lengths){
            $sqtotal += ($mean - $_) ** 2;
        }
        $stdev = ($sqtotal / (@lengths - 1)) ** 0.5;
        $ste = $stdev / sqrt(scalar @lengths);
	return ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq, $stdev, $ste);
}
