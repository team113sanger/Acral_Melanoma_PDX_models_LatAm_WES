#!/usr/bin/env perl

use strict;
use warnings;

my $totbases = 0;
my %counts;
my $totsum = 0;

my @limits = (11, 21, 31, 41, 51, 61, 71, 81, 91, 101, 111, 121, 131, 141, 151);
while (<>) {
	chomp;
	$totbases++;
	my @col = split(/\t/, $_);
	my $dp = $col[2];
	$totsum += $dp;
	foreach my $lim (reverse(@limits)) {
		if ($dp >= $lim) {
			$counts{$lim}++;
#			last;
		}
	}
}

# Summary
my @printvals;
foreach my $lim (@limits) {
	my $percent = sprintf ("%.2f", 100 * $counts{$lim} / $totbases);
	push @printvals, $percent;	
}

my $mean = sprintf("%.2f", $totsum / $totbases);
 
print join("+\t", @limits) . "+\tCov_Mean\n";
print join("\t", @printvals) . "\t$mean\n";

