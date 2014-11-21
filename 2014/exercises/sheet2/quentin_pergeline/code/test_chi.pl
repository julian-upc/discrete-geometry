#!/usr/bin/perl

# perl script to execute the c++ program with all pair of test in test/input

use 5.010;
use warnings;
use strict;

my $dir = "test/input";
my $output = "test/output";

if (! -e "./sheet2") {say "Error: sheet2 does not exist. You should probably execute make before."; exit;}


opendir DIR, $dir || die $!;
my @inputfiles = grep {$_ =~ m/^matrix_[a-zA-Z]+\w*\.txt$/} readdir DIR; #goes in test/input and list all the text files of the form: matrix_*.txt
foreach my $i (@inputfiles) {
	my ($id1) = $i =~ /_(.*?)\.txt/;
	foreach my $j (@inputfiles) {
		my ($id2) = $j =~ /_(.*?)\.txt/;

		say "sheet2.exe  ". $dir."/".$i." ".$dir."/".$j." > ".$output."/result_".$id1."_".$id2.".txt";
		my @args = ("./sheet2", $dir."/".$i, $dir."/".$j, $output."/result_".$id1."_".$id2.".txt"); #execute the c++ program with all pair of well-titled test files available in test/input
		system(@args) == 0 or say "system @args failed: $?"
	}
}


