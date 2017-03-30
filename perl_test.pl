#!/usr/bin/perl

use List::Util first;

%hash = ("first" => 1, "second" => 2, "third" => 3);

my $str = "testing";

print %hash;
print "\n";

if ($str =~ m/ing/)
{
  print "bar\n";
}

if (!$hash{"fourth"})
{
  print "foo\n";
}





# use Bio::Root::Version;
# use Bio::SeqIO;
# use Bio::SearchIO;
# use Bio::Seq;
# use Bio::SeqFeature::Generic;
#
#
# my $input_seqs = Bio::SeqIO->new(-file=>$ARGV[0], -format=>'fasta');
#
# while (my $curr_seq = $input_seqs->next_seq())
# {
#
#   print $curr_seq->subseq(1660593,1661435)."\n".length $curr_seq->seq();
#   print"\n";
#
# }

#
# #first occurrence in a list (t/f)
# $var1 = "two";
# @arr = ("zero", "one", "two", "three");
#
# if ( first {$_ eq $var1} @arr )
# {
# print "foo\n";
# }
#
#
# #user input array (scalar then split to array)
# my $temp_str;
# chomp($temp_str = <>);
# my @arr = split(',', $temp_str);
#
# print "arr: @arr\n";
#
#
# #array deep copy
# my @arr1 = ("zero", "one", "two", "three");
# my @arr2 = @arr1;
#
# print "arr1: @arr1\narr2: @arr2\n\n";
#
# $arr1[1] = 2;
#
# print "arr1: @arr1\narr2: @arr2\n\n";
#
#
# #if-elsif-else
# my $valid_input = ""; # "" and 0 evaluate to false
#
# if (not($valid_input))
# {
# print "false";
# }
# else
# {
# print "true";
# }
#
#
# #Scope
# my $test = 1;
# print $test;
# test();
# print $test;
#
# sub test
# {
# $test = 2;
# }
#
#
# #Index and Substring
# my $STR = "Hello world!";
#
# print substr($STR,index($STR,"o"),4,"oops!");
# print "\n";
# print $STR;
# print "\n";
#
#
# #String concatenation
# my $STR = "foo";
# $STR .= ".bar";
#
# print $STR."\n"; #prints foo.bar
#
#
# #array length
# my @arr = (1,2,3,4,5);
# print "@arr length:".scalar @arr."\n";
#
# #array length with user input
# my $str;
#
# trim($str = <>);
# @arr = split(" ",$str);
#
# print scalar @arr."\n";
#
#
# #integer division
# my @arr = (0,1,2,3,4,5,6);
#
# my $temp = int(scalar @arr/2);
#
# print $temp."\n";
# #alternatively "use integer;" at top of file or in block
#
#
# #getting filenames from paths
# my $path = "/usr/local/bin/file.txt";
#
# $path =~ m/^.*\//;
#
# print $`."\n".$&."\n".$'."\n";
#
#
# #make a new directory
# mkdir("test"); #this function does error checking (no error if dir aleady exists)
#
#
# #Regex substitution
# my $str = "hello world";
# print $str."\n";
# $str =~ s/ /\//;
# print $str."\n";
#
#
# #File I/O (write mode, '<' for read mode, others available)
# open(OUTFILE,">test.txt") or die "Couldn't open file file.txt, $!\n";
# print OUTFILE "TEST!!!!\n" || die "Couldn't close file properly\n";
# close(OUTFILE);
