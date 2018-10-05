#!/usr/bin/perl
#To work at CINECA

use lib '/cineca/prod/compilers/perl/5.x/none/lib/perl5';

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED

#VarGenius - Variant Discovery and Annotation Tool
#Copyright (C) <2017>  <Francesco Musacchia>

use strict;
use warnings;
use Getopt::Long;

binmode STDOUT, ":utf8";
use utf8;
 
use JSON;

my $file = $ARGV[0];

my $json;
{
  local $/; #Enable 'slurp' mode
  open my $fh, "<", $file;
  $json = <$fh>;
  close $fh;
}

my $data = decode_json($json);
my $DemuxResults = $data->{"ConversionResults"}->{'DemuxResults'};
#print $response;

my $studies_hash;
#############################CODE TO USE FOR NIG JSON OUTPUT
foreach my $SampleId (@$DemuxResults){
	print "\n\n$SampleId";
	#print "\n\n".$arr_el->{'attributes'}->{'accession'}." ".$arr_el->{'attributes'}->{'name'}."\n"; 
	#$studies_hash->{$arr_el->{'attributes'}->{'name'}} = $arr_el->{'attributes'}->{'accession'};
	#foreach my $key (keys %$arr_el){
	 #print $key."\n"; 
	#}
}
