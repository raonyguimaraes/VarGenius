#!/usr/bin/perl

#Libraries
use strict;
use warnings;

#Get the file with the parameters for the module
my $command = $ENV{COMMAND};

 
unless ( (system $command) == 0) {
 if ($? == -1) {
			print "failed to execute: $!\n";
	}
	elsif ($? & 127) {
			printf "child died with signal %d, %s coredump\n",
					($? & 127),  ($? & 128) ? 'with' : 'without';
	}
	else {
			printf "child exited with value %d\n", $? >> 8;
			
	}	
}



