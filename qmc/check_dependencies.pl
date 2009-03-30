#!/usr/bin/perl
# Written by Frank Petruzielo
# 3/26/09
# If an object is used, make certain that it is declared in an object needed statement
# It is necessary to use the use only syntax for any external objects.
# make off the beginning of the objects with the comment !objects (including those from other modules)
# see for example backflow_mod.f90

use diagnostics;
use strict;
use warnings;
use 5.010;

die "Requires a single file\n" unless (@ARGV == 1);
die "Requires a .f90 file \n" unless ($ARGV[0] =~ /\.f90/);
die "$ARGV[0] does not exist\n" unless -e $ARGV[0];

chomp(my @file_contents = <>);
my $in_bld_routine = 0;
my @bld_routines;
for (@file_contents) {
    /^\s*subroutine.*_bld/ && ($in_bld_routine = 1);
    if ($in_bld_routine && !/^\s*!/) {
	push @bld_routines, $_ ;
    }
    /^\s*end\s+subroutine.*_bld/ && ($in_bld_routine = 0);
}


#which objects are available to the the bld routines
my @obj_available;
my $start=0;
for (@file_contents) {
    /^\s*contains\s*$/ && last;
    /!\s*objects/ && ($start=1);
    if ($start) {
	s/^.*::// && push @obj_available, /\w+/g;
	s/^.*only\s*:// && push @obj_available, /\w+/g;
    }
}

#determine objects called with object_needed and object_provide_in_node
my $obj_needed_in_sub_href={};
my $obj_provide_in_node_in_sub_href={};
for (@bld_routines) {
    /^\s*subroutine\s+(.*)_bld/ && ($in_bld_routine = $1);
    /object_needed\s*\(\s*\'(.*)\'\s*\)/ && ($obj_needed_in_sub_href -> {$in_bld_routine} -> {$1} = 1 );
    /object_provide_in_node\s*\(\s*lhere\s*,\s*\'(.*)\'\s*\)/ && ($obj_provide_in_node_in_sub_href -> {$in_bld_routine} -> {$1} = 1 );
}

#which objects are used in each node
my $obj_used_in_sub_href={};
for (@bld_routines) {
    /^\s*subroutine\s+(.*)_bld/ && ($in_bld_routine = $1);
    for my $temp_obj (@obj_available) {
	unless ($temp_obj eq $in_bld_routine) {
	    /\b$temp_obj\b/ && ($obj_used_in_sub_href -> {$in_bld_routine} -> {$temp_obj} = 1 );
	}
    }
}

#check is using objects that aren't asked for with need or provide_in_node
for my $sub_name (sort keys %$obj_used_in_sub_href) {
    for my $obj (sort keys %{$obj_used_in_sub_href -> {$sub_name}}) {
	if (! ($obj_needed_in_sub_href -> {$sub_name} -> {$obj}  || $obj_provide_in_node_in_sub_href -> {$sub_name} -> {$obj})) {
	    say "BUG:In building of $sub_name, $obj is used but not asked for!";
	}
#	elsif ( ! $obj_needed_in_sub_href -> {$sub_name} -> {$obj} ) {
#	    say "WARNING:In building of $sub_name, $obj is requested by object_provide_in_node!";
#	}
    }
}

