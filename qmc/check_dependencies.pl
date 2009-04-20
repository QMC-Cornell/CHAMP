#!/usr/bin/perl
# Written by Frank Petruzielo
# 3/26/09
# updated on 4/15/09 to account for the possibility of multiple objects created in one build routine
# and to allow for the use of multiple files
# If an object is used, make certain that it is declared in an object needed statement
# It is necessary to use the use only syntax for any external objects.
# mark off the beginning of the objects with the comment !begin objects (including those from other modules)
# mark off the end of the objects with the comment !end objects
# see for example backflow_mod.f90

use diagnostics;
use strict;
use warnings;
use 5.010;

#don't want subroutines to have access to variables so change scope with brackets
{
    die "Usage: ./check_dependencies file1.f90 file2.f90 ... \n" unless (@ARGV > 0);

    #only want .f90 files
    my @file_names = grep { /.f90$/ } @ARGV;
    die "Requires at least one .f90 file \n" unless (@file_names > 0);
    for (@file_names) { 
	say "Checking $_...";
	#extract contents from all build routines
	my @file_contents = get_bld_routine_contents($_);	
	#what objects are available to the build routines?
	my @obj_available = get_objects_available($_);
        #which objects are requested with object needed in each bld routine (return as reference)
	my $obj_needed_in_sub_href = {}; 
	$obj_needed_in_sub_href = get_objects_needed(@file_contents);
	#which objects are requested with object provide in node in each bld routine (return as reference)
	my $obj_provide_in_node_in_sub_href = {};
	$obj_provide_in_node_in_sub_href = get_objects_provide_in_node(@file_contents);	
	#which objects are built with object create in each bld routine (return as reference)
	my $obj_create_in_sub_href = {};
	$obj_create_in_sub_href = get_objects_create(@file_contents);	
	#which objects of those available are used in each build routine (ignoring those created in the build routine)?
	my $obj_used_in_sub_href = {};
	$obj_used_in_sub_href = get_objects_used(\@file_contents,\@obj_available,$obj_create_in_sub_href);	
	#are any objects used that are not requested with object needed or object provide in node?
	find_bugs($obj_used_in_sub_href,$obj_needed_in_sub_href,$obj_provide_in_node_in_sub_href);
    }
}

#==============================================================
sub get_bld_routine_contents {
#==============================================================
# Frank Petruzielo 4/15/09
# Extract contents from all of the build routines in a .f90 file
# Expects a .f90 file name as input
#==============================================================
    my $file_name = pop @_;
    open my $file, $file_name or die "Can't open $file_name: $!\n";
    chomp(my @file_contents = <$file>);
    my $in_bld_routine = 0;
    my @bld_routines_contents;
    # ignore blank lines
    @file_contents = grep {!/^\s*$/} @file_contents;
    for (@file_contents) {
	#in build routine
	if (/^\s*subroutine.*_bld/) {
	    $in_bld_routine = 1;
	}
	#add lines in bld routines that are not comment lines
	if ($in_bld_routine and !/^\s*!/) {
	    push @bld_routines_contents, $_ ;
	}
	#leaving build routine
	if (/^\s*end\s+subroutine.*_bld/) {
	    $in_bld_routine = 0;
	}
    }
    return @bld_routines_contents;
}
#==============================================================

sub get_objects_available {
#==============================================================
# Frank Petruzielo 4/15/09
# Determine which objects are available to the build routines
# by looking at the global variables declared and asked for with the module use only syntax
# Expects a .f90 file name as input
#==============================================================
    my $file_name = pop @_;
    open my $file, $file_name or die "Can't open $file_name: $!\n";
    chomp(my @file_contents = <$file>);
    my $start=0;
    my @obj_available;
    for (@file_contents) {
	if (/^\s*contains\s*$/) {
	    #no more objects available after contains statement
	    last;
	}
	if (/!\s*begin\s+objects/) {
	    #this is the flag we use to designate the beginning of the objects
	    $start = 1;
	}
	if ($start) {
	    if (s/^.*::// or s/^.*only\s*://) {
		#remove declaration statement or use only statment
		push @obj_available, /\w+/g; #don't want dimensions (:,:) or commas so only take "perl words"
	    }
	}
    }
    return @obj_available;
}
#==============================================================

sub get_objects_needed {
#==============================================================
# Frank Petruzielo 4/15/09
# Determine which objects are requested with object needed each bld routine (return as reference)
# Expects the contents of all bld routines of a particular file as input
#==============================================================
    my @bld_routines_contents = @_;
    my $obj_needed_in_sub_href={};
    my $in_bld_routine;
    for (@bld_routines_contents) {
	if (/^\s*subroutine\s+(.*)_bld/) {
	    # in new build routine so store the name
	    $in_bld_routine = $1;
	}	
	if (/object_needed\s*\(\s*\'(.*)\'\s*\)/) {
	    #found an object needed request in this routine so mark this object in this routine with a 1 (placeholder)
	    $obj_needed_in_sub_href -> {$in_bld_routine} -> {$1} = 1;
	}
    }
    return $obj_needed_in_sub_href;
}
#==============================================================

#==============================================================

sub get_objects_provide_in_node {
#==============================================================
# Frank Petruzielo 4/15/09
# Determine which objects are requested with object provide in node in each bld routine (return as reference)
# Expects the contents of all bld routines of a particular file as input
#==============================================================
    my @bld_routines_contents = @_;
    my $obj_provide_in_node_in_sub_href={};
    my $in_bld_routine;
    for (@bld_routines_contents) {
	if (/^\s*subroutine\s+(.*)_bld/) {
	    # in new build routine so store the name
	    $in_bld_routine = $1;
	}	
	if (/object_provide_in_node\s*\(\s*lhere\s*,\s*\'(.*)\'\s*\)/) {
	    #found an object provide in node request request in this routine so mark this object in this routine with a 1 (placeholder)
	    $obj_provide_in_node_in_sub_href -> {$in_bld_routine} -> {$1} = 1;
	}
    }
    return $obj_provide_in_node_in_sub_href;
}
#==============================================================

sub get_objects_create {
#==============================================================
# Frank Petruzielo 4/15/09
# Determine which objects are built with object create in each bld routine (return as reference)
# Expects the contents of all bld routines of a particular file as input
#==============================================================
    my @bld_routines_contents = @_;
    my $obj_create_in_sub_href={};
    my $in_bld_routine;
    for (@bld_routines_contents) {
	if (/^\s*subroutine\s+(.*)_bld/) {
	    # in new build routine so store the name
	    $in_bld_routine = $1;
	}	
	if (/object_create\s*\(\s*\'(.*)\'\s*\)/) {
	    #found an object create request in this routine so mark this object in this routine with a 1 (placeholder)
	    $obj_create_in_sub_href -> {$in_bld_routine} -> {$1} = 1;
	}
    }
    return $obj_create_in_sub_href;
}
#==============================================================

sub get_objects_used {
#==============================================================
# Frank Petruzielo 4/15/09
# Determine which objects of those available are used in each bld routine (return as reference)
# Ignores those created in the build routine as it should
# Expects the contents of all bld routines of a particular file, the available objects, and those objects created in the nodes as input (must use references)
#==============================================================
    my $bld_routines_contents_ref = shift;
    my $obj_available_ref = shift;
    my $obj_create_in_sub_href = shift;
    my @bld_routines_contents = @$bld_routines_contents_ref;
    my @obj_available = @$obj_available_ref;
    my $obj_used_in_sub_href={};
    my $in_bld_routine;
    for (@bld_routines_contents) {
	if (/^\s*subroutine\s+(.*)_bld/) {
	    # in new build routine so store the name
	    $in_bld_routine = $1;
	}	
	#loop over available objects checking if they are used
	for my $temp_obj (@obj_available) {
	    #ignore use of objects created in the node
	    unless ($obj_create_in_sub_href -> {$in_bld_routine} -> {$temp_obj}) {
		#check for object
		/\b$temp_obj\b/ && ($obj_used_in_sub_href -> {$in_bld_routine} -> {$temp_obj} = 1 );
	    }
	}
    }
    return $obj_used_in_sub_href;
}
#==============================================================

sub find_bugs {
#==============================================================
# Frank Petruzielo 4/15/09
# Determine which objects are used but not requested with object needed or object provide_in_node, and print error messages
# Expects those objects used in the nodes, those objects requested with object needed, and those objects requested with object_provide_in_node as input
#==============================================================
    my $obj_used_in_sub_href = shift;
    my $obj_needed_in_sub_href = shift;
    my $obj_provide_in_node_in_sub_href = shift;
    #check if using objects that aren't asked for with need or provide_in_node
    for my $sub_name (sort keys %$obj_used_in_sub_href) {
	for my $obj (sort keys %{$obj_used_in_sub_href -> {$sub_name}}) {
	    unless ($obj_needed_in_sub_href -> {$sub_name} -> {$obj}) {
		#object not requested with object needed
		unless ($obj_provide_in_node_in_sub_href -> {$sub_name} -> {$obj}) {
		    #object requested with object provide in node
		    say "BUG:In building of $sub_name, $obj is used but not asked for!";
		}
		#else {
		#    object requested with object provide in node
		#    say "WARNING:In building of $sub_name, $obj is requested by object_provide_in_node!";
		#}
	    }
	}
    }
}


