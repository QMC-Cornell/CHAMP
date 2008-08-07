#!/usr/bin/perl

use File::Basename;
use File::Copy;

# Documentation
if ($#ARGV < 0){
 print "Usage: ",basename($0)," list_of_champ_output_files\n\n";
 print "Description: reduce size of champ output files\n";
 print "by removing printings of iterations bewteen 'BEGINNING OF' and 'END OF'\n";
exit;
}

# Arguments
@allfiles=@ARGV;

# loop over all files
foreach $file (@allfiles){
 print "examining $file\n";

# get timestamp info 
 $mtime = (stat($file))[9];

 open INP, "<$file";

# temporary files
 $file_temp = $file.".temp";
 open OUT, ">$file_temp";

 $i = 0;
 $between_beginning_and_end = 0;

 while (<INP>){
   $i = $i + 1;

   $line_previous_10 = $line_previous_9;
   $line_previous_9 = $line_previous_8;
   $line_previous_8 = $line_previous_7;
   $line_previous_7 = $line_previous_6;
   $line_previous_6 = $line_previous_5;
   $line_previous_5 = $line_previous_4;
   $line_previous_4 = $line_previous_3;
   $line_previous_3 = $line_previous_2;
   $line_previous_2 = $line_previous;
   $line_previous = $line_current;
   $line_current = $_;

# detect 'BEGINNING OF'
  if ($between_beginning_and_end == 0){
   if (/BEGINNING OF EQUILIBRATION/ || /Beginning of equilibration/){
    $between_beginning_and_end = 1;
    $line_beginning_index = $i;
   }
  }

# print lines if outside 'BEGINNING and END' or 10 lines after 'BEGINNING OF'
  if ($between_beginning_and_end == 0 or $i < $line_beginning_index + 10){
    print OUT $_;
  }

# print last lines when encountering 'Final write:' for DMC
  if ($between_beginning_and_end == 1){
   if (/Final write:/){
    $between_beginning_and_end = 0;
    print OUT "............ lines cut ............\n";
    print OUT $line_previous_4;  
    print OUT $line_previous_3;  
    print OUT $line_previous_2;  
    print OUT $line_previous;  
    print OUT $line_current;  
   }  
  }

# print last lines when encountering 'END       OF all' for VMC
  if ($between_beginning_and_end == 1){
   if (/END       OF all/ || /End       of accumulation/ || /End       of  accumulation/){
    $between_beginning_and_end = 0;
    print OUT "............ lines cut ............\n";
    print OUT $line_previous_3;  
    print OUT $line_previous_2;  
    print OUT $line_previous;  
    print OUT $line_current;  
   }  
  }
   
  } # end of loop over lines of file
  
# print last lines if exit of file without encountering 'END       OF all' or 'Final write:'
  if ($between_beginning_and_end == 1){
   print OUT "............ lines cut ............\n";
   print OUT $line_previous_10;  
   print OUT $line_previous_9;  
   print OUT $line_previous_8;  
   print OUT $line_previous_7;  
   print OUT $line_previous_6;  
   print OUT $line_previous_5;  
   print OUT $line_previous_4;  
   print OUT $line_previous_3;  
   print OUT $line_previous_2;  
   print OUT $line_previous;  
   print OUT $line_current;  
  }

 close INP;
 close OUT;

# replace original file with new reduced one
 rename $file_temp, $file;

# reset original timestamp
  utime $mtime, $mtime, $file or die "Error setting timestamp for $file: $!\n";

 } # end of loop over files

exit;

