#! /usr/bin/perl
#
# Usage: makemake {<program name> {<F95 compiler or fc or f77 or cc or c>}}
#
# Generate a Makefile from the sources in the current directory.  The source
# files may be in either C, FORTRAN 77, Fortran 90 or some combination of
# these languages.  If the F95 compiler specified is cray or parasoft, then
# the Makefile generated will conform to the conventions of these compilers.
# To run makemake, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Written by Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
use File::Find;
use List::MoreUtils;
#
# HARDCODE FOLDERS TO SCAN FOR SRCS FILES
# (the orders sometimes matters)
#
@scan = ('commons','MED_tools','fit','vmc','dmc','dmc/dmc_elec','fit/MPI','vmc/MPI','dmc/dmc_elec/MPI','dmc/dmc_elec/MPI_global_pop','dmc/dmc_elec/MPI_global_pop_big');
#
# Header
#
print STDOUT "Producing 'Makefile'\n";
open(MAKEFILE, "> Makefile");
print MAKEFILE "# Warning: this 'Makefile' is automatically generated by makemake.pl by typing 'make make'.\n";
print MAKEFILE "# Warning: DO NOT MODIFY THIS 'Makefile' BY HAND!\n";
print MAKEFILE "# Warning: However, you may need to modify the file 'makefile.inc' for machine-specific options.\n\n";
#
# Global variable
#
print MAKEFILE "include ../makefile.inc\n\n";
print MAKEFILE "OBJDIR=compiler_files\n";
print MAKEFILE "MODDIR=compiler_files\n";
print MAKEFILE "INCDIR=.\n\n";
#
# Source folder listing
# (gather all folder containing respectively .f90, .f and .c files)
#
find(\&dir,'.');
sub dir {
  unless (grep{$_ eq $File::Find::dir} @srcdir1) {
    push(@srcdir1,$File::Find::dir) if ($_ =~ /\.f90$/);
  }
  unless (grep{$_ eq $File::Find::dir} @srcdir2) {
    push(@srcdir2,$File::Find::dir) if ($_ =~ /\.f$/)  ;
  }
  unless (grep{$_ eq $File::Find::dir} @srcdir3) {
    push(@srcdir3,$File::Find::dir) if ($_ =~ /\.c$/)  ;
  }
}
foreach (@srcdir1,@srcdir2,@srcdir3){ s/^\.//g };
foreach (@srcdir1,@srcdir2,@srcdir3){ s/^\///g };
foreach $file (@srcdir1,@srcdir2,@srcdir3){ $file =~ s/$/\//g unless $file eq '' }; 
#
# Source listing
#
@srcs=<*.f90 *.c *.f>;
foreach $dir (@scan) {
  push(@srcs,<$dir/*.f>);
  push(@srcs,<$dir/*.f90>);
  push(@srcs,<$dir/*.c>);
}
print MAKEFILE "SRCS =\t";
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
foreach (@objs) { s/^/\$(OBJDIR)\// };
&PrintWords(8, 0, @objs);
print MAKEFILE "\n\n";
#
# Define common macros
#
print MAKEFILE "# PHONY is a Make keyword that tells it to execute the rule regardless of\n";
print MAKEFILE "# the dependencies.  We use it here for 2 reasons:\n";
print MAKEFILE "# 1) prevents it from getting confused if there happens to be a file named\n";
print MAKEFILE "#    'clean' or 'clean_all' in the directory.\n";
print MAKEFILE "# 2) Also, for the libraries, linpack etc. it does the make even though there\n";
print MAKEFILE "#    are directories named linpack etc.\n";
print MAKEFILE ".PHONY: clean_local clean clean_all clean_all_lib make\n\n";
#
# Default
#
print MAKEFILE "\$(EXE): dir \$(LIBS_MAKE) \$(OBJS) revision_and_date.inc\n";
print MAKEFILE "\t\$(LD) \$(LD_FLAGS) -o \$@ \$(OBJS) \$(LIBS) \$(LD_END)\n\n";
#
# Make the "compiler_files" directories
#
print MAKEFILE "dir:\n";
print MAKEFILE "\t@\$(call mk_fortran_dir)\n\n";
print MAKEFILE "define mk_fortran_dir\n";
print MAKEFILE "\tmkdir -p \$(MODDIR)\n";
foreach $dir (sort(@srcdir1,@srcdir2,@srcdir3)) {
  print MAKEFILE "\tmkdir -p \$(OBJDIR)/$dir\n";
};
print MAKEFILE "endef\n\n";

#
# revision_and_date
#
print MAKEFILE "revision_and_date.inc: \$(SRCS)\n";
print MAKEFILE "\tbash ../tools/revision_and_date.sh\n\n";
#
# Makefile
#
print MAKEFILE "make:\n";
print MAKEFILE "\tperl makemake.pl $ARGV[0]\n";
#print MAKEFILE "\t/usr/bin/ctags *.f *.f90\n";
#print MAKEFILE "\tfind . -iregex '.*\\.f\\(90\\)?' | etags -\n";
print MAKEFILE "\n";
#
# Cleans
#
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -fr *exe *.dif *.lst *.sav compiler_files\n";
print MAKEFILE "\n";
print MAKEFILE "clean_all_lib:\n";
print MAKEFILE "\tmake clean_all\n";
print MAKEFILE "\tcd ../lib ; make clean_all\n";
print MAKEFILE "\n";
#
# debug
#
print MAKEFILE "debug:\n";
print MAKEFILE "\tmake \"EXE=\$\{DEBUG_EXE\}\" \"F95_FLAGS=\$\{F95_DEBUG_FLAGS\}\" \"F77_FLAGS=\$\{F77_DEBUG_FLAGS\}\" \"CC_FLAGS=\$\{CC_DEBUG_FLAGS\}\" \"LD_FLAGS=\$\{LD_DEBUG_FLAGS\}\"\n\n";
#
# prof
#
print MAKEFILE "prof:\n";
print MAKEFILE "\tmake \"EXE=\$\{PROF_EXE\}\" \"F95_FLAGS=\$\{F95_PROF_FLAGS\}\" \"F77_FLAGS=\$\{F77_PROF_FLAGS\}\" \"CC_FLAGS=\$\{CC_PROF_FLAGS\}\"  \"LD_FLAGS=\$\{LD_PROF_FLAGS\}\"\n\n";
#
# mpi
#
print MAKEFILE "mpi:\n";
print MAKEFILE "\tmake \"EXE=\$\{MPI_EXE\}\" \"F95=\$\{F95_MPI\}\" \"F77=\$\{F77_MPI\}\" \"CC=\$\{CC_MPI\}\" \"F95_FLAGS=\$\{F95_MPI_FLAGS\}\" \"F77_FLAGS=\$\{F77_MPI_FLAGS\}\" \"CC_FLAGS=\$\{CC_MPI_FLAGS\}\" \"LD=\$\{LD_MPI\}\" \"LD_END=\$\{LD_END_MPI\}\"\n\n";
#
# debug_mpi
#
print MAKEFILE "debug_mpi:\n";
print MAKEFILE "\tmake \"EXE=\$\{MPI_EXE\}\" \"F95=\$\{F95_MPI\}\" \"F77=\$\{F77_MPI\}\" \"CC=\$\{CC_MPI\}\" \"F95_FLAGS=\$\{F95_DEBUG_MPI_FLAGS\}\" \"F77_FLAGS=\$\{F77_DEBUG_MPI_FLAGS\}\" \"CC_FLAGS=\$\{CC_MPI_FLAGS\}\" \"LD=\$\{LD_MPI\}\" \"LD_END=\$\{LD_END_MPI\}\"\n\n";
#
# Libraries
#
#print MAKEFILE "ifdef EINSPLINE\n";
#print MAKEFILE "LIBS = ../lib/lib/libcyrus.a ../lib/lib2/blas/libblas.a ../lib/lib2/lapack/liblapack.a ../lib/lib2/linpack/liblinpack.a ../lib/lib2/einspline/lib/libeinspline.a ../lib/lib2/pspline/pspline/libpspline.a ../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a ../lib/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a\n";
#print MAKEFILE "else\n";
#print MAKEFILE "LIBS = ../lib/lib/libcyrus.a ../lib/lib2/blas/libblas.a ../lib/lib2/lapack/liblapack.a ../lib/lib2/linpack/liblinpack.a ../lib/lib2/pspline/pspline/libpspline.a ../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a ../lib/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a\n";
#print MAKEFILE "endif\n\n";
print MAKEFILE "../lib/lib/libcyrus.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
print MAKEFILE "../lib/lib2/blas/libblas.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
print MAKEFILE "../lib/lib2/lapack/liblapack.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
print MAKEFILE "../lib/lib2/linpack/liblinpack.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
print MAKEFILE "../lib/lib2/einspline/lib/libeinspline.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
#print MAKEFILE "../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a: ../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a(main/anneal.o)\n";
print MAKEFILE "../lib/SimulatedAnnealing/quench_anneal/lib/libquench.a:\n";
print MAKEFILE "\tcd ../lib ; make\n";
print MAKEFILE "\n";
#
# Suffixes
#
print MAKEFILE ".SUFFIXES: \n\n";
print MAKEFILE ".SUFFIXES: .f90 .o .f .c\n\n";
#
# Rules for all objects
#
foreach $dir (@srcdir1) {
  print MAKEFILE "\$(OBJDIR)/%.o:$dir%.f90\n";
  print MAKEFILE "\t\$(F95) \$(F95_FLAGS) \$(INCCMD)\$(INCDIR) \$(MODCMD)\$(MODDIR) -o \$@ -c \$<\n";
}
foreach $dir (@srcdir2) {
  print MAKEFILE "\$(OBJDIR)/$dir%.o:$dir%.f\n";
  print MAKEFILE "\t\$(F77) \$(F77_FLAGS) \$(INCCMD)\$(INCDIR) \$(MODCMD)\$(MODDIR) -o \$@ -c \$<\n";
}
foreach $dir (@srcdir3) {
  print MAKEFILE "\$(OBJDIR)/$dir%.o:$dir%.c\n";
  print MAKEFILE "\t\$(CC) \$(CC_FLAGS) \$(INCCMD)\$(INCDIR) \$(MODCMD)\$(MODDIR) -o \$@ -c \$<\n";
}
print MAKEFILE "\n\n";
#
# Dependency listings
#
&MakeDependsf95($ARGV[1]);
&MakeDepends("*.f", '^\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.c", '^\s*#\s*include\s+["\']([^"\']+)["\']');

#
# SUB used in the scripts ==========================================================
#

#
# &PrintWords(current output column, extra tab?, word list); --- print words nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }
#
# &LanguageCompiler(compiler, sources); --- determine the correct language compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "F77"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "F95";
         }
      }
   else {
      CASE: {
         grep(/\.(f90|f95)$/, @srcs)   && do { $compiler = "F95"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "F77";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }
#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }
#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }
#
# &MakeDepends(language pattern, include file sed pattern); --- dependency maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
      if (@incs) {
         $file =~ s/\.[^.]+$/.o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE "\n";
         undef @incs;
         }
      }
   }
#
# &MakeDependsf95(f95 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf95 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   #foreach $file (<*.f90 commons/*.f90 MED_tools/*.f90 tools/*.f tools/*.f90 tools/*.c initialization/*.f90 dmc/*.f90 >) {
   foreach $file (@srcs) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.(f90|f95)$/.o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   #foreach $file (<*.f90 commons/*.f90 MED_tools/*.f90 tools/*.f tools/*.f90 tools/*.c initialization/*.f90 dmc/*.f90>) {
   foreach $file (@srcs) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*include \'mpif\.h\'/i && pop(@incs);  # JT: remove mpif.h from the list of dependancies
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (@incs || @modules) {
         ($objfile = $file) =~ s/\.(f90|f95)$/.o/;
         #$objfile =~ s/.*\///;
         $objfile =~ s/^/\$(OBJDIR)\//;
         print MAKEFILE "$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         foreach $f (@dependencies) {
           #$f =~ s/.*\///;
           $f =~ s/^/\$(OBJDIR)\//;
         }
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         #
         # Cray F95 compiler
         #
         if ($compiler eq "cray") {
            print MAKEFILE "\t\$(F95) \$(F95_FLAGS) -c ";
            foreach $depend (@dependencies) {
               push(@modules, "-p", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         #
         # ParaSoft F95 compiler
         #
         if ($compiler eq "parasoft") {
            print MAKEFILE "\t\$(F95) \$(F95_FLAGS) -c ";
            foreach $depend (@dependencies) {
               $depend =~ s/\.o$/.(f90|f95)/;
               push(@modules, "-module", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         }
      }
   }
