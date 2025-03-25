#! /usr/local/bin/perl
# Stripped version of makemake.perl (originally: Michael Wester)
#
# Adapted to write dependencies for Fortran90 which can be included
# in Makefile.am in a libtool based configuration
#
# Thomas Reiter <thomasr@nikhef.nl> 2010
# ---------------------------------------------------------------------------
#
# Usage: makemake {<program name> {<F90 compiler or fc or f77 or cc or c>}}
#
# Generate a Makefile from the sources in the current directory.  The source
# files may be in either C, FORTRAN 77, Fortran 90 or some combination of
# these languages.  If the F90 compiler specified is cray or parasoft, then
# the Makefile generated will conform to the conventions of these compilers.
# To run makemake, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Written by Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
# ---------------------------------------------------------------------------
open(MAKEFILE, "> Makefile.dep");

#
# Source listing
#
@srcs = <*.f90>;

#
# Dependency listings
#
&MakeDependsf90($ARGV[1]);
#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
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
# &MakeDependsf90(); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local(@dependencies);
   local(%filename);
   local(%lfilename);
   local(@incs);
   local(@odeps);
   local(@modules);
   local($objfile);
   #local($lobjfile);
   local($prefix);
   local($mname);
   local($fname);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$//;
      }
   }

   #print MAKEFILE "# Correspondance between .mod files and .o files\n";
   #while ( ($mname, $fname) = each %filename )
   #{
   #   print MAKEFILE "$mname.mod: $fname\n";
   #}

   print MAKEFILE "# Module dependencies\n";

   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90>) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
      }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.f90$//;

         undef @dependencies;
         foreach $module (@modules) {
            if ($filename{$module}) {
               push(@dependencies, $filename{$module});
            }
         }
         @dependencies = &uniq(sort(@dependencies));

         if (scalar(@dependencies) > 0) {
            foreach $ext (".o", ".lo", ".obj") {
               @odeps = map { $_ . $ext } @dependencies;

               print MAKEFILE "$objfile$ext: ";
               &PrintWords(length($objfile) + 2, 0,
                   @odeps, &uniq(sort(@incs)));
               print MAKEFILE "\n";
            }
         }

         undef @incs;
         undef @modules;
      }
      close(FILE)
   }
}
