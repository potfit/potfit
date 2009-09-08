#!/usr/bin/perl
#
# 09.05.2003 --- U. Koschella
#
# Creates directory STATIC and a static
# version of the html manual therein
# (i.e. parses the server side includes)
#
# NOTE: This is a very simple parser.
#       Server side include commands may
#       not exceed one line!
#       They may also not be recursive
#       for the moment (i.e. an included
#       file may not contain SSIs).

use File::Copy;

if ( ! -e "STATIC" ) { mkdir "STATIC" ; }
if ( ! -d "STATIC" ) { die "File STATIC exists but is no directory! Exiting"; }

opendir D, "." or die "Cannot open \".\". Exiting";
@all = readdir D;
closedir D;

@html = grep /^[^\.].*\.html/, @all;
@css  = grep /^[^\.].*\.css/, @all;
@gif  = grep /.*\.gif/, @all;
@png  = grep /.*\.png/, @all;
@jpg  = grep /.*\.jpg/, @all;
@jpeg = grep /.*\.jpeg/, @all;

@cp    = ( @css, @gif, @png, @jpg, @jpeg );
@parse = ( @html );

foreach $i ( @cp ) {
  copy $i, "STATIC/$i"  or die "Cannot copy $i to STATIC/. Exiting";
}
`cp -r potentials STATIC`;

foreach $i ( @parse ) {
  open IN,  "<$i"         or die "Cannot open $i for reading. Exiting";
  open OUT, ">STATIC/$i"  or die "Cannot open STATIC/$i for writing. Exiting";
  while (<IN>) {
    if ($_ =~ /<!--#include .*-->/) {
      $before = $after = $command = $_;
      $before  =~ s/<!--#include .*-->.*\n//;
      $after   =~ s/.*<!--#include .*-->//;
      $command =~ s/.*<!--#include (.*)-->.*\n/$1/;
      print OUT $before;
      if ($command =~ /\s*virtual=".*"\s*/) {
        $file = $command;
        $file =~ s/\s*virtual="(.*)"\s*/$1/;
        if (open IN2, "<$file") {
          while (<IN2>) {
            print OUT $_;
          }
          close IN2;
        } else {
          print STDERR "WARNING: Cannot open $file for reading. Ignoring SSI command in $i.\n"; 
        }
      } else {
        print STDERR "WARNING: Ignoring unknown server side include command $command in $i.\n";
      }
      print OUT $after;
    } else {
      print OUT $_;
    }
  }
  close OUT;
  close IN;
}

