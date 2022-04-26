#!/usr/bin/env perl

# This script in part relies on code from the import scripts of Krona
# https://github.com/marbl/Krona and makes use of a perl module from 
# this project.
# See also https://github.com/marbl/Krona/blob/master/KronaTools/LICENSE.txt


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';
use lib dirname(abs_path($0)).
"/KronaTools/lib";
use KronaTools;


my $out = 'krona.html';
my $noQuantity = 0;
my $color = undef;
my $hueLow = 0;
my $hueHigh = 120;
my $help = 0;

GetOptions(
  'output|o=s' => \$out,
  'no-quantity|q' => \$noQuantity,
  'color|c=s' => \$color,
  'hue-low=i' => \$hueLow,
  'hue-high=i' => \$hueHigh,
  "help|h|?" => \$help,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage("$0: No input files.\n") if !@ARGV;

setOption('out', $out);

my $tree = newTree();
my @datasetNames;
my $set = 0;

foreach my $input(@ARGV) {

  my($fileName, $magFile, $name) = parseDataset($input);

  push @datasetNames, $name;

  open INFILE, "<$fileName" or die $!;

  while ( <INFILE> ) {
    if (/^#/) {
      next;
    }

    chomp;

    my @lineage = split /\t/ ;
    my $quant;
    my $val;

    if ($noQuantity) {
      $quant = 1;
    } else {
      $quant = shift @lineage;
    }

    if (defined($color)) {
      $val = shift @lineage;
      if (!length $val) {
        $val = undef;
      }
    } else {
      $val = undef;
    }

    addByLineage($tree, $set, \@lineage, undef, $quant, $val);
  }

  $set++;
  close INFILE;
}

my @attributeNames =
  (
    'magnitude',
    'magnitudeUnassigned'
  );

my @attributeDisplayNames =
  (
    'Total',
    'Unassigned'
  );

if (defined($color)) {
  push(@attributeNames, 'score');
  push(@attributeDisplayNames, $color);
} else {
	$hueLow = $hueHigh = undef;
}

writeTree
  (
    $tree, 
	\@attributeNames,
	\@attributeDisplayNames,
	\@datasetNames,
    $hueLow,
    $hueHigh
  );

__END__

=head1 NAME

my-prog.pl - Customized mysqldump utility

=head1 SYNOPSIS

        krona_import.pl [OPTIONS] file1 [file2...]

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exits.

=item B<--name>

Sets the Krona chart name

=item B<--output>

Sets the Krona HTML output file

=item B<--no-quantity>

No first numeric column with abundance (quantity, magnitude) information

=item B<--color>

Name of a color variable, whose numeric values are assumed to be in the
second column (or first if --no-quantity).

=item B<--hue-low>

Hue for lowest values

=item B<--hue-high>

Hue for highest values

=back

=cut