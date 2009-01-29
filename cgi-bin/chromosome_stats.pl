#!/usr/local/bin/perl -wT
#########
# Author:     Magnus Manske (mm6@sanger.ac.uk)
# Group:      Team 112
#

BEGIN { push @INC, "."; }

use strict ;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use DBI;
use settings ;

my $cgi = new CGI;
my %chromosomes ;

# Read the sequence of a chromosome from a FASTA file
sub get_chromosomes_from_genome_file {
	my $genome_file = shift ;
	$genome_file = "$datapath/$genome_file" ;
	open GENOME , $genome_file or die "NOT FOUND" ;
	my $s = '' ;
	my $lchr = '' ;
	while ( <GENOME> ) { # For each line of input file
		if ( $_ =~ /^>/ ) { # If it starts with an ">", it indicates the end of a chromosome and the beginnig of a new one...
			$chromosomes{$lchr} = $s if $lchr ne '' ;
			$s = "" ;
			chomp ;
			$lchr = substr ( $_ , 1 ) ;
			$lchr =~ s/^\s+// ;
			$lchr =~ s/\s+$// ;
			$lchr =~ /^([\S+]+)/ ;
			$lchr = $1 ;
			$lchr =~ s/[^A-Za-z0-9]/_/g ;
		} else { # Otherwise, DNA
			chomp ;
			$s .= uc $_ ;
		}
	}
	$chromosomes{$lchr} = $s if $lchr ne '' ;
}

get_chromosomes_from_genome_file ( $genome_file ) ;

print $cgi->header(-type=>'text/plain',-expires=>'-1s');
foreach ( sort keys %chromosomes ) {
	print "$_\t" . length ( $chromosomes{$_} ) . "\n" ;
}
