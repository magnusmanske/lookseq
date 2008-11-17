#!/usr/bin/perl

use strict ;
use warnings ;

my ( $reference , $map , $fragment , $variance , $from , $to ) = @ARGV ; # $from, $to are optional
my %chrs ;

# Read reference genome
get_reference_from_file ( $reference ) ;

$from = 1 unless defined $from ;
unless ( defined $to ) {
	$to = $from ;
	foreach ( keys %chrs ) {
		$to = length ( $chrs{$_} ) if $to < length ( $chrs{$_} ) ;
	}
}
$from = int ( $from ) ;
$to = int ( $to ) ;

# Create database tables
&create_tables ;

my $seqlen ;
my $lastname = '' ;
#my $lastchr = '' ;
my @lastdata = () ;
#open MAPVIEW , "sort $map |" ;
open MAPVIEW , $map ;
my $cnt = 0 ;
while ( <MAPVIEW> ) {
	chomp ;
	my ( $basename , @b ) = split "\t" , $_ ;
	last unless $basename ;
	
	$basename =~ s|/[12]$|| ;
	
#	last if $cnt++ > 10000 ; # TESTING
	
	my @a ; # Patching for mapview
	$a[13] = $b[3] ; # sequence
	$a[0] = $b[1] ; # chromosome
	$a[1] = int ( $b[2] ) + 1 ; # position; ASSUMES ZERO-BASED POSITION, WHICH IS DEFAULT IN BOWTIE!
	$a[2] = $b[0] ; # strand
	
	
	if ( $seqlen ) {
		die "Length mismatch!\n" if $seqlen != length $a[13] ;
	} else {
		$seqlen = length $a[13] ;
	}
	
	if ( $lastname eq $basename ) {
		add_pair ( $basename , \@lastdata , \@a ) ;
		$lastname = '' ;
	} elsif ( $lastname ) {
		add_single ( $lastname , \@lastdata ) ; # This does not seem to occur. Strange...
		$lastname = $basename ;
		@lastdata = @a ;
	} else {
		$lastname = $basename ;
		@lastdata = @a ;
	}
}
close MAPVIEW ;

&finish_db ;

# Add a read pair
sub add_pair {
	my ( $name , $d1 , $d2 ) = @_ ;
	return if $d1->[0] ne $d2->[0] ; # Not on the same chromosome, can't handle this...
	my $chr = $d1->[0] ;
	my $pos1 = $d1->[1] ;
	my $pos2 = $d2->[1] ;
	return unless defined $chrs{$chr} ;
	return if $pos1 >= length $chrs{$chr} ;
	return if $pos2 >= length $chrs{$chr} ;
	
	return if $pos1 > $to or $pos2 < $from ;
	
	my $seq1 = get_match ( $d1 ) ;
	my $seq2 = get_match ( $d2 ) ;
	my $inv = $d1->[2] eq $d2->[2] ;
	
	if ( $pos1 > $pos2 ) {
		( $pos1 , $pos2 ) = ( $pos2 , $pos1 ) ;
		( $seq1 , $seq2 ) = ( $seq2 , $seq1 ) ;
	}
	
	my $p = sprintf "%s %12d" , $chr , $pos1 ;
	
	if ( $inv ) { # INVERSIONS
		print "/*I$p*/ INSERT INTO $chr"."_inversions (read_name,pos1,seq1,pos2,seq2) VALUES (\"$name\",$pos1,\"$seq1\",$pos2,\"$seq2\") ;\n" ;
	} elsif ( $seq1 or $seq2 ) { # SNP MATCH
		print "/*S$p*/ INSERT INTO $chr"."_snp_match (read_name,pos1,seq1,pos2,seq2) VALUES (\"$name\",$pos1,\"$seq1\",$pos2,\"$seq2\") ;\n" ;
	} else { # PERFECT MATCH
		print "/*P$p*/ INSERT INTO $chr"."_perfect_match (read_name,pos1,pos2) VALUES (\"$name\",$pos1,$pos2) ;\n" ;
	}
}

# Add a single read
sub add_single {
	my ( $name , $d1 ) = @_ ;
	my $chr = $d1->[0] ;
	my $pos = $d1->[1] ;
	return unless defined $chrs{$chr} ;
	return if $pos >= length $chrs{$chr} ;
	return if $pos > $to or $pos < $from ;

	my $seq = get_match ( $d1 ) ;
	my $p = sprintf "%s %12d" , $chr , $pos ;
	print "/*X$p*/ INSERT INTO $chr"."_single (read_name,pos1,seq1) VALUES (\"$name\",$pos,\"$seq\") ;\n" ;
}

sub get_match {
	my ( $data ) = @_ ;
	my $chr = $data->[0] ;
	my $pos = $data->[1] ;
	my $str = $data->[2] ;
	my $seq = uc $data->[13] ;
	my $ref = substr $chrs{$chr} , $pos-1 , length $seq ;
	return '' if $seq eq $ref ;
=cut	
	if ( $str eq '-' ) {
		$seq =~ tr/ACGT/TGCA/ ;
		$seq = reverse $seq ;
	}
=cut
	return uc $ref if length ( $seq ) != length ( $ref ) ;

	my $ret = '' ;
	foreach ( 0 .. $seqlen - 1 ) {
		my $s = substr ( $seq , $_ , 1 ) ;
		if ( $s eq substr ( $ref , $_ , 1 ) ) {
			$ret .= lc $s ;
		} else {
			$ret .= $s ;
		}
	}
	return $ret ;
}


# Read the sequence of a chromosome from a FASTA file
sub get_reference_from_file {
	my $genome_file = shift ;
	open GENOME , $genome_file or return "" ;
	my $s = "" ;
	my $lchr ;
	while ( <GENOME> ) { # For each line of input file
		chomp ;
		if ( $_ =~ /^>/ ) { # If it starts with an ">", it indicates the end of a chromosome and the beginnig of a new one...
			$chrs{$lchr} = $s if $lchr ;
			$s = "" ;
			$lchr = substr ( $_ , 1 ) ;
			$lchr =~ s/^\s+// ;
			$lchr =~ s/\s+$// ;
		} else { # Otherwise, DNA
			$s .= uc $_ ;
		}
	}
	$chrs{$lchr} = $s if $lchr ;
}

sub create_tables {
	my $dbinit = '/*000*/' ;
	print "/*0*/ BEGIN EXCLUSIVE TRANSACTION;\n" ;
	print "$dbinit CREATE TABLE meta ( key VARCHAR[64] , value VARCHAR[256] );\n" ;
	print "$dbinit CREATE TABLE chromosomes ( name VARCHAR[256] , size INTEGER );\n" ;
	foreach my $chr ( sort keys %chrs ) {
		print "$dbinit CREATE TABLE $chr" . "_perfect_match ( read_name VARCHAR[32], pos1 INTEGER, pos2 INTEGER );\n" ;
		print "$dbinit CREATE TABLE $chr" . "_inversions ( read_name VARCHAR[32], pos1 INTEGER, seq1 VARCHAR[64] , pos2 INTEGER , seq2 VARCHAR[64] );\n" ;
		print "$dbinit CREATE TABLE $chr" . "_single ( read_name VARCHAR[32], pos1 INTEGER, seq1 VARCHAR[64] );\n" ;
		print "$dbinit CREATE TABLE $chr" . "_snp_match ( read_name VARCHAR[32], pos1 INTEGER, seq1 VARCHAR[64] , pos2 INTEGER , seq2 VARCHAR[64] );\n" ;
	}
}

sub finish_db {
	my $dbend = '/*ZZZ*/' ;
	my $dbend2 = '/*ZZZZ*/' ;
	print "$dbend INSERT INTO meta ( key , value ) VALUES ( \"fragment\",\"$fragment\" ) ;\n" ;
	print "$dbend INSERT INTO meta ( key , value ) VALUES ( \"variance\",\"$variance\" ) ;\n" ;
	print "$dbend INSERT INTO meta ( key , value ) VALUES ( \"read_length\",\"$seqlen\" ) ;\n" ;
	print "$dbend INSERT INTO meta ( key , value ) VALUES ( \"name_prefix\",\"\" ) ;\n" ;
	print "$dbend INSERT INTO meta ( key , value ) VALUES ( \"dbversion\",\"2\" ) ;\n" ;
	
	print "/*ZZ*/ COMMIT;\n" ;
	print "$dbend BEGIN EXCLUSIVE TRANSACTION;\n" ;
	foreach my $chr ( sort keys %chrs ) {
		print "$dbend2 INSERT INTO chromosomes ( name , size ) VALUES ( \"$chr\"  , " . length ( $chrs{$chr} ) . " ) ;\n" ;
		print "$dbend2 CREATE INDEX $chr"."_sin_index ON $chr"."_single ( pos1 );\n" ;
		print "$dbend2 CREATE INDEX $chr"."_snp_index ON $chr"."_snp_match ( pos1 , pos2 );\n" ;
		print "$dbend2 CREATE INDEX $chr"."_inv_index ON $chr"."_inversions ( pos1 , pos2 );\n" ;
		print "$dbend2 CREATE INDEX $chr"."_perfect_index ON $chr"."_perfect_match ( pos1 , pos2 );\n" ;

	}
	print "/*ZZZZZ*/ COMMIT;\n" ;
}
