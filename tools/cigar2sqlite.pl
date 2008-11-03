#!/usr/bin/perl

use strict ;
use warnings ;

my %chromosomes ;
my %data ;

while ( <> ) {
	chomp ;
	$_ =~ s/^cigar:\s*\S+\s*//i ;
	#next unless $_ =~ /^(\S+)\s(\S+)\s(\d+)\s(\d+)\s([+-])\s(\S+)\s(\d+)\s(\d+)\s([+-])\s(\d+)\s(.+)$/ ;
	next unless $_ =~ /^(\S+)\s(\d+)\s(\d+)\s([+-])\s(\S+)\s(\d+)\s(\d+)\s([+-])\s(\d+)\s(.+)$/ ;
	#my ( $strain , $qid , $qstart , $qend , $qstrand , $tid , $tstart , $tend , $tstrand , $score , $rest ) = ( $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 ) ;
	my ( $qid , $qstart , $qend , $qstrand , $tid , $tstart , $tend , $tstrand , $score , $rest ) = ( $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 ) ;
	$rest =~ s/\s*([MDIN]) /|$1/g ;
	$rest = substr $rest , 1 ; # Remove leading "|"
	
	$chromosomes{$tid} = 1 ;
	
	my $subid = 1 ;
	my $key ;
	do {
		$key = sprintf "%s|%12d|%d" , $tid , $tstart , $subid++ ;
	} while ( defined $data{$key} ) ;
	my $value = qq { INSERT INTO cigar ( name , chromosome , start , stop , orientation , data ) VALUES ( "$qid" , "$tid" , $tstart , $tend , "$qstrand" , "$rest" ) } ;
	$data{$key} = $value ;
}

#print scalar ( keys %data ) . "\n" ;

print "BEGIN EXCLUSIVE TRANSACTION;\n" ;
print "CREATE TABLE cigar ( name VARCHAR[128] , chromosome VARCHAR[32] , start INTEGER , stop INTEGER , orientation CHAR , data TINYTEXT ) ;\n" ;
print "CREATE TABLE chromosomes ( name VARCHAR[256] , size INTEGER );\n" ;
print "CREATE TABLE meta ( key VARCHAR[64] , value VARCHAR[256] );\n" ;

foreach ( sort keys %chromosomes ) {
	print "INSERT INTO chromosomes ( name , size ) VALUES ( \"$_\" , 0 ) ;\n" ; # Dummy
}

foreach ( sort keys %data ) {
	print $data{$_} . ";\n" ;
}

print "CREATE INDEX cigar_idx ON cigar ( chromosome , start , stop ) ;\n" ;
print "COMMIT ;\n" ;
