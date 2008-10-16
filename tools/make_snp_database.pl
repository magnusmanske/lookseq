#!/usr/local/bin/perl

use strict ;
use warnings ;

my %chrs ;

#print "CREATE TABLE meta ( key VARCHAR[64] , value VARCHAR[256] ) ;\n" ;
print "BEGIN EXCLUSIVE TRANSACTION;\n" ;
while ( <> ) {
	chomp ;
	my ( $chr , $pos , $ref , $alt ) = split /\s/ , $_ ;
	next unless defined $alt ;
	
	$chr .= "_snps" ;

	unless ( defined $chrs{$chr} ) {
		$chrs{$chr} = 1 ;
		print "COMMIT;\n" ;
		print "BEGIN EXCLUSIVE TRANSACTION;\n" ;
		print "CREATE TABLE $chr ( pos INTEGER PRIMARY KEY, ref CHAR, alt CHAR );\n" ;
	}
	
	print "INSERT INTO $chr ( pos , ref , alt ) VALUES ( $pos , \"$ref\" , \"$alt\");\n" ;
}
print "COMMIT;\n" ;
