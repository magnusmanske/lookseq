#!/usr/local/bin/perl

use strict ;
use warnings ;

my ( $features , $tags ) = @ARGV ;
die "Usage: make_feature_list_from_embl.pl FEATURES_FILE TAGS_FILE\n" unless defined $tags ;

my $dummy ;
my @out ;
my %chromosomes ;

open FEATURES , $features ;
$dummy = <FEATURES> ;
while ( <FEATURES> ) {
	chomp ;
	next unless $_ = /^(\d+)\s(\S+)\s(\S+)\s(\d+)\s(\d+)\s([+-])$/ ;
	$chromosomes{$2} = 1 ;
	push @out , "INSERT INTO chr_$2 ( type , start , end , strand , id ) VALUES ( \"$3\" , $4 , $5 , \"$6\" , $1 ) ;\n" ;
}
close FEATURES ;

open TAGS , $tags ;
$dummy = <TAGS> ;
while ( <TAGS> ) {
	chomp ;
	next unless $_ =~ /^(\d+)\s+([a-z_]+)\s+(\S+)$/ ;
	push @out , "INSERT INTO tags ( id , key , value ) VALUES ( $1 , \"$2\" , \"$3\") ;\n" ;
}
close TAGS ;

print "BEGIN EXCLUSIVE TRANSACTION;\n" ;
print "CREATE TABLE chromosomes ( name VARCHAR[64] ) ;\n" ;
foreach ( sort keys %chromosomes ) {
	print "CREATE TABLE chr_$_ ( type VARCHAR[64] , start INTEGER , end INTEGER , strand CHAR , id INTEGER ) ;\n" ;
	print "INSERT INTO chromosomes ( name ) VALUES ( \"$_\" ) ;\n" ;
}
print "CREATE TABLE tags ( id INTEGER , key VARCHAR[64] , value VARCHAR[256] ) ;\n" ;

foreach ( @out ) {
	print $_ ;
}
foreach ( sort keys %chromosomes ) {
	print "CREATE INDEX idx_$_ ON chr_$_ ( start , end ) ;\n" ;
}
print "CREATE INDEX idx_tags ON tags ( id ) ;\n" ;
print "COMMIT;\n" ;
