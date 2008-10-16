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
use GD ;
use settings ;

my $cgi = new CGI;
my $search = $cgi->param('search') ;
print $cgi->header(-type=>'text/plain',-expires=>'-1s');

sub me_die {
	my ( $out ) = @_ ;
	print "$out\n$search" ;
	exit ;
}

my $dbfile = "$datapath/$annotation_file";
my $dbh = DBI->connect(
	"dbi:SQLite:dbname=$dbfile", # DSN: dbi, driver, database file
	"",                          # no user
	"",                          # no password
	{ RaiseError => 1 },         # complain if something goes wrong
) or me_die "ERROR\n".$DBI::errstr;

me_die "ERROR\nNo search phrase\n" if $search eq '' ;

my $tags = $dbh->selectall_arrayref ( qq { SELECT id,key,value FROM tags WHERE value LIKE \"%$search%\" } ) ;
me_die "ERROR\nNo results\n" if scalar ( @{$tags} ) == 0 ;
me_die "ERROR\nToo many results\n" if scalar ( @{$tags} ) > 50 ;

my %ids ;
foreach ( @{$tags} ) {
	$ids{$_->[0]} = 1 ;
}

my $id_list = join "," , sort keys %ids ;
my $chromosomes = $dbh->selectall_arrayref ( qq { SELECT name FROM chromosomes } ) ;
my @features ;

foreach ( @{$chromosomes} ) {
	my $chr = $_->[0] ;
	my $sql ="SELECT type,start,end,strand,id FROM chr_$chr WHERE id IN ( $id_list )" ;
	my $ref = $dbh->selectall_arrayref ( $sql ) ;
	push @features , $ref ;
}

my %out ;
foreach my $chrid ( 0 .. $#features ) {
	my @a = @{$features[$chrid]} ;
	my $chr = $chromosomes->[$chrid] ;
	foreach ( @a ) {
		my $key = $chr->[0] . "|" . $_->[1] . "|" . $_->[2] . "|" . $_->[0] ;
		next if defined $out{$key} ;
		
		$_->[5] = $chr->[0] ;
		foreach my $tag ( @{$tags} ) {
			next unless $tag->[0] eq $_->[4] ;
			$_->[6] = $tag->[2] ;
			last ;
		}
		
		$out{$key} = $_ ;
		print join ( "\t" , @{$_} ) . "\n" ;
	}
}

me_die "ERROR\nNo results\n" if scalar ( keys %out ) == 0 ;
