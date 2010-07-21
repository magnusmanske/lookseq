#!/usr/bin/perl -w
use strict;

#use Bio::Tools::GFF;

my $file = '/nfs/team112/annotation/plasmodium/falciparum/Pfalciparum_PlasmoDB-5.5.gff' ;

#=cut

my @vt = ( 'gene' , 'exon' , 'CDS' , 'mRNA' , 'rRNA' ) ;
my %valid_type ;
$valid_type{$_} = 1 foreach ( @vt ) ;


my @out ;
my %chrs ;
my $id = 0 ;
open GFF, $file or die "Can't open GFF file $file\n" ;
while(<GFF>){
	next if $_ =~ m/^#/ ;
	last if $_ =~ m/^>/ ;
	chomp ;
	my ($seqid, undef, $type, $start, $end, undef, $strand, undef, $attrs) = split;
	next unless defined $valid_type{$type} ;
	my $chr = $seqid ;
	$chr =~ s/^.*\|// ;
	next unless $chr =~ m/^MAL/ or $chr =~ m/^MITO/ or $chr =~ m/^PLAST/ ;
	$chrs{$chr} = 1 ;
	
	my @d = split '\|' , $attrs ;
	@d = split ';' , ( $d[1] || '' ) ;
	if ( 1 < scalar @d ) {
		$id++ ;
		push @out , "INSERT INTO chr_$chr (type,start,end,strand,id) VALUES (\"$type\",$start,$end,\"$strand\",$id);" ;
		$d[0] = "ID=" . $d[0] ;
		foreach ( @d ) {
			next unless $_ =~ m/^(\w+)=(.+)$/ ;
			my ( $k , $v ) = ( $1 , $2 ) ;
			$v =~ tr/+/ / ;
			push @out , "INSERT INTO tags (id,key,value) VALUES ($id,\"$k\",\"$v\");" ;
		}
	}
	
#	push @{$gff{$chr}}, [$start, $end];
}
close GFF ;

print "CREATE TABLE chr_$_ ( type VARCHAR[64] , start INTEGER , end INTEGER , strand CHAR , id INTEGER );\n" foreach ( keys %chrs ) ;
print "CREATE TABLE chromosomes ( name VARCHAR[64] );\n" ;
print "CREATE TABLE tags ( id INTEGER , key VARCHAR[64] , value VARCHAR[256] );\n" ;
print "BEGIN EXCLUSIVE;\n" ;
print "INSERT INTO chromosomes ( name ) VALUES ( \"$_\");\n" foreach ( keys %chrs ) ;
print join ( "\n" , @out ) ;
print "\nCOMMIT;\n" ;
print "CREATE INDEX idx_$_ ON chr_$_ ( start , end );\n" foreach ( keys %chrs ) ;
print "CREATE INDEX idx_tags ON tags ( id );\n" ;

#print join "\n" , keys %gff ;

=cut

my $x = read_exons_from_gff ( $file ) ;
print join "\n" , keys %{$x} ;

sub read_exons_from_gff {
	my ( $file ) = @_ ;
	my %ret ;
	open FILE , $file ;
	my $gffio = Bio::Tools::GFF->new(-fh => \*FILE, -gff_version => 3);
	while(my $feat = $gffio->next_feature()) {
		my $loc = $feat->seq_id() ;
		$loc =~ s/^.*(MAL\d+).*$/$1/ ;
		next unless 'exon' eq lc $feat->primary_tag ;
		push @{$ret{$loc}} , [ $feat->start , $feat->end ] ;
	}
	$gffio->close();
	return \%ret ;
}
=cut