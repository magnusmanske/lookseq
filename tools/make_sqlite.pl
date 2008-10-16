#!/usr/bin/perl

my ( $lane , $pair , $fragment , $variance , $group ) = @ARGV ;

`./findknownsnps --noshortcuts --mspi=10 --genome=../3D7_pm.fa --fastq=/nfs/malaria/data2/solexa/lane_data/$lane/$lane.fastq --pair=$pair --fragment=$fragment --variance=$variance --snps=MERGED_CAPILLARY_MAQ.tab.subset --foum --sqlite=$lane.sqlite.text` ;
`sqlite3 $lane.sqlite < $lane.sqlite.text` ;

if ( defined $group ) {
	`sqlite3 $lane.sqlite "INSERT INTO meta ( key , value ) VALUES ( 'group' , '$group' );"` ;
}

`rm $lane.sqlite.text` ;
`mv $lane.sqlite /nfs/WWWdev/INTWEB_docs/data/teams/team112/solexa/` ;
