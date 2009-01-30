#!/usr/local/bin/perl -T
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
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval nanosleep );
use settings ;

my $cgi = new CGI;
my $time0 = [gettimeofday()];
my $debug_output = 0 ;
$debug_output = $cgi->param('debug') if defined $cgi->param('debug') ;
if ( $debug_output ) {
	use warnings ;
}

# Parameters
my $output = $cgi->param('output') || 'text' ;
my $view = $cgi->param('view') || '' ;
my $database = $cgi->param('lane') ;
my $from = $cgi->param('from') ;
my $to = $cgi->param('to') ;
my $chromosome = $cgi->param('chr') ;
my $width = $cgi->param('width') || 1024 ;
my $height = $cgi->param('height') || 512 ;
my $display = $cgi->param('display') || '|perfect|snps|inversions|' ;


# Sanger cache - others, ignore
my ( $sanger_web , $sanger_cache_db , $sanger_cache_hours , $sanger_cache_key ) ;
if ( $use_sanger_cache and ( $to - $from > 500000 ) ) { # Only for large displays
	use SangerPaths qw(core);
	use SangerWeb;
	$sanger_web = SangerWeb->new();
	$sanger_cache_db = $sanger_web->dbstore();
	$sanger_cache_hours = 24 * 10 ;
	$sanger_cache_key = "LookSeq/" . CGI::url(-query=>1) ;

	my $data = $sanger_cache_db->get ( $sanger_cache_key ) ;
	
	if ( 0 < length $data ) {
		if ( $output eq 'text' ) {
			print $sanger_web->cgi()->header(-type=>'text/plain',-expires=>'-1s');
		} else {
			print $sanger_web->cgi()->header(-type=>'image/png',-expires=>'-1s');
			binmode STDOUT;
		}
		print $data ;
		exit ;
	}

}

my $display_perfect = $display =~ m/\|perfect\|/ ;
my $display_snps = $display =~ m/\|snps\|/ ;
my $display_inversions = $display =~ m/\|inversions\|/ ;
my $display_inversions_ext = $display =~ m/\|inversions_ext\|/ ;
my $display_single = $display =~ m/\|single\|/ ;
my $display_pair_links = $display =~ m/\|pairlinks\|/ ;
my $display_pot_snps = $display =~ m/\|potsnps\|/ ;
my $display_noscale = $display =~ m/\|noscale\|/ ;
my $display_uniqueness = $display =~ m/\|uniqueness\|/ ;
my $display_inline_annotation = 1 ; # FIXME

$display_inversions_ext = 0 unless $display_inversions ;

# Parameter paranoia
$from =~ s/\D//g ;
$to =~ s/\D//g ;
$chromosome =~ s/[^A-Za-z0-9\-_\.]//g ;
my $orig_chr = $chromosome ;
$chromosome = "chr_$chromosome" if $chromosome =~ /^\d/ ;

unless ( $view eq 'annotation' or $view eq 'gc' ) {
	die "ERROR\nNot a valid lane\n" if $database eq '' ;
}
die "ERROR\nTO not larger than FROM\n" if $from >= $to ;

my @all_smr ;
my @all_pmr ;
my @all_sin ;
my @all_inv ;
my @all_meta ;
my @all_cig ;
my @all_ann ;
my @snpcache ;
my @dotcache ;
my @anno_sum ;
my @anno_count ;
my $total_reads = 0 ;
my $show_chars ; # For pileup; global var easier that passing it through lots'o'methods
my $do_coverage = $view eq 'coverage' ? 1 : 0 ;
my @coverage ;
my $refseq ;
my %meta ;
my @known_snps ;
my $known_snp_color ;
my $number_of_databases ;
my %cur_meta ;
my @inversion_left ;
my @inversion_right ;
my @inversion_middle ;
my $avg_frag_size ;
my $scale_height = $display_noscale ? 0 : 20 ;

my %iupac = (
	'A'=>'A',
	'C'=>'C',
	'G'=>'G',
	'T'=>'T',
	'R'=>'AG',
	'Y'=>'CT',
	'M'=>'AC',
	'K'=>'GT',
	'W'=>'AT',
	'S'=>'CG',
	'B'=>'CGT',
	'D'=>'AGT',
	'H'=>'ACT',
	'V'=>'ACG',
	'N'=>'ACGT'
) ;
my %iupac_reverse ;
foreach ( keys %iupac ) {
	$iupac_reverse{$iupac{$_}} = $_ ;
}

if ( $display_pot_snps and defined  $snp_file ) {
	my $dbfile = "$datapath/$snp_file";
	my $dbh = DBI->connect(
	    "dbi:SQLite:dbname=$dbfile", # DSN: dbi, driver, database file
	    "",                          # no user
	    "",                          # no password
	    { RaiseError => 1 },         # complain if something goes wrong
	) or die "ERROR\n".$DBI::errstr;

	my $table = $chromosome . '_snps' ;
	my $snps = [] ;
	my $where ;
	$where = "WHERE pos BETWEEN $from AND $to" if $to - $from < 500000 ;
	eval { $snps = $dbh->selectall_arrayref ( "SELECT pos,ref,alt FROM $table $where" ) } ;
	foreach ( @{$snps} ) {
		next if $_->[0] < $from ;
		$known_snps[$_->[0]-$from] = $_->[1] . $_->[2] ;
	}
}

my @databases = split ( ',' , $database ) ;
foreach ( @databases ) {
	my $dbh ;

	if ( $use_mysql ) {
		next unless in_array ( \@mysql_dbs , $_ ) ;
		$dbh = DBI->connect("DBI:mysql:database=$_;host=$mysql_server;port=$mysql_port",
                   $mysql_user, $mysql_password,
                   {'RaiseError' => 1});

	} else {
		$_ =~ s/[^a-zA-Z0-9_\-\.]//g ;
		next if "$_.sqlite" eq $snp_file ;
		next if "$_.sqlite" eq $annotation_file ;
	
		my $dbfile = "$datapath/$_.sqlite";
		next unless -e $dbfile ;
		$dbh = DBI->connect(
		    "dbi:SQLite:dbname=$dbfile", # DSN: dbi, driver, database file
		    "",                          # no user
		    "",                          # no password
		    { RaiseError => 1 },         # complain if something goes wrong
		) or die "ERROR\n".$DBI::errstr;
#		print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output
#		print $dbfile ;
	}

	# Check chromosomes
	my $chrsr = $dbh->selectall_arrayref ( "SELECT name FROM chromosomes" ) ;
	my $chr_is_valid = 0 ;

	foreach ( @{$chrsr} ) {
		$chr_is_valid = 1 if $_->[0] eq $chromosome ;
	}

	next unless $chr_is_valid ;
	
	$number_of_databases++ ;

	# Read metadata
	my $keyname = $use_mysql ? 'the_key' : 'key' ;
	my $metar = $dbh->selectall_arrayref ( "SELECT $keyname,value FROM meta" ) ;
	my %meta = {} ;
	foreach ( @{$metar} ) {
		$meta{$_->[0]} = $_->[1] ;
	}
	$meta{'dbversion'} = 1 unless defined $meta{'dbversion'} ;

	my $from2 = $from - $meta{'read_length'} ;
	my $to2 = $to + $meta{'read_length'} ;
	$from2 = 1 if $from2 < 1 ;

	my $cpm = $chromosome . '_perfect_match' ;
	my $csm = $chromosome . '_snp_match' ;

	my ( $smr , $pmr , $sin , $inv , $cig , $ann , $rl_ref ) ;
	
	# Turn off index search over 25K of sequence
	my ( $pos1 , $pos2 , $irl ) ;
	if ( $use_mysql ) {
		$pos1 = 'pos1' ;
		$pos2 = 'pos2' ;
		$irl = ",read_length_index" ;
		my $dummy = $dbh->selectall_arrayref ( "SELECT id,len1,len2 FROM read_length" ) ;
		foreach ( @{$dummy} ) {
			my $k = $_->[0] ;
			$rl_ref->{$k}->{'len1'} = $_->[1] ;
			$rl_ref->{$k}->{'len2'} = $_->[2] ;
		}
	} else {
		$pos1 = ( $to - $from ) > $use_db_index_cutoff ? '+pos1' : 'pos1' ;
		$pos2 = ( $to - $from ) > $use_db_index_cutoff ? 'pos2' : 'pos2' ;
		$irl = '' ;
	}
	
	# Get SNP matches
	$smr = [] ;
	eval {
		$smr = $display_snps ? $dbh->selectall_arrayref ( "SELECT read_name,pos1,pos2,seq1,seq2 $irl FROM $csm WHERE ( $pos1 BETWEEN $from2 AND $to2 ) OR ( $pos2 BETWEEN $from2 AND $to2 )" ) : [] ;
	} ;
	
	# Get perfect matches
	$pmr = [] ;
	eval {
		$pmr = $display_perfect ? $dbh->selectall_arrayref ( "SELECT read_name,pos1,pos2 $irl FROM $cpm WHERE ( $pos1 BETWEEN $from2 AND $to2 ) OR ( $pos2 BETWEEN $from2 AND $to2 )" ) : [] ;
	} ;
	
	# Get CIGAR data, if any
	$cig = [] ;
	eval { 
		$cig = $dbh->selectall_arrayref ( "SELECT name,start,stop,orientation,data $irl FROM cigar WHERE chromosome=\"$chromosome\" AND start <= $to AND stop >= $from" ) || [] ;
	} ;
	
	if ( $meta{'dbversion'} > 1 ) {
		my $sin_table = $chromosome . '_single' ;
		my $inv_table = $chromosome . '_inversions' ;
		$sin = $display_single ? $dbh->selectall_arrayref ( "SELECT read_name,pos1,seq1 $irl FROM $sin_table WHERE $pos1 BETWEEN $from2 AND $to2" ) : [] ;
		$inv = [] ;
		if ( $meta{'dbversion'} > 2 ) {
			eval {
				$inv = $display_inversions ? $dbh->selectall_arrayref ( "SELECT read_name,pos1,pos2,seq1,seq2,dir $irl FROM $inv_table WHERE ( $pos1 BETWEEN $from2 AND $to2 ) OR ( $pos2 BETWEEN $from2 AND $to2 )" ) : [] ;
			} ;
		} else {
			eval {
				$inv = $display_inversions ? $dbh->selectall_arrayref ( "SELECT read_name,pos1,pos2,seq1,seq2 $irl FROM $inv_table WHERE ( $pos1 BETWEEN $from2 AND $to2 ) OR ( $pos2 BETWEEN $from2 AND $to2 )" ) : [] ;
			}
		}
	} else {
		$sin = [] ;
		$inv = [] ;
	}
	
	$ann = [] ;
	if ( $meta{'dbversion'} > 3 ) {
		eval { 
			$ann = $dbh->selectall_arrayref ( "SELECT text,pos,height,length,colorcode FROM ${chromosome}_annotation WHERE pos <= $to AND pos+length >= $from" ) || [] ;
		} ;
	}
	
	$total_reads += scalar ( @{$smr} ) ;
	$total_reads += scalar ( @{$pmr} ) ;
	$total_reads += scalar ( @{$sin} ) ;
	$total_reads += scalar ( @{$inv} ) ;
	$total_reads += scalar ( @{$cig} ) ;
	
	if ( $use_mysql ) {
		fix_index_read_length_mysql ( $smr , $rl_ref , 1 ) ;
		fix_index_read_length_mysql ( $pmr , $rl_ref , 1 ) ;
		fix_index_read_length_mysql ( $inv , $rl_ref , 1 ) ;
		fix_index_read_length_mysql ( $sin , $rl_ref , 0 ) ;
	} else {
		fix_index_read_length_sqlite ( $smr , $meta{'read_length'} , 1 ) ;
		fix_index_read_length_sqlite ( $pmr , $meta{'read_length'} , 1 ) ;
		fix_index_read_length_sqlite ( $inv , $meta{'read_length'} , 1 ) ;
		fix_index_read_length_sqlite ( $sin , $meta{'read_length'} , 0 ) ;
	}
	
	push @all_smr , $smr ;
	push @all_pmr , $pmr ;
	push @all_sin , $sin ;
	push @all_inv , $inv ;
	push @all_cig , $cig ;
	push @all_ann , $ann ;
	push @all_meta , \%meta ;
}


sub fix_index_read_length_mysql {
	my ( $ref , $rl_ref , $pair ) = @_ ;
	foreach ( 0 .. scalar(@{$ref})-1 ) {
		my $k = pop @{$ref->[$_]} ;
		my $rl1 = $rl_ref->{$k}->{'len1'} ;
		my $rl2 = $rl_ref->{$k}->{'len2'} ;
		push @{$ref->[$_]} , $rl1 ;
		push @{$ref->[$_]} , $rl2 if $pair ;
	}
}

sub fix_index_read_length_sqlite {
	my ( $ref , $rl , $pair ) = @_ ;
	foreach ( 0 .. scalar(@{$ref})-1 ) {
		push @{$ref->[$_]} , $rl ;
		push @{$ref->[$_]} , $rl if $pair ;
	}
}


#_________________________________________


# Read the sequence of a chromosome from a FASTA file
sub get_chromosome_from_genome_file {
  my $genome_file = shift ;
  $genome_file = "$datapath/$genome_file" ;
  my $chr = shift ;
  open GENOME , $genome_file or return "" ;
  my $s = "" ;
  my $lchr = " " ;
  while ( <GENOME> ) { # For each line of input file
	if ( $_ =~ /^>/ ) { # If it starts with an ">", it indicates the end of a chromosome and the beginnig of a new one...
		return $s if $lchr eq $chr ;

		$s = "" ;
		chomp ;
		$lchr = substr ( $_ , 1 ) ;
		$lchr =~ s/^\s+// ;
		$lchr =~ s/\s+$// ;
		$lchr =~ /^([\S+]+)/ ;
		$lchr = $1 ;
		$lchr =~ s/[^A-Za-z0-9]/_/g ;
		$chr = $lchr if $chr eq "" ;

	} else { # Otherwise, DNA
		chomp ;
		$s .= uc $_ ;
	}
  }
  
  return $s if $lchr eq $chr ;
  return "" ;
}

# Cut out the interesting part
sub get_chromosome_part {
	my ( $file , $chromosome , $from , $to ) = @_ ;
	my $s = get_chromosome_from_genome_file ( $file , $chromosome ) ;
	return $s if $s eq '' ;
	$s = substr $s , $from - 1 , $to - $from + 1 ;
	return $s ;
}


#_________________________________________

sub dump_text {
	print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output

	my %meta = %{$all_meta[0]} ;
	my $smr = @all_smr[0] ;
	my $pmr = @all_pmr[0] ;

	# Print meta data
	print scalar ( keys %meta ) . "\n" ; # Number of meta key/value pairs
	foreach ( keys %meta ) {
		print "$_\t" . $meta{$_} . "\n" ;
	}

	print scalar ( @{$smr} ) . "\n" ; # Number of perfectly matching reads
	foreach ( @{$smr} ) {
		print join "\t" , @{$_} ;
		print "\n" ;
	}

	print scalar ( @{$pmr} ) . "\n" ; # Number of perfectly matching reads
	foreach ( @{$pmr} ) {
		print join "\t" , @{$_} ;
		print "\n" ;
	}

	# Reference sequence
	print get_chromosome_part ( $genome_file , $chromosome , $from , $to ) . "\n" ;
}

sub pile2bands_cigar {
	my ( $data , $bands , $btype , $mode ) = @_ ;
	foreach ( @{$data} ) { # name,start,stop,orientation,data
		add_pileup_cigar ( $_->[1] , $_->[4] , $bands , $btype , $mode ) ;
	}
}

sub add_pileup_cigar {
	my ( $start , $cigar , $bands , $btype , $mode ) = @_ ;
	my $seq = decompile_cigar ( $cigar ) ;
	$seq =~ s/I//g ; # Can't show insertions in pileup - yet!
	if ( $start < $from ) {
		$seq = substr $seq , $from - $start ;
		$start = $from ;
	}
	my $out = '' ;
	my $pos = $start - $from ;
	foreach ( 0 .. length ( $seq ) - 1 ) {
		my $ch = substr $seq , $_ , 1 ;
		if ( $ch eq 'D' ) {
			$out .= 'D' if $pos >= 0 ;
		} elsif ( $ch eq 'M' ) {
			$out .= lc ( substr $refseq , $pos , 1 ) if $pos >= 0 ;
		} else {
			$out .= '?' ; # This should never happen
		}
		$pos++ ;
	}
	
	add_pileup ( $start , $out , length $out , $bands , $btype , $mode , 0 ) ;
}

sub dump_image_pileupview {
	my ( $paired_pileup ) = @_ ;
	my @bands ;
	my @btype ;
	my $ft = $to - $from + 1 ;

	$refseq = get_chromosome_part ( $genome_file , $orig_chr , $from , $to ) ;
	
#	print $cgi->header(-type=>'text/plain',-expires=>'-1s') ;
	foreach my $current_db ( 0 .. $#all_meta ) {
		my %meta = %{$all_meta[0]} ;
#		my $rl = $meta{'read_length'} ;
		$show_chars = $width / $ft > 4 ? 1 : 0 ;
		my $sin = $all_sin[$current_db] ;

		foreach ( 0 .. @$sin ) {
			next unless defined $sin->[$_] ;
			my $rl = pop @{$sin->[$_]} ;
			$sin->[$_]->[3] = $sin->[$_]->[2] ;
			$sin->[$_]->[2] = -1 ; 
			push @{$sin->[$_]} , undef ;
			push @{$sin->[$_]} , $rl ;
			push @{$sin->[$_]} , $rl ;
#			print join ", " , @{$sin->[$_]} ;
#			print "\n" ;
		}

		my $pmr = $all_pmr[$current_db] ;
		foreach ( 0 .. @$pmr ) {
			$pmr->[$_]->[5] = $pmr->[$_]->[3] ;
			$pmr->[$_]->[6] = $pmr->[$_]->[4] ;
			$pmr->[$_]->[3] = undef ;
			$pmr->[$_]->[4] = undef ;
		}
	}
#exit;

#	print "___1\n" ;
	foreach my $current_db ( 0 .. $#all_meta ) {
		%cur_meta = %{$all_meta[$current_db]} ;
		pile2bands_cigar ( $all_cig[$current_db] , \@bands , \@btype , 3 ) ;
	}
#	print "___2\n" ;
	foreach my $current_db ( 0 .. $#all_meta ) {
		%cur_meta = %{$all_meta[$current_db]} ;
		pile2bands ( $all_smr[$current_db] , \@bands , \@btype , $all_meta[$current_db]->{'read_length'} , 0 , $paired_pileup , $all_meta[$current_db]->{'fragment'} ) ;
	}
#	print "___3\n" ;
	foreach my $current_db ( 0 .. $#all_meta ) {
		%cur_meta = %{$all_meta[$current_db]} ;
		pile2bands ( $all_pmr[$current_db] , \@bands , \@btype , $all_meta[$current_db]->{'read_length'} , 0 , $paired_pileup , $all_meta[$current_db]->{'fragment'} ) ;
	}
#	print "___4\n" ;
	foreach my $current_db ( 0 .. $#all_meta ) {
		%cur_meta = %{$all_meta[$current_db]} ;
		pile2bands ( $all_inv[$current_db] , \@bands , \@btype , $all_meta[$current_db]->{'read_length'} , 2 , $paired_pileup , $all_meta[$current_db]->{'fragment'} ) ;
	}
#	print "___5\n" ;
	foreach my $current_db ( 0 .. $#all_meta ) {
		%cur_meta = %{$all_meta[$current_db]} ;
		pile2bands ( $all_sin[$current_db] , \@bands , \@btype , $all_meta[$current_db]->{'read_length'} , 1 , 0 , $all_meta[$current_db]->{'fragment'} ) ;
	}
	
#	exit ;
	
	if ( $output eq 'text' ) {
		print $cgi->header(-type=>'text/plain',-expires=>'-1s');
		print "$chromosome : $from - $to\n$refseq\n" ;
		foreach ( @bands ) {
			print "$_\n" ;
		}
		return ;
	}

	my $charheight = 10 ;
	
	if ( $show_chars ) {
		$height = ( $#bands + 3 ) * $charheight ;
	} else {
		$height = $#bands + 5 ;
	}

	my $im = new GD::Image ( $width , $height + 20 ) ;
	my $white = $im->colorAllocate ( 255 , 255 , 255 ) ;
	my $black = $im->colorAllocate ( 0 , 0 , 0 ) ;
	my $blue = $im->colorAllocate ( 0 , 0 , 255 ) ;
	my $single_color = $im->colorAllocate ( 0 , 128 , 0 ) ;
	my $inversion_color = $im->colorAllocate ( 210 , 105 , 30 ) ;
	my $cigar_color = $im->colorAllocate ( 0x44 , 0xB4 , 0xD5 ) ;
	my $red = $im->colorAllocate ( 255 , 0 , 0 ) ;
	my $ltgrey = $im->colorAllocate ( 220 , 220 , 220 ) ;
	
	my $nob = $#bands ;
	
	if ( $show_chars ) {
		foreach my $pos ( 0 .. length ( $refseq ) - 1 ) {
			my $x1 = int ( $width * $pos / $ft ) ;
			$im->string ( gdSmallFont , $x1 , int ( $height - $charheight - 1 ) , substr ( $refseq , $pos , 1 ) , $black ) ;
		}
		foreach my $band ( 0 .. $nob ) {
			$bands[$band] =~ s/\s+$// ; # Trim trailing space
			my $col ;
			my $lb = length ( $bands[$band] ) - 1 ;
			foreach my $pos ( 0 .. $lb ) {
				my $ch = substr ( $bands[$band] , $pos , 1 ) ;
				next if $ch eq ' ' ;
				my $type = ord ( substr ( $btype[$band] , $pos , 1 ) ) - 65 ;
				my $dist_marker = $type >= 100 ;
				$type -= 100 if $dist_marker ;
				if ( $type == 1 ) {
					$col = ( $ch eq lc $ch ) ? $single_color : $red ;
				} elsif ( $type == 2 ) {
					$col = ( $ch eq lc $ch ) ? $inversion_color : $red ;
				} elsif ( $type == 3 ) {
					$col = ( $ch eq lc $ch ) ? $cigar_color : $red ;
				} else {
					$col = ( $ch eq lc $ch ) ? $blue : $red ;
				}
				my $x1 = int ( $width * $pos / $ft ) ;
				my $y = int ( $height - ( $band + 3 ) * $charheight - 1 ) ;
				if ( $dist_marker ) {
					$im->filledRectangle ( $x1 , $y+2 , $x1+6 , $y+11 , $ltgrey ) ;
				}
				if ( $ch eq '_' ) {
					next unless $display_pair_links ;
					my $x2 = int ( $width * ( $pos + 1 ) / $ft - 1 ) ;
					$im->line ( $x1 , $y + 6 , $x2 , $y + 6 , $ltgrey ) ;
				} else {
					$im->string ( gdSmallFont , $x1 , $y , $ch , $col ) ;
				}
			}
		}
	} else {
		foreach my $band ( 0 .. $nob ) {
			$bands[$band] =~ s/\s+$// ; # Trim trailing space
			my @theband = split // , $bands[$band] ;
			my @thetype = split // , $btype[$band] ;
			$bands[$band] = '' ;
			my $pos = -1 ;
			my $lastx = -1 ;
			my $y = $height - $band - 1 ;
			my $col ;

			foreach my $ch ( @theband ) {
				$pos++ ;
				next if $ch eq ' ' ;
				my $type = ord ( $thetype[$pos] ) - 65 ;
				$type -= 100 if $type >= 100 ;
				if ( $type == 1 ) {
					$col = $single_color ;
				} elsif ( $type == 2 ) {
					$col = $inversion_color ;
				} elsif ( $type == 3 ) {
					$col = $cigar_color ;
				} elsif ( $ch eq '_' ) {
					next unless $display_pair_links ;
					$col = $ltgrey ;
				} else {
					$col = $blue ;
				}

				my $x1 = int ( $width * $pos / $ft ) ;
				next if $x1 <= $lastx ;
				$lastx = $x1 ;
				my $x2 = int ( $width * ( $pos + 1 ) / $ft - 1 ) ;
				if ( $x1 < $x2 ) {
					$im->line ( $x1 , $y , $x2 , $y , $col ) ;
				} else {
					$im->setPixel ( $x1 , $y , $col ) ;
				}
			}

			$pos = -1 ;
			$lastx = -1 ;
			foreach my $ch ( @theband ) {
				$pos++ ;
				next unless ( $ch ge 'A' and $ch le 'Z' ) or ( $ch lt ' ' ) ;
				$ch = chr ( ord ( $ch ) + 64 ) if $ch lt ' ' ;
				my $x1 = int ( $width * $pos / $ft ) ;
				next if $x1 <= $lastx ;
				$lastx = $x1 ;
				my $x2 = int ( $width * ( $pos + 1 ) / $ft - 1 ) ;
				if ( $x1 < $x2 ) {
					$im->line ( $x1 , $y , $x2 , $y , $red ) ;
				} else {
					$im->setPixel ( $x1 , $y , $red ) ;
				}
			}
		}
	}
	
	$height += 20 ;
	draw_h_axis ( $im , $black ) ;

	my $elapsed = tv_interval ( $time0 );
	$im->string ( gdSmallFont , 5 , 0 , "Rendering time : $elapsed, $total_reads total reads" , $black ) if $debug_output ;

	write_png ( $im ) ;
}

sub pile2bands {
	my ( $r_reads , $r_bands , $r_btype , $rl , $mode , $paired_pileup , $fragment ) = @_ ;

	if ( $paired_pileup == 1 ) {
		foreach my $read ( @{$r_reads} ) {
			my ( $rl1 , $rl2 ) = ( $read->[scalar(@{$read})-2] , $read->[scalar(@{$read})-1] ) ;
			my ( $from1 , $from2  ) = ( $read->[1] , $read->[2] ) ;
			next if $from2 < 0 ; # No single reads...
			next if $from1 + $rl < $from and $from2 > $to ;
			my $dist = $from2 - $from1 + 1 + $rl2 ;
			my $between = $from2 - $from1 - $rl1 ;
			my $seq = pileup_my_refseq_part ( $from1 , $read->[3] , $rl1 ) ;
			$seq .= '_' x $between ;
			$seq .= pileup_my_refseq_part ( $from2 , $read->[4] , $rl2 ) ;
			$seq =~ tr/ /_/ ;
			next if $seq =~ /^_*$/ ;
			add_pileup ( $from1 , $seq , length ( $seq ) , $r_bands , $r_btype , $mode , $dist ) ;
		}
	} else {
		foreach my $read ( @{$r_reads} ) {
			my ( $rl1 , $rl2 ) = ( $read->[scalar(@{$read})-2] , $read->[scalar(@{$read})-1] ) ;
			my ( $from1 , $from2  ) = ( $read->[1] , $read->[2] ) ;
			my $dist = $from2 > 0 ? $from2 - $from1 + 1 + $rl2 : $fragment ;
			
#			print ":" . join ( ',' , @{$read} ) . "\n" ;
			
			if ( $from1 <= $to and $from1 >= $from - $rl1 ) {
				add_pileup ( $from1 , $read->[3] , $rl1 , $r_bands , $r_btype , $mode , $dist ) ;
			}
			if ( $from2 <= $to and $from2 >= $from - $rl2 ) {
				add_pileup ( $from2 , $read->[4] , $rl2 , $r_bands , $r_btype , $mode , $dist ) ;
			}
		}
	}
}

sub pileup_my_refseq_part {
	my ( $rfrom , $seq , $rl ) = @_ ;
	return $seq if $seq ;
	my $start = $rfrom - $from ;
	$seq = '' ;
	while ( $rl > 0 and $start < 0 ) {
		$start++ ;
		$rl-- ;
		$seq .= 'x' ;
	}
	$seq .= substr lc $refseq , $start , $rl ;
	return $seq ;
}

sub add_pileup {
	my ( $start , $read_sequence , $rl , $rbands , $rbtype , $mode , $dist ) = @_ ;
	my $min ;
	my $max ;
	my $out = '' ;
	my $reflen = length $refseq ;
	foreach ( 0 .. $rl-1 ) {
		my $onref = $start + $_ - $from ;
		next if $onref < 0 ;
		last if $onref >= $reflen ;
		$min = $onref unless defined $min ;
		$max = $onref ;
		if ( $read_sequence ) {
			$out .= substr ( $read_sequence , $_ , 1 )  ;
		} else {
			$out .= lc substr ( $refseq , $onref , 1 ) ;
		}
	}
#	print "$read_sequence\t$out\n" ;
	return unless defined $min ;
	
	my $use_band ;
	my $limit = scalar ( @{$rbands} ) - 1 ;
	foreach ( 0 .. $limit ) {
		next unless substr ( $rbands->[$_] , $min , $max - $min + 1 ) =~ m/^ +$/ ;
		$use_band = $_ ;
		last ;
	}
	
	unless ( defined $use_band ) {
		$use_band = scalar ( @{$rbands} ) ;
		push @{$rbands} , ( ' ' x $reflen ) ;
		push @{$rbtype} , ( ' ' x $reflen ) ;
	}
	
	if ( $mode != 1 and $mode != 3 ) {
		$mode += 100 if $dist > $cur_meta{'fragment'} + $cur_meta{'variance'} or $dist < $cur_meta{'fragment'} - $cur_meta{'variance'}  ;
	}
	substr ( $rbands->[$use_band] , $min , length $out ) = $out ;
	substr ( $rbtype->[$use_band] , $min , length $out ) = chr ( $mode + 65 ) x length ( $out ) ;
}


#________________________________________________________________________________________________________________________________________________

sub add_coverage {
	my ( $r_reads , $reflen , $single ) = @_ ;
	
	my ( $start , $end , $rl , $le ) ;
	foreach my $read ( @{$r_reads} ) {
		$le = scalar ( @{$read} ) - 1 ;
		if ( $single ) {
			$rl = $read->[$rl] ;
		} else {
			$rl = $read->[$rl-1] ;
		}
		
		$start = $read->[1] - $from ;
		$end = $start + $rl - 1 ;
		$start = 1 if $start < 1 ;
		$coverage[$_]++ foreach ( $start .. $end ) ;
		
		next if $single ;
		$rl = $read->[$rl] ;
		$start = $read->[2] - $from ;
		$end = $start + $rl - 1 ;
		$start = 1 if $start < 1 ;
		$end = $#coverage if $end > $#coverage ;
		$coverage[$_]++ foreach ( $start .. $end ) ;
	}
}


sub dump_image_coverageview {
	my $ft = $to - $from + 1 ;
	$coverage[$_] = 0 foreach ( 0 .. $ft+1000 ) ;

	foreach my $current_db ( 0 .. $#all_meta ) {
		my %meta = %{$all_meta[$current_db]} ;
		my $rl = $meta{'read_length'} ;
		my $smr = $all_smr[$current_db] ;
		my $pmr = $all_pmr[$current_db] ;
		my $sin = $all_sin[$current_db] ;
		my $inv = $all_inv[$current_db] ;
		
		# CIGAR
		foreach ( @{$all_cig[$current_db]} ) {
			my $start = $_->[1] ;
			my $seq = decompile_cigar ( $_->[4] ) ;
			$seq =~ tr/I// ; # Can't do insertions - yet!
			my $end = $start + length ( $seq ) - 1 ;
			foreach ( $start .. $end ) {
				next if $_ < $from or $_ > $to ;
				next if substr ( $seq , $_ - $start , 1 ) ne 'M' ;
				$coverage[$_-$from]++ ;
			}
		}
		
		# Fix single reads
		$sin->[$_]->[2] = -1000 foreach ( 0 .. @$sin ) ;
		add_coverage ( $smr , $ft , 0 ) ;
		add_coverage ( $pmr , $ft , 0 ) ;
		add_coverage ( $sin , $ft , 1 ) ;
		add_coverage ( $inv , $ft , 0 ) ;
	}

	$height = 0 ;
	foreach ( @coverage ) {
		$height = $_ if $_ > $height ;
	}
#	$scale_height = $display_noscale ? 0 : 20 ;
	$height += $scale_height ; # Scale

	my $im = new GD::Image ( $width , $height ) ;
	my $white = $im->colorAllocate ( 255 , 255 , 255 ) ;
	my $black = $im->colorAllocate ( 0 , 0 , 0 ) ;
	my $blue = $im->colorAllocate ( 0 , 0 , 255 ) ;

	my @max_per_pixel ;
	$max_per_pixel[$_] = 0 foreach ( 0 .. $width ) ;
	foreach my $pos ( 0 .. $ft-1 ) {
		my $x1 = int ( $width * $pos / $ft ) ;
		$max_per_pixel[$x1] = $coverage[$pos] if $max_per_pixel[$x1] < $coverage[$pos] ;
	}
	
	if ( $ft >= $width ) {
		foreach ( 0 .. $width-1 ) {
			my $y = $height - $scale_height - $max_per_pixel[$_] ;
			$im->line ( $_ , $y , $_ , $height - $scale_height , $blue ) ;
		}
	} else {
		foreach my $pos ( 0 .. $ft-1 ) {
			my $x1 = int ( $width * $pos / $ft ) ;
			next if $max_per_pixel[$x1] > $coverage[$pos] ;
			my $x2 = int ( $width * ( $pos + 1 ) / $ft - 1 ) ;
			my $y = $height - $scale_height - $coverage[$pos] ;
			if ( $x2 > $x1 ) {
				$im->filledRectangle ( $x1 , $y , $x2 , $height - $scale_height , $blue ) ;
			} else {
				$im->line ( $x1 , $y , $x1 , $height - $scale_height , $blue ) ;
			}
		}
	}

	unless ( $display_noscale ) {
		draw_h_axis ( $im , $black ) ;
	}
	foreach ( 0 .. $height - $scale_height ) {
		next unless $_ % 50 == 0 ;
		my $y = $height - $scale_height - $_ ;
		$im->line ( 0 , $y , 5 , $y , $black ) ;
		$im->string ( gdSmallFont , 8 , $y-7 , $_ , $black ) ;
	}


	show_debugging_output ( $im , $black ) ;

	write_png ( $im ) ;
}


#________________________________________________________________________________________________________________________________________________


sub dump_image_indelview {
	my $ft = $to - $from + 1 ;
	my $text_mode = $width / $ft > 5 ;
	my $max_dist = 1 ; # Dummy value
	
	my $im = new GD::Image ( $width , $height ) ;
	my $white = $im->colorAllocate ( 255 , 255 , 255 ) ;
	my $black = $im->colorAllocate ( 0 , 0 , 0 ) ;
	my $blue = $im->colorAllocate ( 0 , 0 , 255 ) ;
	my $light_blue = $im->colorAllocate ( 0x44 , 0xB4 , 0xD5 ) ;
	my $red = $im->colorAllocate ( 255 , 0 , 0 ) ;
	my $single_color = $im->colorAllocate ( 0 , 128 , 0 ) ;
	my $inversion_color = $im->colorAllocate ( 0xD2 , 0x69 , 0x1E ) ; # D2691E
	my $variance_color = $im->colorAllocate ( 0xFF , 0xF6 , 0x8F ) ; # FFF68F
	my $inversion_left_color = $im->colorAllocate ( 0xB6 , 0xBA , 0x18 ) ; # B6BA18
	my $inversion_right_color = $im->colorAllocate ( 0x7C , 0xEB , 0x98 ) ; # 7CEB98
	my $inversion_middle_color = $im->colorAllocate ( 0x99 , 0xD2 , 0x58 ) ; # 99D258
	my $orange = $im->colorAllocate ( 0xFF , 0xAA , 0x00 ) ;
	my $grey = $im->colorAllocate ( 200 , 200 , 200 ) ;
	
	my @ann_color = ( $black , $blue , $red , $single_color , $orange , $inversion_right_color, $inversion_middle_color , $variance_color  , $inversion_left_color , $inversion_color ) ;
	
	if ( $number_of_databases > 0 ) {
		my %meta = %{$all_meta[0]} ;
		$max_dist = 0 ;
		
		foreach my $current_db ( 0 .. $number_of_databases - 1 ) {
			my $new_max_dist = $cgi->param('maxdist') || $meta{'fragment'} + $meta{'variance'} * 2 ;
			$new_max_dist = 1200 if scalar ( @{$all_cig[$current_db]} ) > 0 ;
			$max_dist = $new_max_dist if $max_dist < $new_max_dist ;
		}
		
		$max_dist = $height if $max_dist < $height ; # No loss in detail, but avoiding ugly rendering artifacts
		$im->filledRectangle (	0 , $height - ( $meta{'fragment'}+$meta{'variance'} ) * $height / $max_dist ,
								$width , $height - ( $meta{'fragment'}-$meta{'variance'} ) * $height / $max_dist , $variance_color ) ;
	}
	
	# Single matches
	if ( $display_single ) {
		my @stack ;
		$stack[$_] = 0 foreach ( 0 .. $width ) ;
		my @stacklen ;
		foreach my $current_db ( 0 .. $number_of_databases - 1 ) {
			my %meta = %{$all_meta[$current_db]} ;
			my $sin = $all_sin[$current_db] ;
			my $upper = $height - ( $meta{'fragment'} + $meta{'variance'} ) ;
			foreach ( @{$sin} ) {
				my $p1 = $_->[1] ;
				my $x1 = int ( ( $p1 - $from ) * $width / $ft ) ;
				#$x1 = 0 if $x1 < 0 ;
				next if $x1 < 0 ;
				#$stack[$_]++ foreach ( $x1 .. $x1 + $rl ) ;
				$stack[$x1]++ ;
				my $rl = $_->[scalar(@{$_})-1] ;
				my $len = int ( $rl * $width / $ft ) ;
				$len = 1 if $len < 1 ;
				$stacklen[$x1] = $len ;
			}
			foreach ( 0 .. $width-1 ) {
				next unless $stack[$_] ;
				$im->line ( $_ , 0 , $_ , $stack[$_] , $single_color ) if  $stacklen[$_]  == 1 ;
				$im->filledRectangle ( $_ , 0 , $_+$stacklen[$_] , $stack[$_] , $single_color ) if $stacklen[$_] > 1 ;
			}
		}
	}

	# Grey lines between link pairs
	if ( $display_pair_links ) {
		foreach my $current_db ( 0 .. $number_of_databases - 1 ) {
			my %meta = %{$all_meta[$current_db]} ;
			my $rl = $meta{'read_length'} ;
			my $len = int ( $rl * $width / $ft ) ;
			$len = 1 if $len < 1 ;
			draw_indel_pair_links ( $im , $all_pmr[$current_db] , $grey , $max_dist , $rl , $ft , $len ) if $display_perfect ;
			draw_indel_pair_links ( $im , $all_smr[$current_db] , $grey , $max_dist , $rl , $ft , $len ) if $display_snps ;
			draw_indel_pair_links ( $im , $all_inv[$current_db] , $grey , $max_dist , $rl , $ft , $len ) if $display_inversions ;
		}
	}

	# Inversion boundaries
	if ( $display_inversions_ext ) {
		foreach my $current_db ( 0 .. $number_of_databases - 1 ) {
			my %meta = %{$all_meta[$current_db]} ;
			my $rl = $meta{'read_length'} ;
			my $len = int ( $rl * $width / $ft ) ;
			$len = 1 if $len < 1 ;
			$avg_frag_size = $meta{'fragment'} || 0 ;
			draw_indel_matches ( $im , $all_inv[$current_db] , $inversion_color , $max_dist , $rl , $ft , $len , 1 , $meta{'dbversion'} > 2 ) ;
		}
		foreach my $x ( 0 .. $width ) {
			my $y = $height - 20 ;
			if ( $inversion_middle[$x] > 0 ) {
				$im->line ( $x , $y , $x , $y - $inversion_middle[$x] , $inversion_middle_color ) ;
#				$y -= $inversion_middle[$x] ;
			}
			if ( $inversion_left[$x] > 0 ) {
				$im->line ( $x , $y , $x , $y - $inversion_left[$x] , $inversion_left_color ) ;
#				$y -= $inversion_left[$x] ;
			}
			if ( $inversion_right[$x] > 0 ) {
				$im->line ( $x , $y , $x , $y - $inversion_right[$x] , $inversion_right_color ) ;
#				$y -= $inversion_right[$x] ;
			}
		}
	}
	
	# Blue/magenta read lines
	foreach my $current_db ( 0 .. $number_of_databases - 1 ) {
		my %meta = %{$all_meta[$current_db]} ;
		my $rl = $meta{'read_length'} ;
		my $len = int ( $rl * $width / $ft ) ;
		$len = 1 if $len < 1 ;
		$avg_frag_size = $meta{'fragment'} || 0 ;
		
		# Perfect and SNP matches
		draw_inline_annotation ( $im , $all_ann[$current_db] , \@ann_color , $max_dist , $rl , $ft ) if $display_inline_annotation ;
		draw_indel_matches ( $im , $all_pmr[$current_db] , $blue , $max_dist , $rl , $ft , $len ) if $display_perfect ;
		draw_indel_matches ( $im , $all_smr[$current_db] , $blue , $max_dist , $rl , $ft , $len ) if $display_snps ;
		draw_indel_matches ( $im , $all_inv[$current_db] , $inversion_color , $max_dist , $rl , $ft , $len , 0 , $meta{'dbversion'} > 2 ) if $display_inversions ;
	}
	
	# Moving annotation average
	if ( 0 < scalar @anno_sum ) {
		foreach my $col ( 0 .. $#anno_sum ) {
			my ( @x , @y ) ;
			foreach my $pos ( 0 .. $width ) {
				if ( defined $anno_sum[$col]->[$pos] ) {
					$y[$pos] = $anno_sum[$col]->[$pos] / $anno_count[$col]->[$pos] ;
					push @x , $pos ;
				}
			}
			next if scalar ( @x ) == 0 ;
			foreach my $p ( 1 .. $#x ) {
				my $dist = $x[$p] - $x[$p-1] ;
				my $diff = $y[$x[$p]] - $y[$x[$p-1]] ;
				foreach ( $x[$p-1]+1 .. $x[$p]-1 ) {
					$y[$_] = $diff * ( $_ - $x[$p-1] ) / $dist + $y[$x[$p-1]] ;
				}
			}
			foreach ( 0 .. $x[0] ) {
				$y[$_] = $y[$x[0]] ;
			}
			foreach ( $x[$#x] .. $width ) {
				$y[$_] = $y[$x[$#x]] ;
			}

			# Smoothing
			foreach my $iter ( 0 .. 100 ) {
				foreach ( 1 .. $width-2 ) {
					$y[$_] = ( $y[$_-1] + $y[$_] + $y[$_+1] ) / 3 ;
				}
			}

			foreach ( 1 .. $width-1 ) {
				$im->line ( $_-1 , int($y[$_-1]) , $_ , int($y[$_]) , $ann_color[$col] ) ;
			}
		}
	}
	
	# CIGAR
	foreach my $current_db ( 0 .. $number_of_databases - 1 ) {
		my $cig = $all_cig[$current_db] ;
		foreach ( @{$cig} ) {
			my ( $name , $start , $stop , $orientation , $data ) = @{$_} ;
			my $x1 = ( $start - $from ) * $width / $ft ;
			my $x2 = ( $stop - $from ) * $width / $ft ;
			my $y = $height - ( $stop - $start + 1 ) * $height / $max_dist ;
			draw_indel_cigar ( $im , $x1 , $x2 , $y , $data , $light_blue , $red , $single_color ) ;
		}
	}
	
	# Uniqueness
	if ( $display_uniqueness ) {
		my $u = get_chromosome_part ( $uniqueness_file , $orig_chr , $from , $to ) ;
#		print $cgi->header(-type=>'text/plain',-expires=>'-1s'); print $u ; exit ;
		my $box = ( $width / $ft ) > 1 ;
		foreach my $p ( 0 .. length ( $u ) - 1 ) {
			my $y = int ( substr ( $u , $p , 1 ) ) ;
			my $x = $width * $p / $ft ;
			$im->line ( $x , $height , $x , $height - $y * 2 , $red ) unless $box ;
			$im->filledRectangle ( $x , $height - $y * 2 , $width * ( $p + 1 ) / $ft , $height , $red ) if $box ;
		}
	}

	# SNPs
	foreach my $current_db ( 0 .. $number_of_databases - 1 ) {
		my $smr = $all_smr[$current_db] ;
		my %meta = %{$all_meta[$current_db]} ;
		my $rl = $meta{'read_length'} ;

		# SNPs for SNP matches
		if ( $display_snps ) {
			foreach ( @{$smr} ) {
				my ( $rl1 , $rl2 ) = ( $_->[scalar(@{$_})-2] , $_->[scalar(@{$_})-1] ) ;
				my $seq1 = $_->[3] ;
				my $seq2 = $_->[4] ;
				next unless $seq1 or $seq2 ;
				my $p1 = $_->[1] ;
				my $p2 = $_->[2] ;
				my $ofs = $p2 - $p1 + $rl2 ; # Observed fragment size
				my $y = $height - $ofs * $height / $max_dist ;
				draw_indel_snps ( $im , $seq1 , $y , $p1 , $ft , $text_mode , $red ) if $seq1 ;
				draw_indel_snps ( $im , $seq2 , $y , $p2 , $ft , $text_mode , $red ) if $seq2 ;
			}
		}
	}
	
	draw_h_axis ( $im , $black ) ;
	
	# Vertical axis
	my $diff = 100 ;
	$diff = 500 if $max_dist > 2500 ;
	$diff = 2500 if $max_dist > 25000 ;
	foreach ( 1 .. $max_dist ) {
		next unless $_ % $diff == 0 ;
		my $y = $height - $_ * $height / $max_dist ;
		$im->line ( 0 , $y , 10 , $y , $black ) ;
		$im->string ( gdSmallFont , 13 , $y > 5 ? $y - 6 : -3 , $_ , $black ) ;
	}

	show_debugging_output ( $im , $red ) ;

	write_png ( $im ) ;
}

sub draw_inline_annotation {
	my ( $im , $dbref , $colref , $max_dist , $rl , $ft ) = @_ ;
	my $show_text = $ft <= 25000 ? 1 : 0 ;
	foreach ( @{$dbref} ) { # text,pos,height,length,colorcode
		my $p1 = $_->[1] ;
		my $p2 = $p1 + $_->[3] ;
		my $ofs = $_->[2] ;
		my $colorcode = $_->[4] ;
		my $y = $height - $ofs * $height / $max_dist ;
		my $x1 = ( $p1 - $from ) * $width / $ft ;
		my $x2 = ( $p2- $from ) * $width / $ft ;
		my $draw = $colorcode > 0 ;
		$x1 = 0 if $x1 < 0 ;
		next if $x2 < $x1 ;
		if ( $draw ) {
			$im->line ( $x1 , $y , $x2 , $y , $colref->[$colorcode] ) if $x2 > $x1 ;
			$im->setPixel ( $x1 , $y , $colref->[$colorcode] ) if $x2 == $x1 ;
		}
		foreach ( $x1 .. $x2 ) {
			$anno_sum[$colorcode]->[$_] += $y ;
			$anno_count[$colorcode]->[$_]++ ;
		}
#		$im->string ( gdSmallFont , $x1 , $y , $_->[0] , $draw ? $colref->[$colorcode] : $colref->[0] ) if $show_text ;
	}
}

sub draw_indel_cigar {
	my ( $im , $x1 , $x2 , $y , $data , $perfect_col , $del_col , $ins_col ) = @_ ;
	my $ft = $to - $from + 1 ;
	my $w = $x2 - $x1 + 1 ;
	my $l = 0 ;
	my $show_text = $ft <= 25000 ? 1 : 0 ;
	my @i = split /\|/ , $data ;
	foreach ( @i ) {
		$_ =~ m/^(.)(\d+)$/ ;
		$l += int ( $2 ) unless $1 eq 'I' ;
	}
	my $abs_p = 0 ;
	foreach ( @i ) {
		$_ =~ m/^(.)(\d+)$/ ;
		my ( $mode , $num ) = ( $1 , int ( $2 ) ) ;
		my $xa = int ( $x1 + $abs_p * $w / $l ) ;
		my $xb = int ( $x1 + ($abs_p+$num) * $w / $l - 1 ) ;
		my $xw = $num * $width / $ft ;
		if ( $mode eq 'M' ) {
			$abs_p += $num ;
			$im->line ( $xa , $y , $xb , $y , $perfect_col ) ;
		} elsif ( $mode eq 'D' ) {
			$abs_p += $num ;
			$im->line ( $xa , $y , int ( ( $xb + $xa ) / 2 ) , $y - $num , $del_col ) ;
			$im->line ( $xb , $y , int ( ( $xb + $xa ) / 2 ) , $y - $num , $del_col ) ;
			$im->string ( gdSmallFont , int ( $xa + $xw / 2 - 3 ) , $y - $num - 12 , $num , $del_col ) if $show_text ;
		} elsif ( $mode eq 'I' ) {
			$im->line ( $xa , $y , $xa , $y + $num , $ins_col ) ;
			$im->line ( $xa , $y , $xa + $xw , $y + $num , $ins_col ) ;
			$im->string ( gdSmallFont , $xa , $y + $num , $num , $ins_col ) if $show_text ;
		}
	}
}

sub decompile_cigar {
	my ( $in ) = @_ ;
	my @i = split /\|/ , $in ;
	my $ret = '' ;
	foreach ( @i ) {
		$_ =~ m/^(.)(\d+)$/ ;
		$ret .= $1 x int ( $2 ) ;
	}
	return $ret ;
}

sub show_debugging_output {
	return unless $debug_output ;
	my ( $im , $red ) = @_ ;
	my $elapsed = tv_interval ( $time0 );
	$im->string ( gdSmallFont , 5 , 0 , "Rendering time : $elapsed, $total_reads total reads" , $red ) ;
	$im->string ( gdSmallFont , 5 , 10 , "Databases : $database" , $red ) ;
}

sub draw_indel_pair_links {
	my ( $im , $arr_ref , $color , $max_dist , $rl , $ft , $len ) = @_ ;
	return if scalar ( @{$arr_ref} ) > 100000 ; # Hard cutoff # $len == 1 and 
	foreach ( @{$arr_ref} ) {
		my ( $rl2 ) = ( $_->[scalar(@{$_})-1] ) ;
		my $p1 = $_->[1] ;
		my $p2 = $_->[2] ;
		my $ofs = $p2 - $p1 + $rl2 ; # Observed fragment size
		my $y = $height - $ofs * $height / $max_dist ;
		my $x1 = ( $p1 - $from ) * $width / $ft ;
		my $x2 = ( $p2 - $from ) * $width / $ft ;
		next if $x1 + 1 >= $x2 ;
		if ( $x1 == $x2 ) {
			$im->setPixel ( $x1 , $y , $color ) ;
		} else {
			$im->line ( $x1 , $y , $x2 , $y , $color ) ;
		}
	}
}

sub draw_indel_matches {
	my ( $im , $arr_ref , $color , $max_dist , $rl , $ft , $len , $dont_draw , $inversion ) = @_ ;
	my $inversion_height = 15 ; # pixel

	unless ( $dont_draw ) {
		foreach ( @{$arr_ref} ) {
			my ( $rl1 , $rl2 ) = ( $_->[scalar(@{$_})-2] , $_->[scalar(@{$_})-1] ) ;
			my $p1 = $_->[1] ;
			my $p2 = $_->[2] ;
			my $ofs = $p2 - $p1 + $rl2 ; # Observed fragment size
			my $y = $height - $ofs * $height / $max_dist ;
			next if $y < 0 ;
			my $x1 = ( $p1 - $from ) * $width / $ft ;
			my $x2 = ( $p2 - $from ) * $width / $ft ;
			
			my $len1 = int ( $rl1 * $width / $ft ) ;
			my $len2 = int ( $rl2 * $width / $ft ) ;
			
			if ( $len1 > 1 or $len2 > 1 ) {
				$im->line ( $x1 , $y , $x1 + $len1 , $y , $color ) ;
				$im->line ( $x2 , $y , $x2 + $len2 , $y , $color ) ;
			} elsif ( $len == 1 ) {
				unless ( $x1 < 0 or defined $dotcache[$x1]->[$y] ) {
					$im->setPixel ( $x1 , $y , $color ) ;
					$dotcache[$x1]->[$y] = 1 ;
				}
				unless ( $x2 < 0 or defined $dotcache[$x2]->[$y] ) {
					$im->setPixel ( $x2 , $y , $color ) ;
					$dotcache[$x2]->[$y] = 1 ;
				}
			}
		}
	}

	return unless $display_inversions_ext ;

	# BEGIN EXTENDED INVERSION DISPLAY
	foreach ( @{$arr_ref} ) {
		my ( $rl2 ) = ( $_->[scalar(@{$_})-1] ) ;
		my $p1 = $_->[1] ;
		my $p2 = $_->[2] ;
		my $ofs = $p2 - $p1 + $rl ; # Observed fragment size
		my $y = $height - $ofs * $height / $max_dist ;
		my $x1 = ( $p1 - $from ) * $width / $ft ;
		my $x2 = ( $p2 - $from ) * $width / $ft ;

		my $dir = $_->[5] ;
		my @dir = split '' , $dir ;
		my $xm = int ( ( $x1 + $x2 ) / 2 ) ;
		my ( $xf , $xt ) ;
		my $smaller = $p1 > $p2 ? $p2 : $p1 ;
		my $larger = $p1 < $p2 ? $p2 : $p1 ;
		my ( $int_var , $middle ) ;
		if ( $dir[1] eq '+' ) {
			$xf = int ( ( $smaller - $from ) * $width / $ft ) ;
			$xt = int ( ( $smaller - $from + $avg_frag_size ) * $width / $ft + $len ) ;
			$middle = $smaller + $ofs*3/4 ;
			$int_var = $ofs/4 ;
		}  elsif ( $dir[1] eq '-' ) {
			$xt = int ( ( $smaller - $from + $ofs ) * $width / $ft ) ;
			$xf = int ( ( $smaller - $from + $ofs - $avg_frag_size ) * $width / $ft + $len ) ;
			$middle = $larger - $ofs*3/4 ;
			$int_var = $ofs/4 ;
		} else { next } ;
		$xf = 0 if $xf < 0 ;
		if ( '+' eq $dir[1] ) {
			$im->arc ( $xm + $len , $y , ( $x2 - $xm ) * 2 , $inversion_height , 180 , 360 , $color ) unless $dont_draw ;
			foreach my $p ( $xf .. $xt ) { $inversion_left[$p]++ ; }
		} elsif ( '-' eq $dir[1] ) {
			$im->arc ( $xm , $y , ( $x2 - $xm ) * 2 , $inversion_height , 0 , 180 , $color-1 ) unless $dont_draw ;
			foreach my $p ( $xf .. $xt ) { $inversion_right[$p]++ ; }
		}
		
		$xf = int ( ( $middle - $int_var - $from ) * $width / $ft ) ;
		$xt = int ( ( $middle + $int_var - $from ) * $width / $ft ) ;
		$xf = 0 if $xf < 0 ;
		foreach my $p ( $xf .. $xt ) { $inversion_middle[$p]++ ; }
		
	}
	# END EXTENDED INVERSION DISPLAY
}

sub draw_indel_snps {
	my ( $im , $seq , $y , $p1 , $ft , $text_mode , $red ) = @_ ;
	return if $y < 0 ;

	if ( $text_mode ) {
		while ( $seq =~ m/[A-Z]/g ) {
			my $a = pos ( $seq ) - 1 ;
			my $x1 = int ( ( $p1 - $from + $a ) * $width / $ft ) ;
			$im->line ( $x1 , $y-2 , $x1 , $y , $red ) ;
			$im->string ( gdSmallFont , $x1 , $y-1 , $& , $red ) ;
		}
	} else {
		my $ww = int ( $width / $ft ) - 1 ;
		if ( $ww > 0 ) {
			while ( $seq =~ m/[A-Z]/g ) {
				my $a = pos ( $seq ) - 1 ;
				my $x1 = int ( ( $p1 - $from + $a ) * $width / $ft ) ;
				$im->filledRectangle ( $x1 , $y-2 , $x1 + $ww , $y+2 , $red ) ;
			}
		} else {
			return if $y < 0 ;
			while ( $seq =~ m/[A-Z]/g ) {
				my $a = pos ( $seq ) - 1 ;
				my $x1 = int ( ( $p1 - $from + $a ) * $width / $ft ) ;
				next if $x1 < 0 or defined $snpcache[$x1]->[$y] ;
				$im->line ( $x1 , $y-2 , $x1 , $y+2 , $red ) ;
				$snpcache[$x1]->[$y] = 1 ;
			}
		}
	}
}

sub draw_h_axis {
	my ( $im , $color ) = @_ ;
	my $ft = $to - $from + 1 ;

	# Show known SNPs
	if ( 0 < scalar @known_snps ) {
		$known_snp_color = $im->colorAllocate ( 150 , 150 , 150 ) unless defined $known_snp_color ;
		my $show_text = $ft < 300 ? 1 : 0 ;
		my $last_x = -1 ;
		foreach ( 0 .. $#known_snps ) {
			next unless defined $known_snps[$_] ;
			my $x = $_ * $width / $ft + 2 ;
			unless ( $show_text ) {
				next if $x == $last_x ;
				$im->line ( $x , $height - 19 , $x , $height - 12 , $known_snp_color ) ;
				$last_x = $x ;
				next ;
			}
			my $ch = substr ( $known_snps[$_] , 1 , 1 ) ;
			if ( defined $refseq ) {
				$ch = $iupac{$ch} ;
				my $orig = substr $refseq , $_ , 1 ;
				$ch =~ s/$orig//g ;
				$ch = $iupac_reverse{$ch} ;
			}
			$im->string ( gdSmallFont , $x - 2 , $height - 21 , $ch , $known_snp_color ) ;
		}
	}

	# Show horizontal axis
	my $step = 1 ;
	while ( $ft / $step > 10 ) {
		$step *= 5 if ( $ft / $step > 10 ) ;
		$step *= 2 if ( $ft / $step > 10 ) ;
	}
	for ( my $pos = $from - $from % $step ; $pos < $to ; $pos += $step ) {
		next if $pos < 2 ;
		my $x = ( $pos - $from ) * $width / $ft ;
		$im->line ( $x , $height - 12 , $x , $height - 17 , $color ) ;
		$im->string ( gdSmallFont , $x - 20 , $height - 12 , $pos , $color ) ;
	}
	
}

sub dump_image_annotation {
	my ( $return_url ) = @_ ;
	my $posx = $cgi->param('x') ;
	my $posy = $cgi->param('y') ;
	my $found_id ;
	
	my $dbfile = "$datapath/$annotation_file";
	my $dbh ;
	if ( defined $annotation_file ) {
		$dbh = DBI->connect(
		    "dbi:SQLite:dbname=$dbfile", # DSN: dbi, driver, database file
	    	"",                          # no user
		    "",                          # no password
		    { RaiseError => 1 },         # complain if something goes wrong
		) or die "ERROR\n".$DBI::errstr;
	}

	my $ft = $to - $from + 1 ;
	$height = 20 ;
	my $im = new GD::Image ( $width , $height ) ;
	my $white = $im->colorAllocate ( 255 , 255 , 255 ) ;
	my $black = $im->colorAllocate ( 0 , 0 , 0 ) ;
	my $blue = $im->colorAllocate ( 0 , 0 , 255 ) ;
	my $red = $im->colorAllocate ( 255 , 0 , 0 ) ;
	my $green = $im->colorAllocate ( 0 , 255 , 0 ) ;
	my $gray = $im->colorAllocate ( 200 , 200 , 200 ) ;

	#print $cgi->header(-type=>'text/plain',-expires=>'-1s');
	my $features = [] ;
	my $tags_arr = [] ;
	
	if ( defined $annotation_file ) {
		eval { $features = $dbh->selectall_arrayref ( "SELECT type,start,end,strand,id FROM chr_$chromosome WHERE start <= $to AND end >= $from" ) || [] } ;
	}
	
	my %fids ;
	my ( %ktags , %vtags ) ;
	foreach ( @{$features} ) {
		$fids{$_->[4]} = 1 ;
	}
	
	my $many_features = 0 ; #scalar ( @{$features} ) > 1500 ;
	unless ( $many_features ) {
		my $fid = join ',' , sort keys %fids ;
		eval { $tags_arr = $dbh->selectall_arrayref ( "SELECT id,key,value FROM tags WHERE id IN ( $fid )" ) || [] } ;
		foreach ( @{$tags_arr} ) {
			push @{$ktags{$_->[0]}} , $_->[1] ;
			push @{$vtags{$_->[0]}} , $_->[2] ;
		}
	}
	
	my %shown ;
	my $arrow_len = 4 ;
	my @annotation ;
	foreach ( @{$features} ) {
		my ( $type , $start , $end , $strand , $id ) = @{$_} ;
		next if $type eq 'source' ;
		my $col = $gray ;
		my $off = 1 ;
		if ( $type =~ m/^cds/i ) {
			$col = $blue  ;
			$off = 0 ;
		} elsif ( $type =~ m/^repeat/i ) {
			$col = $red ;
			$off = 2 ;
		} elsif ( $type =~ m/^centromere$/i ) {
			$col = $green ;
			$off = 3 ;
		}
		my $x1 = int ( ( $start - $from ) * $width / $ft ) ;
		my $x2 = int ( ( $end - $from ) * $width / $ft ) ;
		my $y1 = int ( $height * $off / 5 ) ;
		my $y2 = int ( $height * ( $off + 1 ) / 5 ) - 1 ;
		$im->filledRectangle ( $x1 , $y1 , $x2 , $y2 , $col ) ;
		
		if ( $return_url and $col != $gray and $posx >= $x1 and $posx <= $x2 ) { # IGNORES MISC_FEATURE!!!
			if ( ( $posy >= $y1 and $posy <= $y2 ) or ( not defined $found_id ) ) {
				$found_id = $id ;
			}
		}
		
		next if $many_features ;

		my $midx = int ( ( $x1 + $x2 ) / 2 ) ;
		$midx = $width - $arrow_len*2 if $midx + $arrow_len*2 >= $width ;
		$midx = $arrow_len*2 if $midx < $arrow_len*2 ;
		my $is_small = $midx - $arrow_len * 3 < $x1 or $midx + $arrow_len * 3 > $x2 ;
		$is_small = 1 if $x2 - $x1 < 10 ;
		unless ( $is_small ) {
			my $midy = int ( ( $y1 + $y2 ) / 2 + 1 ) ;
			$im->line ( $midx - $arrow_len , $midy , $midx + $arrow_len, $midy , $white ) ;
			if ( $strand == '+' ) {
				$im->line ( $midx + $arrow_len-2 , $y1 , $midx + $arrow_len , $midy , $white ) ;
				$im->line ( $midx + $arrow_len-2 , $y2 , $midx + $arrow_len , $midy , $white ) ;
			} else {
				$im->line ( $midx - $arrow_len+2 , $y1 , $midx - $arrow_len , $midy , $white ) ;
				$im->line ( $midx - $arrow_len+2 , $y2 , $midx - $arrow_len , $midy , $white ) ;
			}
		}
		
		# Now show annotation
		next if $x2 - $x1 < 15 ;
		next unless defined $ktags{$id} ;
		my ( $ann , $sysid ) ;
		my @karr = @{$ktags{$id}} ;
		foreach ( 0 .. $#karr ) {
			my $v = $vtags{$id}->[$_] ;
			if ( $karr[$_] eq 'systematic_id' ) {
				$sysid = $v ;
			} elsif ( $karr[$_] eq 'primary_name' ) {
				$ann = $v ;
			}
		}
		
		if ( defined $sysid ) {
			if ( defined $ann ) {
				$ann .= " ($sysid)" ;
			} else {
				$ann = $sysid ;
			}
		}
		
		next unless $ann ;
		next if defined $shown{$ann} and $is_small ;
		$shown{$ann} = 1 ;
		my $annx = $x1 ;
		$annx = 1 if $x1 < 1 ;
		push @annotation , "$annx\t$y2\t$ann" ;
#		$im->string ( gdSmallFont , $annx , $y2 , $ann , $black ) ;
	}
	
	foreach ( @annotation ) {
		my ( $x , $y , $text ) = split "\t" , $_ ;
		$im->string ( gdSmallFont , $x , $y , $text , $black ) ;
	}
	
	if ( $return_url ) {
		print $cgi->header(-type=>'text/plain',-expires=>'-1s');
		if ( defined $found_id ) {
			my @out ;
			my @keys = @{$ktags{$found_id}} ;
			my @values = @{$vtags{$found_id}} ;
			foreach ( 0 .. $#keys ) {
				my $k = $keys[$_] ;
				my $v = $values[$_] ;
				if ( $k eq 'systematic_id' or $k eq 'previous_systematic_id' or $k eq 'synonym' ) {
					my $url = "http://plasmodb.org/plasmo/showRecord.do?name=GeneRecordClasses.GeneRecordClass&primary_key=$v&project_id=PlasmoDB" ;
					my $kn = 'Systematic ID' ;
					$kn = 'Previous ID' if $k eq 'previous_systematic_id';
					$kn = 'Synonym' if $k eq 'synonym';
					push @out , "$kn : <a href=\"$url\" target='_blank'>$v</a>" ;
				} elsif ( $k eq 'db_xref' ) {
					if ( $v =~ m/^PMID:(.+)$/ ) {
						my $url = "http://www.ncbi.nlm.nih.gov/pubmed/$1" ;
						push @out , "PubMed : <a href=\"$url\" target='_blank'>$1</a>" ;
					} else {
						push @out , "Database reference : $v" ;
					}
				} elsif ( $k eq 'primary_name' ) {
					unshift @out , "<b>Name : $v</b>" ;
				} elsif ( $k eq 'ID' ) {
					push @out , "ID : $v" ;
				} else {
					push @out , "$k : $v" ;
				}
			}
			print join "\n" , @out ;
		}
		return ;
	}

	write_png ( $im ) ;
}

sub dump_image_gc {
	my $range = 1+11*2 ;
	my $from2 = $from - $range ;
	my $to2 = $to + $range ;
	$from2 = 1 if $from2 < 1 ;
	$refseq = get_chromosome_part ( $genome_file , $orig_chr , $from2 , $to2 ) ;

	my $ft = $to - $from + 1 ;

	my @gc ;
	foreach my $start ( $from2 .. $to2 - $range ) {
		my $real = int ( $start + $range / 2 ) ;
		next unless $real >= $from and $real <= $to ;
#		my $s = substr $refseq , $start - $from2 , $range ;
		my $cnt = ( substr ( $refseq , $start - $from2 , $range ) =~ tr/CG// ) ;
#		my $percent = int ( $cnt * 100 / $range ) ;
		$gc[$real-$from] = int ( $cnt * 100 / $range ) ;
	}

	$height = 50 ;
	my $im = new GD::Image ( $width , $height ) ;
	my $white = $im->colorAllocate ( 255 , 255 , 255 ) ;
	my $black = $im->colorAllocate ( 0 , 0 , 0 ) ;
	my $blue = $im->colorAllocate ( 0 , 0 , 255 ) ;
	my $red = $im->colorAllocate ( 255 , 0 , 0 ) ;

	$im->line ( 0 , $height/2 , $width , $height/2 , $red ) ;
	
	$im->line ( 0 , $height-1 , 5 , $height-1 , $black ) ;
	$im->string ( gdSmallFont , 7 , $height - 12 , '0%' , $black ) ;
	$im->line ( 0 , $height/2 , 5 , $height/2 , $black ) ;
	$im->string ( gdSmallFont , 7 , $height/2 - 6 , '50%' , $black ) ;
	$im->line ( 0 , 0 , 5 , 0 , $black ) ;
	$im->string ( gdSmallFont , 7 , -2 , '100%' , $black ) ;
	
	$im->string ( gdSmallFont , $width - 30 , -2 , '[' . int ( $range / 2 ) . ']' , $black ) ;

	my ( $lastx , $lasty ) ;
	
	if ( $ft > $width ) {
		my $min = 100 ;
		my $max = 0 ;
		foreach my $pos ( 0 .. $#gc ) {
			my $x = int ( $pos * $width / $ft ) ;
			my $y = $height - int ( $gc[$pos] * $height / 100 ) ;
			
			if ( $x == $lastx ) {
				$min = $y if $y < $min ;
				$max = $y if $y > $max ;
			} else {
				$im->line ( $lastx , $min , $lastx , $max , $blue ) if $max >= $min ;
				$lastx = $x ;
				$min = $y ;
				$max = $y ;
			}
		}
		$im->line ( $lastx , $min , $lastx , $max , $blue ) if $max > $min ;
	} else {
		foreach my $pos ( 0 .. $#gc ) {
			my $x = int ( $pos * $width / $ft ) ;
			my $y = $height - int ( $gc[$pos] * $height / 100 ) ;
			$im->line ( $lastx , $lasty , $x , $y , $blue ) if defined $lastx ;
			( $lastx , $lasty ) = ( $x , $y ) ;
		}
	}

	my $elapsed = tv_interval ( $time0 );
	$im->string ( gdSmallFont , 5 , 0 , "Rendering time : $elapsed, $total_reads total reads" , $red ) if $debug_output ;

	write_png ( $im ) ;
}

sub in_array {
    my ($arr,$search_for) = @_;
    foreach my $value (@$arr) {
        return 1 if $value eq $search_for;
    }
    return 0;
}

sub write_png {
	my ( $im ) = @_ ;
	print $cgi->header(-type=>'image/png',-expires=>'-1s');
	binmode STDOUT;
	my $png = $im->png () ;
	print $png ;
	if ( $use_sanger_cache ) {
		$sanger_cache_db->set ( $png , $sanger_cache_key , $sanger_cache_hours ) ;
	}
}

# if ( $debug_output ) {
	# print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output
	# print tv_interval ( $time0 ); ;
	# exit ;
# }


if ( $output eq 'image' ) {
	if ( $view eq 'indel' ) {
		&dump_image_indelview ;
	} elsif ( $view eq 'pileup' ) {
		dump_image_pileupview ( 0 ) ;
	} elsif ( $view eq 'paired_pileup' ) {
		dump_image_pileupview ( 1 ) ;
	} elsif ( $view eq 'coverage' ) {
		&dump_image_coverageview ;
	} elsif ( $view eq 'annotation' ) {
		dump_image_annotation ( 0 ) ;
	} elsif ( $view eq 'gc' ) {
		&dump_image_gc ;
	}
} elsif ( $output eq 'url' ) {
	if ( $view eq 'annotation' ) {
		dump_image_annotation ( 1 ) ;
	}
} elsif ( $output eq 'text' ) {
	if ( $view eq 'pileup' ) {
		&dump_image_pileupview ;
	} else {
		&dump_text ;
	}
}
