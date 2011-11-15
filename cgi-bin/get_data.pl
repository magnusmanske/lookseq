#!/usr/local/bin/perl -T
#########
# Author:     Magnus Manske (mm6@sanger.ac.uk)  Team 112
#             Petr Danecek pd3@sanger.ac.uk     Team 145
#

my $DEBUG = 0;
if ( $DEBUG )
{
    # Run e.g. as get_data.pl from=95503225 to=95503325 chr=2 output=image width=700 lane=129S1_SvImJ.bam view=indel display='|perfect|snps|inversions|pairlinks|potsnps|uniqueness|gc|coverage|'
    BEGIN { push @INC, "/software/varinf/lib/", '.'; }
}
else
{
    BEGIN { push @INC, "."; }
}

use strict ;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use DBI;
use GD ;
use Time::HiRes qw(gettimeofday tv_interval);
use settings ;
use Data::Dumper ;
use LWP::Simple;
use Digest::MD5 qw(md5);
use File::Temp qw/ tempfile tempdir /;

my $cgi = new CGI;
my $time0 = [gettimeofday];
my $debug_output = 0 ;
$debug_output = $cgi->param('debug') if defined $cgi->param('debug') ;
if ( $debug_output ) {
	use warnings ;
}

#print STDERR "get_data.pl: " . CGI::url(-query=>1) ."\n";

$ENV{PATH}= '/usr/local/bin';   # solely to stop taint from barfing

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
my $mapq_cutoff = $cgi->param('mapq') ? $cgi->param('mapq') : 0;
my $sam_show_read_arrows = $display =~ m/\|orientation\|/ ;
if ( !($mapq_cutoff=~/^\d+$/) ) { $mapq_cutoff=0; }
my $max_insert_size = $cgi->param('maxdist') ? $cgi->param('maxdist') : 'Auto';
$max_insert_size =~ /(\d+)/; $max_insert_size=$1;

$from =~ /(\d+)/ ; $from = $1 ;
$to =~ /(\d+)/ ; $to = $1 ;


my $display_perfect = $display =~ m/\|perfect\|/ ;
my $display_snps = $display =~ m/\|snps\|/ ;
my $display_inversions = $display =~ m/\|inversions\|/ ;
my $display_inversions_ext = $display =~ m/\|inversions_ext\|/ ;
my $display_single = $display =~ m/\|single\|/ ;
my $display_pair_links = $display =~ m/\|pairlinks\|/ ;
my $display_pot_snps = $display =~ m/\|potsnps\|/ ;
my $display_noscale = $display =~ m/\|noscale\|/ ;
my $display_uniqueness = $display =~ m/\|uniqueness\|/ ;
my $display_faceaway = $display =~ m/\|faceaway\|/ ;
my $display_inline_annotation = 1 ; # FIXME

my $MODE_SINGLE       = 1;
my $MODE_DELETION     = 4;
my $MODE_INSERTION    = 5;
my $MODE_MAPQUAL_ZERO = 6;
my $height_nonzero_isize   = 0;
my $height_zero_isize  = 0;
my $max_pileup_rows = 100;   # In the pileup view, print only this many rows.

$display_inversions_ext = 0 unless $display_inversions ;

# Parameter paranoia
$from =~ s/\D//g ;
$to =~ s/\D//g ;
$chromosome =~ s/[^A-Za-z0-9\-_\.]//g ;
my $orig_chr = $chromosome ;
#$chromosome = "chr_$chromosome" if $chromosome =~ /^\d/ ;

unless ( $view eq 'annotation' or $view eq 'gc' ) {
	die "ERROR\nNot a valid lane\n" if $database eq '' ;
}
if ( $from >= $to )
{
    die "ERROR\nTO not larger than FROM .. from=" .$cgi->param('from'). " to=" .$cgi->param('to'). " \n" ;
}

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
my $reflength ;
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
my $max_dist ;
my $ft ;

# SAM/BAM variables
my $using_bam = 0 ;
my $cmd_samtools = "$execpath/samtools" ;
#$cmd_samtools = "/nfs/sf8/G1K/bin/samtools" unless !$DEBUG;
my $sam_show_read_quality = $display =~ m/\|readqual\|/ ;
my %sam_reads ;
my ( @sam_single , @sam_perfect , @sam_snps , @sam_inversions , @sam_capillary, @sam_mapqual_zero, @sam_mapqual_zero_pair, @sam_perfect_singles ) ;
my $sam_max_found_fragment = 0 ;
my $im ;
my ( $sam_white , $sam_black , $sam_col_single_read , $sam_col_mismatch , $sam_col_matching_read , $sam_col_inversion , $sam_col_read_pair_connection ) ;
my ( $sam_col_read_pair_quality , $sam_col_single_read_quality, $sam_col_mapqual_zero, $col_insertion, $col_deletion ) ;

if ( $view eq 'gc' ) 
{
    # Do not call samtools and parse the region when only the reference sequence is used for this.
    &dump_image_gc;
    exit;
}




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
#print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output

if ( defined $render_image_command and 1 == scalar ( @databases ) and -e $render_image_command ) { # Use C code instead
	exit ( 0 ) if render_bam_file () ;
}


foreach ( @databases ) {
	if ( $_ =~ /.bam$/ ) {
		$_ =~ /([a-zA-Z0-9_\-\.]+)/ ;
        my $file = $1;
        my $cwd;
        if ( $bam_ftp )
        {
            # Samtools require chdir to read bam .bai index files.
            $cwd = `/bin/pwd`;
            chomp($cwd);
            $cwd =~ /(.*)/;
            $cwd = $1;
            $datapath =~ /(.*)/;
            $datapath = $1;
            chdir($datapath) or die "Could not chdir $datapath: $!";
            $file = "$bam_ftp/$file";
        }
        else
        {
            $file = "$datapath/$file";
        }
        sam_read_data($file);

        if ( $bam_ftp )
        {
            chdir($cwd) or die "Could not chdir $cwd: $!";
        }
		$using_bam = 1 ;
		next ;
	}

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

if ( $using_bam ) 
{
	$refseq = get_chromosome_part( $genome_file , $orig_chr , $from , $to );
	&sam_bin ;
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
	
	if ( $file eq $reference_fa ) 
    {
		$chromosome =~ m/^([\w\s\d_]+)/;
		my $ch = $1 ;
        
        my $cwd;
        if ( $reference_ftp )
        {
            $cwd = `/bin/pwd`;
            chomp($cwd);
            $cwd =~ /(.*)/;
            $cwd = $1;
            $datapath =~ /(.*)/;
            $datapath = $1;
            chdir($datapath) or die "Could not chdir $datapath: $!";
            $file = "$reference_ftp/$file";
        }
        else
        {
            $file = "$datapath/$reference_fa";
        }

		my $cmd = "$cmd_samtools faidx $file $ch:$from-$to |" ;
		open FILE , $cmd ;
		my $s ;
		while ( <FILE> ) {
			next if $_ =~ /^>/ ;
			chomp ;
			$s .= uc $_ ;
		}
		close FILE ;

        if ( $reference_ftp )
        {
            chdir($cwd) or die "Could not chdir $cwd: $!";
        }
		return $s ;
	}
	
	# Default / fallback
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

    if ( !$refseq )
    {
        # Make these things locally - dump_gc_image needs different range
	    $refseq = get_chromosome_part( $genome_file , $orig_chr , $from , $to );
    }
	
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
	
	if ( $using_bam ) {
		$show_chars = $width / $ft > 4 ? 1 : 0 ;
        if ( !$show_chars ) { &dump_image_indelview; }
#	print $cgi->header(-type=>'text/plain',-expires=>'-1s'); 
		pile2bands_sam ( \@bands , \@btype , 0 , $paired_pileup , \@sam_perfect ) ;
		pile2bands_sam ( \@bands , \@btype , 0 , $paired_pileup , \@sam_perfect_singles ) ;
		pile2bands_sam ( \@bands , \@btype , 0 , $paired_pileup , \@sam_snps ) ;
		pile2bands_sam ( \@bands , \@btype , 2 , $paired_pileup , \@sam_inversions ) ;
		pile2bands_sam ( \@bands , \@btype , $MODE_SINGLE , $paired_pileup , \@sam_single ) ;
		pile2bands_sam ( \@bands , \@btype , $MODE_MAPQUAL_ZERO , $paired_pileup , \@sam_mapqual_zero ) ;
		pile2bands_sam ( \@bands , \@btype , $MODE_MAPQUAL_ZERO , $paired_pileup , \@sam_mapqual_zero_pair ) ;
#	exit ;
	}

	
	if ( $output eq 'text' ) {
		my $potseq = ' ' x length ( $refseq ) ;
		my $altseq = $refseq ;
		
		if ( $display_pot_snps ) {
			foreach ( 0 .. ( $to - $from ) ) {
				next unless defined $known_snps[$_] ;
				substr ( $potseq , $_ , 1 ) = substr $known_snps[$_] , 1 , 1 ;
				substr ( $altseq , $_ , 1 ) = 'N' ;
			}
		}
	
		print $cgi->header(-type=>'text/plain',-expires=>'-1s');
		print "$chromosome : $from - $to\n" ;
		print "$potseq\n" ;
		print "$refseq\n" ;
		print "$altseq\n" ;
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
	my $white = $im->colorAllocate ( @{$lscolor{'white'}} ) ;
	my $black = $im->colorAllocate ( @{$lscolor{'black'}} ) ;
	my $blue = $im->colorAllocate ( @{$lscolor{'blue'}} ) ;
	my $single_color = $im->colorAllocate ( @{$lscolor{'green2'}} ) ;
	my $inversion_color = $im->colorAllocate ( @{$lscolor{'brown'}} ) ;
	my $cigar_color = $im->colorAllocate ( @{$lscolor{'cyan'}} ) ;
	my $red = $im->colorAllocate ( @{$lscolor{'red'}} ) ;
	my $ltgrey = $im->colorAllocate ( @{$lscolor{'ltgrey'}} ) ;
    my $col_mapqual_zero = $im->colorAllocate( @{$lscolor{'mapqual_zero'}} );
	
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
				my $x1 = int ( $width * $pos / $ft ) ;
				my $y = int ( $height - ( $band + 3 ) * $charheight - 1 ) ;
				if ( $type == $MODE_SINGLE ) 
                {
					$col = ( $ch eq lc $ch ) ? $single_color : $red ;
				} 
                elsif ( $type == 2 ) {
					$col = ( $ch eq lc $ch ) ? $inversion_color : $red ;
				} elsif ( $type == 3 ) {
					$col = ( $ch eq lc $ch ) ? $cigar_color : $red ;
                }
                elsif ( $type==$MODE_DELETION )
                {
                    $col = $white;
                    $im->filledRectangle($x1-1,$y+2,$x1+gdSmallFont->width-1,$y+$charheight+1,$red);
				} 
                elsif ( $type==$MODE_INSERTION )
                {
                    $col = $white;
                    $im->filledRectangle($x1-1,$y+2,$x1+gdSmallFont->width-1,$y+$charheight+1,$red);
                }
				elsif ( $type == $MODE_MAPQUAL_ZERO ) 
                {
					$col = ( $ch eq lc $ch ) ? $col_mapqual_zero : $red ;
				} 
				else 
                {
					$col = ( $ch eq lc $ch ) ? $blue : $red ;
				}
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
        if ( $nob>$max_pileup_rows )
        {
            my $msg = "Output truncated, too many sequences in this view.";
            my $w = length($msg)*gdMediumBoldFont->width;
            my $h = gdMediumBoldFont->height;
            my $x = ($width - $w)*0.5;
            my $y = 1;
            $im->filledRectangle($x,$y,$x+$w+1,$y+$h+1,$red);
            $im->string(gdMediumBoldFont, $x,$y, $msg,$white);
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



# The returned position is relative to $from, starting from 0.
# Only the part of the sequence visible in the current window is returned.
#
sub trim_sam_reads
{
    my ($seq,$cigar,$mode,$start) = @_;

    my $seq_len = length($seq);
    my $ref_len = $to - $from + 1;

    # If not visible, return immediately
    if ( $start > $to ) { return ($to-$from+1,'',''); }
    if ( $start + $seq_len <= $from ) { return (0,'',''); }

    my ($out_offset,$out_seq,$out_btype);

    my $ninserts = 0;
    my $ndels = 0;
    my $iseq  = 0;
    while ($cigar)
    {
        if ( !($cigar=~/^(\d+)(\D)/) ) { last }

        my $count = $1;
        my $type  = $2;
        $cigar    = $';

        if ( $type eq 'M' )
        {
            my $iref = $start + $iseq - $ninserts + $ndels - $from;
            for (my $i=0; $i<$count; $i++)
            {
                if ( $iref >= $ref_len ) { last; }

                my $snp = substr($seq,$iseq,1);
                my $ref = substr($refseq,$iref,1);
                
                if ( $snp eq $ref )
                {
                    $out_seq   .= lc($snp);
                    $out_btype .= chr($mode + 65);
                }
                else
                {
                    $out_seq   .= uc($snp);
                    $out_btype .= chr($mode + 65);
                }

                $iref++;
                $iseq++;
            }
            if ( $iref> $ref_len ) { last; }
        }
        elsif ( $type eq 'D' )
        {
            $out_seq   .= '*' x $count;
            $out_btype .= chr($MODE_DELETION + 65) x $count;
            $ndels += $count;
        }
        elsif ( $type eq 'I' )
        {
            if ( $iseq+$start<$from )
            {
                # In case the insertion is not visible (hidden to the left), do not display it.
                $iseq += $count;
                $ninserts += $count;
            }
            else
            {
                $out_seq   .= uc substr($seq,$iseq,$count);
                $out_btype .= chr($MODE_INSERTION + 65) x $count;
                $iseq += $count;
                $ninserts += $count;
            }
        }
        elsif ( $type eq 'H' ) 
        {
            # Ignore hard clips
            next;
        }
        elsif ( $type eq 'S' ) 
        {
            # Ignore soft clips
            $iseq  += $count;
            next;
        }
        else { die "Could not parse the cigar $seq .. $cigar.\n" }
    }
    if ( $start < $from )
    {
        substr($out_seq, 0, $from-$start,'');
        substr($out_btype, 0, $from-$start,'');
        $start = $from;
    }
    if ( $start + length($out_seq) - 1 > $to ) 
    {
        $out_seq = substr($out_seq,0,$to-$start+1);
    }
    return ($start-$from,$out_seq,$out_btype);
}


sub pile2bands_sam 
{
	my ( $r_bands , $r_btype , $mode , $paired_pileup , $data ) = @_ ;

	foreach ( @{$data} ) 
    {
		my $r = $sam_reads{$_};
        
        my ($start,$seq,$btype) = trim_sam_reads($r->[0]->[8],$r->[0]->[4],$mode,$r->[0]->[2]);
		my $seq_len = length $seq;

        if ( scalar @$r == 1 )
        {
            sam_add_pileup($start, $seq, $seq_len, $r_bands, $r_btype, $btype);
            next;
        }

        my ($start2,$seq2,$btype2) = trim_sam_reads($r->[1]->[8],$r->[1]->[4],$mode,$r->[1]->[2]);
		my $seq_len2 = length $seq2;

        if ( !$seq_len && !$seq_len2 ) { next; }

        my $between = $start2 - $start - $seq_len;

        if ( $between<0 )
        {
            # The sequences overlap - the overlapping positions will be overwritten
            #   by $seq2
            substr($seq,$between,-$between) = substr($seq2,0,-$between);
            substr($seq2,0,-$between,'');

            substr($btype,$between,-$between) = substr($btype2,0,-$between);
            substr($btype2,0,-$between,'');
            $between = 0;
        }

        $seq .= '_' x $between ;
        $seq .= $seq2;

        $btype .= chr ($mode + 65) x $between ;
        $btype .= $btype2;

        sam_add_pileup ($start, $seq , length $seq, $r_bands , $r_btype , $btype);
	}
}


sub sam_add_pileup 
{
	my ($seq_start, $seq, $seq_len, $rbands, $rbtype, $btype) = @_ ;

    my $ref_len = $to - $from + 1;

	my $use_band ;
	my $limit = scalar ( @{$rbands} ) - 1 ;

    if ( $limit>$max_pileup_rows ) { return; }

	foreach my $band ( 0 .. $limit ) 
    {
        my $str;
        if ( $seq_start>0 ) 
        { 
            $str = substr($rbands->[$band], $seq_start-1, $seq_len+2); 
        }
        else
        {
            $str = substr($rbands->[$band], $seq_start, $seq_len+1); 
        }

        if ( $str=~/^ +$/ ) 
        { 
		    $use_band = $band;
		    last;
        }
	}
	
	unless ( defined $use_band ) {
		$use_band = scalar ( @{$rbands} ) ;
		push @{$rbands} , ( ' ' x $ref_len ) ;
		push @{$rbtype} , ( ' ' x $ref_len ) ;
	}
	
	substr ( $rbands->[$use_band] , $seq_start, $seq_len ) = $seq ;
	substr ( $rbtype->[$use_band] , $seq_start, $seq_len ) = $btype;
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
	my ( $start , $read_sequence , $rl , $rbands , $rbtype , $mode , $dist , $sam ) = @_ ;
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
			if ( defined $sam ) {
				my $ch = substr ( $read_sequence , $_ , 1 )  ;
				if ( $ch eq '_' or $ch eq lc substr ( $refseq , $onref , 1 ) ) {
					$out .= $ch ;
				} else {
					$out .= uc substr ( $read_sequence , $_ , 1 )  ;
				}
			} else {
				$out .= substr ( $read_sequence , $_ , 1 )  ;
			}
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
	
	if ( $mode != 1 and $mode != 3 and not defined $sam ) {
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

sub add_coverage_sam_single 
{
	my ($coverage, $from, $to, $r) = @_ ;

	my $start = $r->[2];
	my $end   = $start + length $r->[8];
	if ( $end > $to ) { $end = $to; }

	my $ifrom = $start - $from;
	if ( $ifrom < 0 ) { $ifrom = 0; }

	my $ito = $end - $from;
	if ( $ito > $to-$from ) { $ito=$to-$from; }

	for (my $i=$ifrom; $i<=$ito; $i++) 
    {
		$$coverage[$i]++;
	}
}

sub add_coverage_sam 
{
	my ($coverage, $from, $to, $data) = @_;
	foreach ( @{$data} ) 
    {
		my $d = $sam_reads{$_};
		add_coverage_sam_single($coverage, $from, $to, $d->[0]);
		add_coverage_sam_single($coverage, $from, $to, $d->[1]) unless !$d->[1];
	}
}


sub dump_image_coverageview 
{
	if ( !$using_bam ) 
    {
        die "FIXME: the old code not included\n";
    }

    my (@coverage,@coverage_tot);
    for my $i (0 .. ($to-$from))
    {
        $coverage[$i] = 0;
    }
	
    add_coverage_sam(\@coverage, $from, $to, \@sam_perfect);
    add_coverage_sam(\@coverage, $from, $to, \@sam_perfect_singles);
    add_coverage_sam(\@coverage, $from, $to, \@sam_snps);
    add_coverage_sam(\@coverage, $from, $to, \@sam_inversions);
    add_coverage_sam(\@coverage, $from, $to, \@sam_single);

    for my $i (0 .. ($to-$from))
    {
        $coverage_tot[$i] = $coverage[$i];
    }
    add_coverage_sam(\@coverage_tot, $from, $to, \@sam_mapqual_zero);
    add_coverage_sam(\@coverage_tot, $from, $to, \@sam_mapqual_zero_pair);

    my $plot  = Plot->new({width=>$width,height=>80,fgcolor=>$lscolor{'mapqual_zero'},bgcolor=>$lscolor{'white'}});
    $plot->scaled_polygon(\@coverage_tot);
    $plot->set({fgcolor=>$lscolor{'blue'}});
    $plot->scaled_polygon(\@coverage);
    $plot->set({fgcolor=>$lscolor{'black'}});
    $plot->draw_yscalebar(\@coverage);
    write_png($$plot{image});

    return;
}


sub dump_image_coverageview_ori {
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
	
	if ( $using_bam ) {
		add_coverage_sam ( $from , $to , \@sam_perfect ) ;
		add_coverage_sam ( $from , $to , \@sam_snps ) ;
		add_coverage_sam ( $from , $to , \@sam_inversions ) ;
		add_coverage_sam ( $from , $to , \@sam_single , 1 ) ;
	}

	$height = 0 ;
	foreach ( @coverage ) {
		$height = $_ if $_ > $height ;
	}

	my $second = $cgi->param ( 'second' ) ;
	my @mpp2 ;
	if ( $second ) {
		my $u = $cgi->url() ;
		my @p = $cgi->param ;
		my @p2 ;
		foreach my $k ( @p ) {
			next if $k eq 'lane' ;
			my $v = $cgi->param($k) ;
			if ( $k eq 'second' ) {
				$k = 'lane' ;
				$v = $cgi->param('second') ;
			}
			$v = 'text' if $k eq 'output' ;
			push @p2 , "$k=$v" ;
		}
		$u .= "?" . join ( '&' , @p2 ) ;
		my $c = get ( $u ) ;
		@mpp2 = split "\n" , $c ;
		foreach ( @mpp2 ) {
			$height = $_ if $_ > $height ;
		}
		
#		print $cgi->header(-type=>'text/plain',-expires=>'-1s');
#		print $u ;
#		exit ;
	}

	$height += 1 unless $height ;
	$height += $scale_height ; # Scale
	
	my $im = new GD::Image ( $width , $height ) ;
	my $white = $im->colorAllocate ( @{$lscolor{'white'}} ) ;
	my $black = $im->colorAllocate ( @{$lscolor{'black'}} ) ;
	
	my ( $col , $col2 , $col3 ) ;
	if ( $second ) {
		$col = $im->colorAllocate ( @{$lscolor{'ltblue'}} ) ;
		$col2 = $im->colorAllocate ( @{$lscolor{'ltgreen'}} ) ;
		$col3 = $im->colorAllocate ( @{$lscolor{'yellow'}} ) ;
	} else {	
		$col = $im->colorAllocate ( @{$lscolor{'blue'}}  ) ;
	}

	my @max_per_pixel ;
	$max_per_pixel[$_] = 0 foreach ( 0 .. $width ) ;
	foreach my $pos ( 0 .. $ft-1 ) {
		my $x1 = int ( $width * $pos / $ft ) ;
		$max_per_pixel[$x1] = $coverage[$pos] if $max_per_pixel[$x1] < $coverage[$pos] ;
	}

	# Return data as text if requested
	if ( $output eq 'text' ) {
		print $cgi->header(-type=>'text/plain',-expires=>'-1s');
		if ( $ft >= $width ) {
			print join "\n" , @max_per_pixel ;
		} else {
			print join "\n" , @coverage ;
		}
		return ;
	}
	
	if ( $ft >= $width ) {
		foreach ( 0 .. $width-1 ) {
			
#			$im->line ( $_ , $y , $_ , $height - $scale_height , $col ) ;
			if ( $second ) {
				my $y1 = $height - $scale_height - $max_per_pixel[$_] ;
				my $y2 = $height - $scale_height - $mpp2[$_] ;
				$im->line ( $_ , $y2 , $_ , $height - $scale_height , $col2 ) ;
				$im->line ( $_ , $y1 , $_ , $height - $scale_height , $col ) ;
				$y2 = $y1 if $y1 > $y2 ;
				$im->line ( $_ , $y2 , $_ , $height - $scale_height , $col3 ) ;
			} else {
				my $y = $height - $scale_height - $max_per_pixel[$_] ;
				$im->line ( $_ , $y , $_ , $height - $scale_height , $col ) ;
			}
		}
	} else {
		foreach my $pos ( 0 .. $ft-1 ) {
			my $x1 = int ( $width * $pos / $ft ) ;
#			next if $max_per_pixel[$x1] > $coverage[$pos] ;
			my $x2 = int ( $width * ( $pos + 1 ) / $ft - 1 ) ;
			my $y = $height - $scale_height - $coverage[$pos] ;
			if ( $second ) {
				my $y2 = $height - $scale_height - $mpp2[$pos] ;
				if ( $x2 > $x1 ) {
					$im->filledRectangle ( $x1 , $y2 , $x2 , $height - $scale_height , $col2 ) ;
				} else {
					$im->line ( $x1 , $y2 , $x1 , $height - $scale_height , $col2 ) ;
				}
			}
			if ( $x2 > $x1 ) {
				$im->filledRectangle ( $x1 , $y , $x2 , $height - $scale_height , $col ) ;
			} else {
				$im->line ( $x1 , $y , $x1 , $height - $scale_height , $col ) ;
			}
			if ( $second ) {
				my $y2 = $height - $scale_height - $mpp2[$pos] ;
				$y2 = $y if $y2 < $y ;
				if ( $x2 > $x1 ) {
					$im->filledRectangle ( $x1 , $y2 , $x2 , $height - $scale_height , $col3 ) ;
				} else {
					$im->line ( $x1 , $y2 , $x1 , $height - $scale_height , $col3 ) ;
				}
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


# The paired_reads display splits the screen into two parts, where one positions
#   the reads vertically according to their insert size, while in the other are
#   the reads positiioned randomly.
sub sam_n_zero_isize
{
    return scalar @sam_single + scalar @sam_mapqual_zero;
}
sub sam_n_nonzero_isize
{
    return scalar @sam_perfect + scalar @sam_inversions + scalar @sam_mapqual_zero_pair + scalar @sam_perfect_singles;
}


sub dump_image_indelview {
	$ft = $to - $from + 1 ;
	my $text_mode = $width / $ft > 5 ;
	$max_dist = 1 ; # Dummy value

	$im = new GD::Image ( $width , $height ) ;
	my $white = $im->colorAllocate ( @{$lscolor{'white'}} ) ;
	my $black = $im->colorAllocate ( @{$lscolor{'black'}} ) ;
	my $blue = $im->colorAllocate ( @{$lscolor{'blue'}} ) ;
	my $light_blue = $im->colorAllocate ( @{$lscolor{'ltblue'}} ) ;
	my $red = $im->colorAllocate ( @{$lscolor{'red'}} ) ;
	my $single_color = $im->colorAllocate ( @{$lscolor{'green2'}} ) ;
	my $inversion_color = $im->colorAllocate ( @{$lscolor{'brown'}} ) ;
	my $variance_color = $im->colorAllocate ( 0xFF , 0xF6 , 0x8F ) ; # FFF68F
	my $inversion_left_color = $im->colorAllocate ( 0xB6 , 0xBA , 0x18 ) ; # B6BA18
	my $inversion_right_color = $im->colorAllocate ( 0x7C , 0xEB , 0x98 ) ; # 7CEB98
	my $inversion_middle_color = $im->colorAllocate ( 0x99 , 0xD2 , 0x58 ) ; # 99D258
	my $orange = $im->colorAllocate ( 0xFF , 0xAA , 0x00 ) ;
	my $grey = $im->colorAllocate ( @{$lscolor{'grey'}} ) ;
    my $col_mapqual_zero = $im->colorAllocate( @{$lscolor{'mapqual_zero'}} );
    $col_insertion = $im->colorAllocate( 0x66, 0xdd, 0x66 );
    $col_deletion  = $im->colorAllocate( 0xdd, 0x66, 0x66  );
	
	my @ann_color = ( $black , $blue , $red , $single_color , $orange , $inversion_right_color, $inversion_middle_color , $variance_color  , $inversion_left_color , $inversion_color ) ;
	
	$max_dist = $cgi->param('maxdist') ;
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

	if ( $using_bam ) {
		$sam_white = $white ;
		$sam_black = $black ;
		$sam_col_single_read = $single_color ;
		$sam_col_mismatch = $red ;
		$sam_col_matching_read = $blue ;
		$sam_col_inversion = $inversion_color ;
		$sam_col_read_pair_connection = $grey ;
		$sam_col_read_pair_quality = $orange ;
		$sam_col_single_read_quality = $orange ;
		$sam_col_mapqual_zero = $col_mapqual_zero;
        # Uh, what is this?
        #
		#   if ( $max_dist <= 1 ) {
		#   	$max_dist = $sam_max_found_fragment ;
		#   	$max_dist = $height - $scale_height if $max_dist < $height - $scale_height ;
		#   	$max_dist = $height ; # AAARGH DUMMY FIXME !!!!!
		#   }
        #&sam_paint ;
        my $nnonzeros = sam_n_nonzero_isize();
        my $nzeros    = sam_n_zero_isize();
        my $margin    = gdSmallFont->height*2;
        $height_nonzero_isize  = ($nnonzeros+$nzeros) ? int(($height-$margin)*$nnonzeros/($nnonzeros+$nzeros)) : 1;
        if ( !$height_nonzero_isize ) { $height_nonzero_isize=1; }    # to prevent division by 0
		$height_nonzero_isize = $height - 20 unless $display_single ;
        $height_zero_isize = $height - $margin - $height_nonzero_isize;

        $max_dist = $max_insert_size;
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
		if ($using_bam ) {
			sam_paint_single_short_reads ( \@sam_single , $sam_col_single_read ) ;
			sam_paint_single_short_reads ( \@sam_mapqual_zero, $sam_col_mapqual_zero ) ;
			sam_paint_short_single_reads_quality ( \@sam_single , $sam_col_single_read_quality ) if $sam_show_read_quality ;
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
		if ( $using_bam ) {
			sam_paint_short_read_pair_connections ( \@sam_perfect , $sam_col_read_pair_connection ) ;
			sam_paint_short_read_pair_connections ( \@sam_snps , $sam_col_read_pair_connection ) ;
			sam_paint_short_read_pair_connections ( \@sam_inversions , $sam_col_read_pair_connection ) ;
			sam_paint_short_read_pair_connections ( \@sam_mapqual_zero_pair, $sam_col_read_pair_connection ) ;
		}
	}

	# Read quality (BAM only)
	if ( $using_bam and $sam_show_read_quality ) {
		sam_paint_short_read_pairs_quality ( \@sam_perfect , $sam_col_read_pair_quality ) ;
		sam_paint_short_read_pairs_quality ( \@sam_snps , $sam_col_read_pair_quality ) ;
		sam_paint_short_read_pairs_quality ( \@sam_inversions , $sam_col_read_pair_quality ) ;
		sam_paint_short_read_pairs_quality ( \@sam_mapqual_zero_pair, $sam_col_mapqual_zero ) ;
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
	if ( $using_bam ) {
		sam_paint_short_read_pairs ( \@sam_perfect , $sam_col_matching_read , 1 ) ;
		sam_paint_short_read_pairs ( \@sam_snps , $sam_col_matching_read , 1 ) ;
		sam_paint_short_read_pairs ( \@sam_inversions , $sam_col_inversion , 1 ) ;
		sam_paint_short_read_pairs ( \@sam_mapqual_zero_pair, $sam_col_mapqual_zero, 0 ) ;
        sam_paint_short_read_pairs ( \@sam_perfect_singles, $sam_col_matching_read, 0 ) ;
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
	my $diff = 100;
	$diff = 500 if $max_dist > 2500 ;
	$diff = 2500 if $max_dist > 25000 ;

    my ($w,$h) = (gdSmallFont->width,gdSmallFont->height);
    if ( $using_bam )
    {
        my $min_diff = $max_dist*$h/$height_nonzero_isize;
        $diff = 1;
        my $i=0;
        while ( $diff < 2*$min_diff )
        {
            $diff *= $i%2 ? 2 : 5;
            $i++;
        }
    }

    # Why to make so many unnecessary iterations?
    #
	#   foreach ( 1 .. $max_dist ) {
	#   	next unless $_ % $diff == 0 ;
	#   	my $y = $height - $_ * $height / $max_dist ;
    #       push @positions,$y;
	#   	$im->line ( 0 , $y , 10 , $y , $black ) ;
	#   	$im->string ( gdSmallFont , 13 , $y > 5 ? $y - 6 : -3 , $_ , $black ) ;
	#   }
    #
    # .. do this instead:
    #
    # Area covered by the y-tick label
    my @positions = (); # Remember the positions of ticks, so that we can draw the label
    my $pos = 0;
    while ( $pos<=$max_dist )
    {
        my $y = $using_bam ? sam_get_y($pos) : $height - $pos * $height / $max_dist;    
        $im->line ( 0 , $y , 10 , $y , $black ) ;
        $im->string ( gdSmallFont , 13 , $y > 5 ? $y - 6 : -3 , $pos , $black ) ;
        push @positions,$y;
        $pos += $diff;
    }

    # Too bad, there are no margins in the image for the axis labels - place the label
    #   in the middle between two ticks, close to the half of the image.
    #
    my $label = 'Insert Size';
    $w *= length $label;
    my $label_y = $height*0.5;
    if ( @positions ) 
    {
        my $i  = int((scalar @positions)*0.5);
        my $i1 = $i+1;
        if ( $i1 >= scalar @positions ) { $i1=$i; }
        $label_y = ($positions[$i1] + $positions[$i] + $w)*0.5;
    }
    $im->filledRectangle(0,$label_y-$w-1,$h-2,$label_y-2, $white);
    $im->stringUp(gdSmallFont, -2, $label_y, $label, $black);

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
			$middle = $smaller + int ( ( $ofs - $rl + $avg_frag_size ) / 2 )  ;
			$int_var = abs ( $avg_frag_size - 2 * $rl ) ;
		}  elsif ( $dir[1] eq '-' ) {
			$xt = int ( ( $smaller - $from + $ofs ) * $width / $ft ) ;
			$xf = int ( ( $smaller - $from + $ofs - $avg_frag_size ) * $width / $ft + $len ) ;
			$middle = $smaller + int ( ( $ofs - $rl + $avg_frag_size ) / 2 )  ;
			$int_var = abs ( $avg_frag_size - 2 * $rl ) ;
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
	my $white = $im->colorAllocate ( @{$lscolor{'white'}} ) ;
	my $black = $im->colorAllocate ( @{$lscolor{'black'}} ) ;
	my $blue = $im->colorAllocate ( @{$lscolor{'blue'}} ) ;
	my $red = $im->colorAllocate ( @{$lscolor{'red'}} ) ;
	my $green = $im->colorAllocate ( @{$lscolor{'green'}} ) ;
	my $gray = $im->colorAllocate ( @{$lscolor{'grey'}} ) ;

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


sub dump_image_deletions {
	$height = 50 ;
	my $im = new GD::Image ( $width , $height ) ;
	my $white = $im->colorAllocate ( @{$lscolor{'white'}} ) ;
	my $black = $im->colorAllocate ( @{$lscolor{'black'}} ) ;
	my $blue = $im->colorAllocate ( @{$lscolor{'ltblue'}} ) ;
	my $yellow = $im->colorAllocate ( @{$lscolor{'yellow'}} ) ;
	my $green = $im->colorAllocate ( @{$lscolor{'ltgreen'}} ) ;
	my $grey = $im->colorAllocate ( @{$lscolor{'grey'}} ) ;
	
	my $ft = $to - $from + 1 ;
	my $max_del_width = 10000 ;
	my ( @left , @right , @center ) ;
	
	if ( $number_of_databases == 0 ) {
		write_png ( $im ) ;
		return ;
	}
	
	foreach my $current_db ( 0 .. $number_of_databases - 1 ) {
		my %meta = %{$all_meta[$current_db]} ;
		next if ( $meta{'fragment'} || 0 ) < 1 ;
		my $upper = $meta{'fragment'}+$meta{'variance'} ;
		my $lower = $meta{'fragment'}-$meta{'variance'} ;
		my $fs_expected = $meta{'fragment'} ;
		my $rl = $meta{'read_length'} ;
#		$avg_frag_size = $meta{'fragment'} || 0 ;
#		next if $avg_frag_size == 0 ;
		foreach ( @{$all_pmr[$current_db]} , @{$all_smr[$current_db]} ) {
			my ( $rl1 , $rl2 ) = ( $_->[scalar(@{$_})-2] , $_->[scalar(@{$_})-1] ) ;
			my $p1 = $_->[1] ;
			my $p2 = $_->[2] ;
			my $ofs = $p2 - $p1 + $rl2 ; # Observed fragment size
			next if $ofs <= $upper + 20 ;#and $ofs >= $lower ; # INSERTIONS MISSING!!
			next if $ofs > $max_del_width ;
			my $est_var = int ( ( $ofs - $fs_expected ) / 2 ) ;
			my $c = int ( $p1 + $ofs/2 ) ;
			foreach ( $c - $est_var .. $c + $est_var ) {
				$center[$_]++ ;
			}
			foreach ( $p1 .. $p1 + $rl1 ) {
				$left[$_]++ ;
			}
			foreach ( $p2 .. $p2 + $rl2 ) {
				$right[$_]++ ;
			}
#			$center[int(($p1+$p2+$rl2)/2)]++ ;
		}
	}
	
	my $max = 0 ;
	foreach my $pos ( $from ... $to ) {
		$max = $center[$pos] if ( ($center[$pos]||0)) > $max ;
	}
	if ( $max == 0 ) {
		write_png ( $im ) ;
		return ;
	}
	my $factor = ($height-1) / $max ;
#	$factor = 1 ;
	foreach my $pos ( $from ... $to ) {
		my $x = int ( ( $pos - $from ) * $width / $ft ) ;
		my $y = $height - 1 ;
		$im->line ( $x , $y , $x , $y - $center[$pos]*$factor , $green ) if $center[$pos] ;
		$im->line ( $x , $y , $x , $y - $left[$pos]*$factor , $blue ) if $left[$pos] ;
		$im->line ( $x , $y , $x , $y - $right[$pos]*$factor , $yellow ) if $right[$pos] ;
	}
	
	my $legend_color = $black ;
	my @steps = ( 1000 , 500 , 250 , 100 , 50 , 25 , 10 , 5 ) ;
	shift @steps while 1 < scalar @steps and $steps[0] * 5 > $max ;
	for ( my $a = $steps[0] ; $a <= $max ; $a += $steps[0] ) {
		my $y = $height - $a * $factor ;
		$im->line ( 0 , $y , 5 , $y , $legend_color ) ;
		$im->string ( gdSmallFont , 10 , $y-5 , $a , $legend_color ) ;
	}
	$im->string ( gdSmallFont , 40 , 0 , "[reads]" , $legend_color ) ;

	write_png ( $im ) ;
}

sub dump_image_gc
{
    my $win   = 5;
    my $from2 = $from - $win;
    my $to2   = $to   + $win;
    if ( $from2 < 1 ) { $from2 = 1; }
	my $refseq = get_chromosome_part($genome_file, $orig_chr, $from2, $to2);

    my @gc_data = ();
    my $len = length($refseq);
    my $cnt = 0;
    for (my $i=0; $i<2*$win; $i++) 
    { 
        my $base = substr($refseq,$i,1);
        if ( $base eq 'C' or $base eq 'G' ) { $cnt++; }
    }
    for (my $i=$win; $i<$len-$win; $i++)
    {
        my $old = substr($refseq,$i-$win,1);
        my $new = substr($refseq,$i+$win,1);
        if ( $new eq 'C' or $new eq 'G' ) { $cnt++; }
        $gc_data[$i-$win] = 100.*$cnt/(2.*$win+1);
        if ( $old eq 'C' or $old eq 'G' ) { $cnt--; }

        #printf STDERR $i-$win." .. %.1f .. ".substr($refseq,$i,1)."  [$old $new]\n",$gc_data[$i-$win];
    }

    my $plot  = Plot->new({width=>$width,height=>60,ymin=>0,ymax=>100,baseline=>50,fgcolor=>$lscolor{'green2'},bgcolor=>$lscolor{'white'},units=>'%'});
    $plot->scaled_line(\@gc_data);
    $plot->set({fgcolor=>$lscolor{'mapqual_zero'}});
    $plot->baseline();
    $plot->set({fgcolor=>$lscolor{'black'}});
    $plot->draw_yscalebar(\@gc_data);
    $plot->legend('[win '. ($win*2+1) . ']');
    write_png($$plot{image});
}


sub error_exit
{
    my (@msg) = @_;
    my $msg = join('',@msg);
    my $plot = Plot->new({width=>$width,height=>60,fgcolor=>$lscolor{'black'},bgcolor=>$lscolor{'white'}});
    my $img  = $$plot{image};
    $img->string(gdLargeFont, 0, 0, $msg, $$plot{fgcolor});
    write_png($$plot{image});
    exit;
}


sub dump_image_gc_ori {
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
	my $white = $im->colorAllocate ( @{$lscolor{'white'}} ) ;
	my $black = $im->colorAllocate ( @{$lscolor{'black'}} ) ;
	my $blue = $im->colorAllocate ( @{$lscolor{'blue'}} ) ;
	my $red = $im->colorAllocate ( @{$lscolor{'red'}} ) ;

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
	print $cgi->header(-type=>'image/png',-expires=>'-1s') unless $DEBUG;
	binmode STDOUT;
	my $png = $im->png () ;
	print $png ;
}



#______________________________________________________________________________________


# METHODS

# Paint
sub sam_paint {
#	$| = 1;	print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output
#	print "So far...\n" ;
	if ( $display_pair_links ) {
		sam_paint_short_read_pair_connections ( \@sam_perfect , $sam_col_read_pair_connection ) ;
		sam_paint_short_read_pair_connections ( \@sam_snps , $sam_col_read_pair_connection ) ;
		sam_paint_short_read_pair_connections ( \@sam_inversions , $sam_col_read_pair_connection ) ;
	}

	if ( $sam_show_read_quality ) {
		sam_paint_short_read_pairs_quality ( \@sam_perfect , $sam_col_read_pair_quality ) ;
		sam_paint_short_read_pairs_quality ( \@sam_snps , $sam_col_read_pair_quality ) ;
		sam_paint_short_read_pairs_quality ( \@sam_inversions , $sam_col_read_pair_quality ) ;
		sam_paint_short_single_reads_quality ( \@sam_single , $sam_col_single_read_quality ) if $display_single ;
	}

	sam_paint_short_read_pairs ( \@sam_perfect , $sam_col_matching_read , 0 ) ;
	sam_paint_short_read_pairs ( \@sam_snps , $sam_col_matching_read , 1 ) ;
	sam_paint_short_read_pairs ( \@sam_inversions , $sam_col_inversion , 1 ) ;
	sam_paint_single_short_reads ( \@sam_single , $sam_col_single_read ) if $display_single ;
	sam_paint_single_short_reads ( \@sam_mapqual_zero, $sam_col_mapqual_zero ) if $display_single ;
#	print "DONE!" ; exit ;
}

# Read data in
sub sam_read_data {
	my ( $file_bam ) = @_ ;
	
    my $delta = $max_insert_size ? $max_insert_size : 550;
	my $from2 = $from - $delta;
	my $to2   = $to   + $delta;
	$from2 = 1 if $from2 < 1 ;
	
	$chromosome =~ /([a-zA-Z0-9_.]+)/ ;
	my $cmd = "cd $datapath ; $cmd_samtools view $file_bam $1:$from2-$to2" ;
	
	# De-tainting $ENV{'PATH'}
	$ENV{'PATH'} =~ /(.*)/;
	$ENV{'PATH'} = $1;
	
#	print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output
    my $max = 0;
    my %isize_hist  = ();
	my $nlines=0;
    my $max_nlines = 5000000;
	open PIPE , "$cmd 2>&1 |" or die "$cmd |: $!";
	while ( <PIPE> ) 
    {
		next if $_ =~ m/^\[/ ;
		$_ =~ /^(\S+)\s/ ;
		my $id = $1;
		my @a = split "\t" , $';

        if ( 0x0004 & $a[0] ) { next; } # This read is unmapped
        if ( $mapq_cutoff && $a[3]<$mapq_cutoff ) { next }   # This read has lower MAPQ
		push @{$sam_reads{$id}} , \@a;

		$nlines++;
        if ( $nlines>$max_nlines ) 
        { 
            error_exit("Sorry, there are too many (>$max_nlines) reads in this region."); 
        }

        if ( $a[4]=~/^(\d+)S/ )
        {
            $a[4] = $';
            $a[8] = substr($a[8],$1);
        }

        if ( $max_insert_size ) { next }
        my $isize = abs($a[7]);
        if ( !$isize ) { next; }

        # To determine the maximum insert size, we must first know if the read 
        #   will be actually visible in the window. 
        #
        my $pos_from = $a[2];
        my $pos_to   = $pos_from+length($a[8])-1;
        if ( $pos_from>$to ) { next }
        if ( $pos_to<$from ) { next }

        if ( $isize>$max ) { $max=$isize; }
        $isize_hist{$isize}++;
    }
	close PIPE ;
	if ( defined $ls_max_lines and $nlines > $ls_max_lines ) { $display_snps = 0; }

    if ( !$max_insert_size ) { $max_insert_size= $max ? 1.1*$max : 500; }
    if ( $nlines>10000 && $view eq 'paired_pileup' )
    {
        $view = 'indel';
        #print STDERR "get_data.pl: too many lines, switching to indelview ..\n";
    }
	
	if ( $nlines>2000 && $view eq 'indel' ) { # Do not show quality for too many reads
		$sam_show_read_quality = 0 ;
	}
	
    #print STDERR "get_data.pl: nlines=$nlines\n";
}


# Separate into bins
sub sam_bin 
{
    foreach my $read ( keys %sam_reads )
    {
        my $r = $sam_reads{$read};
        my $nreads = scalar @$r;

        if ( $nreads==2 and $r->[0]->[2] == $r->[1]->[2] and $r->[0]->[1] == $r->[1]->[1])
        {
            # There are two reads, but they are mapped to the same position.
            # Can this happen with unmapped reads filtered out? Just delete
            #   the redundant one.
            delete($$r[1]);
            $nreads = 1;
        }

        if ( $nreads==1 )
        {
            if ( $r->[0]->[0] & 0x0008 )
            {
                # The mate is unmapped
                push @sam_single, $read;
                next;
            }

            if ( !$r->[0]->[3] )
            {
                # The mapping quality is zero - marks non-unique mapping
                if ( $r->[0]->[7] )
                {
                    # Although we do not have the mate, the read has a non-zero insert size.
                    push @sam_mapqual_zero_pair, $read;
                }
                else
                {
                    push @sam_mapqual_zero, $read;
                }
                next;
            }

            # These don't need to be exactly "perfect". For example, zero-insert
            #   size reads with mate mapped to a different chromosome fall in this
            #   category as well.
            push @sam_perfect_singles, $read;
            next;
        }

        if ( ( $r->[0]->[0] & 0x0010 ) == ( $r->[1]->[0] & 0x0010 ) )
        {
            # Both reads are mapped to the same strand - inversion
            push @sam_inversions, $read;
            next;
        }

        if ( !$r->[0]->[3] )
        {
            # The mapping quality is zero - marks non-unique mapping
            push @sam_mapqual_zero_pair, $read;
            next;
        }

        push @sam_perfect, $read;
    }
}



sub sam_check4snps {
#    return 0;   # this is treated properly in trim_sam_reads
	if ( !$display_snps ) { return 0; }

	my ($r) = @_;
	my $seq_start = $r->[2];
	my $seq       = $r->[8];
	my $seq_len   = length($seq);
	my $seq_end   = $seq_start + $seq_len - 1;

	my $iseq_offset = 0;
	my $iref_offset = 0;
	if ( $seq_start < $from ) { $iseq_offset = $from - $seq_start; }
	elsif ( $seq_start > $from ) { $iref_offset = $seq_start - $from; }

	my $ref_len = $to - $from + 1;

	my $iseq = $iseq_offset;
	my $iref = $iref_offset;
	
#	print substr ( $refseq , $iref , $seq_len ) . "1\n" ;
#	print substr ( $seq , $iseq , $seq_len ) . "2\n\n" ;

	my $mismatches = 0 ;
	while ( $iseq<$seq_len && $iref<$ref_len ) {
		my $ref = substr($refseq,$iref,1);
		my $snp = substr($seq,$iseq,1);

		if ( $ref ne $snp ) {
			substr($r->[8], $iseq, 1) = lc($snp);
			$mismatches++;
		}

		$iseq++;
		$iref++;
	}

	return $mismatches ;
}

=cut
sub sam_check4snps_old {
	my ( $r ) = @_ ;
	my $seq = $r->[8] ;
	my $pos = $r->[2] ;
	if ( $pos >= $from ) {
		if ( $seq eq substr $refseq , $pos - $from , length $seq ) {
			return 0 ;
		}
	}

	my $mismatches = 0 ;

	my $rpos = $pos - $from - 1 ;
	my $start = 0 ;
	if ( $rpos < -1 ) {
		$start = -$rpos ;
		$rpos = 0 ;
	}
	my $end = length ( $seq ) - 1 ;
	$end-- while ( $rpos + $end - $start >= $reflength ) ;
	foreach my $spos ( $start .. $end ) {
		$rpos++ ;
		if ( substr ( $seq , $spos , 1 )  ne substr ( $refseq , $rpos , 1 ) ) {
			substr ( $r->[8] , $spos , 1 ) = lc substr ( $seq , $spos , 1 )  ;
			$mismatches++ ;
		}
	}

	return $mismatches ;
}
=cut

sub sam_get_y {
	my ( $r ) = @_ ;

    my $y = !ref($r) ? $r : $$r[7];
    return $height_nonzero_isize - $height_nonzero_isize * abs($y)/$max_insert_size;

	#   my $y = abs ( $r->[7] ) ;#- length $r->[8] ;
	#   $y = ( $height  ) * $y / $max_dist ;
	#   #$y = $height - $y - $scale_height ;
	#   $y = $height - $y ; #ori
	#   return $y ;
}

sub sam_paint_short_read_pair_connections {
	my ( $r , $col ) = @_ ;
	foreach ( @{$r} ) {
		my $read = $sam_reads{$_} ;
		my $y = sam_get_y ( $read->[0] ) ;
#		my $y = abs ( $read->[0]->[7] ) ;
#		$y = ( $height - $scale_height ) * $y / $max_dist ;
#		$y = $height - $scale_height - $y ;
		my $x1 = $read->[0]->[2] - $from ;
		my $x2 = $read->[1]->[2] - $from ;
		$x1 = int ( $x1 * $width / $ft ) ;
		$x2 = int ( $x2 * $width / $ft ) ;
#		print "$x1/$y -> $x2/$y | $height | $scale_height | $max_dist\n" ;
        #if ( !($x1==289 && $x2==957) ) { next; }
		$im->line ( $x1 , $y , $x2 , $y , $col ) ;
	}
}

####

sub sam_get_single_read_height {
	my ( $r ) = @_ ;
    
    return $height_nonzero_isize + abs($$r[2] - $from)%$height_zero_isize;

	# return abs ( ( $r->[2] - $from ) % 50 ) ;
	# my $y = hex ( substr md5_hex ( $r->[9] ) , 0 , 8 ) % ( $height - $scale_height ) ;
	# $y = $height - $scale_height - $y ; ori
	# return $y ;
}


sub sam_paint_single_short_reads {
	my ( $r , $col ) = @_ ;
#	print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output
#	my @out ;
	foreach ( @{$r} ) {
		my $r = $sam_reads{$_}->[0] ;
		my $y = sam_get_single_read_height ( $r ) ;
#		push @out , $r->[2] . "\t$from\n" ;
		sam_paint_single_short_read ( $r , $col , $y , 1 ) ;
	}
#	print join '' , sort @out ;
#	exit ;
}

sub sam_paint_short_single_reads_quality {
	my ( $r , $col ) = @_ ;
	foreach ( @{$r} ) {
		my $r = $sam_reads{$_}->[0] ;
		my $y = sam_get_single_read_height ( $r ) ;
		sam_paint_single_short_read_quality ( $r , $y , $col ) ;
	}
}

####

sub sam_paint_short_read_pairs_quality {
	my ( $r , $col ) = @_ ;
	foreach ( @{$r} ) {
		sam_paint_short_read_pair_quality ( $sam_reads{$_} , $col ) ;
	}
}

sub sam_paint_short_read_pair_quality {
	my ( $r , $col ) = @_ ;
	my $y = sam_get_y ( $r->[0] ) ;
#	my $y = abs ( $r->[0]->[7] ) ;
#	$y = ( $height - $scale_height ) * $y / $max_dist ;
#	$y = $height - $scale_height - $y ;
	sam_paint_single_short_read_quality ( $r->[0] , $y , $col ) ;
	sam_paint_single_short_read_quality ( $r->[1] , $y , $col ) ;
}

sub sam_paint_single_short_read_quality {
	my ( $r , $y , $col ) = @_ ;
	my $qs = $r->[9] ;
	chomp $qs ;
	my $l = length ( $qs ) - 1 ;
	my ( $lastx , $lasty ) ;
	
	my $lq = 99999 ;
	my ( $lqx , $lqy ) ;
	
	foreach ( 0 .. $l ) {
		my $ch = substr $qs , $_ , 1 ;
		my $x = int ( ( $r->[2] - $from + $_ ) * $width / $ft ) ;
		my $q = ord ( $ch ) - 33 ;
		my $qh = 30 - $q ;
		if ( $q < $lq ) {
			$lq = $q ;
			$lqx = $x ;
			$lqy = $y + $qh ;
		}
		$im->line ( $lastx , $lasty , $x , $y+$qh , $col ) if defined $lastx ;
		$lastx = $x ;
		$lasty = $y + $qh ;
	}
	
	$im->string(gdSmallFont,$lqx,$lqy,$lq, $col) ;
}

#______________________________________________________________________________________


####

sub sam_paint_short_read_pairs {
	my ( $r , $col , $draw_snps ) = @_ ;
	foreach ( @{$r} ) {
		sam_paint_short_read_pair ( $sam_reads{$_} , $col , $draw_snps ) ;
	}
}

sub sam_paint_short_read_pair {
	my ( $r , $col , $draw_snps ) = @_ ;
	my $y = sam_get_y ( $r->[0] ) ;
#	my $y = abs ( $r->[0]->[7] ) ;
#	$y = ( $height - $scale_height ) * $y / $max_dist ;
#	$y = $height - $scale_height - $y ;
	sam_paint_single_short_read ( $r->[0] , $col , $y , $draw_snps ) ;
	sam_paint_single_short_read ( $r->[1] , $col , $y , $draw_snps ) ;
}

sub sam_paint_single_short_read 
{
	my ( $r , $col , $y , $draw_snps ) = @_ ;

    if ( $y < 0 ) { return; }
    if ( $y > $height ) { return; }

    # Draw each cigar segment separately
    my $xfrom = $r->[2] - $from;
    my $cigar = $r->[4];
    my @coords;
    while ($cigar=~/^(\d+)(\D)/)
    {
        my $nbases = $1;
        my $type   = $2;
        $cigar = $';

        my $x1 = $xfrom;
        my $x2 = $xfrom + $nbases;

        $xfrom += $nbases;
        if ( $x1 > $to-$from or $x2 < 0 ) { next; }

        $x1 = int( $x1*$width/$ft);
        $x2 = int( $x2*$width/$ft);

        if ( $x1<0 ) { $x1 = 0; }
        if ( $x2>=$width ) { $x2 = $width-1; }
        if ( $x1 == $x2 ) { next; }

        if ( $type eq 'M' ) 
        { 
            $im->line($x1, $y, $x2, $y, $col);

            # Save the coordinates to draw the arrows
            if ( !@coords or !$cigar ) { push @coords, $x1,$x2; }

        }
        elsif ( $type eq 'I' )
        {
            $im->line($x1, $y-1, $x2, $y-1, $col_insertion);
            $xfrom -= $nbases;
        }
        elsif ( $type eq 'D' )
        {
            $im->line($x1, $y, $x2, $y, $col_deletion);
        }
    }

    if ( $sam_show_read_arrows && @coords )
    {
        if ( $r->[0] & 0x0010 ) 
        { 
            # reverse
            $im->line($coords[0], $y , $coords[0]+2 , $y+2, $col );
        } 
        else 
        {
            $im->line($coords[-1], $y , $coords[-1]-2 , $y-2, $col );
        }
    }

    if ( $draw_snps )
    {
        # while ( $r->[8] =~ m/[a-z]/g ) 
        # {
        #     $x1 = int ( ( $r->[2] - $from + pos ( $r->[8] ) - 1  ) * $width / $ft ) ;
        #     $im->line ( $x1 , $y-2 , $x1 , $y+2 , $sam_col_mismatch ) ;
        # }
    }
	#   return unless $draw_snps ;
	#   return if 0 == sam_check4snps ( $r ) ;
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
	} elsif ( $view eq 'deletions' ) {
		&dump_image_deletions ;
	}
} elsif ( $output eq 'url' ) {
	if ( $view eq 'annotation' ) {
		dump_image_annotation ( 1 ) ;
	}
} elsif ( $output eq 'text' ) {
	if ( $view eq 'pileup' ) {
		&dump_image_pileupview ;
	} elsif ( $view eq 'coverage' ) {
		&dump_image_coverageview ;
	} else {
		&dump_text ;
	}
}



# Call external renderer
sub render_bam_file {
	return 0 if $view ne 'indel' and $view ne 'coverage' and $view ne 'pileup' ;
	return 0 if 1 != scalar @databases ;
	my $file = $databases[0] ;
	return 0 unless $file =~ m/\.bam$/ ;
	return 0 unless defined $reference_fa ;
	

	# BAM file
	my $cwd;
	if ( $bam_ftp )
	{
		# Samtools require chdir to read bam .bai index files.
		$cwd = `/bin/pwd`;
		chomp($cwd);
		$cwd =~ /(.*)/;
		$cwd = $1;
		$datapath =~ /(.*)/;
		$datapath = $1;
		chdir($datapath) or die "Could not chdir $datapath: $!";
		$file = "$bam_ftp/$file";
	}
	else
	{
		$file = "$datapath/$file";
	}

	# Reference file
	my $ref_file ;
	if ( $reference_ftp )
	{
		$cwd = `/bin/pwd`;
		chomp($cwd);
		$cwd =~ /(.*)/;
		$cwd = $1;
		$datapath =~ /(.*)/;
		$datapath = $1;
		chdir($datapath) or die "Could not chdir $datapath: $!";
		$ref_file = "$reference_ftp/$reference_fa";
	}
	else
	{
		$ref_file = "$datapath/$reference_fa";
	}
	
	# PNG output
	my ( $tmp_fh , $tmp_fn ) = tempfile ( "LookSeqXXXXXX" , '/tmp' ) ;
	close $tmp_fh ;
	
	# Taint vodoo
	$file =~ /([a-zA-Z0-9_\-\.\/\:]+)/ ;
	$file = $1 ;
	$chromosome =~ /([a-zA-Z0-9_\-\.\/\:]+)/ ;
	$chromosome = $1 ;
	$width =~ /(\d+)/ ;
	$width = $1 ;
	$height =~ /(\d+)/ ;
	$height = $1 ;
	$view =~ /(\w+)/ ;
	$view = $1 ;

	# Construct command
	my $cmd = $render_image_command ;
	$cmd .= " --bam=\"$file\"" ;
	$cmd .= " --ref=$ref_file" ;
	$cmd .= " --region=\"$chromosome:$from-$to\"" ;
	$cmd .= " --png=$tmp_fn" ;
	$cmd .= " --width=$width" ;
	$cmd .= " --height=$height" ;
	$cmd .= " --view=$view" ;
	$cmd .= " --vmax=$max_insert_size" if $max_insert_size =~ m/^\d+$/ ;
	
	my @options ;
	push @options , 'pairs' if $display_perfect ;
	push @options , 'snps' if $display_snps and $to - $from < 2000000 ; # HARD SNP CUTOFF; will only show red instead of blue anyway
	push @options , 'inversions' if $display_inversions ;
	push @options , 'single' if $display_single ;
	push @options , 'linkpairs' if $display_pair_links ;
	push @options , 'noscale' if $display_noscale ;
	push @options , 'arrows' if $sam_show_read_arrows ;
	push @options , 'readqual' if $sam_show_read_quality ;
	push @options , 'text' if $output eq 'text' and $view eq 'pileup' ;
	push @options , 'faceaway' if $display_faceaway ;
	push @options , 'colordepth' ; # FIXME always on
# my $display_pot_snps = $display =~ m/\|potsnps\|/ ; # FIXME
	
	$cmd .= " --options=" . join ( ',' , @options ) ;
	$cmd .= " --mapq=$mapq_cutoff" if $mapq_cutoff > 0 ;
	
	# Run command and print output
#	print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output
#	print "$cmd\n" ; exit ( 0 ) ;
	`$cmd` ;
#	print -s $tmp_fn ;
	
	if ( $output eq 'text' and $view eq 'pileup' ) {
		print $cgi->header(-type=>'text/plain',-expires=>'-1s');
		open FILE , $tmp_fn ;
		while ( <FILE> ) {
			print $_ ;
		}
		close FILE ;
	} else {
		print $cgi->header(-type=>'image/png',-expires=>'-1s');
		open PNG , $tmp_fn ;
		binmode STDOUT ;
		binmode PNG ;
		my $buff ;
		while (read(PNG, $buff, 8 * 2**10)) {
			print STDOUT $buff;
		}
		close PNG ;
	}
	
	
	# Cleanup
	unlink $tmp_fn ;
	
	
	
	## --bam=/nfs/users/nfs_m/mm6/ftp/ag/bam/AC0001-C.bam 
	## --options=snps,pairs,arrows,single,faceaway,inversions,linkpairs,colordepth 
	## --ref=/nfs/users/nfs_m/mm6/ftp/ag/Anopheles_gambiae.clean.fa 
	# --region="2L:1-200000" 
	## --png=2L.a.png
	
#	print "$file\n" ;
#			sam_read_data($file);
	return 1 ;
}




exit;


#----------- Plot ----------------------

package Plot;

use GD;

sub Plot::new
{
    my ($class,$args) = @_;

    my $self = $args ? $args : {};
    
    if ( !$$self{width} ) { die "Expected width parameter.\n"; }
    if ( !$$self{height} ) { die "Expected height parameter.\n"; }
    if ( !$$self{bgcolor} ) { die "Expected bgcolor parameter\n"; }
    if ( !$$self{fgcolor} ) { die "Expected fgcolor parameter\n"; }

    if ( !exists($$self{margin_left}) ) { $$self{margin_left}=0; }
    if ( !exists($$self{margin_right}) ) { $$self{margin_right}=0; }
    if ( !exists($$self{margin_top}) ) { $$self{margin_top}=10; }
    if ( !exists($$self{margin_bottom}) ) { $$self{margin_bottom}=10; }

    $$self{image} = new GD::Image ($$self{width}, $$self{height});

    $$self{fgcolor} = $$self{image}->colorAllocate(@{$$self{fgcolor}});
    $$self{bgcolor} = $$self{image}->colorAllocate(@{$$self{bgcolor}});

    $$self{image}->filledRectangle(0,0,$width,$height,$$self{bgcolor});

    # Image dimensions without the margins
    $$self{area_width}  = $$self{width} - $$self{margin_left} - $$self{margin_right};
    $$self{area_height} = $$self{height} - $$self{margin_top} - $$self{margin_bottom};
    $$self{baseline} = exists($$self{baseline}) ? $$self{baseline} : 0;

    bless $self, ref($class) || $class;
    return $self;
}

sub Plot::set
{
    my ($self,$args) = @_;
    while (my ($key,$value)=each %$args)
    {
        if ( $key eq 'fgcolor' ) { $$self{fgcolor}=$$self{image}->colorAllocate(@$value); }
        if ( $key eq 'bgcolor' ) { $$self{bgcolor}=$$self{image}->colorAllocate(@$value); }
    }
}

sub Plot::set_scale_factor
{
    my ($self,$data) = @_;

    my ($ymax,$ymin);

    if ( !exists($$self{ymin}) || !exists($$self{ymax}) )
    {
        $ymax = $$data[0];
        $ymin = $$data[0];
        for my $value (@$data)
        {
            if ( $value>$ymax ) { $ymax=$value; }
            if ( $value<$ymax ) { $ymin=$value; }
        }
    }
    $ymin = exists($$self{ymin}) ? $$self{ymin} : $ymin;
    $ymax = exists($$self{ymax}) ? $$self{ymax} : $ymax;

    $$self{xscale} = $$self{area_width}/(scalar @$data-1);
    $$self{yscale} = $ymax ? $$self{area_height}/$ymax : 1;
}

sub Plot::polygon
{
    my ($self,$data,$poly) = @_;

    if ( !$$self{xscale} ) { $self->set_scale_factor($data); }

    my $xscale = $$self{xscale};
    my $yscale = $$self{yscale};

    my $ndata  = scalar @$data;
    my $prev_x = -1;
    my $prev_y = 0;

    if ( !$poly ) { $poly = new GD::Polygon; }
    for (my $i=0; $i<$ndata; $i++)
    {
        my $x = int($i*$xscale) + $$self{margin_left};
        my $y = $$self{area_height} - int($$data[$i]*$yscale) + $$self{margin_top};

        if ( $prev_x!=-1 && $prev_x==$x && $prev_y<$y ) { next; }
        $poly->addPt($x,$y);

        $prev_x = $x;
        $prev_y = $y;
    }
    return $poly;
}

sub Plot::scaled_polygon
{
    my ($self,$data) = @_;

    if ( !$$self{xscale} ) { $self->set_scale_factor($data); }

    my $poly = new GD::Polygon;
    $poly->addPt($$self{margin_left} + 0,$$self{margin_top} + $$self{area_height} - int($$self{baseline}*$$self{yscale}));
    $self->polygon($data,$poly);
    $poly->addPt($$self{width} - $$self{margin_right}-1,$$self{margin_top} + $$self{area_height} - int($$self{baseline}*$$self{yscale}));
    $$self{image}->setAntiAliased($$self{fgcolor});
    $$self{image}->filledPolygon($poly,$$self{fgcolor});
}

sub Plot::scaled_line
{
    my ($self,$data) = @_;
    my $poly = $self->polygon($data);
    $$self{image}->setAntiAliased($$self{fgcolor});
    $$self{image}->unclosedPolygon($poly,$$self{fgcolor});
}

sub Plot::baseline
{
    my ($self) = @_;
    my $y = $$self{margin_top} + $$self{area_height} - int($$self{baseline}*$$self{yscale});
    $$self{image}->line($$self{margin_left},$y,$$self{width}-$$self{margin_right}-1,$y,$$self{fgcolor});
}

sub Plot::legend
{
    my ($self,$legend) = @_;

    # How can one determine text dimensions of the rendered text in GD?
    #   ImageMagick seems to be much better in this.
    my $len = length $legend;

    my $y = 0;
    my $x = $$self{width} - $$self{margin_right} - gdSmallFont->width*$len;

    $$self{image}->string(gdSmallFont,$x,$y,$legend, $$self{fgcolor});
}

sub Plot::draw_yscalebar
{
    my ($self,$data) = @_;

    if ( !$$self{yscale} ) { $self->set_scale_factor($data); }

    my $max_value = $$self{area_height} / $$self{yscale};
    my $min_diff = 1.2 * gdSmallFont->height / $$self{yscale};
    my $units = $$self{units} ? $$self{units} : '';

    my $ticks_len = 10;
    my $pos = 0;

    my $i = 0;
    my $diff = 1;
    while ( $diff < $min_diff )
    {
        #print STDERR "diff=$diff\n";
        $diff *= $i%2 ? 2 : 5;
        $i++;
    }
    if ( $pos+$diff>$max_value )
    {
        # If we are here, there is only one tick at 0. Try to do multiples of 5 instead.
        $i = 0;
        $diff = 1;
        while ( $diff < $min_diff )
        {
            #print STDERR "diff=$diff\n";
            $diff *= 5;
            $i++;
        }
    }
    if ( $pos+$diff>$max_value && ($max_value-0)>$min_diff )
    {
        # Still no tick? If there should be only one tick (at 0), print at least the maximum.
        $diff = $max_value;
    }

    #print STDERR "diff=$diff max_val=$max_value min_diff=$min_diff font_height=",gdSmallFont->height," yscale=$$self{yscale}\n";

    while ( $pos<=$max_value )
    {
        my $y = $$self{area_height} - int($pos*$$self{yscale}) + $$self{margin_top};
        #print STDERR "$pos $y\n";

        $$self{image}->line($$self{margin_left}, $y, $$self{margin_left}+$ticks_len, $y, $$self{fgcolor});
        $$self{image}->string(gdSmallFont, $$self{margin_left} + $ticks_len*1.5, $y - gdSmallFont->height*0.5, $pos.$units, $$self{fgcolor});
        $pos += $diff;
    }
}

1;

