#!/usr/bin/perl
#########
# Author:     Magnus Manske (mm6@sanger.ac.uk)
# Group:      Team 112
#

use strict ;
use warnings ;
use Cwd ;
use CGI;
#use CGI::Carp qw/fatalsToBrowser warningsToBrowser/ ;
use URI::Escape;
use Digest::SHA1  qw(sha1 sha1_hex sha1_base64);
use JSON;
use Data::Dumper ;

my $cgi = new CGI ;
my $debug = $cgi->param ( 'debug' ) ;
my $action = $cgi->param ( 'action' ) || '' ;
my $callback = $cgi->param ( 'callback' ) ;

print $cgi->header(-type=>'text/plain',-expires=>'-1s') if $debug ;

# Read data
my $config = '' ;
open FILE , 'config.json' ;
while ( <FILE> ) {
    $config .= $_ ;
}
close FILE ;
my $data = from_json ( $config , {utf8 => 1} ) ;

my $render_image_command = $data->{'misc'}->{'render_command'} ;

#print Dumper ( $data ) . "\n\n" if $debug ;

my %samples = %{$data->{'samples'}} ;
my %species = %{$data->{'species'}} ;
my %groups ;
foreach my $sample ( keys %samples ) {
    my $g = $samples{$sample}->{'group'} || 'No group' ;
    push @{$groups{$g}} , $sample ;
}

#print Dumper ( %groups ) . "\n\n" if $debug ;

# Process action
my $json ;    

if ( $action eq 'login' ) { # No login
    $json->{'error'} = 'OK' ;
} elsif ( $action eq 'logout' ) { # No logout
    $json->{'error'} = 'OK' ;
} elsif ( $action eq 'test_login' ) { # Everybody's welcome!
    $json->{'error'} = 'OK' ;
    $json->{'username'} = "user" ;
} elsif ( $action eq 'render_image' ) {

	my $sample = $cgi->param('sample') ;
	my $species = $samples{$sample}->{'species'} ;
	my $algorithm = $cgi->param('alg') ;
	$algorithm =~ s/[^\w]//g ;
	my $bam_ftp = $samples{$sample}->{'bam'} ;
	render_image ( $species , $bam_ftp , $algorithm ) ;
	exit 0 ;

} elsif ( $action eq 'clear_temp' ) { # Obsolete
    $json->{'error'} = 'OK' ;
} elsif ( $action eq 'get_species_data' ) {

    my %sd ;
	$json->{'error'} = 'OK' ;
	foreach my $species ( keys %species ) {
		$sd{$species}->{'chromosomes'} = readFAI ( $species{$species}->{'fasta'} . '.fai' ) ;
	}
	$json->{'species'} = \%sd ;

} elsif ( $action eq 'mysamples' ) {

    my %algs ;
    foreach my $sample ( keys %samples ) {
        my $g = $samples{$sample}->{'group'} ;
        my $sp = $samples{$sample}->{'species'} ;
        
        my $alg = $samples{$sample}->{'alg'} ;
        $algs{$alg} = 1 ;
        
        unless ( defined $json->{'projects'}->{$g} ) {
            $json->{'projects'}->{$g}->{'code'} = $g ;
            $json->{'projects'}->{$g}->{'desc'} = '' ;
            $json->{'projects'}->{$g}->{'countries'} = {} ;
            $json->{'projects'}->{$g}->{'species'} = {} ;
        }
        
        $json->{'samples'}->{$sample}->{'sid'} = $sample ;
        $json->{'samples'}->{$sample}->{'pid'} = $g ;
        $json->{'samples'}->{$sample}->{'species'} = $sp ;
        $json->{'samples'}->{$sample}->{'alg'}->{$alg} = 1 ;
        $json->{'projects'}->{$g}->{'species'}->{$species{$sp}->{'name'}}++ ;
    }
	my @algs = sort keys %algs ;
    $json->{'error'} = 'OK' ;
	$json->{'algorithms'} = \@algs ;

} elsif ( $action eq 'getannotation' ) {

	my $chr = $cgi->param('chr') ;
	my $species = $cgi->param('species') ;
	if ( defined $species{$species} ) {
		$json->{'annotation'} = read_gff3 ( $species{$species}->{'gff'} , $chr , $cgi->param('query') ) ;
		$json->{'error'} = 'OK' ;
	} else {
		$json->{'error'} = 'Species or annotation not found' ;
	}

} else {
	$json->{'error'} = "Unknown action : $action" ;
}


print $cgi->header(-type=>'application/json',-expires=>'-1s') unless $debug ;
print "$callback(" if defined $callback ;
print to_json ( $json ) unless defined $debug ;
print Dumper ( $json ) if defined $debug ;
print ");" if defined $callback ;


0 ;

sub readFAI {
	my ( $file ) = @_ ;
	my %ret ;
	open FILE , $file ;
	while ( <FILE> ) {
		my @d = split "\t" , $_ ;
		next unless defined $d[1] ;
		$ret{$d[0]} = $d[1] ;
	}
	close FILE ;
	return \%ret ;
}

sub render_image {
	my ( $species , $bam_ftp , $algorithm ) = @_ ;
	
	my $display = $cgi->param('display') || '|perfect|snps|inversions|' ;
	my $output = $cgi->param('output') || 'text' ;
	my $view = $cgi->param('view') || '' ;
	my $database = $cgi->param('sample') ;
	my $from = $cgi->param('from') ;
	my $to = $cgi->param('to') ;
	my $chromosome = $cgi->param('chr') ;
	my $width = $cgi->param('width') || 1024 ;
	my $height = $cgi->param('height') || 512 ;
	my $sam_show_read_arrows = $display =~ m/\|orientation\|/ ;
	my $sam_show_base_quality = $display =~ m/\|basequal\|/ ;
	my $mapq_cutoff = $cgi->param('mapq') || 0;
	
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
	
	my $max_insert_size = $cgi->param('maxdist') ? $cgi->param('maxdist') : 'Auto';
	$max_insert_size =~ /(\d+)/; $max_insert_size=$1;
	
	$from =~ /(\d+)/ ; $from = $1 ;
	$to =~ /(\d+)/ ; $to = $1 ;
	
	my @databases = split ( ',' , $database ) ;
	
	my $reference_fa = $species{$species}->{'fasta'} ; #$refpath . '/' . $sd{$species}->{'reference'} ;
	my $reference_ftp = $reference_fa ;
	
	
	exit 0 if $view ne 'indel' and $view ne 'coverage' and $view ne 'pileup' ;
	exit 0 if 1 != scalar @databases ;

	my $file = $bam_ftp ;
	return 0 unless defined $reference_fa ;
	
	$file =~ m|^(.+)/[^/]+$| ;
	my $d = $1 ;

    my $use_ftp = $file =~ m|^ftp\://| ? 1 : 0 ;
	if ( $use_ftp ) {
		# Working directory
		my $cwd = ( $data->{'misc'}->{'tmp'} || "/tmp" ) . "/$algorithm" ;
		mkdir $cwd unless -d $cwd ;
		$d = $cwd ;
		#chdir ( $cwd ) ;
		
		# Remove .bai files older >1h
		`find /lookseq/data/tmp/ -name "*.bai" -ctime 0.04 -exec rm {} \;` ;
		
    }

	if ( $cgi->param('debug') eq '1' ) {
		print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output
    }
    
    print "cd $d\n" if ( $cgi->param('debug') eq '1' ) ;
	chdir $d ;
	
	# Reference file
	my $ref_file = $reference_fa ;
	
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
	$cmd .= " --png=-" ;
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
	push @options , 'readqual' if $sam_show_base_quality ;
	push @options , 'text' if $output eq 'text' and $view eq 'pileup' ;
	push @options , 'faceaway' if $display_faceaway ;
	push @options , 'colordepth' ; # FIXME always on
	# my $display_pot_snps = $display =~ m/\|potsnps\|/ ; # FIXME
	
	$cmd .= " --options=" . join ( ',' , @options ) ;
	$cmd .= " --mapq=$mapq_cutoff" if $mapq_cutoff > 0 ;
	
	# Run command and print output
	if ( $cgi->param('debug') eq '1' ) {
		print "$cmd\n" ; exit ( 0 ) ;
	}
#	`$cmd` ;
	#	print -s $tmp_fn ;
	
	if ( $output eq 'text' and $view eq 'pileup' ) {
		print $cgi->header(-type=>'text/plain',-expires=>'-1s');
		open FILE , "$cmd |" ;#$tmp_fn ;
		while ( <FILE> ) {
			print $_ ;
		}
		close FILE ;
	} else {
		print $cgi->header(-type=>'image/png',-expires=>'-1s');
		open PNG , "$cmd |" ;#$tmp_fn ;
		binmode STDOUT ;
		binmode PNG ;
		my $buff ;
		while (read(PNG, $buff, 8 * 2**10)) {
			print STDOUT $buff;
		}
		close PNG ;
	}
	
	
	# Cleanup
#	unlink $tmp_fn ;
}

sub read_gff3 {
	my ( $file , $chromosome , $query ) = @_ ;
	my %ret ;
	my $af = $file ;
	return \%ret unless -f $af ;
	
	my $cmd = "grep -i '$chromosome' $af |" ;
	$cmd = "cat $af |" unless defined $chromosome ;

	my @matches ;
	if ( defined $query ) {
		$query =~ s/[^\w.:-]//g ; # Paranoia
		my @q = split /\s+/ , $query ;
		foreach ( @q ) {
			next if $_ eq '' ;
			push @matches , $_ ;
			$cmd .= " grep -i '$_' |" ;
		}
	}
	
	open FILE , $cmd ;
	while ( <FILE> ) {
		next if $_ =~ m/^#/ ;
		chomp ;
		my ( $chr , $db , $type , $from , $to , $x1 , $strand , $x2 , $l ) = split "\t" , $_ ;
		next unless defined $l ;
		$chr =~ s/^.*\|// ;
		next unless not defined $chromosome or lc $chr eq lc $chromosome ;
		next if $type eq 'supercontig' ;
		next if $type eq 'mRNA' or $type eq 'CDS' ; # HACK deemed unneccessary ("exon" covers it all...)

		my $next = 0 ;
		foreach my $q ( @matches ) {
			$next = 1 unless $l =~ m/$q/i ;
		}

		my %d ;
		$d{'from'} = $from ;
		$d{'to'} = $to ;
		$d{'strand'} = $strand unless $strand eq '' ;
		$d{'chr'} = $chr unless defined $chromosome ;
		
		if ( $l =~ m/^ID=([^;]+)/ ) {
			my $name = $1 ;
			$name =~ s/^.+\|// ;
			$d{'name'} = $name if $name ne '' ;
		} else {
			$l =~ m/;name=([^;]+);{0,1}/i ;
			$d{'name'} = $1 if defined $1 and $1 ne lc $type ;
		}

		push @{$ret{$type}} , \%d ;
	}
	close FILE ;
	return \%ret ;
}
