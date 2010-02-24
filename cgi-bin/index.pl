#!/usr/local/bin/perl -T
#########
# Author:     Magnus Manske (mm6@sanger.ac.uk)
#             Petr Danecek (pd3@sanger.ac.uk), Team 145
# Group:      Team 112
#

BEGIN { push @INC, "."; }

use strict;
use warnings ;
use DBI ;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Basename;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval nanosleep );
use settings ;

my ( $cgi , $myscript , $language ) ;
my %i18n ;

sub prepare_myscript {
	$myscript = '' ;
	my $test = $cgi->param('test') || 0 ;
	
	$myscript .= "var current_lane = '' ;\n" ;
	$myscript .= "var lanes = new Array () ;\n" ;
	
	my %lanes ;
	if ( $use_mysql ) {
		foreach my $db ( @mysql_dbs ) {
			my $dbh = DBI->connect("DBI:mysql:database=$db;host=$mysql_server;port=$mysql_port", $mysql_user, $mysql_password, {'PrintError'=>0});
			if ( $dbh ) {
				my $sth = $dbh->prepare("SELECT * FROM meta WHERE the_key='name'");
				$sth->execute();
				my $name ;
				while (my $ref = $sth->fetchrow_hashref()) {
					$name = $ref->{'value'} ;
				}
				if ( defined $name ) {
					$lanes{$name} = "lanes.push ( \"$db|$name\" ) ;\n" ;
					#$myscript .= "lanes.push ( \"$db|$name\" ) ;\n" ;
				} else {
					$lanes{$db} = "lanes.push ( \"$db\" ) ;\n" ;
					#$myscript .= "lanes.push ( \"$db\" ) ;\n" ;
				}
			} else {
				$lanes{$db} = "lanes.push ( \"$db\" ) ;\n" ;
				#$myscript .= "lanes.push ( \"$db\" ) ;\n" ;
			}
		}
	} else {
		opendir(DIR, $datapath) || die "can't opendir $datapath: $!";
	    my @dbs = grep { /\.sqlite$/ } readdir(DIR);
	    closedir DIR;
		opendir(DIR, $datapath) || die "can't opendir $datapath: $!";
	    my @dbs_bam = grep { /\.bam$/ } readdir(DIR);
	    closedir DIR;

		if ( 0 == scalar @dbs_bam ) {
			opendir(DIR, $datapath) || die "can't opendir $datapath: $!";
			@dbs_bam = grep { /\.bai$/ } readdir(DIR);
			foreach ( 0 .. $#dbs_bam ) {
				$dbs_bam[$_] =~ s/\.bai$// ;
			}
			closedir DIR;
		}
		
		push @dbs , @dbs_bam ;
		
#		print $cgi->header(-type=>'text/plain',-expires=>'-1s'); # For debugging output
#		print join "\n" , @dbs ;
		
		foreach ( sort @dbs ) {
			next if (defined $annotation_file && $_ eq $annotation_file);
			next if (defined $snp_file && $_ eq $snp_file);
			$lanes{$_} = "lanes.push ( \"$_\" ) ;\n" ;
			#$myscript .= "lanes.push ( \"$_\" ) ;\n" ;
		}
	}
	
	if ( $paranoia_mode ) {
		my $show_samples = $cgi->param('showSamples') || '' ;
		$myscript .= "paranoia_mode = true ;\nshow_samples = \"$show_samples\" ;\n" ;
		my @keep = split "," ,  $show_samples ;
		my $lane = $cgi->param('lane') ;
		push @keep , $lane if defined $lane ;
		my %n ;
		foreach ( @keep ) {
			$n{$_} = $lanes{$_} if defined $lanes{$_} ;
		}
		%lanes = %n ;
	} else {
		$myscript .= "paranoia_mode = false ;\n" ;
	}
	
	foreach ( sort keys %lanes ) {
#		print $lanes{$_} ;
		$myscript .= $lanes{$_} ;
	}

    if ( -e "$datapath/data.info" )
    {
        open(my $fh,'<',"$datapath/data.info") or die "$datapath/data.info: $!"; 
        my @lines = <$fh>;
        close($fh);
        $myscript .= join('',@lines);
    }

	my $skin = $cgi->param('skin') || 'lookseq' ;
	$skin =~ s/[^A-Za-z]//g ;
	$myscript .= "var skin = \"$skin\" ;\n" ;
	
	my $width = $cgi->param('width') || 1024 ;
	$myscript .= "var img_width = $width ;\n" ;
	$myscript .= "var test = $test ;\n" ;
	$myscript .= "var display_init = \"" . ( $cgi->param('display') || '' ) . "\" ;\n" ;
	$myscript .= "var cgi_path = '$cgi_path' ;\n" ;
	$myscript .= "var language = '$language' ;\n" ;
	$myscript .= "var indel_zoom = '" . ( $cgi->param('maxdist') || 'auto' ) . "' ;\n" ;
	
	$myscript .= "\ni18n = new Array() ;\n" ;
	foreach my $k ( keys %i18n ) {
		my $v = $i18n{$k} ;
		$myscript .= "i18n['$k'] = \"$v\" ;\n" ;
	}
	$myscript .= "\n" ;

	my $show = $cgi->param('show') || '' ;
	my $lane = $cgi->param('lane') || '' ;
    my $win  = $cgi->param('win')  || '';
	if ( $lane ne '' and $show =~ m/^([\w\d\._]+):(\d+)-(\d+),([a-z_]+)$/ ) 
    {
        my $chrom = $1;
        my $from  = $2;
        my $to    = $3;
        my $mode  = $4;
        if ( $win ) 
        { 
            $to = $from+$win;
            $myscript .= "var init_max_window=$win;\n";
        }
        else 
        { 
            $myscript .= "var init_max_window=0;\n";
        }
		$myscript .= "var use_init = true ;\nvar init_chr = \"$chrom\" ;\nvar init_from = $from ;\nvar init_to = $to ;\nvar init_mode = \"$mode\" ;\nvar init_lane = \"$lane\" ;\n" ;
	} else {
		$myscript .= "var use_init = false ;\n" ;
		$myscript .= "var init_max_window=0;\n";
	}
}

sub load_i18n_data {
	my $fn = shift ;
	my $interface_file = "$htmlpath/$fn.$language" ;
	$interface_file = "$htmlpath/$fn.en" unless -e $interface_file ; # Default to en
	return unless -e $interface_file ; # Even that's not there!
	open INTERFACE , $interface_file or die "$interface_file: $!";
	while ( <INTERFACE> ) {
		chomp ;
		my ( $key , $value ) = split "\t" , $_ , 2 ;
		next if (!$key || $key eq '');
		$key =~ tr/ /_/ ;
		$i18n{$key} = $value ;
	}
	close INTERFACE ;
}

sub debug
{
    my (@msg) = @_;
    my $msg = join('',@msg);
    $msg =~ s{\n}{_}g;
    print STDERR "index.pl: $msg\n";
}

sub main {
	my $sw ;

	$cgi = new CGI ;

	
	$language = lc ( $cgi->param('display') || 'en' ) ;
	$language =~ s/[^a-z]//g ;
	
	load_i18n_data ( 'interface' ) ;
	load_i18n_data ( 'custom' ) ;

	&prepare_myscript ;
	
	my $skin = $cgi->param('skin') || 'lookseq' ;
	$skin =~ s/[^A-Za-z]//g ;
	
	my @js_files = ( $myscript ,
					{ -type => 'text/javascript', -src => "$webroot/$skin.js" }
					) ;
	unshift @js_files , { -type => 'text/javascript', -src => "$webroot/custom.js" } if -e "$htmlpath/custom.js" ;

	print $cgi->header;
	#print $cgi->start_html ( -script => \@js_files,
	my $x = $cgi->start_html ( -script => \@js_files,
							 -onLoad => "init()" ,
							 -title => $tooltitle,
							 -meta => {
								'robots' => $robots_flag
							 }
							 );
	$x =~ s{<!DOCTYPE[^>]+>}{<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">};
	$x =~ s{<html[^>]+>}{<html xmlns="http://www.w3.org/1999/xhtml">};
	print $x;
	
	my $out = '' ;
    my $layout = $cgi->param('xxx') ? 'test.html' : "$skin.html";
	open HTML , "$htmlpath/$layout" or die "$htmlpath/$layout: $!";
	while ( <HTML> ) {
		$_ =~ s/__HTMLPATH__/$webroot/g ;
		$out .= $_ ;
	}
	close HTML ;
	
	my $mapq_param = $cgi->param('mapq') || '' ;
	my $mapq = '' ;
	foreach ( 1 .. 50 ) {
		my $checked = $mapq_param eq $_ ? ' selected' : '' ;
		$mapq .= "<option value='$_'$checked>$_</option>\n" ;
	}
	$out =~ s/MAPQ_SELECT/$mapq/ ;
	
	foreach ( keys %i18n ) {
		my $key = "%$_%" ;
		my $value = $i18n{$_} ;
		$out =~ s/$key/$value/gi ;
	}
	
	print "$out" ;
#	print "!$htmlpath!" ;# print $sw->footer() ; exit ( 0 ) ;

	print $cgi->end_html() ;
}

&main;
0;
