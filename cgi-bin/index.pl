#!/usr/local/bin/perl -T
#########
# Author:     Magnus Manske (mm6@sanger.ac.uk)
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

if ( $use_sanger_layout ) {
	use SangerPaths qw(core);
	use SangerWeb;
}

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
		foreach ( sort @dbs ) {
			next if $_ eq $annotation_file ;
			next if $_ eq $snp_file ;
			$lanes{$_} = "lanes.push ( \"$_\" ) ;\n" ;
			#$myscript .= "lanes.push ( \"$_\" ) ;\n" ;
		}
	}
	foreach ( sort keys %lanes ) {
		$myscript .= $lanes{$_} ;
	}
	
	my $width = $cgi->param('width') || 1024 ;
	$myscript .= "var img_width = $width ;\n" ;
	$myscript .= "var test = $test ;\n" ;
	$myscript .= "var sanger_layout = $use_sanger_layout ;\n" ;
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
	if ( $lane ne '' and $show =~ m/^([\w\d\._]+):(\d+)-(\d+),([a-z]+)$/ ) {
		$myscript .= "var use_init = true ;\nvar init_chr = \"$1\" ;\nvar init_from = $2 ;\nvar init_to = $3 ;\nvar init_mode = \"$4\" ;\nvar init_lane = \"$lane\" ;\n" ;
	} else {
		$myscript .= "var use_init = false ;\n" ;
	}
}

sub load_i18n_data {
	my $fn = shift ;
	my $interface_file = "$htmlpath/$fn.$language" ;
	$interface_file = "$htmlpath/$fn.en" unless -e $interface_file ; # Default to en
	return unless -e $interface_file ; # Even that's not there!
	open INTERFACE , $interface_file ;
	while ( <INTERFACE> ) {
		chomp ;
		my ( $key , $value ) = split "\t" , $_ , 2 ;
		next if $key eq '' ;
		$key =~ tr/ /_/ ;
		$i18n{$key} = $value ;
	}
	close INTERFACE ;
}

sub main {
	my $sw ;

	if ( $use_sanger_layout ) {
		my @js_files = ( "$webroot/lookseq.js" ) ;
		unshift @js_files , "$webroot/custom.js" if -e "$htmlpath/custom.js" ;
	    $sw  = SangerWeb->new({
		    'title'   => $tooltitle,
		    'banner'  => $tooltitle,
		    'author'  => q(mm6),
		    'inifile' => "$docroot/Info/header.ini",
		    'jsfile' => \@js_files,
		    'onload' => "init()",
		});
	    $cgi = $sw->cgi();
	} else {
		$cgi = new CGI ;
	}

	
	$language = lc ( $cgi->param('display') || 'en' ) ;
	$language =~ s/[^a-z]//g ;
	
	load_i18n_data ( 'interface' ) ;
	load_i18n_data ( 'custom' ) ;

	&prepare_myscript ;

	
	if ( $use_sanger_layout ) {
		my $dummy = $sw->header( { 'script' => $myscript } );
		$dummy =~ s/<script/<meta name="robots" content="robots_flag" \/><script/ ; # Hacking around stupid Sanger cgi object
		print $dummy ;
	} else {
		my @js_files = ( $myscript ,
						{ -type => 'text/javascript', -src => "$webroot/lookseq.js" }
						) ;
		unshift @js_files , { -type => 'text/javascript', -src => "$webroot/custom.js" } if -e "$htmlpath/custom.js" ;

		print $cgi->header ;
		print $cgi->start_html ( -script => \@js_files,
								 -onLoad => "init()" ,
								 -title => $tooltitle,
								 -meta => {
									'robots' => $robots_flag
								 }
								 ) ;
	}
	

	my $out = '' ;
	open HTML , "$htmlpath/lookseq.html" ;
	while ( <HTML> ) {
		$_ =~ s/__HTMLPATH__/$webroot/g ;
		$out .= $_ ;
	}
	close HTML ;
	
	foreach ( keys %i18n ) {
		my $key = "%$_%" ;
		my $value = $i18n{$_} ;
		$out =~ s/$key/$value/gi ;
	}
	
	print "$out" ;
#	print "!$htmlpath!" ;# print $sw->footer() ; exit ( 0 ) ;

	if ( $use_sanger_layout ) {
		print $sw->footer() ;
	} else {
		print $cgi->end_html() ;
	}
}

&main;
0;
