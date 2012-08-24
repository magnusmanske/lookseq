if ( typeof ( console ) == 'undefined' ) {
	console = {
		log : function ( msg ) { } // Ignore
	} ;
}


String.prototype.ucFirst = function()
{
    return this.charAt(0).toUpperCase() + this.substring(1);
};


//_____________________________________________________________________________________________________________________________________________________________________
// BEGIN LOOKSEQ OBJECT

var lookseq = {

//_____________________________________________________________________________________________________________________________________________________________________
// CLASSES

Track : function () {
	this.id = '' ;
	this.img_id = '' ;
	this.type = '' ;
	this.sample = '' ;
	this.alg = '' ;
	this.note = '' ;
	this.fixed_top = '' ;
	this.track_id = '' ;
	this.image_height = 100 ;

	this.height = function () { 
		if ( this.type == 'scroll' ) return 20 ;
		if ( this.type == 'position' ) return 20 ;
		if ( this.type == 'annotation' ) return 40 ;
		if ( lookseq.mode == 'indel' ) return 200 ;
		if ( this.img_id == '' || $('#'+this.img_id).length == 0 ) return this.image_height ;
		return $('#'+this.img_id).height() ;
	} ;
	
	this.width = function () {
		return lookseq.getMainWidth() ;
	} ;
	
	this.getImageURL = function ( output ) {
		var url = lookseq.api_url ;
		url += '?action=render_image' ;
		url += '&alg=' + this.alg ;
		url += '&from=' + lookseq.getFrom() ;
		url += '&to=' + lookseq.getTo() ;
		url += '&chr=' + lookseq.getChr() ;
		url += '&sample=' + this.sample ;
		url += '&width=' + this.width() ;
		url += '&height=' + this.height() ;
		url += '&maxdist=' + $('#select_maxdist').val() ;
		url += '&view=' + lookseq.mode ;
		if ( output == 'text' ) url += '&output=text' ;
		else url += '&output=image' ;
		url += '&display=|noscale|' ;
		
		$.each ( lookseq.display_options , function ( k , v ) {
			if ( $('#display_'+v).is(':checked') ) url += v + '|' ;
		} ) ;
		
		url += '&debug=0' ;
		return url ;
	} ;
	
	this.paint = function () {
		if ( lookseq.initializing ) return ;

		if ( this.type == 'sample' ) this.paint_image();
		if ( this.type == 'position' ) this.paint_position();
		if ( this.type == 'annotation' ) this.paint_annotation();
		if ( this.type == 'scroll' ) this.paint_scroll();
	} ;
	
	this.paint_scroll = function () {
		var mw = lookseq.getMainWidth() ;
		var cl = lookseq.current_chromosome_length ;
		var f = lookseq.getFrom() ;
		var t = lookseq.getTo() ;
		var hw = Math.floor((t-f)*mw/cl) ;
		if ( hw < 5 ) hw = 5 ;

		$('#'+this.div).css ( { height : this.height , width : mw } ) ;
		
		$('#'+this.div).html('<div id="slider"></div>') ;
		$('#slider').slider ( {
			max : cl - ( t - f ) ,
			value : f ,
			slide: function(event, ui) {
				var f = $('#slider').slider('value') ;
				var t = f + parseInt(lookseq.to) - parseInt(lookseq.from) ;
				$('#from').val(f) ;
				$('#to').val(t) ;
			} ,
			change: function(event, ui) {
				lookseq.updatePositions();
			} ,
		} );
		$('#slider .ui-slider-handle').css ( 'width' , hw+'px' ) ;
//		$('#slider').slider ( 'option' , 'max' , cl ) ;
//		$('#slider').slider ( 'value' , f ) ;
	} ;
	
	this.paint_position = function () {
		$('#'+this.div).css ( { height : this.height , width : lookseq.getMainWidth() } ) ;
		var tw = this.width() ;

		var h = '' ;
		if ( this.note != '' ) h += "<div class='tracknote'>" + this.note + "</div>" ;
		h += "<div class='trackimage' id='" + this.img_id + "' style='width:" + tw + "px'>" ;
		
		var w = lookseq.to - lookseq.from ;
		var step = 10 ;
		while ( w / ( step * 1 ) > 10 ) step *= 10 ;
		
		var start = lookseq.from - w ;
		if ( start < 1 ) start = 1 ;
		start = start - ( start % step ) ;
		var maxp = lookseq.to + w ;
		if ( maxp > lookseq.current_chromosome_length ) maxp = lookseq.current_chromosome_length ;
		for ( var p = start ; p <= maxp ; p += step ) {
			var p2 = p == 0 ? 1 : p ;
			var ox = ( p - lookseq.from ) * tw / w ;
			h += "<div style='left:" + Math.floor(ox) + "px' class='position_marker'>" + lookseq.formatPosition ( p2 ) + "</div>" ;
		}
		
		h += "</div>" ;
		$('#'+this.div).html ( h ) ;
	} ;
	
	this.paint_annotation = function () {
		$('#'+this.div).css ( { height : this.height , width : lookseq.getMainWidth() } ) ;
		var tw = this.width() ;

		var h = '' ;
		if ( this.note != '' ) h += "<div class='tracknote'>" + this.note + "</div>" ;
		h += "<div class='trackimage' id='" + this.img_id + "' style='width:" + tw + "px'>" ;
		
		var w = lookseq.to - lookseq.from ;
		
		var start = lookseq.from - w ;
		if ( start < 1 ) start = 1 ;
		var stop = start + w * 3 ;
		
		var chr = lookseq.getChr();
		var rh = Math.floor ( this.height() / 3 ) - 1 ;
		if ( typeof lookseq.annotation === "undefined" ) lookseq.annotation = [] ;
		$.each ( lookseq.annotation , function ( thetype , typedata ) {
			var y = 2 ;
			var col = '#79FC4E' ;
			if ( thetype == 'gene' ) { y = 0 ; col = '#2F74D0' ; }
			if ( thetype == 'exon' ) { y = 1 ; col = '#FF7575' ; }
			y = y + y * rh ;
			
			$.each ( typedata , function ( k , v ) {
				if ( v.from > stop ) return ;
				if ( v.to < start ) return ;
				var x1 = ( v.from - lookseq.from ) * tw / w ;
				var x2 = ( v.to - lookseq.from ) * tw / w ;
				var xw = x2 - x1 ;
				var name = v.name ;
				if ( name === undefined ) name = '' ;
				var hover = name ;
				if ( hover != '' ) hover += ": \n" ;
				hover += thetype + "; \n" + v.from + '-' + v.to ;
				if ( v.strand !== undefined ) hover += ' (' + v.strand + ')' ;
				var f = parseInt(v.from) - 50 ;
				var t = parseInt(v.to) + 50 ;
				
				if ( v.strand == '+' ) name += ' &rarr;' ;
				if ( v.strand == '-' ) name += ' &larr;' ;
				
				h += "<div style='left:" + Math.floor(x1) + "px;width:" + Math.floor(xw) + "px;top:" + y + "px;height:" + rh + "px;background-color:" + col + "' " ;
				h += "title=\"" + hover.replace('"', '\'') + "\" " ;
				h += "ondblclick='lookseq.showRegion(\""+chr+"\","+f+","+t+")' " ;
				h += "class='annotation_marker'>" + name + "</div>" ;
			} ) ;
			
		} ) ;
		
		h += "</div>" ;
		$('#'+this.div).html ( h ) ;
	} ;
	
	this.paint_image = function () {
		$('#'+this.div).css ( { height : this.height , width : lookseq.getMainWidth() } ) ;
		this.img_id = this.id+'_img' ;

		var h = '' ;
		if ( this.note != '' ) h += "<div class='tracknote'>" + this.note + "</div>" ;
		
		if ( lookseq.mode == 'pileup' ) {
			h += "<div class='pileuptextlink'><a href='" ;
			h += this.getImageURL('text') ;
			h += "' target='_blank'>as text</a></div>" ;
		}
		
		h += "<div class='trackimage'>" ;
		h += "<img track_id='" + this.track_id ;
		h += "' id='" + this.img_id + "' " ;
		h += "width='" + this.width() + "px' " ;
		if ( lookseq.mode == 'indel' ) h += "height='" + this.height() + "px' " ;
		h += " class='trackimage_img' " ;
		if ( lookseq.mode != 'indel' ) h += " onload='lookseq.recalculateTrackLayout()' " ;
		h += "/></div>" ;
		$('#'+this.div).html ( h ) ;
		var reload_counter = 10 ;
		$('#'+this.img_id).error ( function () {
			if ( reload_counter <= 0 ) return ;
			reload_counter-- ;
			var src = $(this).attr('src') ;
			$(this).removeAttr('src') ;
			$(this).attr('src',src) ;
		} ) . attr ( 'src' , this.getImageURL() ) ;
		$('#'+this.img_id).dblclick ( function (e) { lookseq.trackimageDblClick(e,this); return false } ) ;
		$('#'+this.img_id).mousemove ( function (e) { lookseq.trackimageMouseMove(e,this);} ) ;
		$('#'+this.img_id).draggable ( {
			axis: 'x',
			start : function(event, ui) {
				$('#position_on_chromosome').hide() ;
			} ,
			drag: function(event, ui) {
				var x = ui.position.left ;
				lookseq.showNewPosition ( x ) ;
				$.each ( lookseq.tracks , function ( k , v ) {
					$('#'+v.img_id).css({left:x+'px'});
				} ) ;
			} ,
			stop: function(event, ui) {
				$('#position_on_chromosome').show() ;
				lookseq.updatePositions();
			}
		} ) ;
	}
} ,

//_____________________________________________________________________________________________________________________________________________________________________
// VARIABLES

mybaseurl : './index.html' ,
api_url : '' ,
logged_in : false ,
username : '' ,
email : '' ,
algorithms : new Object ,
species : new Object ,
samples : new Array() ,
projects : new Array() ,
got_samples : false ,
selected_samples : new Object ,
current_chromosome_length : 0 ,
tracks : new Array () ,
track_id : 0 ,
scrollbar_width : 0 ,
from : 0 ,
to : 0 ,
chr : '' ,
current_species : '',
min_size : 200 ,
hard_chrlen_limit : 5000000 , // 5mbp
max_pileup : 60000 ,
display_options : ['perfect','snps','single','inversions','pairlinks','orientation','faceaway','basequal'] ,
first_init : true ,
loaded_from_params : false ,
initializing : false ,


//_____________________________________________________________________________________________________________________________________________________________________
// INIT

init : function () {
	$('#select_maxdist').val('500') ;
	$('#select_maxdist').change(function(){lookseq.updateCurrentURL();lookseq.repaintTracks()}) ;
	lookseq.scrollbarWidth() ;
	lookseq.clear() ;
	lookseq.getSpeciesData() ;
} ,

clear : function () {
	lookseq.species = new Object ;
	lookseq.samples = new Array() ;
	lookseq.projects = new Array() ;
	lookseq.username = '' ;
	lookseq.got_samples = false ;
	lookseq.selected_samples = new Object ;
	lookseq.current_chromosome_length = 0 ;
	lookseq.tracks = new Array () ;
	lookseq.current_species = '' ;
	lookseq.mode = 'indel' ;
	lookseq.setMode ( 'indel' ) ;
	$('#species').html() ;
	$('.headerdetails').css({display:'none'}) ;
	$('#lookseq_header').css({display:''}) ;
	
	if ( !lookseq.first_init ) return ;
	lookseq.first_init = false ;
	
	$(window).resize(function() { lookseq.windowResizeHandler() ; } ) ;
	$('#main').mouseout ( function () { $('#position_on_chromosome').html('&nbsp;'); } ) ;
	$("#search_query").keyup(function(event){  if(event.keyCode == 13){  $("#search_button").click(); } });
} ,

scrollbarWidth : function () { 
	if ( lookseq.scrollbar_width != 0 ) return ; // Already done
    var div = $('<div style="width:50px;height:50px;overflow:hidden;position:absolute;top:-200px;left:-200px;"><div style="height:100px;"></div>'); 

    $('#lookseq_container').append(div); 
    var w1 = $('div', div).innerWidth(); 
    div.css('overflow-y', 'scroll'); 
    var w2 = $('div', div).innerWidth(); 
    $(div).remove(); 
    lookseq.scrollbar_width = w1 - w2 ; 
} ,


//_____________________________________________________________________________________________________________________________________________________________________
// EVENT HANDLERS

initializeFromParams : function () {
	if ( lookseq.loaded_from_params ) return ;

	lookseq.loaded_from_params = true ;

	var vars = lookseq.getUrlVars() ;
	if ( vars['samples'] === undefined ) return ; // No samples, no joy
	
	lookseq.initializing = true ;
	
	// Options
	if ( vars['options'] === undefined ) vars['options'] = '' ;
	
	$.each ( lookseq.display_options , function ( k , v ) {
		var id = 'display_' + v ;
		$('#'+id).attr('checked', false );
	} ) ;
	
	$.each ( vars['options'].split(',') , function ( k , v ) {
		var opt = v.split('|')[0] ;
		var val = v.split('|')[1] ;
		var id = 'display_' + opt ;
		$('#'+id).attr('checked', val==1?true:false );
	} ) ;
	
	if ( vars['maxdist'] !== undefined ) {
		$('#select_maxdist').val(vars['maxdist']) ;
	}
	
	// Samples
	var accid = '' ;
	$.each ( vars['samples'].split(',') , function ( k , v ) {
		var sample = v.split('|')[0] ;
		var alg = v.split('|')[1] ;
		var id = "sample_cb_"+sample+"_"+alg ;
		accid = parseInt ( $('#'+id).attr('accid') ) - 1 ;
		$('#'+id).attr('checked', true);
		lookseq.cb_changed ( $('#'+id) ) ;
	} ) ;
	
	// Position
	if ( vars['chr'] !== undefined && vars['from'] !== undefined && vars['to'] !== undefined ) {
		lookseq.showRegion ( vars['chr'] , vars['from'] , vars['to'] ) ;
	}
	
	lookseq.initializing = false ;
	if ( vars['mode'] !== undefined ) lookseq.setMode ( vars['mode'] ) ;
	lookseq.updatePositions ( true ) ;

	lookseq.load_annotation () ;
	if ( accid != '' ) {
		setTimeout ( function () { $('#sample_list').accordion ( 'activate' , accid ) ; } , 0 ) ;
	}
} ,

searchQuery : function () {
	var q = $('#search_query').val() ;
	$.getJSON ( lookseq.api_url , { action:'getannotation' , species:lookseq.current_species , query:q } , function ( data ) {
		var dh = '' ;
		$.each ( data.annotation , function ( typename , list ) {
			$.each ( list , function ( k , v ) {
				var f = parseInt(v.from) - 50 ;
				var t = parseInt(v.to) + 50 ;
				dh += "<tr>" ;
				if ( undefined === v.name ) dh += "<th/>" ;
				else dh += "<th nowrap>" + v.name + "</th>" ;
				dh += "<td>" + typename + "</td>" ;
				dh += "<td>" + v.chr + "</td>" ;
				dh += "<td>" + v.from + "</td>" ;
				dh += "<td>" + v.to + "</td>" ;
				dh += "<td style='text-align:center'>" + v.strand + "</td>" ;
				dh += "<td><a href='#' onclick='lookseq.showRegion(\"" + v.chr + "\","+f+","+t+",1);return false'>Show</a></td>" ;
				dh += "</tr>" ;
			} ) ;
		} ) ;

		if ( dh == '' ) {
			alert ( "Your search for \"" + q + "\" returned no results." ) ;
			return ;
		}
		
		dh = "<table border='1' cellspacing='0' cellpadding='1'><tr><th>Name</th><th>Type</th><th>Chr</th><th>From</th><th>To</th><th>Strand</th><th>Go</th></tr>" + dh + "</table>" ;
		dh = "<div style='max-height:500;overflow:auto;padding-right:20px'>" + dh + "</div>" ;
		$('#search_dialog').html ( dh ) ;
		
		$('#search_dialog').dialog( { 
			modal: true, 
			 width:'auto'
		});


	} ) ;
} ,

windowResizeHandler : function () {
	lookseq.adjustMainSize() ;
} ,

cb_changed : function ( o ) {
	var ox = $(o).attr('ox') ;
	var alg = $(o).attr('alg') ;
	var state = $(o).is(':checked') ;
	
	if ( state ) {
		if ( !lookseq.confirmReferenceChange ( ox ) ) {
			$(o).attr('checked', false);
			return ;
		}
		$(o).attr('checked', true);
	}
	
	if ( undefined === lookseq.selected_samples[ox] ) lookseq.selected_samples[ox] = new Object ;
	lookseq.selected_samples[ox][alg] = state ;
	
	if ( state ) {
		if ( lookseq.tracks.length == 0 ) {
			lookseq.appendScrollTrack() ;
			lookseq.appendPositionTrack() ;
			lookseq.appendAnnotationTrack() ;
		}
		lookseq.appendSampleTrack ( ox , alg ) ;
	} else {
		lookseq.removeSampleTrack ( ox , alg ) ;
	}
	
} ,

appendTypedTrack : function ( type ) {
	var t = new lookseq.Track ;
	lookseq.track_id++ ;
	t.id = lookseq.track_id ;
	t.div = 'track_' + t.id ;
	t.img_id = t.div + '_img' ;
	t.type = type ;
	lookseq.tracks.push ( t ) ;
	var h = "<div class='track' id='" + t.div + "'></div>" ;
	$('#trackstop').append ( h ) ;
	lookseq.recalculateTrackLayout() ;
	t.paint() ;
} ,

appendScrollTrack : function () {
	lookseq.appendTypedTrack ( 'scroll' ) ;
	$.each ( lookseq.tracks , function ( k , v ) {
		if ( v.type != 'scroll' ) return ;
		var id = v.div ;
//		console.log(id);
//		$('#'+id).html('<div id="slider"></div>') ;
//		$('#slider').slider ();
//		$('#slider .ui-slider-handle').css ( 'width' , '500px' ) ;
//		$('#slider').tinyscrollbar({ axis: 'x'});
	} ) ;
} ,

appendPositionTrack : function () {
	lookseq.appendTypedTrack ( 'position' ) ;
} ,

appendAnnotationTrack : function () {
	lookseq.appendTypedTrack ( 'annotation' ) ;
} ,

appendSampleTrack : function ( ox , alg ) {
	var t = new lookseq.Track ;
	lookseq.track_id++ ;
	t.id = lookseq.track_id ;
	t.div = 'track_' + t.id ;
	t.type = 'sample' ;
	t.sample = ox ;
	t.alg = alg ;
	t.note = t.sample + " [" + t.alg + "]" ;
	lookseq.tracks.push ( t ) ;
	var h = "<div class='track' id='" + t.div + "'></div>" ;
	$('#tracks').append ( h ) ;
	lookseq.recalculateTrackLayout() ;
	t.paint() ;
} ,

removeSampleTrack : function ( ox , alg ) {
	$.each ( lookseq.tracks , function ( k , v ) {
		if ( v === undefined ) return ;
		if ( v.type != 'sample' ) return ;
		if ( v.sample != ox ) return ;
		if ( v.alg != alg ) return ;
		$('#'+v.div).remove() ;
		lookseq.tracks.splice(k,1) ;
	} ) ;
	lookseq.recalculateTrackLayout() ;
} ,

recalculateTrackLayout : function () {
	// Pileup display fix
	$('#si_pileup').fadeTo('slow', ( lookseq.getTo() - lookseq.getFrom() > lookseq.max_pileup ? 0.1 : 1 ) ) ;

	lookseq.updateCurrentURL() ;

	for ( var r = 0 ; r <= 1 ; r++ ) {
		var h = 0 ;
		$.each ( lookseq.tracks , function ( k , v ) {
			if ( r == 0 && v.type == 'sample' ) return ;
			if ( r == 1 && v.type != 'sample' ) return ;
			$('#'+v.div).css ( { top : h+'px' , height : v.height() } ) ;
			v.fixed_top = h ;
			v.track_id = k ;
			h += v.height() + 1 ;
		} ) ;
		if ( r == 0 ) $('#trackstop').height(h) ;
	}

	var tt = $('#trackstop').height()+1 ;
	var nh2 = $('#main').height() - $('#trackstop').height() - 2 ;
	$('#trackstop').css ( { top : '0px' , left : '0px' } ) ;
	$('#tracks').css ( { top : tt+'px' , left : '0px' , height : nh2+'px' } ) ;
} ,

confirmReferenceChange : function ( ox ) {
	var nspecies = lookseq.samples[ox].species.toLowerCase() ;
	var cspecies = lookseq.getCurrentSpecies() ;
	
	// First time species is set
	if ( cspecies == '' ) {
		lookseq.showSpeciesData ( nspecies ) ;
		return true ;
	}
	
	// Same species
	if ( lookseq.species[cspecies].reference == lookseq.species[nspecies].reference ) return true ;

	// New reference?
	var ret = confirm ( "Selecting this sample will change the current reference sequence, and de-select all other samples. Continue?" ) ;

	if ( ret ) { // New reference!
		lookseq.mode = 'indel' ;
		lookseq.setMode ( 'indel' ) ;
		lookseq.selected_samples = new Object ;
		lookseq.tracks = new Array() ;
		lookseq.recalculateTrackLayout() ;
		$('#trackstop').html('') ;
		$('#tracks').html('') ;
		$('#sample_list input:checkbox').attr('checked',false);
		lookseq.showSpeciesData ( nspecies ) ;
	}
	
	return ret ;
} ,

onChromosomeChanged : function () {
	if ( lookseq.current_species == '' ) return ;
	if ( lookseq.species[lookseq.current_species] === undefined ) return ;
	
	var nchr = $('#chromosome').val() ;
	var ncl = lookseq.species[lookseq.current_species].chromosomes[nchr] ;
	if ( ncl > lookseq.hard_chrlen_limit ) ncl = lookseq.hard_chrlen_limit ;
	$('#from').val(1) ;
	$('#to').val(ncl) ;
	lookseq.updatePositions();
	lookseq.load_annotation () ;
} ,

updatePositions : function ( force ) {
	if ( lookseq.current_species == '' ) return ;
	if ( lookseq.species[lookseq.current_species] === undefined ) return ;

	var nfrom = parseInt ( $('#from').val() ) ;
	var nto = parseInt ( $('#to').val() ) ;
	var nchr = $('#chromosome').val() ;
//	console.log ( nchr ) ;
	var ncl = lookseq.species[lookseq.current_species].chromosomes[nchr] ;
	
	// Sanitize positions
	if ( nfrom < 1 ) nfrom = 1 ;
	if ( nfrom > ncl ) nfrom = ncl ;
	if ( nto < 1 ) nto = 1 ;
	if ( nto > ncl ) nto = ncl ;
	
	while ( nfrom + lookseq.min_size > nto ) {
		if ( nto == ncl ) nfrom-- ;
		else nto++ ;
	}
	
	$('#from').val(nfrom) ;
	$('#to').val(nto) ;
	
	if ( nchr == lookseq.chr && nfrom == lookseq.from && nto == lookseq.to && force === undefined ) {
		$.each ( lookseq.tracks , function ( k , v ) {
			if ( '' == v.img_id ) return ;
			$('#'+v.img_id).animate({left:'0px'},200);
		} ) ;
		return ; // No change
	}
	
	if ( nto - nfrom > lookseq.max_pileup && lookseq.mode == 'pileup' ) {
		lookseq.mode = 'indel' ;
		lookseq.setMode ( 'indel' ) ;
	}
	
	lookseq.chr = nchr ;
	lookseq.from = nfrom ;
	lookseq.to = nto ;
	lookseq.current_chromosome_length = ncl ;
	
	lookseq.recalculateTrackLayout() ;
	
	$.each ( lookseq.tracks , function ( k , v ) {
//		console.log ( k ) ;
		v.paint() ;
	} ) ;
	
} ,

onMoveLeft : function () {
	var w = lookseq.to - lookseq.from ;
	var nfrom = Math.floor ( lookseq.from - w / 4 ) ;
	if ( nfrom < 1 ) nfrom = 1 ;
	var nto = nfrom + w ;
	if ( nto > lookseq.current_chromosome_length ) return ;
	$('#from').val(nfrom) ;
	$('#to').val(nto) ;
	lookseq.updatePositions();
} ,

onMoveRight : function () {
	var w = lookseq.to - lookseq.from ;
	var nfrom = Math.floor ( lookseq.from + w / 4 ) ;
	if ( nfrom < 1 ) nfrom = 1 ;
	var nto = nfrom + w ;
	if ( nto > lookseq.current_chromosome_length ) return ;
	$('#from').val(nfrom) ;
	$('#to').val(nto) ;
	lookseq.updatePositions();
} ,

onZoomSize : function ( size ) {
	var middle = ( Math.floor(lookseq.from) + Math.floor(lookseq.to) ) / 2 ;
	var nfrom = Math.floor ( middle - size / 2 ) ;
	var nto = nfrom + Math.floor ( size ) ;
	$('#from').val(nfrom) ;
	$('#to').val(nto) ;
	lookseq.updatePositions();
} ,

onZoomIn : function () {
	var w = lookseq.to - lookseq.from ;
	lookseq.onZoomSize ( w / 2 ) ;
} ,

onZoomOut : function () {
	var w = lookseq.to - lookseq.from ;
	lookseq.onZoomSize ( w * 2 ) ;
} ,

onZoom121 : function () {
	lookseq.onZoomSize ( Math.floor ( lookseq.getMainWidth() / 8 ) ) ;
} ,

onZoomFull : function () {
	$('#from').val(1) ;
	$('#to').val(lookseq.current_chromosome_length) ;
	lookseq.updatePositions();
} ,

event2pos : function ( e , o ) {
	var x = e.pageX - $(o).offset().left ;
	var tid = $(o).attr('track_id') ;
	if ( lookseq.tracks[tid] === undefined ) return 0 ;
	var pos = Math.floor ( lookseq.from + ( lookseq.to - lookseq.from ) * x / lookseq.tracks[tid].width() ) ;
	return pos ;
} ,

trackimageDblClick : function ( e , o ) {
	var pos = lookseq.event2pos ( e , o ) ;
	var w = lookseq.to - lookseq.from ;
	var nw = w / 10 ;
	var nfrom = Math.floor ( pos - nw ) ;
	var nto = Math.floor ( nfrom + nw * 2 ) ;
	$('#from').val(nfrom) ;
	$('#to').val(nto) ;
	lookseq.updatePositions();
	return false ;
} ,

trackimageMouseMove : function ( e , o ) {
	var pos = lookseq.event2pos ( e , o ) ;
	$('#position_on_chromosome').html ( lookseq.formatPosition ( pos ) ) ;
	return false ;
} ,

showRegion : function ( nchr , nfrom , nto , extra ) {
	if ( extra == 1 ) $('#search_dialog').dialog('close');
	
	if ( nchr != lookseq.chr ) {
		$('#chromosome').val(nchr) ;
		lookseq.onChromosomeChanged() ;
	}
	$('#from').val(nfrom) ;
	$('#to').val(nto) ;
	lookseq.updatePositions();
	return false ;
} ,

formatPosition : function ( p ) {
	p += '' ;
	var ret = '' ;
	while ( p.length > 3 ) {
		ret = ',' + p.substr(-3) + ret ;
		p = p.substr(0,p.length-3);
	}
	ret = p + ret ;
	return ret ;
} ,

updateCurrentURL : function () {
	var url = lookseq.getCurrentURL() ;
	$('#current_url').attr ( 'href' , url ) ;
} ,

//_____________________________________________________________________________________________________________________________________________________________________
// DATA RETRIEVAL

getCurrentURL : function () {
	var url = lookseq.mybaseurl+'?' ;
	var parts = new Array() ;

	// Misc
	parts.push ( "chr=" + escape(lookseq.getChr()) ) ;
	parts.push ( "from=" + escape(lookseq.getFrom()) ) ;
	parts.push ( "to=" + escape(lookseq.getTo()) ) ;
	parts.push ( "mode=" + escape(lookseq.mode) ) ;
	
	// Samples
	var samples = new Array() ;
	$('#sample_list input[type=checkbox]:checked').each ( function ( k , v ) {
		samples.push ( $(v).attr('ox') + "|" + $(v).attr('alg') ) ;
	} ) ;
	parts.push ( "samples=" + escape(samples.sort().join(',')) ) ;
	
	// Display options
	var options = new Array() ;
	$.each ( lookseq.display_options , function ( k , v ) {
		options.push ( v + '|' + ( $('#display_'+v).is(':checked') ? 1 : 0 ) ) ;
	} ) ;
	parts.push ( "options=" + escape(options.join(',')) ) ;
	parts.push ( "maxdist=" + $('#select_maxdist').val() ) ;

	url += parts.join ( '&' ) ;
	
	return url ;
} ,

getUrlVars : function () {
    var vars = [], hash;
    var hashes = window.location.href.slice(window.location.href.indexOf('?') + 1).split('&');
    for(var i = 0; i < hashes.length; i++)
    {
        hash = hashes[i].split('=');
        vars.push(hash[0]);
        vars[hash[0]] = unescape ( hash[1] ) ;
    }
    return vars;
} ,

getMainWidth : function () {
	return $('#tracks').innerWidth() - lookseq.scrollbar_width ;
} ,

getFrom : function () {
	var ret = lookseq.from ;
//	var ret = $('#from').val() ;
//	if ( ret === undefined ) return 0 ;
	return parseInt(ret) ;
} ,

getTo : function () {
	var ret = lookseq.to ;
//	var ret = $('#to').val() ;
//	if ( ret === undefined ) return 0 ;
	return parseInt(ret) ;
} ,

getChr : function () {
	var ret = $('#chromosome').val() ;
	if ( ret === undefined ) return 0 ;
	return ret ;
} ,

getCurrentSpecies : function () {
	var ret = '' ;
	$.each ( lookseq.selected_samples , function ( o , algs ) {
		if ( ret != '' ) return ;
		$.each ( algs , function ( a , x ) {
			if ( !x ) return ;
			ret = lookseq.samples[o].species.toLowerCase() ;
		} ) ;
	} ) ;
	return ret ;
} ,

getSpeciesData : function () {
	$.getJSON ( lookseq.api_url , { action : 'get_species_data' } , function ( data ) {
		lookseq.species = data.species ;
		lookseq.adjustMainSize() ;
		lookseq.test_login() ;
	} ) ;
} ,

getMySamples : function ( show_later ) {
	if ( lookseq.got_samples ) return ;

	$('#sample_list').html("<i>loading samples...</i>") ;
	$.getJSON ( lookseq.api_url , { action : 'mysamples' } , function ( data ) {
		if ( data.error == 'OK' ) {
			lookseq.samples = data.samples ;
			lookseq.projects = data.projects ;
			lookseq.algorithms = data.algorithms ;
			lookseq.got_samples = true ;
			if ( show_later ) lookseq.showSamples() ;
		}
	} ) ;
} ,

formatMBP : function ( v ) {
	if ( v > 100000 ) {
		v = Math.floor ( v / 100000 ) / 10 ;
		return v + "mbp" ;
	}
	if ( v > 100 ) {
		v = Math.floor ( v / 100 ) / 10 ;
		return v + "kbp" ;
	}
	return v + "bp" ;
} ,


//_____________________________________________________________________________________________________________________________________________________________________
// DISPLAY

repaintTracks : function () {
	$.each ( lookseq.tracks , function ( k ,v ) {
//		if ( v.type != 'sample' ) return ;
		v.paint() ;
	} ) ;
} ,

setMode : function ( mode ) {
	if ( mode == 'pileup' && lookseq.max_pileup < lookseq.getTo() - lookseq.getFrom() ) return ;
	
	$.each ( ['indel','pileup','coverage'] , function ( k , v ) {
		$('#si_'+v).css ( { 'border-style' : ( mode == v ? 'groove' : 'none' ) , padding : ( mode == v ? '0px' : '3px' ) } ) ;
	} ) ;

	if ( mode == 'pileup' ) $('#display_basequal').removeAttr('disabled') ;
	else $('#display_basequal').attr('disabled','disabled');

	if ( mode == 'indel' ) {
		$('#display_pairlinks').removeAttr('disabled') ;
		$('#display_orientation').removeAttr('disabled') ;
		$('#select_maxdist').removeAttr('disabled') ;
	} else {
		$('#display_pairlinks').attr('disabled','disabled');
		$('#display_orientation').attr('disabled','disabled');
		$('#select_maxdist').attr('disabled','disabled');
	}

	if ( mode == lookseq.mode ) return ;
	
	lookseq.mode = mode ;
	lookseq.recalculateTrackLayout() ;
	$.each ( lookseq.tracks , function ( k , v ) {
		v.paint();
	} ) ;
} ,

showNewPosition : function ( x ) {
	var w = lookseq.getMainWidth();
	var posw = lookseq.to - lookseq.from ;
	var nfrom = lookseq.from - Math.floor ( x * posw / w ) ;
	if ( nfrom < 1 ) nfrom = 1 ;
	var nto = nfrom + posw ;
	while ( nto > lookseq.current_chromosome_length && nfrom > 1 ) {
		nfrom-- ;
		nto-- ;
	}

	$('#from').val(nfrom) ;
	$('#to').val(nto) ;
} ,

showSpeciesData : function ( spec ) {
	var h = "<i>" + spec.ucFirst() + "</i><br/>" ;
	
	h += "<table cellspacing=0 cellpadding=0>" ;
	h += "<tr><td>Chr</td><td><select id='chromosome' onchange='lookseq.onChromosomeChanged()'>" ;
	var first = true ;
	var ck = new Array() ;
	$.each ( lookseq.species[spec].chromosomes , function ( k , v ) {
		ck.push ( k ) ;
	} ) ;
	ck.sort() ;
	$.each ( ck , function ( x , chr ) {
		var v = lookseq.species[spec].chromosomes[chr] ;
		if ( v < 1000 ) return ; // HARD CUTOFF MINIMUM CHROMOSOME LENGTH
		h += "<option value='" + chr + "'" ;
		if ( first ) {
			first = false ;
			lookseq.current_chromosome_length = v ;
			lookseq.chr = chr ;
			h += ' selected' ;
		}
		h += ">" + chr + " (" + lookseq.formatMBP ( v ) + ")</option>" ;
	} ) ;
	h += "</select></td></tr>" ;
	h += "<tr><td>From</td><td><input id='from' type='number' value='1' size='15' /></td></tr>" ;

	var ncl = lookseq.current_chromosome_length ;
	if ( ncl > lookseq.hard_chrlen_limit ) ncl = lookseq.hard_chrlen_limit ;

	h += "<tr><td>To</td><td><input id='to' type='number' value='" + ncl + "' size='15' /></td></tr>" ;
	h += "<tr><td/><td><input type='button' value='Update' id='updatePositions' onclick='lookseq.updatePositions()' />" 
	h += "</td></tr>" ;
	h += "</table>" ;
	
	$('#species').html ( h ) ;
	$('.headerdetails').css({display:''}) ;
	$('#lookseq_header').css({display:'none'}) ;

	$("#from").keyup(function(event){  if(event.keyCode == 13){  $("#updatePositions").click(); } });
	$("#to").keyup(function(event){  if(event.keyCode == 13){  $("#updatePositions").click(); } });
	
	lookseq.from = 1 ;
	lookseq.to = ncl ;
	lookseq.current_species = spec ;
	
	lookseq.load_annotation () ;
} ,

load_annotation : function () {
	var spec = lookseq.current_species ;
	$('#load_anno').html('<i>Loading annotation...</i>') ;
	$('.annotation_marker').remove();
	lookseq.annotation = new Object ;
	
	$.getJSON ( lookseq.api_url , { action : 'getannotation' , species : spec , chr : lookseq.chr } , function ( data ) {
		$('#load_anno').html('') ;
		lookseq.annotation = data.annotation ;
		$.each ( lookseq.tracks , function ( k , v ) {
			if ( v.type == 'annotation' ) v.paint() ;
		} ) ;
	} ) ;
} ,

showSamples : function () {
	if ( !lookseq.got_samples ) {
		lookseq.getMySamples(true);
		return ;
	}

	var sl = '' ;

	// Sorting projects by code
	var pk = new Array() ;
	var pc2pid = new Array() ;
	$.each ( lookseq.projects , function (pid,pd) {
		var code = pd.code ;
		if ( code == '' ) code = '#' + pid ;
		pk.push ( code ) ;
		pc2pid[code] = pid ;
	} ) ;
	pk = pk.sort();
	
	var tabid = 1 ;
	$.each ( pk , function ( k1 , v1 ) {
		var pid = pc2pid[v1] ;
		var pd = lookseq.projects[pid] ;
		var code = pd.code ;
		if ( code == '' ) code = '#' + pid ;
		var has_sample = false ;
		
		var sb = '' ;
		sb += "<h3><a href='#'>" + code + "</a></h3><div>" ;
		
		var countries = '' ;
		var species = '' ;
		$.each ( pd.countries , function (k,v) { if ( countries != '' ) countries += ', ' ; countries += k ; } ) ;
		$.each ( pd.species , function (k,v) { if ( species != '' ) species += ', ' ; species += k ; } ) ;
		sb += "<div style='font-size:7pt;text-align:center;border-bottom:1px solid #DDDDDD' id='acc_tab_"+k1+"'><i>" + species + "</i><br/>" + countries + "</div>" ;

		if ( pd.desc != '' && lookseq.logged_in ) sb += "<div style='font-size:7pt'>" + pd.desc + "</div>" ;
		
		sb += "<table style='font-size:8pt' cellpadding=1 cellspacing=0><tr><td/>" ;
		$.each ( lookseq.algorithms , function ( k2 , v2 ) {
			sb += "<td class='sample_cell'>" + v2 + "</td>" ;
		} ) ;
		sb += "</tr>" ;
		var sk = new Array() ;
		$.each ( lookseq.samples , function (ox,sd) {
			sk.push ( ox ) ;
		} ) ;
		sk.sort();
		$.each ( sk , function (dummy,ox) {
			if ( lookseq.samples[ox].pid != pid ) return ;
			if ( lookseq.samples[ox].species === undefined ) return ;
			if ( undefined === lookseq.samples[ox].alg ) return ;
			var title = '' ;
			if ( undefined !== lookseq.samples[ox].country ) title = " title='" + lookseq.samples[ox].country + "'" ;
			sb += "<tr><td nowrap " + title + " style='padding-right:3px;'>" + ox + "</td>" ;
			$.each ( lookseq.algorithms , function ( k2 , v2 ) {
				if ( undefined === lookseq.samples[ox].alg[v2] ) sb += "<td class='sample_cell'>&mdash;</td>" ;
				else sb += "<td class='sample_cell'><input accid='"+tabid+"' type='checkbox' id='sample_cb_"+ox+"_"+v2+"' ox='"+ox+"' alg='"+v2+"' onclick='lookseq.cb_changed(this)' /></td>" ;
			} ) ;
//			sb += "<td style='text-align:center'><input type='checkbox' ox='"+ox+"' alg='som' onclick='lookseq.cb_changed(this)' /></td>" ;
			sb += "</tr>" ;
			has_sample = true ;
		} ) ;
		sb += "</table></div>" ;
		
		if ( has_sample ) {
			sl += sb ;
			tabid++ ;
		}
	} ) ;
//alert(sl);
//	alert(JSON.stringify(lookseq.samples));

	$('#sample_list').html(sl) ;
	$('#sample_list').accordion({autoHeight:false,collapsible: true,active:false}) ;
	$('#sample_list').css({fontSize:'8pt',maxHeight:$('#sidebar').height()*2});
	$('#sample_list .ui-accordion-content').css({padding:'1px'});

	lookseq.initializeFromParams() ;
} ,

adjustMainSize : function () {
	var sbw = $('#sidebar').width() ;
	var foh = $('#footer').height() ;
	if ( !$('#sidebar').is(':visible') ) sbw = 0 ; else sbw += 30 ;
	if ( !$('#footer').is(':visible') ) foh = 20 ; else foh += 30 ;

	var nw = $('#lookseq_container').width() - sbw ;
	var nh = $('#lookseq_container').height() - 10*2 - foh - 105 ;
	
	$('#sidebar').css({maxHeight:nh,overflow:'auto'});
	
	if ( $('#main').width() == nw && $('#main').height() == nh && ( $('#sidebar').height() == nh || !$('#sidebar').is(':visible') ) ) return ;
	
	var t = 500 ;
	if ( $('#main').width() == 0 ) t = 0 ;
	
	$('#main').animate({
			width : nw,
            height: nh
        }, t);
	
	var nh2 = nh - $('#trackstop').height() ;
	
	$('#trackstop').css ( { width : nw , maxWidth : nw } ) ;
	$('#tracks').css ( { width : nw , maxWidth : nw , height : nh2 , maxHeight : nh2 } ) ;
	$('#main').css ( { maxWidth : nw , maxHeight : nh } ) ;
	
	if ( $('#sidebar').is(':visible') ) $('#sidebar').animate ( { height : nh } , t ) ;
} ,

//_____________________________________________________________________________________________________________________________________________________________________
// LOGIN-RELATED

showLoginDialog : function () {
	var vars = lookseq.getUrlVars() ;
	if ( vars['skip_login'] == 1 ) { // Public demo mode
		lookseq.showSamples();
		return ;
	}

	$("#username").keyup(function(event){  if(event.keyCode == 13){  $("#login").click(); } });
	$("#email").keyup(function(event){  if(event.keyCode == 13){  $("#login").click(); } });
	$("#password").keyup(function(event){  if(event.keyCode == 13){  $("#login").click(); } });
	
	$('#login_dialog').dialog( { 
		modal: true, 
		width:'auto',
		open:function(){
			$('#email').focus();
		},
		close:function() {
			if ( !lookseq.logged_in ) {
				lookseq.showSamples();
			}
		} 
	});
} ,


showLoginStatus : function () {
	if ( lookseq.logged_in ) {
		$('#login_status').html('Wellcome, '+lookseq.username+'!<br/><a href="/">Main page</a> | <a href="#" onclick="lookseq.logout();return false">Log out</a>');
	} else {
		$('#login_status').html('Not logged in<br/><a href="/">Main page</a> | <a href="#" onclick="lookseq.showLoginDialog()">Log in</a>');
	}
} ,

test_login : function () {
	if ( lookseq.logged_in ) return ;
	$.getJSON ( lookseq.api_url , { action : 'test_login' } , function ( data ) {
		if ( data.error == 'OK' ) {
			lookseq.showSamples();
			lookseq.logged_in = true ;
			lookseq.username = data.username ;
		} else {
			lookseq.showLoginDialog () ;
		}
		lookseq.showLoginStatus();
	} ) ;
} ,

login : function () {
	lookseq.email = $('#email').val() ;
	var password = $('#password').val() ;
	$.getJSON ( lookseq.api_url , { action : 'login' , email : lookseq.email , password : password } , function ( data ) {
		if ( data.error == 'OK' ) {
			lookseq.username = data.username ;
			lookseq.got_samples = false ;
			lookseq.logged_in = true ;
			lookseq.showLoginStatus();
			lookseq.showSamples() ;
			$('#login_dialog').dialog('close');
		} else {
			$('#login_warning').html(data.error).show();
		}
	} ) ;
	return false ;
} ,

logout : function () {
	$.getJSON ( lookseq.api_url , { action : 'logout' } , function ( data ) {
		lookseq.clear() ;
		lookseq.logged_in = false ;
		lookseq.test_login() ;
	} ) ;
	return false ;
} ,

//_____________________________________________________________________________________________________________________________________________________________________
// FIN


	the_end_of_the_object : '' // To avoid worrying about the final comma...

} ; // END LOOKSEQ OBJECT

