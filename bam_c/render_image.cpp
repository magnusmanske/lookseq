#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>

//int global_test_flag = 0 ;

#define PNG_DEBUG 3
#include <png.h>

#include "sam.h"  
#include "faidx.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std ;

#define ARROW_LENGTH 3
#define SATURATION_THRESHOLD 64

#define NUMBER_WIDTH 4
#define NUMBER_HEIGHT 6
#define NUMBER_ALIGN_HCENTER 1
#define NUMBER_ALIGN_RIGHT 2
#define NUMBER_ALIGN_VCENTER 4

#define PC_EMPTY 0
#define PC_REFERENCE 1
#define PC_SNP 2
#define PC_SINGLE 4
#define PC_INSERTION 8
#define PC_DELETION 16

typedef int64_t postype ;

typedef struct {  
	int beg, end;  
	samfile_t *in;  
} tmpstruct_t;  

typedef vector <postype> TVI ;

vector <TVI> lcd_chars ;
char bam2char[255] ;

class Tbam2png ;

class Tpileupchar {
	public :
	Tpileupchar ( char _c = ' ' , char _type = PC_EMPTY , char _q = -1 ) { c = _c ; type = _type ; quality = _q ; }
	char c ;
	char type ;
	char quality ;
} ;

typedef vector <Tpileupchar> VPC ;

class Tbam_draw {
	public :
	Tbam_draw ( Tbam2png *_base ) ;
	virtual void set_range () ;
	void set_vertical_range ( int a , int b ) ;
	virtual void draw_single_read ( const bam1_t *b, void *data ) {} ;
	virtual void draw_paired_read ( const bam1_t *b, void *data) {} ;
	virtual void merge_into_png ( unsigned char *p , int red , int green , int blue ) ;
	virtual void merge_all () {} ;
	virtual postype get_start () { return start ; }
	virtual int get_neccessary_height () ;

	postype single_offset ;

	protected :
	virtual void render_char ( char ch , int x , int y ) ;
	virtual void fill_rect ( int x1 , int y1 , int x2 , int y2 , int r , int g , int b ) ;
	virtual void render_number ( int number , int x , int y , int align = 0 ) ;
	virtual void draw_axis () {} ;
//	virtual postype pos2mem ( postype x , postype y , bool outer_right = false ) { return -1 ; } ;
	virtual void paint_single_read ( unsigned char *bucket , const bam1_t *b, void *data , int y ) {} ;
	virtual void paint_single_read_cigar ( unsigned char *bucket , const bam1_t *b, void *data , int y ) {} ;
	virtual void set_pixel ( int x , int y ) ;
	virtual int opt_axis ( int from , int to , int max_steps ) ;
	unsigned char *alloc_p () ;
	
	Tbam2png *base ;
	postype max_mem ;
	postype left , right , top , bottom , w , h ;
	postype start , end , size ;
	postype vstart , vend , vsize ;
	int basecol_red , basecol_green , basecol_blue ;
} ;

class Tbam_draw_paired : public Tbam_draw {
	public :
	Tbam_draw_paired ( Tbam2png *_base ) ;
	virtual void draw_single_read ( const bam1_t *b, void *data) ;
	virtual void draw_paired_read ( const bam1_t *b, void *data) ;
	virtual void merge_all () ;

	protected :
	virtual void draw_axis () ;
	virtual void paint_single_read ( unsigned char *bucket , const bam1_t *b, void *data , int y ) ;
	virtual void paint_single_read_cigar ( unsigned char *bucket , const bam1_t *b, void *data , int y ) ;
	virtual void hline ( unsigned char *bucket , postype from , postype to , postype y , bool count_in_coverage = true ) ;
	virtual void paint_arrow ( unsigned char *bucket , postype from , postype to , postype y , bool reverse ) ;
	inline postype pos2mem ( postype x , postype y , bool outer_right = false ) ;
	unsigned char *psingle , *ppaired , *psnps , *pconn , *pfaceaway1 , *pfaceaway2 , *pinversions1 , *pinversions2 , *plowq , *p_cigar_insertion , *p_cigar_deletion ;
} ;

class Tbam_draw_coverage : public Tbam_draw_paired {
	public :
	Tbam_draw_coverage ( Tbam2png *_base ) ;
	virtual void set_range () ;
	virtual int get_neccessary_height () ;
	virtual void merge_all () ;
	
	protected:
	virtual void draw_axis () ;
	virtual void hline ( unsigned char *bucket , postype from , postype to , postype y , bool count_in_coverage = true ) ;
	
	int *cov ;
} ;

class Tbam_draw_pileup : public Tbam_draw_paired {
	public :
	Tbam_draw_pileup ( Tbam2png *_base ) ;
	virtual void set_range () ;
	virtual int get_neccessary_height () ;
	virtual void merge_all () ;
	virtual string get_text_rendering () ;
	
	protected:
	virtual void draw_axis () ;
	virtual void hline ( unsigned char *bucket , postype from , postype to , postype y , bool count_in_coverage = true ) ;
	virtual void paint_single_read ( unsigned char *bucket , const bam1_t *b, void *data , int y ) ;
	void get_pileup_read_start ( const bam1_t *b , int &row , int &pos , bool has_insertion = false ) ;
	
	vector <VPC> pile ;
	bool render_as_text ;
} ;

class Tbam2png {
	public :
	Tbam2png () {} ;
	void init ( string _bamfile , string _region , string _png_out , int _mapq = 0 ) ;
	void set_options ( string options ) ;
	
	// BAM methods
	static int fetch_func(const bam1_t *b, void *data) ;
	void read_bam_file () ;

	// PNG methods
	void create_png () ;
	void write_png_file(char* file_name) ;
	void process_file(void) ;
	void set_image_size ( postype w , postype h ) ;
	inline postype get_width() { return width ; }
	inline postype get_height() { return height ; }
	inline png_bytep * get_row_pointers () { return row_pointers ; }
	inline postype get_start () { return tmp.beg ; }
	inline postype get_end () { return tmp.end ; }

	char *refseq ;
	string refseq_file ;
	
	Tbam_draw *draw ;
	bool o_single , o_pairs , o_arrows , o_snps , o_faceaway , o_inversions , o_linkpairs , o_colordepth ;
	bool o_noscale , o_readqual , o_text ;
	int total_snps , total_reads ;
	int highlight_from , highlight_to ;
	bool use_highlight ;
	
	private :

	void abort_(const char * s, ...) ;

	// BAM variables
	string bam_file , png_file , region ;
	tmpstruct_t tmp;
	int mapq ;
	
	// PNG variables
	postype width, height;
	png_byte color_type;
	png_byte bit_depth;
	png_structp png_ptr;
	png_infop info_ptr;
	int number_of_passes;
	png_bytep * row_pointers;
	
} ;

Tbam2png *b2p ; // Neccessary for hack around BAM needing static function


//////////////////////////////////////////////////////////////////////////////////////
// Tbam_draw
Tbam_draw::Tbam_draw ( Tbam2png *_base ) {
	base = _base ;
	base->draw = this ;
	max_mem = 0 ;
	vstart = 0 ;
	vend = 1000 ;
	single_offset = 0 ;
	basecol_red = basecol_green = basecol_blue = 0 ;
}

void Tbam_draw::set_vertical_range ( int a , int b ) {
	vstart = a ;
	vend = b ;
}

int Tbam_draw::get_neccessary_height () {
	return base->get_height() ;
}


int Tbam_draw::opt_axis ( int from , int to , int max_steps ) {
	int step = 1 ;
	if ( max_steps <= 0 ) return step ;
	while ( 1 ) {
//		printf ( "%d\t%d\t%d\n" , step , max_steps , to - from ) ;
		if ( ( step * 5 ) * max_steps > to - from ) break ;
		step *= 5 ;
		if ( ( step * 2 ) * max_steps > to - from ) break ;
		step *= 2 ;
	}
	return step ;
}

void Tbam_draw::fill_rect ( int x1 , int y1 , int x2 , int y2 , int r , int g , int b ) {
	if ( x1 > x2 ) { int x = x1 ; x1 = x2 ; x2 = x ; }
	if ( y1 > y2 ) { int y = y1 ; y1 = y2 ; y2 = y ; }
	
	if ( x1 < 0 ) x1 = 0 ;
	if ( y1 < 0 ) y1 = 0 ;
	if ( x2 >= base->get_width()) x2 = base->get_width() - 1 ;
	if ( y2 >= base->get_height() ) y2 = base->get_height() - 1 ;
	if ( x1 > x2 || y1 > y2 ) return ;

	int r2 = basecol_red ;
	int g2 = basecol_green ;
	int b2 = basecol_blue ;
	
	basecol_red = r ;
	basecol_green = g ;
	basecol_blue = b ;
	
	for ( int x = x1 ; x <= x2 ; x++ ) {
		for ( int y = y1 ; y <= y2 ; y++ ) {
			set_pixel ( x , y ) ;
		}
	}
	
	basecol_red = r2 ;
	basecol_green = r2 ;
	basecol_blue = r2 ;
}

void Tbam_draw::render_char ( char ch , int x , int y ) {
//	if ( digit < 0 || digit > 9 ) return ;
	for ( int a = 0 ; a < lcd_chars[ch].size() ; a++ ) {
		int b = lcd_chars[ch][a] ;
		int ox = x , oy = y ;
		switch ( b ) {
			case 0 : break ;
			case 1 : break ;
			case 2 : ox += NUMBER_WIDTH ; break ;
			case 3 : oy += NUMBER_HEIGHT/2 ; break ;
			case 4 : oy += NUMBER_HEIGHT/2 ; break ;
			case 5 : ox += NUMBER_WIDTH ; oy += NUMBER_HEIGHT/2 ; break ;
			case 6 : oy += NUMBER_HEIGHT ; break ;
			case 7 : ox += NUMBER_WIDTH/2 ; b = 1 ; break ;
			case 8 : ox += NUMBER_WIDTH/2 ; oy += NUMBER_HEIGHT/2 ; b = 1 ; break ;
		}
		if ( b == 0 || b == 3 || b == 6 ) {
			for ( int a = 1 ; a < NUMBER_WIDTH ; a++ ) {
				set_pixel ( ox + a , oy ) ;
			}
		} else {
			for ( int a = 1 ; a < NUMBER_HEIGHT / 2 ; a++ ) {
				set_pixel ( ox , oy + a ) ;
			}
		}
	}
}

void Tbam_draw::render_number ( int number , int x , int y , int align ) {
	char num[100] ;
	sprintf ( num , "%d" , number ) ;
	
	int dw = NUMBER_WIDTH + 2 ;
	int len = strlen ( num ) ;
	
	if ( align & NUMBER_ALIGN_HCENTER ) x -= len * dw / 2 ;
	else if ( align & NUMBER_ALIGN_RIGHT ) x -= len * dw ;
	
	if ( align & NUMBER_ALIGN_VCENTER ) y -= NUMBER_HEIGHT / 2 ;
	
	for ( char *c = num ; *c ; c++ ) {
		render_char ( *c , x , y ) ;
		x += dw ;
	}
}

void Tbam_draw::set_pixel ( int x , int y ) {
	int bw = base->get_width() ;
	if ( x < 0 || y < 0 || x >= bw || y >= base->get_height() ) return ;
	png_bytep * row_pointers = base->get_row_pointers () ;
	png_byte* row = row_pointers[y];
	png_byte* ptr = &(row[x*4]);
	ptr[0] = basecol_red ;
	ptr[1] = basecol_green ;
	ptr[2] = basecol_blue ;
	ptr[3] = 255 ;
}

void Tbam_draw::set_range () {
	start = base->get_start() + 1 ;
	end = base->get_end() + 1 ;
	size = end - start ;
	vsize = vend - vstart ;
}

unsigned char *Tbam_draw::alloc_p () {
	unsigned char *p = (unsigned char*) malloc ( max_mem ) ;
	for ( int a = 0 ; a < max_mem ; a++ ) *(p+a) = (unsigned char) 0 ;
	return p ;
}

void Tbam_draw::merge_into_png ( unsigned char *p , int red , int green , int blue ) {
	png_bytep * row_pointers = base->get_row_pointers () ;
	
	unsigned char max = 0 ;
	if ( base->o_colordepth ) {
		for ( int y = 0 ; y < h ; y++) {
			for ( int x = 0 ; x < w ; x++) {
				unsigned char c = *(p+x+y*w) ;
				if ( max < c ) max = c ;
			}
		}
	}
	
	for ( int y = 0 ; y < h ; y++) {
		png_byte* row = row_pointers[y+top];
		for ( int x = 0 ; x < w ; x++) {
			unsigned char c = *(p+x+y*w) ;
			if ( c == 0 ) continue ;
			png_byte* ptr = &(row[(x+left)*4]);
			
			if ( base->o_colordepth ) c = (int) 128 - c * 128 / max ;

			ptr[0] = red ? red : c ;
			ptr[1] = green ? green : c ;
			ptr[2] = blue ? blue : c ;
			ptr[3] = 255 ;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////
// Tbam_draw_paired

Tbam_draw_paired::Tbam_draw_paired ( Tbam2png *_base ) : Tbam_draw ( _base ) {
	left = 0 ;
	right = 0 ;
	top = 0 ;
	bottom = 30 ;
	w = base->get_width() - left - right ;
	h = base->get_height() - top - bottom ;
	max_mem = w * h ;
	psingle = alloc_p () ;
	ppaired = alloc_p () ;
	plowq = alloc_p () ;
	psnps = alloc_p () ;
	p_cigar_insertion = alloc_p () ;
	p_cigar_deletion = alloc_p () ;
	pconn = alloc_p () ;
	pfaceaway1 = alloc_p () ;
	pfaceaway2 = alloc_p () ;
	pinversions1 = alloc_p () ;
	pinversions2 = alloc_p () ;
}

void Tbam_draw_paired::draw_single_read ( const bam1_t *b, void *data ) {
	paint_single_read ( psingle , b , data , b->core.pos % single_offset ) ;
}

void Tbam_draw_paired::draw_paired_read ( const bam1_t *b, void *data) {
	postype isize = b->core.isize ;
	unsigned char *bucket = ppaired ;
	

	// Faceaway
	bool is_faceaway = false ;
	bool needs_to_become_inversion = false ;
	if ( isize > 0 && b->core.flag & BAM_FREVERSE ) is_faceaway = true ;
	else if ( isize < 0 && 0 == ( b->core.flag & BAM_FREVERSE ) ) is_faceaway = true ;
	if ( is_faceaway ) {
		if ( base->o_faceaway ) {
			if ( isize < 0 ) bucket = pfaceaway1 ;
			else bucket = pfaceaway2 ;
		} else if ( base->o_inversions ) {
			is_faceaway = false ;
			needs_to_become_inversion = true ;
		} else return ;
	}

	// Inversion
	bool is_inversion = false ;
	if ( ( b->core.flag & BAM_FREVERSE ) == ( b->core.flag & BAM_FMREVERSE ) ) {
//		cout << "Inversion!\n" ;
		is_inversion = true ;
		if ( base->o_inversions ) {
			if ( isize < 0 ) bucket = pinversions1 ;
			else bucket = pinversions2 ;
		} else return ;
	} else if ( needs_to_become_inversion ) {
		return ;
	}
	
	if ( !is_inversion && !is_faceaway && !base->o_pairs ) return ; // No drawing normal pairs if unwanted
	
	if ( bucket == ppaired && b->core.qual == 0 ) bucket = plowq ;

	// Paint read
	paint_single_read ( bucket , b , data , abs(isize) ) ;

	if ( !base->o_linkpairs ) return ;

	// Paint linkage
	if ( isize > 0 ) {
		postype x = b->core.pos ;
		hline ( pconn , x , x + isize , abs(isize) , false ) ;
	} else {
		postype x = b->core.pos + b->core.l_qseq - 1 ;
		hline ( pconn , x + isize , x , abs(isize) , false ) ;
	}
}

void Tbam_draw_paired::paint_single_read ( unsigned char *bucket , const bam1_t *b, void *data , int y ) {
	if ( y < vstart || y >= vend ) return ;

	if ( b->core.n_cigar > 1 ) {
		paint_single_read_cigar ( bucket , b , data , y ) ;
		return ;
	}
	
	postype from = b->core.pos ;
	postype to = from + b->core.l_qseq - 1 ;
	

	hline ( bucket , from , to , y ) ;
	
	// SNPs
	if ( base->o_snps && base->refseq ) {
		uint8_t *s = bam1_seq ( b ) ;
		postype p = b->core.pos + 1 ;
		for ( postype a = 0 ; a < b->core.l_qseq ; a++ , p++ ) {
			if ( p < start ) continue ;
			if ( p >= end ) break ;

			if ( *(base->refseq+p-start) != bam2char[bam1_seqi(s,a)] ) {
				hline ( psnps , p , p , y , false ) ;
			}
		}
	}

	// Arrows
	if ( base->o_arrows ) {
		paint_arrow ( bucket , from , to , y , b->core.flag & BAM_FREVERSE ) ;
	}
}

void Tbam_draw_paired::paint_single_read_cigar ( unsigned char *bucket , const bam1_t *b, void *data , int y ) {
	uint32_t *cigdata = bam1_cigar(b) ;
	postype p = b->core.pos ;
	int rp = 0 ;
	uint8_t *s = bam1_seq ( b ) ;
	
	for ( int cigcnt = 0 ; cigcnt < b->core.n_cigar ; cigcnt++ ) {
		uint32_t ciglen = cigdata[cigcnt] >> 4 ;
		uint32_t cigtype = cigdata[cigcnt] & 15 ;
		
		if ( cigtype == BAM_CMATCH ) { // OK
			for ( int b = 0 ; b < ciglen ; b++ , p++ , rp++ ) {
				if ( *(base->refseq+p-start+1) != bam2char[bam1_seqi(s,rp)] ) {
					hline ( psnps , p , p , y ) ;
				} else {
					hline ( bucket , p , p , y ) ;
				}
			}
		} else if ( cigtype == BAM_CINS ) { // OK
			for ( int b = 0 ; b < ciglen ; b++ ) {
				int yo = b > ciglen/2 ? ciglen - b : b ;
				int xo = b - ciglen/2 ;
				yo = ciglen/2 - yo ;
				hline ( p_cigar_insertion , p + xo  , p + xo + 1 , y+yo*2 , false ) ;
				hline ( p_cigar_insertion , p + xo  , p + xo + 1 , y+yo*2-1 , false ) ;
			}
			rp += ciglen ;
		} else if ( cigtype == BAM_CDEL ) { // OK
			for ( int b = 0 ; b < ciglen ; b++ , p++ ) {
				int yo = b > ciglen/2 ? ciglen - b : b ;
				hline ( p_cigar_deletion , p , p+1 , y-yo*2 , false ) ;
				hline ( p_cigar_deletion , p , p+1 , y-yo*2+1 , false ) ;
			}
		} else if ( cigtype == BAM_CREF_SKIP ) { // UNTESTED
			rp += ciglen ;
//			p += ciglen ;
		} else if ( cigtype == BAM_CSOFT_CLIP ) { // UNTESTED
			rp += ciglen ;
		} else if ( cigtype == BAM_CHARD_CLIP ) { // UNTESTED
			p += ciglen ;
		} else if ( cigtype == BAM_CPAD ) { // UNTESTED
			p += ciglen ;
		} else {
//			cerr << cigtype << endl ;
		}
		
	}
	
	// Arrows
	if ( base->o_arrows ) {
		paint_arrow ( bucket , b->core.pos , p , y , b->core.flag & BAM_FREVERSE ) ;
	}
}


void Tbam_draw_paired::paint_arrow ( unsigned char *bucket , postype from , postype to , postype y , bool reverse ) {
	if ( reverse ) {
		to = from + ARROW_LENGTH ;
	} else {
		from = to - ARROW_LENGTH ;
		y -= ARROW_LENGTH ;
	}

	if ( y < vstart || y >= vend ) return ;
	
	if ( from < start ) from = start ;
	if ( to < start ) to = start ;
	if ( from >= end ) from = end-1 ;
	if ( to >= end ) to = end-1 ;
	
	for ( int a = 0 ; a < ARROW_LENGTH ; a++ ) {
		int m = pos2mem ( from , y ) ;
		if ( m >= 0 && m < max_mem ) {
			if ( *(bucket+m) < SATURATION_THRESHOLD ) (*(bucket+m))++ ;
		}
		from++ ;
		y++ ;
	}
}

void Tbam_draw_paired::hline ( unsigned char *bucket , postype from , postype to , postype y , bool count_in_coverage ) {
	if ( y < vstart || y >= vend ) return ;
	
	if ( from < start ) from = start ;
	if ( to < start ) to = start ;
	if ( from >= end ) from = end-1 ;
	if ( to >= end ) to = end-1 ;
	
	int mp1 = pos2mem ( from , y , false ) ;

	if ( mp1 < 0  || mp1 >= max_mem ) {
//		cerr << "OUT-OF-BOUNDS HLINE!" << endl ;
		return ;
	}

	
//	int mp2 = pos2mem ( to , y , true ) ;
	int mp2 = from == to && size > w * 2 ? mp1 : pos2mem ( to , y , true ) ;

	if ( mp1 == mp2 ) {
		if ( *(bucket+mp1) < SATURATION_THRESHOLD ) (*(bucket+mp1))++ ;
		return ;
	}

	if ( mp2 < 0 || mp2 >= max_mem ) {
//		cerr << "OUT-OF-BOUNDS HLINE!" << endl ;
		return ;
	}
	

	for ( unsigned char *b = bucket + mp1 ; mp1 <= mp2 ; mp1++ , b++ ) {
		if ( *b < SATURATION_THRESHOLD ) (*b)++ ;
	}
}

postype Tbam_draw_paired::pos2mem ( postype x , postype y , bool outer_right ) {
	if ( x < start || x > end ) return -1 ;
	y = ( y - vstart ) * h / vsize ;
	if ( y < 0 || y >= h ) return -1 ;
	
	x = ( x - start ) * w ;

	if ( outer_right ) {
		x = ( x+w ) / size - 1 < x / size ? x / size : ( x+w ) / size - 1 ;
		if ( x >= w ) x = w-1 ;
	} else {
		x /= size ;
	}
	
	return ( h - y - 1 ) * w + x ;
}

void Tbam_draw_paired::draw_axis () {
	postype a ;
	postype bot = base->get_height() - bottom ;
	postype rig = base->get_width() - right ;

	if ( left ) {
		for ( a = top ; a <= bot ; a++ ) {
			set_pixel ( left-1 , top+a ) ;
		}
	}
	for ( a = left ; a <= rig ; a++ ) {
		set_pixel ( a , bot ) ;
	}
	
	postype step_h = opt_axis ( start , end , w / 200 ) ;
	for ( a = start / step_h ; a * step_h <= end ; a++ ) {
		if ( a * step_h < start ) continue ;
		postype x = ( ( a * step_h ) - start ) * w / size + left ;
		for ( postype b = 0 ; b < 5 ; b++ ) set_pixel ( x , bot + b ) ;
		render_number ( a * step_h , x , bot + 10 , NUMBER_ALIGN_HCENTER ) ;
	}
	
	postype vstart2 = vstart + single_offset ;
	postype step_v = opt_axis ( vstart2 , vend , h / 100 ) ;
	for ( a = vstart / step_v ; a * step_v <= vend ; a++ ) {
		if ( a * step_v < vstart ) continue ;
		if ( a * step_v == 0 && single_offset == 0 ) continue ;
		postype y = bot - ( (a * step_v) - vstart + single_offset ) * h / vsize ;
		if ( left ) {
			for ( postype b = 0 ; b < 5 ; b++ ) set_pixel ( left - b , y ) ;
			render_number ( a * step_v , left-7 , y , NUMBER_ALIGN_RIGHT|NUMBER_ALIGN_VCENTER ) ;
		} else {
			for ( postype b = 0 ; b < 5 ; b++ ) set_pixel ( b , y ) ;
			render_number ( a * step_v , 7 , y , NUMBER_ALIGN_VCENTER ) ;
		}
	}
}

void Tbam_draw_paired::merge_all () {
	// Clear canvas
	png_bytep * row_pointers = base->get_row_pointers () ;
	for ( int y = 0 ; y < base->get_height() ; y++) {
		png_byte* row = row_pointers[y];
		for ( int x = 0 ; x < base->get_width() ; x++) {
			png_byte* ptr = &(row[x*4]);
			ptr[0] = 255 ;
			ptr[1] = 255 ;
			ptr[2] = 255 ;
			ptr[3] = 255 ;
		}
	}

	if ( base->use_highlight ) {
		postype bot = base->get_height() - bottom ;
		postype rig = base->get_width() - right ;
		basecol_red = 0xFF ;
		basecol_green = 0xF3 ;
		basecol_blue = 0x80 ;
		postype y_upper = bot - ( base->highlight_to - vstart + single_offset ) * h / vsize ;
		postype y_lower = bot - ( base->highlight_from - vstart + single_offset ) * h / vsize ;
		for ( int x = left+1 ; x < rig ; x++ ) {
			for ( int y = y_upper ; y <= y_lower ; y++ ) {
				set_pixel ( x , y ) ;
			}
		}
		basecol_red = 0 ;
		basecol_green = 0 ;
		basecol_blue = 0 ;
	}
	
	if ( base->o_linkpairs ) merge_into_png ( pconn , 220 , 220 , 220 ) ;
	if ( base->o_single ) merge_into_png ( psingle , 0 , 255 , 0 ) ;
	if ( base->o_pairs ) {
		merge_into_png ( plowq , 0xFF , 0x80 , 0x40 ) ;
		merge_into_png ( ppaired , 0 , 0 , 255 ) ;
	}
	if ( base->o_faceaway ) {
		merge_into_png ( pfaceaway1 , 0 , 128 , 128 ) ;
		merge_into_png ( pfaceaway2 , 128 , 128 , 0 ) ;
	}
	if ( base->o_inversions ) {
		merge_into_png ( pinversions1 , 0x80 , 0x80 , 0 ) ;
		merge_into_png ( pinversions2 , 0x40 , 0x80 , 0x80 ) ;
	}
	
	if ( base->o_snps ) {
		merge_into_png ( p_cigar_insertion , 0 , 200 , 0 ) ;
		merge_into_png ( p_cigar_deletion  , 200 , 0 , 0 ) ;
		merge_into_png ( psnps , 255 , 0 , 0 ) ;
	}
	if ( !base->o_noscale ) draw_axis () ;
//	render_number ( base->total_reads , 200 , 50 ) ;
}


//////////////////////////////////////////////////////////////////////////////////////
// Tbam_draw_coverage

Tbam_draw_coverage::Tbam_draw_coverage ( Tbam2png *_base ) : Tbam_draw_paired ( _base ) {
//	left = right = top = bottom = 0 ;
}

void Tbam_draw_coverage::set_range () {
	Tbam_draw::set_range () ;
	cov = new int[size] ;
	for ( int i = 0 ; i < size ; i++ ) cov[i] = 0 ;
}

void Tbam_draw_coverage::hline ( unsigned char *bucket , postype from , postype to , postype y , bool count_in_coverage ) {
	if ( !count_in_coverage ) return ;
	if ( from < start ) from = start ;
	if ( to >= end ) to = end - 1 ;
	for ( int i = from ; i <= to ; i++ ) cov[i-start]++ ;
}

int Tbam_draw_coverage::get_neccessary_height () {
	int max = 1 ;

	int i ;
	int width = base->get_width() ;
	int upper = size > width ? width : size ;
	int *p1 = new int[upper] ;
	int *p2 = new int[upper] ;
	for ( i = 0 ; i < upper ; i++ ) p1[i] = p2[i] = 0 ;

	for ( i = 0 ; i < size ; i++ ) {
		int x = upper * i / size ;
		p1[x] += cov[i] ;
		p2[x]++ ;
	}
	for ( i = 0 ; i < upper ; i++ ) p1[i] /= ( p2[i] == 0 ? 1 : p2[i] ) ;
	for ( i = 0 ; i < upper ; i++ ) {
		if ( max < p1[i] ) max = p1[i] ;
	}
	delete [] p1 ;
	delete [] p2 ;
//	for ( i = 0 ; i < size ; i++ ) {
//		if ( max < cov[i] ) max = cov[i] ;
//	}
	return max+1 ;
}

void Tbam_draw_coverage::draw_axis () {
	postype a ;
	postype bot = base->get_height() ;
	
	postype step_v = opt_axis ( 0 , bot , bot/20 ) ;
	if ( step_v <= 1 ) return ;
	for ( a = 1 ; a * step_v <= bot ; a++ ) {
		postype y = bot - a * step_v ;
		for ( postype b = 0 ; b < 5 ; b++ ) set_pixel ( b , y ) ;
		render_number ( a * step_v , 7 , y , NUMBER_ALIGN_VCENTER ) ;
	}
}

void Tbam_draw_coverage::merge_all () {
	// Clear canvas
	png_bytep * row_pointers = base->get_row_pointers () ;
	for ( int y = 0 ; y < base->get_height() ; y++) {
		png_byte* row = row_pointers[y];
		for ( int x = 0 ; x < base->get_width() ; x++) {
			png_byte* ptr = &(row[x*4]);
			ptr[0] = 255 ;
			ptr[1] = 255 ;
			ptr[2] = 255 ;
			ptr[3] = 255 ;
		}
	}

	int width = base->get_width() ;
	int height = base->get_height() ;
	int i ;
	int *p1 = new int[width] ;
	int *p2 = new int[width] ;
	for ( i = 0 ; i < width ; i++ ) p1[i] = p2[i] = 0 ;
	
	if ( size >= width ) {
		for ( i = 0 ; i < size ; i++ ) {
			int x = width * i / size ;
			p1[x] += cov[i] ;
			p2[x]++ ;
		}
	} else {
		for ( i = 0 ; i < size ; i++ ) {
			for ( int x = width * i / size ; x < width * ( i + 1 ) / size ; x++ ) {
				p1[x] += cov[i] ;
				p2[x]++ ;
			}
		}
	}
	
	for ( i = 0 ; i < width ; i++ ) p1[i] /= ( p2[i] == 0 ? 1 : p2[i] ) ;
	
	for ( i = 0 ; i < width ; i++ ) {
		if ( p1[i] == 0 ) continue ;
		for ( int j = 0 ; j < p1[i] ; j++ ) {
			int y = height - j - 1 ;
			png_byte* row = row_pointers[y];
			png_byte* ptr = &(row[i*4]);
			ptr[0] = 0 ;
			ptr[1] = 0 ;
			ptr[2] = 255 ;
			ptr[3] = 255 ;
		}
	}

	//if ( !base->o_noscale ) 
	draw_axis () ;
}


//////////////////////////////////////////////////////////////////////////////////////
// Tbam2png
Tbam_draw_pileup::Tbam_draw_pileup ( Tbam2png *_base ) : Tbam_draw_paired ( _base ) {
}

void Tbam_draw_pileup::set_range () {
	Tbam_draw_paired::set_range() ;
}

int Tbam_draw_pileup::get_neccessary_height () {
	render_as_text = base->get_width() >= NUMBER_WIDTH * size ;
	if ( render_as_text ) return (pile.size()+2)*10 + 2 ;
	return (pile.size()+2) ;
}

void Tbam_draw_pileup::draw_axis () {
}

void Tbam_draw_pileup::get_pileup_read_start ( const bam1_t *b , int &row , int &pos , bool has_insertion ) {
	pos = b->core.pos + 1 - start ;
	for ( row = 0 ; row < pile.size() ; row++ ) {
		bool occupied = false ;
		postype p = b->core.pos + 1 ;
		for ( postype a = -1 ; a < b->core.l_qseq + 1 ; a++ , p++ ) {
			if ( p < start ) continue ;
			if ( p >= end ) break ;
			if ( pile[row][pos+a].type == PC_EMPTY ) continue ;
			if ( has_insertion && ( pile.size() < row && pile[row+1][pos+a].type == PC_EMPTY ) ) continue ;
			occupied = true ;
			break ;
		}
		if ( occupied ) continue ;
		return ;
	}
	
	row = pile.size() ;
	pile.push_back ( VPC ( size , Tpileupchar() ) ) ;
//	if ( has_insertion ) pile.push_back ( VPC ( size , Tpileupchar() ) ) ;
}

void Tbam_draw_pileup::paint_single_read ( unsigned char *bucket , const bam1_t *b, void *data , int y ) {
	postype from , to ;
	bool has_insertion = false ;
	uint32_t *cigdata = bam1_cigar(b) ;

	// TODO use bam_cigar2qlen, maybe
	if ( b->core.n_cigar == 1 ) {
		from = b->core.pos ;
		to = from + b->core.l_qseq - 1 ;
	} else { // cigar
		from = b->core.pos ;
		to = from + b->core.l_qseq - 1 ;
		
		for ( int cigcnt = 0 ; cigcnt < b->core.n_cigar ; cigcnt++ ) {
			if ( BAM_CINS == cigdata[cigcnt] & 15 ) has_insertion = true ;
		}
	}
	
//	if ( y < vstart || y >= vend ) return ;
//	if ( global_test_flag ) cout << "!\n" ;
	
	int pile_row , pile_pos ;
	get_pileup_read_start ( b , pile_row , pile_pos , has_insertion ) ;

	uint8_t *s = bam1_seq ( b ) ;
	postype p = b->core.pos + 1 ;
	char *q = (char*) bam1_qual(b) ;
	
	uint32_t rp = 0 ;
	for ( int cigcnt = 0 ; cigcnt < b->core.n_cigar ; cigcnt++ ) {
		uint32_t ciglen = cigdata[cigcnt] >> 4 ;
		uint32_t cigtype = cigdata[cigcnt] & 15 ;

		if ( cigtype == BAM_CMATCH ) {
			for ( postype a = 0 ; a < ciglen ; a++ , p++ , pile_pos++ , q++ ) {
				if ( p < start ) continue ;
				if ( p >= end ) return ;

				char base_read = bam2char[bam1_seqi(s,a+rp)] ;
				char base_ref = base->refseq ? *(base->refseq+p-start) : base_read ;
				pile[pile_row][pile_pos].c = base_read ;
				pile[pile_row][pile_pos].quality = *q ;
				if ( base_read == base_ref ) {
					pile[pile_row][pile_pos].type |= PC_REFERENCE ;
				} else {
					pile[pile_row][pile_pos].type |= PC_SNP ;
				}
				if ( b->core.flag & BAM_FMUNMAP ) pile[pile_row][pile_pos].type |= PC_SINGLE ;
			}
			rp += ciglen ;
		} else if ( cigtype == BAM_CINS ) {
			for ( int b = 0 ; b < ciglen ; b++ ) {
				while ( pile.size() <= pile_row+1 ) pile.push_back ( VPC ( size , Tpileupchar() ) ) ;
				char base_read = bam2char[bam1_seqi(s,b+rp)] ;
				if ( pile_pos+b > 0 && pile_pos+b < pile[pile_row+1].size() ) {
					pile[pile_row+1][pile_pos+b].type |= PC_SNP ;
					pile[pile_row+1][pile_pos+b].type |= PC_INSERTION ;
					pile[pile_row+1][pile_pos+b].c = base_read ;
				}
			}
			rp += ciglen ;
		} else if ( cigtype == BAM_CDEL ) { // OK
			for ( int b = 0 ; b < ciglen ; b++ , p++ , pile_pos++ ) {
				if ( p < start ) continue ;
				if ( p >= end ) return ;
				pile[pile_row][pile_pos].type = PC_DELETION ;
				pile[pile_row][pile_pos].c = '-' ;
			}
		} else if ( cigtype == BAM_CREF_SKIP ) { // UNTESTED
			rp += ciglen ;
			p += ciglen ;
		} else if ( cigtype == BAM_CSOFT_CLIP ) { // OK
			rp += ciglen ;
		} else if ( cigtype == BAM_CHARD_CLIP ) { // UNTESTED
/*			rp += ciglen ;
			pile_pos += ciglen ;
			p += ciglen ;*/
		} else if ( cigtype == BAM_CPAD ) { // UNTESTED
			pile_pos += ciglen ;
			p += ciglen ;
		} else {
//			cerr << cigtype << endl ;
		}
		if ( rp >= b->core.l_qseq ) return ;
//		if ( cigtype != BAM_CMATCH ) pile[pile_row][pile_pos].c = cigtype ;
	}
}

void Tbam_draw_pileup::hline ( unsigned char *bucket , postype from , postype to , postype y , bool count_in_coverage ) {
}

string Tbam_draw_pileup::get_text_rendering () {
	string ret ;
	
	char dummy[100] ;
	sprintf ( dummy , "%d-%d\n" , (int)start , (int)end ) ;
	ret += dummy ;

	if ( base->refseq ) {
		for ( int col = 0 ; col < size ; col++ ) {
			char c = *(base->refseq+col) ;
			ret += c ;
		}
		ret += "\n\n" ;
	}

	for ( int row = 0 ; row < pile.size() ; row++ ) {
		for ( int col = 0 ; col < pile[row].size() && col < size ; col++ ) {
			if ( pile[row][col].type == PC_EMPTY ) {
				ret += ' ' ;
			} else if ( pile[row][col].type &= PC_INSERTION ) {
				ret += pile[row][col].c - 'A' + 'a' ;
			} else {
				ret += pile[row][col].c ;
			}
		}
		ret += "\n" ;
	}
	return ret ;
}

void Tbam_draw_pileup::merge_all () {
	// Clear canvas
	png_bytep * row_pointers = base->get_row_pointers () ;
	for ( int y = 0 ; y < base->get_height() ; y++) {
		png_byte* row = row_pointers[y];
		for ( int x = 0 ; x < base->get_width() ; x++) {
			png_byte* ptr = &(row[x*4]);
			ptr[0] = 255 ;
			ptr[1] = 255 ;
			ptr[2] = 255 ;
			ptr[3] = 255 ;
		}
	}
	
	int width = base->get_width() ;
	int height = base->get_height() ;
	
	bool use_quality = base->o_readqual ;
	
	if ( render_as_text ) { // Bases as text
		int row_height = 10 ;
		for ( int row = 0 ; row < pile.size() ; row++ ) {
			int y = height - ( row + 3 ) * row_height ;
			for ( int col = 0 ; col < pile[row].size() && col < size ; col++ ) {
				if ( pile[row][col].type == PC_EMPTY ) continue ;
				
				int x = col * width / size ;
				if ( use_quality && pile[row][col].quality != -1 ) {
					int i = pile[row][col].quality * 6 ; // FIXME if quality > 40
					if ( i > 255 ) i = 255 ;
					fill_rect ( x-1 , y-1 , ( col + 1 ) * width / size - 2 , y + row_height - 2 , i , i , i ) ;
				}

				if ( pile[row][col].type & PC_REFERENCE ) { basecol_red = basecol_green = 0 ; basecol_blue = 255 ; }
				if ( pile[row][col].type & PC_SINGLE ) { basecol_red = 0 ; basecol_green = 255 ; basecol_blue = 0 ; }
				if ( pile[row][col].type & PC_SNP ) { basecol_red = 255 ; basecol_green = basecol_blue = 0 ; }

				if ( pile[row][col].type & PC_INSERTION ) {
					basecol_red = 0 ;
					basecol_green = 0 ;
					basecol_blue = 0 ;
					fill_rect ( x , y-1 , ( col + 1 ) * width / size - 2 , y + row_height - 2 , 0 , 200 , 0 ) ;
				}

				if ( pile[row][col].type & PC_DELETION ) {
					basecol_red = 255 ;
					basecol_green = 0 ;
					basecol_blue = 0 ;
				}

				render_char ( pile[row][col].c , x , y ) ;

			}
		}
		
		basecol_red = basecol_green = basecol_blue = 0 ;
		
		if ( base->refseq ) {
			int y = height - 10 ;
			for ( int col = 0 ; col < size ; col++ ) {
				char c = *(base->refseq+col) ;
				int x = col * width / size ;
				render_char ( c , x , y ) ;
			}
		}
		
	} else { // Bases as dots
	
		// Reference
		basecol_red = basecol_green = 0 ; basecol_blue = 255 ;
		for ( int row = 0 ; row < pile.size() ; row++ ) {
			int y = height - ( row + 1 ) ;
			for ( int col = 0 ; col < pile[row].size() && col < size ; col++ ) {
				if ( 0 == ( pile[row][col].type & PC_REFERENCE ) ) continue ;

				if ( use_quality ) {
					int i ;
					if ( pile[row][col].quality != -1 ) {
						i = pile[row][col].quality * 6 ; // FIXME if quality > 40
						if ( i > 255 ) i = 255 ;
					} else i = 255 ;
					basecol_blue = i ;
				}
				
				int x1 = col * width / size ;
				int x2 = ( col + 1 ) * width / size - 1 ;
				if ( x1 > x2 ) x2 = x1 ;
				for ( int x = x1 ; x <= x2 ; x++ ) {
					set_pixel ( x , y ) ;
				}
				
			}
		}

		// InDels
//		basecol_red = 255 ; basecol_green = basecol_blue = 0 ;
		for ( int row = 0 ; row < pile.size() ; row++ ) {
			int y = height - ( row + 1 ) ;
			for ( int col = 0 ; col < pile[row].size() && col < size ; col++ ) {
				bool doit = false ;
				if ( pile[row][col].type & PC_INSERTION ) {
					basecol_red = 0 ; basecol_green = 128 ; basecol_blue = 0 ;
					doit = true ;
				} else if ( pile[row][col].type & PC_DELETION ) {
					basecol_red = 128 ; basecol_green = 0 ; basecol_blue = 0 ;
					doit = true ;
				}
				
				if ( !doit ) continue ;
				
				int x1 = col * width / size ;
				int x2 = ( col + 1 ) * width / size - 1 ;
				if ( x1 > x2 ) x2 = x1 ;
				for ( int x = x1 ; x <= x2 ; x++ ) {
					set_pixel ( x , y ) ;
				}
				
			}
		}

		// SNPs
		basecol_red = 255 ; basecol_green = basecol_blue = 0 ;
		for ( int row = 0 ; row < pile.size() ; row++ ) {
			int y = height - ( row + 1 ) ;
			for ( int col = 0 ; col < pile[row].size() && col < size ; col++ ) {
				if ( pile[row][col].type != PC_SNP ) continue ;
				
				if ( use_quality ) {
					int i ;
					if ( pile[row][col].quality != -1 ) {
						i = pile[row][col].quality * 6 ; // FIXME if quality > 40
						if ( i > 255 ) i = 255 ;
					} else i = 255 ;
					basecol_red = i ;
				}

				int x1 = col * width / size ;
				int x2 = ( col + 1 ) * width / size - 1 ;
				if ( x1 > x2 ) x2 = x1 ;
				for ( int x = x1 ; x <= x2 ; x++ ) {
					set_pixel ( x , y ) ;
				}
				
			}
		}
		
		basecol_red = basecol_green = basecol_blue = 0 ;
	}


	//if ( !base->o_noscale ) 
	draw_axis () ;
}




//////////////////////////////////////////////////////////////////////////////////////
// Tbam2png

void Tbam2png::init ( string _bam_file , string _region , string _png_file , int _mapq ) 
{
	use_highlight = false ;
	bam_file = _bam_file ;
	region = _region ;
	png_file = _png_file ;
	png_file = _png_file ;
	mapq = _mapq ;
	draw = NULL ;
	refseq = NULL ;
	width = 1024 ;
	height = 768 ;
	o_single = o_pairs = o_arrows = o_snps = o_faceaway = o_inversions = o_linkpairs = o_colordepth = false ;
	o_noscale = o_readqual = o_text = false ;
}

void Tbam2png::set_image_size ( postype w , postype h ) {
	width = w ;
	height = h ;
}

void Tbam2png::set_options ( string options ) {
	const char *last = options.c_str() ;
	char *c = (char*) options.c_str() ;
	vector <string> ov ;
	for ( c++ ; *c ; c++ ) {
		if ( *c == ',' ) {
			*c = 0 ;
			ov.push_back ( last ) ;
			last = c + 1 ;
		}
	}
	ov.push_back ( last ) ;
	for ( int a = 0 ; a < ov.size() ; a++ ) {
		if ( ov[a] == "snps" ) o_snps = true ;
		else if ( ov[a] == "pairs" ) o_pairs = true ;
		else if ( ov[a] == "arrows" ) o_arrows = true ;
		else if ( ov[a] == "single" ) o_single = true ;
		else if ( ov[a] == "faceaway" ) o_faceaway = true ;
		else if ( ov[a] == "inversions" ) o_inversions = true ;
		else if ( ov[a] == "linkpairs" ) o_linkpairs = true ;
		else if ( ov[a] == "colordepth" ) o_colordepth = true ;
		else if ( ov[a] == "noscale" ) o_noscale = true ;
		else if ( ov[a] == "readqual" ) o_readqual = true ;
		else if ( ov[a] == "text" ) o_text = true ;
	}
	
	if ( o_single ) draw->single_offset = 50 ;
}


void Tbam2png::read_bam_file () {
	if ( draw == NULL ) abort_ ( "No drawing class instanced!\n" ) ;

	
	tmp.beg = 0 ;
	tmp.end = 0;   
	tmp.in = samopen(bam_file.c_str(), "rb", 0); 
	
//	cout << tmp.in->header->text << endl ;
	for ( char *c = tmp.in->header->text ; *c ; c++ ) {
		if ( *c == '\n' && *(c+1) == '@' && *(c+2) == 'C' && *(c+3) == 'O' ) {
			string s ;
			use_highlight = false ;
			for ( char *d = c+5 ; *d > 13 ; d++ ) {
				if ( *d == ' ' ) {
					if ( s == "HIGHLIGHT" ) use_highlight = true ;
					s = "" ;
				} else s += *d ;
			}
			if ( use_highlight ) {
				highlight_from = atoi ( (char*) s.c_str() ) ;
				const char *d ;
				for ( d = s.c_str() ; *(d-1) != '-' ; d++ ) ;
				highlight_to = atoi ( d ) ;
			}
		}
	}
	
	
	int ref;  
	bam_index_t *idx;  
	bam_plbuf_t *buf;  
	idx = bam_index_load(bam_file.c_str()); // load BAM index  
	if (idx == 0) {  
		fprintf(stderr, "BAM indexing file is not available.\n");  
		return ;  
	}  
	bam_parse_region(tmp.in->header, region.c_str(), &ref,  &tmp.beg, &tmp.end); // parse the region  
	
	if ( !refseq_file.empty() ) {
		faidx_t *fai = fai_load ( refseq_file.c_str() ) ;
		int len = tmp.end - tmp.beg + 2 ;
		refseq = fai_fetch ( fai , region.c_str() , &len ) ;
		if ( len < tmp.end - tmp.beg ) tmp.end = len + tmp.beg ;
/*		if ( tmp.beg == 0 ) { // No range set
			tmp.beg = 1 ;
			tmp.end = len ;
		}*/
		fai_destroy ( fai ) ;
	}
	
	if ( width * 50 < tmp.end - tmp.beg ) o_arrows = false ;
	
	draw->set_range () ;
	total_snps = 0 ;
	total_reads = 0 ;
	bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, NULL, fetch_func);  
	bam_index_destroy(idx);  

//	cout << "Total SNPs : " << total_snps << endl ;
}


int Tbam2png::fetch_func(const bam1_t *b, void *data) {  
	if ( b->core.qual < b2p->mapq ) return 0 ;
	b2p->total_reads++ ;


/*	global_test_flag = 0 ;
	if ( b->core.flag == 83 ) { // TEST FIXME
		char s[1000] ;
		for ( postype a = 0 ; a < b->core.l_qseq ; a++ ) {
			s[a] = (char) bam2char[bam1_seqi(bam1_seq(b),a)] ;
		}
		s[b->core.l_qseq] = 0 ;
		if ( strstr ( s , "GAATTAAACGATTG" ) ) global_test_flag = 1 ;
	}*/

	if ( b->core.flag & BAM_FPROPER_PAIR ) {

		//if ( b2p->o_pairs ) 
		b2p->draw->draw_paired_read ( b , data ) ;
	} else if ( b->core.flag & BAM_FUNMAP ) {
//		cout << "UNMAPPED\t" ;
	} else if ( b->core.flag & BAM_FMUNMAP ) {
		if ( b2p->o_single ) b2p->draw->draw_single_read ( b , data ) ;
	} else if ( b->core.isize != 0 ) {
		b2p->draw->draw_paired_read ( b , data ) ;
	} else {
		if ( b2p->o_single ) b2p->draw->draw_single_read ( b , data ) ;
/*		cout << "DUNNO\t" ;
		bool read_reverse = b->core.flag & BAM_FMUNMAP ;
		bool mate_reverse = b->core.flag & BAM_FMREVERSE ;
		bool read1 = b->core.flag & BAM_FREAD1 ;
		bool read2 = b->core.flag & BAM_FREAD2 ;

		cout << b->core.pos << "\t" << b->core.flag << "\t" << b->core.isize << "\t" << read_reverse << "\t" << mate_reverse ;
		cout << "\t" << read1 << "\t" << read2 << endl ;
*/
	}
	return 0;  
}  



void Tbam2png::create_png () {
	int x, y;
	color_type = 6 ;
	bit_depth = 8 ;
	number_of_passes = 1 ;
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (!png_ptr)
		abort_("[read_png_file] png_create_read_struct failed");

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
		abort_("[read_png_file] png_create_info_struct failed");

	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[read_png_file] Error during init_io");

//
	png_set_sig_bytes(png_ptr, 8);

//	png_read_info(png_ptr, info_ptr);

	info_ptr->width = width;
	info_ptr->height = height;
	info_ptr->color_type = color_type;
	info_ptr->bit_depth = bit_depth;

	number_of_passes = png_set_interlace_handling(png_ptr);
	png_read_update_info(png_ptr, info_ptr);
//
		
	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
	for (y=0; y<height; y++) {
		row_pointers[y] = (png_byte*) malloc(info_ptr->rowbytes);
	}
}

void Tbam2png::write_png_file(char* file_name)
{
	int x, y;

	// create file 
	string fn ( file_name ) ;
	bool is_stdout = fn == "-" ;
	FILE *fp = is_stdout ? freopen(NULL, "wb", stdout) : fopen(file_name, "wb");
	if (!fp)
		abort_("[write_png_file] File %s could not be opened for writing", file_name);


	// initialize stuff 
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	
	if (!png_ptr)
		abort_("[write_png_file] png_create_write_struct failed");

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
		abort_("[write_png_file] png_create_info_struct failed");

	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[write_png_file] Error during init_io");

	png_init_io(png_ptr, fp);


	// write header 
	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[write_png_file] Error during writing header");

	png_set_IHDR(png_ptr, info_ptr, width, height,
		     bit_depth, color_type, PNG_INTERLACE_NONE,
		     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);


	// write bytes 
	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[write_png_file] Error during writing bytes");

	png_write_image(png_ptr, row_pointers);


	// end write 
	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[write_png_file] Error during end of write");

	png_write_end(png_ptr, NULL);

        // cleanup heap allocation 
	for (y=0; y<height; y++)
		free(row_pointers[y]);
	free(row_pointers);

        if ( !is_stdout ) fclose(fp);
}


void Tbam2png::process_file(void) {
	int x , y ;
	if (info_ptr->color_type != PNG_COLOR_TYPE_RGBA)
		abort_("[process_file] color_type of input file must be PNG_COLOR_TYPE_RGBA (is %d)", info_ptr->color_type);

	for (y=0; y<height; y++) {
		png_byte* row = row_pointers[y];
		for (x=0; x<width; x++) {
			png_byte* ptr = &(row[x*4]);
//			printf("Pixel at position [ %d - %d ] has the following RGBA values: %d - %d - %d - %d\n", x, y, ptr[0], ptr[1], ptr[2], ptr[3]);

			if ( x == y || x == 0 || y == 0 ) ptr[0] = ptr[1] = ptr[2] = 255 ;
			else ptr[0] = ptr[1] = ptr[2] =  ( ( x + y ) / 2 ) % 255 ;
			ptr[3] = 255 ;

			// set red value to 0 and green value to the blue one 
//			ptr[0] = 0;
//			ptr[1] = ptr[2];
		}
	}

}

void Tbam2png::abort_(const char * s, ...)
{
	va_list args;
	va_start(args, s);
	vfprintf(stderr, s, args);
	fprintf(stderr, "\n");
	va_end(args);
	abort();
}



// MAIN STUFF

void init_lcd_chars () {
	for ( int a = 0 ; a < 256 ; a++ ) lcd_chars.push_back ( TVI() ) ;
	
	lcd_chars['-'].push_back(3);
	
	
	lcd_chars['0'].push_back(0);
	lcd_chars['0'].push_back(1);
	lcd_chars['0'].push_back(2);
	lcd_chars['0'].push_back(4);
	lcd_chars['0'].push_back(5);
	lcd_chars['0'].push_back(6);
	
	lcd_chars['1'].push_back(2);
	lcd_chars['1'].push_back(5);
	
	lcd_chars['2'].push_back(0);
	lcd_chars['2'].push_back(2);
	lcd_chars['2'].push_back(3);
	lcd_chars['2'].push_back(4);
	lcd_chars['2'].push_back(6);
	
	lcd_chars['3'].push_back(0);
	lcd_chars['3'].push_back(2);
	lcd_chars['3'].push_back(3);
	lcd_chars['3'].push_back(5);
	lcd_chars['3'].push_back(6);
	
	lcd_chars['4'].push_back(1);
	lcd_chars['4'].push_back(2);
	lcd_chars['4'].push_back(3);
	lcd_chars['4'].push_back(5);
	
	lcd_chars['5'].push_back(0);
	lcd_chars['5'].push_back(1);
	lcd_chars['5'].push_back(3);
	lcd_chars['5'].push_back(5);
	lcd_chars['5'].push_back(6);
	
	lcd_chars['6'].push_back(0);
	lcd_chars['6'].push_back(1);
	lcd_chars['6'].push_back(3);
	lcd_chars['6'].push_back(4);
	lcd_chars['6'].push_back(5);
	lcd_chars['6'].push_back(6);
	
	lcd_chars['7'].push_back(0);
	lcd_chars['7'].push_back(2);
	lcd_chars['7'].push_back(5);
	
	lcd_chars['8'].push_back(0);
	lcd_chars['8'].push_back(1);
	lcd_chars['8'].push_back(2);
	lcd_chars['8'].push_back(3);
	lcd_chars['8'].push_back(4);
	lcd_chars['8'].push_back(5);
	lcd_chars['8'].push_back(6);
	
	lcd_chars['9'].push_back(0);
	lcd_chars['9'].push_back(1);
	lcd_chars['9'].push_back(2);
	lcd_chars['9'].push_back(3);
	lcd_chars['9'].push_back(5);
	lcd_chars['9'].push_back(6);
	
	lcd_chars['A'].push_back(0);
	lcd_chars['A'].push_back(1);
	lcd_chars['A'].push_back(2);
	lcd_chars['A'].push_back(3);
	lcd_chars['A'].push_back(4);
	lcd_chars['A'].push_back(5);
	
	lcd_chars['C'].push_back(0);
	lcd_chars['C'].push_back(1);
	lcd_chars['C'].push_back(4);
	lcd_chars['C'].push_back(6);

	lcd_chars['G'].push_back(0);
	lcd_chars['G'].push_back(1);
	lcd_chars['G'].push_back(3);
	lcd_chars['G'].push_back(4);
	lcd_chars['G'].push_back(5);
	lcd_chars['G'].push_back(6);
	
	lcd_chars['T'].push_back(0);
	lcd_chars['T'].push_back(7);
	lcd_chars['T'].push_back(8);
	
	// FIXME
	lcd_chars['N'].push_back(1);
	lcd_chars['N'].push_back(2);
	lcd_chars['N'].push_back(3);
	lcd_chars['N'].push_back(4);
	lcd_chars['N'].push_back(5);
}



int die_usage () {
	cout << "Usage : render_image <options> --bam=FILE --png=FILE --region=\"REGION\"" << endl ;
	cout << "-w --width   NUMBER   Image width (1024)" << endl ;
	cout << "-h --height  NUMBER   Image height (768)" << endl ;
	cout << "-v --vmin    NUMBER   Vertical lower limit (0)" << endl ;
	cout << "-V --vmax    NUMBER   Vertical upper limit (1000)" << endl ;
	cout << "-r --ref     FILE     Reference sequence" << endl ;
	cout << "-i --view    TEXT     coverage or indel (default)" << endl ;
	cout << "-o --options TEXT     arrows,snps,pairs,single,faceaway,inversions,linkpairs,colordepth,noscale" << endl ;
	return 0 ;
}

int main(int argc, char **argv) {
	for ( int a = 0 ; a < 256 ; a++ ) bam2char[a] = '?' ;
	bam2char[1] = 'A' ;
	bam2char[2] = 'C' ;
	bam2char[4] = 'G' ;
	bam2char[8] = 'T' ;
	bam2char[15] = 'N' ;

	string view = "indel" ;
	string bam_file , ref_file , png_file , region , options ;
	int width = 1024 ;
	int height = 768 ;
	int vmin = 0 ;
	int vmax = 1000 ;
	int mapq = 0 ;
	static struct option long_options[] = {
		{ "bam" , optional_argument , 0 , 'b' } ,
		{ "ref" , optional_argument , 0 , 'r' } ,
		{ "region" , optional_argument , 0 , 'R' } ,
		{ "png" , optional_argument , 0 , 'p' } ,
		{ "width" , optional_argument , 0 , 'w' } ,
		{ "height" , optional_argument , 0 , 'h' } ,
		{ "vmin" , optional_argument , 0 , 'v' } ,
		{ "vmax" , optional_argument , 0 , 'V' } ,
		{ "options" , optional_argument , 0 , 'o' } ,
		{ "view" , optional_argument , 0 , 'i' } ,
		{ "mapq" , optional_argument , 0 , 'm' } ,
//		{ "snpsonly" , optional_argument , &snpsonly , true } ,
		{ 0 , 0 , 0 , 0 }
	} ;
	int c = 0 ;
	while ( -1 != ( c = getopt_long (argc, argv, "",long_options, NULL) ) ) {
		switch ( c ) {
			case 0 : break ;
			case 'b' : bam_file = optarg ; break ;
			case 'r' : ref_file = optarg ; break ;
			case 'p' : png_file = optarg ; break ;
			case 'R' : region = optarg ; break ;
			case 'o' : options = optarg ; break ;
			case 'i' : view = optarg ; break ;
			case 'w' : width = atoi ( optarg ) ; break ;
			case 'h' : height = atoi ( optarg ) ; break ;
			case 'v' : vmin = atoi ( optarg ) ; break ;
			case 'V' : vmax = atoi ( optarg ) ; break ;
			case 'm' : mapq = atoi ( optarg ) ; break ;
		}
	}
	
	init_lcd_chars () ;
	
	if ( bam_file.empty() ) return die_usage () ;
	if ( png_file.empty() ) return die_usage () ;
	if ( region.empty() ) return die_usage () ;

	b2p = new Tbam2png () ;
	b2p->refseq_file = ref_file ;
	b2p->init( bam_file , region , png_file , mapq ) ;
	b2p->set_image_size ( width , height ) ;
	Tbam_draw *dp ;
	if ( view == "coverage" ) dp = new Tbam_draw_coverage ( b2p ) ;
	else if ( view == "pileup" ) dp = new Tbam_draw_pileup ( b2p ) ;
	else dp = new Tbam_draw_paired ( b2p ) ;
	dp->set_vertical_range ( vmin , vmax ) ;
	b2p->set_options ( options ) ;
	
	if ( view == "coverage" ) {
		b2p->o_linkpairs = false ;
		b2p->o_snps = false ;
	} else if ( view == "pileup" ) {
		b2p->o_linkpairs = false ;
	}

	
	b2p->read_bam_file () ;
	
	if ( view == "coverage" ) {
		height = dp->get_neccessary_height () ;
		b2p->set_image_size ( width , height ) ;
	} else if ( view == "pileup" ) {
		height = dp->get_neccessary_height () ;
		b2p->set_image_size ( width , height ) ;
	}


	b2p->create_png () ;
	dp->merge_all () ;

	if ( b2p->o_text ) {
		string s = ( (Tbam_draw_pileup*) dp)->get_text_rendering() ;
		if ( png_file == "-" ) {
			cout << s ;
		} else {
			FILE *fp = fopen(png_file.c_str(), "wb");
			fprintf ( fp , "%s" , s.c_str() ) ;
			fclose ( fp ) ;
		}
	} else b2p->write_png_file((char*)png_file.c_str());

	return 0;
}

/*
\rm render_image ; g++ render_image.cpp -O3 -o render_image -lpng -L . -lbam ; time ./render_image --view=pileup --bam=/nfs/disk69/ftp-team112/pf.som/bam/PD0009-01.bam --options=snps,pairs,arrows,single,faceaway,inversions,linkpairs,colordepth --ref=/nfs/team112/refseq/plasmodium/falciparum/3D7_pm.fa --region="MAL1:500605-500805" --png=test.png


\rm render_image ; g++ render_image.cpp -O3 -o render_image -lpng -L . -lbam ; time ./render_image --bam=/nfs/users/nfs_m/mm6/ftp/ag/bam/AC0001-C.bam --options=snps,pairs,arrows,single,faceaway,inversions,linkpairs,colordepth --ref=/nfs/users/nfs_m/mm6/ftp/ag/Anopheles_gambiae.clean.fa --region="2L:1-200000" --png=2L.a.png
\rm render_image ; g++ render_image.cpp -O3 -o render_image -lpng -L . -lbam ; time ./render_image --bam="ftp://ftp.sanger.ac.uk/pub/team112/ag/bam/AC0001-C.bam" --options=pairs,arrows,single,faceaway,inversions,colordepth,snps --ref=/nfs/users/nfs_m/mm6/ftp/ag/Anopheles_gambiae.clean.fa --region="2L" --png=2L.a.png
\rm render_image ; g++ render_image.cpp -O3 -o render_image -lpng -L . -lbam ; cp render_image ~/wwwdev_data_marker3/..

*/
