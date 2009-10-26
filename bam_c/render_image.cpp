#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>

#define PNG_DEBUG 3
#include <png.h>

#include "sam.h"  
#include "faidx.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std ;

#define ARROW_LENGTH 4
#define SATURATION_THRESHOLD 64

#define NUMBER_WIDTH 4
#define NUMBER_HEIGHT 6
#define NUMBER_ALIGN_HCENTER 1
#define NUMBER_ALIGN_RIGHT 2
#define NUMBER_ALIGN_VCENTER 4

typedef long postype ;

typedef struct {  
	int beg, end;  
	samfile_t *in;  
} tmpstruct_t;  

typedef vector <postype> TVI ;

vector <TVI> numbers ;
char bam2char[255] ;

class Tbam2png ;

class Tbam_draw {
	public :
	Tbam_draw ( Tbam2png *_base ) ;
	void set_range () ;
	void set_vertical_range ( int a , int b ) ;
	virtual void draw_single_read ( const bam1_t *b, void *data ) {} ;
	virtual void draw_paired_read ( const bam1_t *b, void *data) {} ;
	virtual void merge_into_png ( unsigned char *p , int red , int green , int blue ) ;
	virtual void merge_all () {} ;
	virtual postype get_start () { return start ; }

	postype single_offset ;

	protected :
	virtual void render_digit ( int digit , int x , int y ) ;
	virtual void render_number ( int number , int x , int y , int align = 0 ) ;
	virtual void draw_axis () {} ;
//	virtual postype pos2mem ( postype x , postype y , bool outer_right = false ) { return -1 ; } ;
	virtual void paint_single_read ( unsigned char *bucket , const bam1_t *b, void *data , int y ) {} ;
	virtual void set_pixel ( int x , int y ) ;
	virtual int opt_axis ( postype from , postype to , postype max_steps ) ;
	unsigned char *alloc_p () ;
	
	Tbam2png *base ;
	postype max_mem ;
	postype left , right , top , bottom , w , h ;
	postype start , end , size ;
	postype vstart , vend , vsize ;
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
	virtual void hline ( unsigned char *bucket , postype from , postype to , postype y ) ;
	virtual void paint_arrow ( unsigned char *bucket , postype from , postype to , postype y , bool reverse ) ;
	inline postype pos2mem ( postype x , postype y , bool outer_right = false ) ;
	unsigned char *psingle , *ppaired , *psnps , *pconn , *pfaceaway1 , *pfaceaway2 , *pinversions1 , *pinversions2 ;
} ;

class Tbam2png {
	public :
	Tbam2png () {} ;
	void init ( string _bamfile , string _region , string _png_out ) ;
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
	bool o_noscale ;
	int total_snps ;
	
	private :

	void abort_(const char * s, ...) ;

	// BAM variables
	string bam_file , png_file , region ;
	tmpstruct_t tmp;
	
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
}

void Tbam_draw::set_vertical_range ( int a , int b ) {
	vstart = a ;
	vend = b ;
}

int Tbam_draw::opt_axis ( postype from , postype to , postype max_steps ) {
	int step = 1 ;
	while ( 1 ) {
		if ( ( step * 5 ) * max_steps > to - from ) break ;
		step *= 5 ;
		if ( ( step * 2 ) * max_steps > to - from ) break ;
		step *= 2 ;
	}
	return step ;
}

void Tbam_draw::render_digit ( int digit , int x , int y ) {
	if ( digit < 0 || digit > 9 ) return ;
//	cout << x << endl ;
//	cout << "Rendering " << digit << "\n" ;
	for ( int a = 0 ; a < numbers[digit].size() ; a++ ) {
		int b = numbers[digit][a] ;
		int ox = x , oy = y ;
		switch ( b ) {
			case 0 : break ;
			case 1 : break ;
			case 2 : ox += NUMBER_WIDTH ; break ;
			case 3 : oy += NUMBER_HEIGHT/2 ; break ;
			case 4 : oy += NUMBER_HEIGHT/2 ; break ;
			case 5 : ox += NUMBER_WIDTH ; oy += NUMBER_HEIGHT/2 ; break ;
			case 6 : oy += NUMBER_HEIGHT ; break ;
		}
		if ( b == 0 || b == 3 || b == 6 ) {
			for ( int a = 1 ; a < NUMBER_WIDTH ; a++ ) {
				set_pixel ( ox + a , oy ) ;
//				cout << (ox+a) << "\t" << oy << endl ;
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
		render_digit ( *c-'0' , x , y ) ;
		x += dw ;
	}
}

void Tbam_draw::set_pixel ( int x , int y ) {
	int bw = base->get_width() ;
	if ( x < 0 || y < 0 || x >= bw || y >= base->get_height() ) return ;
	png_bytep * row_pointers = base->get_row_pointers () ;
	png_byte* row = row_pointers[y];
	png_byte* ptr = &(row[x*4]);
	ptr[0]=ptr[1]=ptr[2]=0 ;
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
	psnps = alloc_p () ;
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
	if ( isize > 0 && b->core.flag & BAM_FREVERSE ) is_faceaway = true ;
	else if ( isize < 0 && 0 == ( b->core.flag & BAM_FREVERSE ) ) is_faceaway = true ;
	if ( is_faceaway ) {
		if ( base->o_faceaway ) {
			if ( isize < 0 ) bucket = pfaceaway1 ;
			else bucket = pfaceaway2 ;
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
	}
	
	if ( !is_inversion && !is_faceaway && !base->o_pairs ) return ; // No drawing normal pairs if unwanted

	// Paint read
	paint_single_read ( bucket , b , data , abs(isize) ) ;

	if ( !base->o_linkpairs ) return ;

	// Paint linkage
	if ( isize > 0 ) {
		postype x = b->core.pos ;
		hline ( pconn , x , x + isize , abs(isize) ) ;
	} else {
		postype x = b->core.pos + b->core.l_qseq - 1 ;
		hline ( pconn , x + isize , x , abs(isize) ) ;
	}
}

void Tbam_draw_paired::paint_single_read ( unsigned char *bucket , const bam1_t *b, void *data , int y ) {
	postype from , to ;
	// TODO use bam_cigar2qlen, maybe
	if ( b->core.n_cigar == 1 ) {
		from = b->core.pos ;
		to = from + b->core.l_qseq - 1 ;
	} else { // FIXME cigar
		from = b->core.pos ;
		to = from + b->core.l_qseq - 1 ;
	}
	
	if ( y < vstart || y >= vend ) return ;

	hline ( bucket , from , to , y ) ;
	
	// SNPs
//	int has_indels = 0 ; //bam_aux2i ( bam_aux_get ( b , "H1" ) ) +  bam_aux2i ( bam_aux_get ( b , "H2" ) ) ;
//	int has_snps = bam_aux2i ( bam_aux_get ( b , "NM" ) ) + has_indels ;
//	cout << bam_aux2i ( bam_aux_get ( b , "MD" ) ) << "!" << endl ;
//	bool has_snps = bam_aux2i ( bam_aux_get ( b , "MD" ) ) != 0 ;
//	bool has_snps = true ;
//	int has_snps = b->core.l_qseq - bam_aux2i ( bam_aux_get ( b , "H0" ) ) ;
	if ( base->o_snps && base->refseq ) {
		uint8_t *s = bam1_seq ( b ) ;
		postype p = b->core.pos + 1 ;
		for ( postype a = 0 ; a < b->core.l_qseq ; a++ , p++ ) {
			if ( p < start ) continue ;
			if ( p >= end ) break ;

			if ( *(base->refseq+p-start) != bam2char[bam1_seqi(s,a)] ) {
//				char c = *(base->refseq+p-start) ;
//				cout << p << "\t" << c << "\t" << (int) *(bam1_qual(b)+a) << endl ;
//				base->total_snps++ ;
				hline ( psnps , p , p , y ) ;
			}
		}
	}

	// Arrows
	if ( base->o_arrows ) {
		paint_arrow ( bucket , from , to , y , b->core.flag & BAM_FREVERSE ) ;
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

void Tbam_draw_paired::hline ( unsigned char *bucket , postype from , postype to , postype y ) {
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
	if ( mp2 < 0 || mp2 >= max_mem ) {
//		cerr << "OUT-OF-BOUNDS HLINE!" << endl ;
		return ;
	}
	
	for ( int m = mp1 ; m <= mp2 ; m++ ) {
		if ( *(bucket+m) < SATURATION_THRESHOLD ) (*(bucket+m))++ ;
	}
}

postype Tbam_draw_paired::pos2mem ( postype x , postype y , bool outer_right ) {
	y = ( y - vstart ) * h / vsize ;
	if ( y < 0 || y >= h ) return -1 ;
	
	if ( x < start || x > end ) return -1 ;

	x = ( x - start ) * w ;

	if ( outer_right ) {
		if ( ( x+w ) / size - 1 < x / size ) {
			x /= size ;
		} else {
			x = ( x+w ) / size - 1 ;
		}
		
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
	
	if ( base->o_linkpairs ) merge_into_png ( pconn , 200 , 200 , 200 ) ;
	if ( base->o_single ) merge_into_png ( psingle , 0 , 255 , 0 ) ;
	if ( base->o_pairs ) merge_into_png ( ppaired , 0 , 0 , 255 ) ;
	if ( base->o_faceaway ) {
		merge_into_png ( pfaceaway1 , 0 , 128 , 128 ) ;
		merge_into_png ( pfaceaway2 , 128 , 128 , 0 ) ;
	}
	if ( base->o_inversions ) {
		merge_into_png ( pinversions1 , 128 , 128 , 255 ) ;
		merge_into_png ( pinversions2 , 128 , 255 , 128 ) ;
	}
	if ( base->o_snps ) merge_into_png ( psnps , 255 , 0 , 0 ) ;
	if ( !base->o_noscale ) draw_axis () ;
//	render_number ( base->total_snps , 200 , 50 ) ; // Number of SNPs
}

//////////////////////////////////////////////////////////////////////////////////////
// Tbam2png

void Tbam2png::init ( string _bam_file , string _region , string _png_file ) 
{
	bam_file = _bam_file ;
	region = _region ;
	png_file = _png_file ;
	png_file = _png_file ;
	draw = NULL ;
	refseq = NULL ;
	width = 1024 ;
	height = 768 ;
	o_single = o_pairs = o_arrows = o_snps = o_faceaway = o_inversions = o_linkpairs = o_colordepth = false ;
	o_noscale = false ;
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
	}
	
	if ( o_single ) draw->single_offset = 50 ;
}


void Tbam2png::read_bam_file () {
	if ( draw == NULL ) abort_ ( "No drawing class instanced!\n" ) ;

	tmp.beg = 0 ;
	tmp.end = 0;   
	tmp.in = samopen(bam_file.c_str(), "rb", 0); 
	
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
	
	if ( width < tmp.end - tmp.beg ) o_arrows = false ;
	
	draw->set_range () ;
	total_snps = 0 ;
	bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, NULL, fetch_func);  
	bam_index_destroy(idx);  
//	cout << "Total SNPs : " << total_snps << endl ;
}


int Tbam2png::fetch_func(const bam1_t *b, void *data) {  
	if ( b->core.flag & BAM_FPROPER_PAIR ) {
		//if ( b2p->o_pairs ) 
		b2p->draw->draw_paired_read ( b , data ) ;
	} else if ( b->core.flag & BAM_FUNMAP ) {
//		cout << "UNMAPPED\t" ;
	} else if ( b->core.flag & BAM_FMUNMAP ) {
//		if ( b2p->o_single ) b2p->draw->draw_single_read ( b , data ) ;
	} else if ( b->core.isize != 0 ) {
		b2p->draw->draw_paired_read ( b , data ) ;
	} else {
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
	FILE *fp = fopen(file_name, "wb");
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

        fclose(fp);
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

void init_numbers () {
	for ( int a = 0 ; a < 10 ; a++ ) numbers.push_back ( TVI() ) ;
	numbers[0].push_back(0);
	numbers[0].push_back(1);
	numbers[0].push_back(2);
	numbers[0].push_back(4);
	numbers[0].push_back(5);
	numbers[0].push_back(6);
	
	numbers[1].push_back(2);
	numbers[1].push_back(5);
	
	numbers[2].push_back(0);
	numbers[2].push_back(2);
	numbers[2].push_back(3);
	numbers[2].push_back(4);
	numbers[2].push_back(6);
	
	numbers[3].push_back(0);
	numbers[3].push_back(2);
	numbers[3].push_back(3);
	numbers[3].push_back(5);
	numbers[3].push_back(6);
	
	numbers[4].push_back(1);
	numbers[4].push_back(2);
	numbers[4].push_back(3);
	numbers[4].push_back(5);
	
	numbers[5].push_back(0);
	numbers[5].push_back(1);
	numbers[5].push_back(3);
	numbers[5].push_back(5);
	numbers[5].push_back(6);
	
	numbers[6].push_back(0);
	numbers[6].push_back(1);
	numbers[6].push_back(3);
	numbers[6].push_back(4);
	numbers[6].push_back(5);
	numbers[6].push_back(6);
	
	numbers[7].push_back(0);
	numbers[7].push_back(2);
	numbers[7].push_back(5);
	
	numbers[8].push_back(0);
	numbers[8].push_back(1);
	numbers[8].push_back(2);
	numbers[8].push_back(3);
	numbers[8].push_back(4);
	numbers[8].push_back(5);
	numbers[8].push_back(6);
	
	numbers[9].push_back(0);
	numbers[9].push_back(1);
	numbers[9].push_back(2);
	numbers[9].push_back(3);
	numbers[9].push_back(5);
	numbers[9].push_back(6);
}



int die_usage () {
	cout << "Usage : render_image <options> --bam=FILE --png=FILE --region=\"REGION\"" << endl ;
	cout << "-w --width   NUMBER   Image width (1024)" << endl ;
	cout << "-h --height  NUMBER   Image height (768)" << endl ;
	cout << "-v --vmin    NUMBER   Vertical lower limit (0)" << endl ;
	cout << "-V --vmax    NUMBER   Vertical upper limit (1000)" << endl ;
	cout << "-r --ref     FILE     Reference sequence" << endl ;
	cout << "-o --options TEXT     arrows,snps,pairs,single,faceaway" << endl ;
	return 0 ;
}


int main(int argc, char **argv) {
	for ( int a = 0 ; a < 256 ; a++ ) bam2char[a] = '?' ;
	bam2char[1] = 'A' ;
	bam2char[2] = 'C' ;
	bam2char[4] = 'G' ;
	bam2char[8] = 'T' ;
	bam2char[15] = 'N' ;

	string bam_file , ref_file , png_file , region , options ;
	int width = 1024 ;
	int height = 768 ;
	int vmin = 0 ;
	int vmax = 1000 ;
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
			case 'w' : width = atoi ( optarg ) ; break ;
			case 'h' : height = atoi ( optarg ) ; break ;
			case 'v' : vmin = atoi ( optarg ) ; break ;
			case 'V' : vmax = atoi ( optarg ) ; break ;
		}
	}
	
	init_numbers () ;
	
	if ( bam_file.empty() ) return die_usage () ;
	if ( png_file.empty() ) return die_usage () ;
	if ( region.empty() ) return die_usage () ;

	b2p = new Tbam2png () ;
	b2p->refseq_file = ref_file ;
	b2p->init( bam_file , region , png_file ) ;
	b2p->set_image_size ( width , height ) ;
	Tbam_draw_paired dp ( b2p ) ;
	dp.set_vertical_range ( vmin , vmax ) ;
	b2p->set_options ( options ) ;
	b2p->read_bam_file () ;

	b2p->create_png () ;
	dp.merge_all () ;

	b2p->write_png_file((char*)png_file.c_str());

	return 0;
}

/*
\rm render_image ; g++ render_image.cpp -O3 -o render_image -lpng -L . -lbam ; time ./render_image --bam=/nfs/users/nfs_m/mm6/ftp/ag/bam/AC0001-C.bam --options=snps,pairs,arrows,single,faceaway,inversions,linkpairs,colordepth --ref=/nfs/users/nfs_m/mm6/ftp/ag/Anopheles_gambiae.clean.fa --region="2L:1-200000" --png=2L.a.png


\rm render_image ; g++ render_image.cpp -O3 -o render_image -lpng -L . -lbam ; time ./render_image --bam="ftp://ftp.sanger.ac.uk/pub/team112/ag/bam/AC0001-C.bam" --options=pairs,arrows,single,faceaway,inversions,colordepth,snps --ref=/nfs/users/nfs_m/mm6/ftp/ag/Anopheles_gambiae.clean.fa --region="2L" --png=2L.a.png

*/
