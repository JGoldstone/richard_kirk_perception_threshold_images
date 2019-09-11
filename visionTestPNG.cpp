#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <parse.h>
#include <time.h>
#include <png.h>

/*
  Make PPM 'T' image for standard monitor 
  for testing colour vision, assuming SonyHD
  monitor with D65 whitepoint.
*/

int main(int argc, char *argv[]) {
  double col[25][13][3], L, dL;
  double da=0.0, db=0.0, dP;

  // confusion points for defective vision
  double xc[3] = {0.747, 1.080, 0.171};
  double yc[3] = {0.253, -0.080, 0.00};
  double white[3] = {0.950173, 1.0, 1.08726}; // XYZ at 1 ft-L D65
  double shift[3], a, b, de, contrast;
  int seed;

  FILE *out;

  int x, y, i, nx, ny, NX, NY, square, border, width, height, n;
  int tristimulus = -1;

  argc = parseArgs(argc, argv);
  if (argc != 2 || argFlag['i']) {
    printf("%s [options] <output_file> makes 3x6 PPM 'T' test image\n"
	   "options:\n"
	   "  L=%%f  luminance        (default 70)\n"
	   "  d=%%f  luminance dither (default 15)\n"
	   "  p=%%f  pixel RGB dither (default 3.0)\n"
	   "  S=%%d  square size      (default 15)\n"
	   "  B=%%d  border size      (default 1)\n"
	   "  r=%%d  random seed      (default uses time)\n"
	   "  -s    make solution image, all patches with large contrast\n"
	   "  -p    output a PNG image (default PPM)\n"
	   "  -i    info message (this is it)\n"
	   , argv[0] );
    exit(1);
  }

  /* Seed random function with time... */
  seed = (int)time(NULL);
  srand( parseInteger('r', seed, seed) );

  out = fopen(argv[1], "w");
  if (out == NULL) {
    fprintf(stderr, "Error: unable to open file %s.\n", argv[1]);
    exit(1);
  }

  L = parseDouble('L', 70.0, -1.0);
  if (L <= 10.0) {
    fprintf(stderr, "Error: illegal or malformed 'L=' argument\n");
    exit(1);
  }
 
  dL = parseDouble('d', 15.0, -1.0);
  if (L <= 0.0) {
    fprintf(stderr, "Error: illegal or malformed 'd=' argument\n");
    exit(1);
  }
  L -= dL/2.0;
  dL /= (double)RAND_MAX;

  dP = parseDouble('p', 3.0, -1.0);
  if (dP <= 0.0) {
    fprintf(stderr, "Error: illegal or malformed 'p=' argument\n");
    exit(1);
  }
  dP /= (double)RAND_MAX;

  square = parseInteger('S', 15, -1);
  if (square < 1) {
    fprintf(stderr, "Error: illegal or malformed 'S=' argument\n");
    exit(1);
  }

  border = parseInteger('B', 1, -1);
  if (border < 0) {
    fprintf(stderr, "Error: illegal or malformed 'B=' argument\n");
    exit(1);
  }

  width = 25*square + 26*border;
  height = 13*square + 14*border;

  /* Set luminance in first channel, zero others... */
  for (x=0; x<25; ++x) {
    for (y=0; y<13; ++y) {
      col[x][y][0] = L + dL*rand();
      col[x][y][1] = 0.0;            
      col[x][y][2] = 0.0;     
    }
  }

  /* Make 'T' symbol in second channel ... */
  for (NX=0; NX<6; ++NX) {
    double del =  argFlag['s'] ? 64.0 : (double)(32 >> NX);
    for (NY=0; NY<3; ++NY) {
      switch (rand() % 4) { // Clockwise from vertical
      case 0:
	for (x=1; x<=3; ++x) col[NX*4+x][NY*4+1][1] = del;
	for (y=1; y<=3; ++y) col[NX*4+2][NY*4+y][1] = del;
	break;
      case 1:
	for (x=1; x<=3; ++x) col[NX*4+x][NY*4+2][1] = del;
	for (y=1; y<=3; ++y) col[NX*4+3][NY*4+y][1] = del;
	break;
      case 2:
	for (x=1; x<=3; ++x) col[NX*4+x][NY*4+3][1] = del;
	for (y=1; y<=3; ++y) col[NX*4+2][NY*4+y][1] = del;
	break;
      case 3:
	for (x=1; x<=3; ++x) col[NX*4+x][NY*4+2][1] = del;
	for (y=1; y<=3; ++y) col[NX*4+1][NY*4+y][1] = del;
	break;
      }
      del = -del;
    }
  }

  /*
    Generate a & b shifts.
    Add a small fraction of the tristimulus confusion point
    to the white point. Calculate the shift in a, and b.
  */
  for (NY=0; NY<3; ++NY) {
    shift[0] = xc[NY];
    shift[1] = yc[NY];
    shift[2] = 1.0 - shift[0] - shift[1];
    for (n=0; n<3; ++n) shift[n] += 20.0*white[n];
    shift[0] /= shift[1];
    shift[2] /= shift[1];
    shift[1] = 1.0;

    for (n=0; n<3; ++n) shift[n] = pow(shift[n]/white[n], 1.0/3.0);

    a = 500.0*(shift[0] - shift[1]);
    b = 200.0*(shift[1] - shift[2]);

    de = sqrt(a*a + b*b);
    da = a/de;
    db = b/de;

    for (y=(4*NY); y<(4*NY+4); ++y) {
      for (x=0; x<25; ++x) {
	col[x][y][2] = col[x][y][1]*db;
	col[x][y][1] *= da;
      }          
    }
  }

  /* Convert Lab to Rec.1886 RGB */
  for (x=0; x<25; ++x) {
    for (y=0; y<13; ++y) {
      double L, a, b, fX, fY, fZ, X, Y, Z, R, G, B;
      L = col[x][y][0];
      a = col[x][y][1];
      b = col[x][y][2];
 
      fY = (L+16.0)/116.0;
      fX = fY + a/500.0;
      fZ = fY - b/200.0;
 
      X = (fX < 24.0/116.0) ? (fX - 16.0/116.0)/7.787 : fX*fX*fX;
      Y = (fY < 24.0/116.0) ? (fY - 16.0/116.0)/7.787 : fY*fY*fY;
      Z = (fZ < 24.0/116.0) ? (fZ - 16.0/116.0)/7.787 : fZ*fZ*fZ;
 
      R = +3.080014*X -1.537119*Y -0.542894*Z;
      G = -0.921307*X +1.876051*Y +0.045256*Z;
      B = 0.0528790*X -0.203989*Y +1.151110*Z;
     
      R = R<0.0 ? 0.0 : pow(R, 1.0/2.6);
      G = G<0.0 ? 0.0 : pow(G, 1.0/2.6);
      B = B<0.0 ? 0.0 : pow(B, 1.0/2.6);

      col[x][y][0] = 255.0*R;
      col[x][y][1] = 255.0*G;
      col[x][y][2] = 255.0*B;
    }
  }

  /* Replicate, dither, then output as 8-bit image */
  if (argFlag['p']) {     /* Output PNG */
    /* Initialize pnglib memory and error handling */
    unsigned char * * graph_;
    graph_ = (png_bytepp)malloc(height * sizeof(png_bytep));
    for (n = 0; n <= (height+10); n++) {
      graph_[n] = (png_bytep)malloc(3*width * sizeof(png_byte));
    }
  
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
      fprintf(stderr, "Error: PNGLIB initialization error.\n");
      exit(1);
    }
    
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
      png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
      fprintf(stderr, "Error: PNGLIB info initialization error.\n");
      exit(1);
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
      png_destroy_write_struct(&png_ptr, &info_ptr);
      fclose(out);
      fprintf(stderr, "Error: Unknown PNGLIB error.\n");
      exit(1);
    }

    png_init_io(png_ptr, out);

    png_set_IHDR(png_ptr, info_ptr, width, height,
        8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    png_set_filter(png_ptr, PNG_FILTER_TYPE_BASE, PNG_ALL_FILTERS);

    time_t gmt;
    png_time mod_time;
    png_text text_ptr[3];
    time(&gmt);
    png_convert_from_time_t(&mod_time, gmt);
    png_set_tIME(png_ptr, info_ptr, &mod_time);
    text_ptr[0].key = strdup("Title");
    text_ptr[0].text = strdup("Vision Test");
    text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[1].key = strdup("Author");
    text_ptr[1].text = strdup("FilmLight Ltd.");
    text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[2].key = strdup("Creation Time");
    text_ptr[2].text = strdup(png_convert_to_rfc1123(png_ptr, &mod_time));
    text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;
    png_set_text(png_ptr, info_ptr, text_ptr, 3);
 
    /* Replicate, dither, then store pixel values */
    int pixelpos=0;
    int width_t3=width*3;
    for (y=0; y<13; ++y) {
      for (n=3*border*width; n; --n) {
        graph_[(int)floor(pixelpos/width_t3)][(pixelpos%width_t3)] = (char) 0;
        pixelpos++;
      }
      for (ny=0; ny<square; ++ny) {
        for (x=0; x<25; ++x) {
          for (n=3*border; n; --n) {
            graph_[(int)floor(pixelpos/width_t3)][(pixelpos%width_t3)] = (char) 0;
            pixelpos++;
          }
          for (nx=0; nx<square; ++nx) {
            for (n=0; n<3; ++n) {
              int v = (int)floor(col[x][y][n] + dP*rand());
              graph_[(int)floor(pixelpos/width_t3)][pixelpos%width_t3] = (char) ( (v<0) ? 0 : (v>255) ? 255 : v );
              pixelpos++;
            }
          }
        }
        for (n=3*border; n; --n) {
          graph_[(int)floor(pixelpos/width_t3)][pixelpos%width_t3] = (char) 0;
          pixelpos++;
        }
      }
    }
    for (n=3*border*width; n; --n) {
      graph_[(int)floor(pixelpos/width_t3)][pixelpos%width_t3] = (char) 0;
      pixelpos++;
    }

    /* Write PNG file */
    png_write_info(png_ptr, info_ptr);
    png_write_image(png_ptr, graph_);
    png_write_end(png_ptr, info_ptr);
    // fclose(out);

    png_destroy_write_struct(&png_ptr, &info_ptr);
  } else {     /* Output PPM */
    fprintf(out, "P6 %d %d 255\n", width, height); 
    for (y=0; y<13; ++y) {
      for (n=3*border*width; n; --n) fputc(0, out);
      for (ny=0; ny<square; ++ny) {
        for (x=0; x<25; ++x) {
          for (n=3*border; n; --n) fputc(0, out);
	      for (nx=0; nx<square; ++nx) {
            for (n=0; n<3; ++n) {
              int v = (int)floor(col[x][y][n] + dP*rand());
              fputc( (v<0) ? 0 : (v>255) ? 255 : v, out);
            }
          }
        }
        for (n=3*border; n; --n) fputc(0, out);
      }
    }
    for (n=3*border*width; n; --n) fputc(0, out);
    fclose(out);
  }
}
