/* xrw2pngmos.c
 *
 * Copyright 2012 David G. Barnes
 *
 * This file is part of S2VOLSURF.
 *
 * S2VOLSURF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * S2VOLSURF is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with S2VOLSURF.  If not, see <http://www.gnu.org/licenses/>. 
 *
 * We would appreciate it if research outcomes using S2VOLSURF would
 * provide the following acknowledgement:
 *
 * "Three-dimensional visualisation was conducted with the S2PLOT
 * progamming library"
 *
 * and a reference to
 *
 * D.G.Barnes, C.J.Fluke, P.D.Bourke & O.T.Parry, 2006, Publications
 * of the Astronomical Society of Australia, 23(2), 82-93.
 *
 * $Id: xrw2pngmos.c 224 2014-04-23 22:54:03Z barnesd $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <png.h>
#include "libxrw.c"

void usage(char *exename) {
  fprintf(stderr, "usage: %s [options] -f xrwfilename -o pngfilename\n", exename);
  fprintf(stderr, "-s s1 s2 s3\t\t set data averaging along each axis (1 1 1)\n");
  fprintf(stderr, "-R n\t\t use n rows in output image\n");
}

#if (0)
// png functions and associated structs learned from
// Ben Bullock's code @ http://www.lemoda.net/c/write-png/, which itself 
// was derived from "Sam"'s code @
// http://stackoverflow.com/questions/1821806/how-to-encode-png-to-buffer-using-libpng
typedef struct {
  uint8_t r,g,b;
} pixel_t;
typedef struct  {
  pixel_t *pix;
  size_t w, h;
} bitmap_t;
static int savePNG (FILE *fp, bitmap_t *bitmap);
#endif

int main(int argc, char *argv[]) {

  #define FNAMELEN 400
  char ifname[FNAMELEN+1];
  int haveifname = 0;
  char ofname[FNAMELEN+1];
  int haveofname = 0;
  int use_rows = 0;

  int stride[3] = {1,1,1}; //{3,3,1};

  if (argc < 5) {
    usage(argv[0]);
    return -1;
  }

  int ic = 1;
  while (ic < argc) {
    if (!strcmp(argv[ic], "-f")) {
      strncpy(ifname, argv[++ic], FNAMELEN);
      haveifname = 1;
    } else if (!strcmp(argv[ic], "-o")) {
      strncpy(ofname, argv[++ic], FNAMELEN);
      haveofname = 1;
    } else if (!strcmp(argv[ic], "-s")) {
      stride[0] = atoi(argv[++ic]);
      stride[1] = atoi(argv[++ic]);
      stride[2] = atoi(argv[++ic]);
      if (stride[0] < 1 || stride[1] < 1 || stride[2] < 1) {
	fprintf(stderr, "stride must be at least 1 in each direction!\n");
	exit(-1);
      }
    } else if (!strcmp(argv[ic], "-R")) {
      use_rows = atoi(argv[++ic]);
    }
    ic++;
  }

  if (!haveifname || !haveofname) {
    usage(argv[0]);
    return -1;
  }
    
  XRAW_STRUCT *xr = loadXraw(ifname);

  if (!xr) {
    fprintf(stderr, "Failed to open or read '%s'.\n", argv[1]);
    return -1;
  }

  showXraw(xr);
  fprintf(stdout, "- - - - - - - - - - - - - - - - - - - - - - - - -\n");
  
  VOL_STRUCT *vol = Xraw2Xvol(xr, stride);
  if (!vol) {
    fprintf(stderr, "Failed to parse data volume.\n");
    return -1;
  }

  showXvol(vol);

  int nx, ny, nz;
  nx = vol->nx;
  ny = vol->ny;
  nz = vol->nz;

  bitmap_t bits;

  if (use_rows == 1) {
    bits.w = nx * nz;
    bits.h = ny;
    bits.pix = (pixel_t *)malloc(nx*ny*nz*sizeof(pixel_t));
    
    int i, j, k;
    float val;
    uint8_t grey;
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
	for (k = 0; k < nz; k++) {
	  val = vol->data[i][j][k];
	  grey = val * 255.0;
	  bits.pix[(nx*nz)*j + (k*nx+i)].r = grey;
	  bits.pix[(nx*nz)*j + (k*nx+i)].g = grey;
	  bits.pix[(nx*nz)*j + (k*nx+i)].b = grey;
	}
      }
    }

  } else {
   
    int Nzx, Nzy;
    if (use_rows > 1) {
      Nzy = use_rows;
      Nzx = nz/Nzy;
      if (nz % Nzy) {
	Nzx++;
      }
    } else {
      Nzx = (int)(floorf(sqrtf((float)nz)));
      if (Nzx * Nzx == nz) {
	Nzy = Nzx;
      } else {
	Nzx++;
	Nzy = (nz / Nzx) + 1;
      }
    }

    fprintf(stderr, "mosaic dimensions are: %d %d\n", Nzx, Nzy);

    bits.w = nx * Nzx;
    bits.h = ny * Nzy;
    bits.pix = (pixel_t *)malloc(nx * ny * Nzx * Nzy * sizeof(pixel_t));

    int i, j, k;
    int izx, izy;
    float val;
    uint8_t grey;
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
	for (k = 0; k < nz; k++) {
	  val = vol->data[i][j][k];
	  grey = val * 255.0;
	  izy = k / Nzx;
	  izx = k % Nzx;
	  bits.pix[(nx*Nzx)*(izy*ny+j) + (izx*nx+i)].r = grey;
	  bits.pix[(nx*Nzx)*(izy*ny+j) + (izx*nx+i)].g = grey;
	  bits.pix[(nx*Nzx)*(izy*ny+j) + (izx*nx+i)].b = grey;
	}
      }
    }
  }

  FILE *fout = fopen(ofname, "w");
  if (fout) {
    savePNG(fout, &bits);
    fclose(fout);
  } else {
    fprintf(stderr, "Failed to open '%s' for writing.\n", ofname);
    return -1;
  }

  return 0;
}

#if (0)
static int savePNG (FILE *fp, bitmap_t *bitmap) {
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;
  int i, j;
  png_byte **row_pointers = NULL;
  int status = -1;
  
  png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) {
    return status;
  }
  info_ptr = png_create_info_struct (png_ptr);
  if (!info_ptr) {
    return status;
  }
  
  /* Set image attributes. */  
  png_set_IHDR(png_ptr,	info_ptr, bitmap->w, bitmap->h,	8,
	       PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  
  /* Initialize rows of PNG. */
  row_pointers = png_malloc (png_ptr, bitmap->h * sizeof (png_byte *));
  for (j = 0; j < bitmap->h; ++j) {
    png_byte *row = png_malloc (png_ptr, 3 * sizeof (uint8_t) * bitmap->w);
    row_pointers[j] = row;
    for (i = 0; i < bitmap->w; ++i) {
      pixel_t *pixel = bitmap->pix + bitmap->w * j + i;
      *row++ = pixel->r;
      *row++ = pixel->g;
      *row++ = pixel->b;
    }
  }
  
  /* Write the image data to "fp". */
  png_init_io (png_ptr, fp);
  png_set_rows (png_ptr, info_ptr, row_pointers);
  png_write_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
  
  status = 0;
  for (j = 0; j < bitmap->h; j++) {
    png_free (png_ptr, row_pointers[j]);
  }
  png_free (png_ptr, row_pointers);
  
  return status;
}

#endif
