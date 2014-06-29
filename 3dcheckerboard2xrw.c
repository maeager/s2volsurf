/* 3dcheckerboard2xrw.c
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
 * $Id: tgastack2xrw.c 70 2012-11-22 23:00:34Z barnesd $
 *
 */

// refer Kim, Wittenbrink & Pang, "Extended Specifications and Test
// Data for Data Level Comparisons of Direct Volume Rendering Algorithms,
// IEEE Trans Viz & Comp Graphics (2001), 7, 299.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "libxrw.c"

void usage(char *exename) {
  fprintf(stderr, "usage: %s [options] -d nx ny nz -o xrwfilename\n", exename);
  fprintf(stderr, " -v low high\t\t set low, high data values (0.0 1.0)\n");
  fprintf(stderr, " -p\t\t set width of pattern (in integer voxels) (3)\n");
  fprintf(stderr, " -w wx wy wz\t\t set world size of pixels (1.0 1.0 1.0)\n");
}

#define FNAMELEN 400

int main(int argc, char *argv[]) {

  char ofname[FNAMELEN+1];
  int haveofname = 0;
  int dims[] = {0, 0, 0};
  int havedims = 0;
  float vals[] = {0., 1.};
  int width = 1;
  float wdelta[] = {1.,1.,1.};

  
  if (argc < 7) {
    usage(argv[0]);
    return -1;
  }

  int ic = 1;
  while (ic < argc) {
    if (!strcmp(argv[ic], "-o")) {
      strncpy(ofname, argv[++ic], FNAMELEN);
      haveofname = 1;
    } else if (!strcmp(argv[ic], "-d")) {
      dims[0] = atoi(argv[++ic]);
      dims[1] = atoi(argv[++ic]);
      dims[2] = atoi(argv[++ic]);
      if (dims[0] < 1 || dims[1] < 1 || dims[2] < 1) {
	fprintf(stderr, "dimensions too small, please >= 1 in every direction!\n");
	exit(-1);
      }
      havedims = 1;
    } else if (!strcmp(argv[ic], "-v")) {
      vals[0] = atoi(argv[++ic]);
      vals[1] = atoi(argv[++ic]);
    } else if (!strcmp(argv[ic], "-p")) {
      width = atoi(argv[++ic]);
    } else if (!strcmp(argv[ic], "-w")) {
      wdelta[0] = atof(argv[++ic]);
      wdelta[1] = atof(argv[++ic]);
      wdelta[2] = atof(argv[++ic]);
    }
    ic++;
  }
  
  if (!haveofname || !havedims) {
    usage(argv[0]);
  }

  VOL_STRUCT *vol = (VOL_STRUCT *)malloc(1 * sizeof(VOL_STRUCT));

  sprintf(vol->filename, "generated_3d_checkerboard");

  vol->nx = dims[0];
  vol->ny = dims[1];
  vol->nz = dims[2];
  
  vol->wdx = wdelta[0];
  vol->wdy = wdelta[1];
  vol->wdz = wdelta[2];
  
  vol->bzero = 0.0;
  vol->bscale = 1.0;

  int i, j, k;
  
  vol->data = (float ***)malloc(vol->nx * sizeof(float **));
  for (i = 0; i < vol->nx; i++) {
    vol->data[i] = (float **)malloc(vol->ny * sizeof(float *));
    for (j = 0; j < vol->ny; j++) {
      vol->data[i][j] = (float *)malloc(vol->nz * sizeof(float));
      for (k = 0; k < vol->nz; k++) {
	if (((i + j + k) / width) % 2) {
	  vol->data[i][j][k] = vals[0];
	} else {
	  vol->data[i][j][k] = vals[1];
	}
      }
    }
  }
  
  int npins = 8;
  float pinpoints[] = {0.0, 29.0, 32.0, 36.0, 127.0, 127.5, 191.25, 255.0};
  COLOUR table[] = {{0.69, 0.68, 0.67},
		    {0.09, 0.10, 0.23},
		    {0.03, 0.09, 0.18},
		    {0.03, 0.09, 0.12},
		    {0.00, 0.06, 0.01},
		    {0.00, 0.05, 0.00},
		    {1.00, 0.03, 0.00},
		    {1.00, 0.00, 0.00}};
  //float alphas[] = {0.5, 0.391, 0.380, 0.365, 0.024, 0.025, 0.263, 0.5};
  
  int pinidx = 0; // index of pin on right of i
  for (i = 0; i < 256; i++) {
    while ((float)i <= pinpoints[pinidx] && pinidx < npins-1) {
      pinidx++;
    }
    if ((pinidx > npins-1) || (pinidx < 1)) {
      fprintf(stderr, "Error finding pin for colortable!\n");
      exit(-1);
    }
    float frac = ((float)i - pinpoints[pinidx-1]) / (pinpoints[pinidx] - pinpoints[pinidx-1]);
    vol->red[i]   = table[pinidx-1].r * (1. - frac) + table[pinidx].r * frac;
    vol->green[i] = table[pinidx-1].g * (1. - frac) + table[pinidx].g * frac;
    vol->blue[i]  = table[pinidx-1].b * (1. - frac) + table[pinidx].b * frac;
  }
    
  showXvol(vol);
  VOL_STRUCT *vol_sm = vol; //rebinXvol(vol, stride);
  normaliseXvol(vol_sm);
  showXvol(vol_sm);

  XRAW_STRUCT *xraw = Xvol2Xraw(vol_sm);
  sprintf(xraw->filename, ofname);
  showXraw(xraw);
  saveXraw(xraw);

  return 0;
}
