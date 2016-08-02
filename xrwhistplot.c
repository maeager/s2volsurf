/* xrwhistoplot.c
 *
 * Michael Eager 2016
 */

/* xrwhisto.c
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
 * $Id: xrwhisto.c 70 2012-11-22 23:00:34Z barnesd $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libxrw.c"

float datastats[4];

void getstats(VOL_STRUCT *xv) {
  fprintf(stdout, "Parsed data (from '%s'):\n", xv->filename);
  fprintf(stdout, "volume dimensions: %6d %6d %6d\n", xv->nx, xv->ny, xv->nz);
  fprintf(stdout, "      voxel sizes: %6.3f %6.3f %6.3f\n", xv->wdx, xv->wdy, xv->wdz);
  fprintf(stdout, "total volume size: %6.1f %6.1f %6.1f\n", xv->wdx*xv->nx, xv->wdy*xv->ny,
	  xv->wdz*xv->nz);
  int i, j, k;
  float v, min, max, mean, meansq;
  min = max = xv->data[0][0][0];
  mean = 0.0;
  meansq = 0.0;
  for (i = 0; i < xv->nx; i++) {
    for (j = 0; j < xv->ny; j++) {
      for (k = 0; k < xv->nz; k++) {
	v = xv->data[i][j][k];
	if (v < min) {
	  min = v;
	}
	if (v > max) {
	  max = v;
	}
	mean += v;
	meansq += v*v;
      }
    }
  }
  meansq /= (float)(xv->nx * xv->ny * xv->nz);
  mean /= (float)(xv->nx * xv->ny * xv->nz);
  float rms = sqrt(meansq - mean*mean);
  fprintf(stdout, "    min, mean(stddev), max: %6.4f %6.4f(%6.4f) %6.4f\n", 
	  min, mean, rms, max);

  datastats[0]=min;
  datastats[1]=max;
  datastats[2]=rms;
  datastats[3]=mean;
}



int main(int argc, char **argv) {
  
  if (argc != 2) {
    fprintf(stderr, "Usage: %s xrwfilename\n", argv[0]);
    return -1;
  }

  XRAW_STRUCT *xr = loadXraw(argv[1]);
  if (!xr) {
    fprintf(stderr, "Failed to open or read '%s'.\n", argv[1]);
    return -1;
  }
  
  int nbin = 256;
  float x[256];
  float data[256];
  int centre = 1;
  
  int i;
  for (i = 0; i < 256; i++) {
    x[i] = ((float)i + 0.5) / 256.0;
    data[i] = 0.;
  }
  
  for (i = 0; i < xr->nx * xr->ny * xr->nz; i++) {
    data[xr->data[i]] += 1.0;
  }

  for (i = 0; i < 256; i++) {
    data[i] = log10f(data[i]+1);
  }

  fprintf(stdout,"Normalised data value\tLog(count)\t%s\n", argv[1]);
  for (i = 0; i < 256; i++) {
    fprintf(stdout,"%4.8f\t%4.8f\n", x[i],data[i]);
  }

  return 0;
}
