/* xrw2xrw.c
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
 * $Id$
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
  fprintf(stderr, "usage: %s [options] -f xrwfilename -o xrwfilename\n", exename);
  fprintf(stderr, "-t t1 t2 t3\t\t set edge trim (centred) along each axis (0 0 0) pre-stride\n");
  fprintf(stderr, "-s s1 s2 s3\t\t set data averaging along each axis (1 1 1)\n");
  fprintf(stderr, "-2         \t\t force dimensions to power-of-two\n");
  fprintf(stderr, "-1G        \t\t ensure image has fewer than 1G pixels\n");
}


int main(int argc, char *argv[]) {

  #define FNAMELEN 400
  char ifname[FNAMELEN+1];
  int haveifname = 0;
  char ofname[FNAMELEN+1];
  int haveofname = 0;

  int trim[3] = {0,0,0};
  int stride[3] = {1,1,1}; //{3,3,1};
  int pow2 = 0;
  int max1g = 0;

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
    } else if (!strcmp(argv[ic], "-t")) {
      trim[0] = atoi(argv[++ic]);
      trim[1] = atoi(argv[++ic]);
      trim[2] = atoi(argv[++ic]);
      if (trim[0] < 0 || trim[1] < 1 || trim[2] < 1) {
	fprintf(stderr, "trim must be at least 0 in each direction!\n");
	exit(-1);
      }
    } else if (!strcmp(argv[ic], "-s")) {
      stride[0] = atoi(argv[++ic]);
      stride[1] = atoi(argv[++ic]);
      stride[2] = atoi(argv[++ic]);
      if (stride[0] < 1 || stride[1] < 1 || stride[2] < 1) {
	fprintf(stderr, "stride must be at least 1 in each direction!\n");
	exit(-1);
      }
    } else if (!strcmp(argv[ic], "-2")) {
      pow2 = 1;
    } else if (!strcmp(argv[ic], "-1G")) {
      max1g = 1;
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

  if (trim[0] + trim[1] + trim[2] > 1) {
    XRAW_STRUCT *xrt = trimXraw(xr, trim);
    if (xrt) {
      deleteXraw(xr);
      xr = xrt;
    }
  }

  if (max1g) {
    while (((long)xr->nx / (long)stride[0] * 
	    (long)xr->ny / (long)stride[1] * 
	    (long)xr->nz / (long)stride[2]) >
	   (long)1024 * (long)1024 * (long)1024) {
      int maxidx = 0;
      int maxn = xr->nx;
      if (xr->ny/stride[1] > xr->nx/stride[0]) {
	maxidx = 1;
	maxn = xr->ny;
      }
      if (xr->nz/stride[2] > maxn/stride[maxidx]) {
	maxidx = 2;
      }
      stride[maxidx] += 1;
    }
  }
  

  fprintf(stderr, "converting...\n");
  VOL_STRUCT *vol = Xraw2Xvol(xr, stride);
  if (!vol) {
    fprintf(stderr, "Failed to parse data volume.\n");
    return -1;
  }

  showXvol(vol);

  if (pow2) {
    fprintf(stderr, "converting to pow2-format...\n");
    VOL_STRUCT *vol2 = makePow2Xvol(vol);
    vol = vol2;
  }


  fprintf(stderr, "converting...\n");
  XRAW_STRUCT *xro = Xvol2Xraw(vol);
  sprintf(xro->filename, ofname);
  showXraw(xro);
  fprintf(stderr, "writing...\n");
  saveXraw(xro);

  return 0;
}
