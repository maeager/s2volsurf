/* ushortraw2xrw.c
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "libxrw.c"

void usage(char *exename) {
  fprintf(stderr, "usage: %s [options] -f ushort_raw_infile -d nx ny nz -o xrwfilename\n", exename);
  fprintf(stderr, " -s s1 s2 s3\t\t set data averaging along each axis (1 1 1)\n");
  fprintf(stderr, " -w wx wy wz\t\t set world size of pixels (1.0 1.0 1.0)\n");
}

#define FNAMELEN 400

int main(int argc, char *argv[]) {

  char ifname[FNAMELEN+1];
  int haveifname = 0;
  char ofname[FNAMELEN+1];
  int haveofname = 0;
  int dims[] = {0, 0, 0};
  int havedims = 0;
  float wdelta[] = {1.,1.,1.};
  int stride[] = {1,1,1};
  
  if (argc < 7) {
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
    } else if (!strcmp(argv[ic], "-d")) {
      dims[0] = atoi(argv[++ic]);
      dims[1] = atoi(argv[++ic]);
      dims[2] = atoi(argv[++ic]);
      if (dims[0] < 1 || dims[1] < 1 || dims[2] < 1) {
	fprintf(stderr, "dimensions too small, please >= 1 in every direction!\n");
	exit(-1);
      }
      havedims = 1;
    } else if (!strcmp(argv[ic], "-s")) {
      stride[0] = atoi(argv[++ic]);
      stride[1] = atoi(argv[++ic]);
      stride[2] = atoi(argv[++ic]);
      if (stride[0] < 1 || stride[1] < 1 || stride[2] < 1) {
	fprintf(stderr, "stride must be at least 1 in each direction!\n");
	exit(-1);
      }
    } else if (!strcmp(argv[ic], "-w")) {
      wdelta[0] = atof(argv[++ic]);
      wdelta[1] = atof(argv[++ic]);
      wdelta[2] = atof(argv[++ic]);
    }
    ic++;
  }
  
  if (!haveifname || !haveofname || !havedims) {
    usage(argv[0]);
  }

  VOL_STRUCT *rawvol = loadUshortRaw(ifname, dims[0], dims[1], dims[2]);
  rawvol->wdx = wdelta[0];
  rawvol->wdy = wdelta[1];
  rawvol->wdz = wdelta[2];

  showXvol(rawvol);
  VOL_STRUCT *rawvol_sm = rebinXvol(rawvol, stride);
  normaliseXvol(rawvol_sm);
  showXvol(rawvol_sm);

  XRAW_STRUCT *xraw = Xvol2Xraw(rawvol_sm);
  sprintf(xraw->filename, ofname);
  showXraw(xraw);
  saveXraw(xraw);

  return 0;
}
