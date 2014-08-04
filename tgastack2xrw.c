/* tgastack2xrw.c
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
  fprintf(stderr, "usage: %s [options] -f infileformat.tga -o xrwfilename\n", exename);
  fprintf(stderr, " * * * (infileformat should include integer specifier e.g. %%04d)\n");
  fprintf(stderr, " -s s1 s2 s3\t\t set data averaging along each axis (1 1 1)\n");
  fprintf(stderr, " -w wx wy wz\t\t set world size of pixels (1.0 1.0 1.0)\n");
}

#define FNAMELEN 400

int main(int argc, char *argv[]) {

  char ifname[FNAMELEN+1];
  int haveifname = 0;
  char ofname[FNAMELEN+1];
  int haveofname = 0;
  float wdelta[] = {1.,1.,1.};
  int stride[] = {1,1,1};
  
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
    } else if (!strcmp(argv[ic], "-w")) {
      wdelta[0] = atof(argv[++ic]);
      wdelta[1] = atof(argv[++ic]);
      wdelta[2] = atof(argv[++ic]);
    }
    ic++;
  }
  
  if (!haveifname || !haveofname) {
    usage(argv[0]);
  }

  // find start/end of filename indices
  fprintf(stderr, "Finding start/end indices for stack: check inputs if this hangs!\n");
  int istart=0, iend;
  char tfname[FNAMELEN+1];
  sprintf(tfname, ifname, istart);
  FILE *tf = fopen(tfname, "r");
  while (!tf) {
    //fclose(tf);
    istart++;
    sprintf(tfname, ifname, istart);
    //fprintf(stderr, "istart = %d\n", istart);
    //fprintf(stderr, "tfname = %s\n", tfname);
    tf = fopen(tfname, "r");
  }
  fclose(tf);
  iend = istart + 1;
  sprintf(tfname, ifname, iend);
  tf = fopen(tfname, "r");
  while (tf) {
    fclose(tf);
    iend++;
    sprintf(tfname, ifname, iend);
    //fprintf(stderr, "iend = %d\n", iend);
    //fprintf(stderr, "tfname = %s\n", tfname);
    tf = fopen(tfname, "r");
  }
  iend--;
  
  fprintf(stderr, "File index appears to go from %d to %d\n", istart, iend);
  
  VOL_STRUCT *tgavol = loadTGAstack(ifname, istart, iend, 1);
  tgavol->wdx = wdelta[0];
  tgavol->wdy = wdelta[1];
  tgavol->wdz = wdelta[2];

  showXvol(tgavol);
  VOL_STRUCT *tgavol_sm;
  if (stride[0] + stride[1] + stride[2] > 3) {
    fprintf(stderr, "rebinning...\n");
    tgavol_sm = rebinXvol(tgavol, stride);
  } else {
    tgavol_sm = tgavol;
  }
  fprintf(stderr, "normalising...\n");
  normaliseXvol(tgavol_sm);
  showXvol(tgavol_sm);

  fprintf(stderr, "converting...\n");
  XRAW_STRUCT *xraw = Xvol2Xraw(tgavol_sm);
  sprintf(xraw->filename, ofname);
  showXraw(xraw);
  fprintf(stderr, "writing...\n");
  saveXraw(xraw);

  return 0;
}
