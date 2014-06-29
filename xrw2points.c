/* xrw2points.c
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
 * $Id: xrw2points.c 74 2012-11-26 02:39:40Z barnesd $
 *
 */

#include "libxrw.c"

void usage(char *exename) {
  fprintf(stderr, "usage: %s options -f xrwfilename -o pointsfilename\n", exename);
  fprintf(stderr, "options:\n");
  fprintf(stderr, "-s s1 s2 s3    set data averaging along each axis (1 1 1)\n");
  fprintf(stderr, "-d d1 d2       set data min, max (default 0.0 1.0)\n");
  fprintf(stderr, "-z             calculate and render volume derivative\n");
}

int main(int argc, char **argv) {

  #define FNAMELEN 400
  char ifname[FNAMELEN+1], ofname[FNAMELEN+1];
  int haveifname = 0, haveofname = 0;
  float dmin = 0, dmax = 1;
  int stride[] = {1,1,1};
  int doDeriv = 0;

  if (argc < 5) {
    usage(argv[0]);
    exit(-1);
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
    } else if (!strcmp(argv[ic], "-d")) {
      dmin = atof(argv[++ic]);
      dmax = atof(argv[++ic]);
    } else if (!strcmp(argv[ic], "-z")) {
      doDeriv = 1;
    }
    ic++; 
  }
  
  if (!haveifname || !haveofname) {
    usage(argv[0]);
    return -1;
  }

  XRAW_STRUCT *xr = loadXraw(ifname);

  VOL_STRUCT *vol = Xraw2Xvol(xr, stride);

  if (doDeriv) {
    derivXvol(vol);
  }

  showXvol(vol);

  int nx, ny, nz;
  nx = vol->nx;
  ny = vol->ny;
  nz = vol->nz;
  
  FILE *fout = fopen(ofname, "w");
  if (!fout) {
    fprintf(stderr, "Failed to open output file '%s'.\n", ofname);
    return(-1);
  }
	   
  int i, j, k;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
	float dval = vol->data[i][j][k];
	if ((dval >= dmin) && (dval <= dmax)) {
	  dval = (dval - dmin) / (dmax - dmin);
	  if (dval < 0.0) {
	    dval = 0.0;
	  } else if (dval > 1.0) {
	    dval = 1.0;
	  }
	  fprintf(fout, "%6.4f,%6.4f,%6.4f,%6.4f\n", 
		  (float)i/(float)(nx-1),
		  (float)j/(float)(ny-1) * vol->wdy*(float)ny / (vol->wdx*(float)nx),
		  (float)k/(float)(nz-1) * vol->wdz*(float)nz / (vol->wdx*(float)nx),
		  dval);
	}
      }
    }
  }
  fclose(fout);

  return 0;
}
