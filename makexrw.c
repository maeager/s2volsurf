/* makexrw.c
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
  fprintf(stderr, "usage: %s [options] -o xrwfilename\n", exename);
  fprintf(stderr, "-b val\t\t set border thickness (0 for no border)\n");
  fprintf(stderr, "-d nx ny nz\t\t set data dimensions (128 128 128)\n");
}


int main(int argc, char *argv[]) {

  #define FNAMELEN 400
  char ofname[FNAMELEN+1];
  int haveofname = 0;

  int border = 0;
  int dim[3] = {128,128,128};

  if (argc < 3) {
    usage(argv[0]);
    return -1;
  }

  int ic = 1;
  while (ic < argc) {
    if (!strcmp(argv[ic], "-o")) {
      strncpy(ofname, argv[++ic], FNAMELEN);
      haveofname = 1;
    } else if (!strcmp(argv[ic], "-d")) {
      dim[0] = atoi(argv[++ic]);
      dim[1] = atoi(argv[++ic]);
      dim[2] = atoi(argv[++ic]);
      if (dim[0] < 1 || dim[1] < 1 || dim[2] < 1) {
	fprintf(stderr, "dim must be at least 1 in each direction!\n");
	exit(-1);
      }
    } else if (!strcmp(argv[ic], "-b")) {
      border = atoi(argv[++ic]);
    }
    ic++;
  }

  if (!haveofname) {
    usage(argv[0]);
    return -1;
  }
    
  XRAW_STRUCT *xr = createXraw(ofname, dim[0], dim[1], dim[2], border);

  if (!xr) {
    fprintf(stderr, "Failed to create new xraw object.\n");
    return -1;
  }

  showXraw(xr);
  fprintf(stdout, "- - - - - - - - - - - - - - - - - - - - - - - - -\n");

  fprintf(stderr, "writing...\n");
  saveXraw(xr);

  return 0;
}
