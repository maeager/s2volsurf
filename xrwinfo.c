/* xrwinfo.c
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
 * $Id: xrwinfo.c 70 2012-11-22 23:00:34Z barnesd $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libxrw.c"

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s xrwfilename\n", argv[0]);
    return -1;
  }

  XRAW_STRUCT *xr = preloadXraw(argv[1]);
  if (!xr) {
    fprintf(stderr, "Failed to open or read '%s'\n", argv[1]);
    return -1;
  }

  showXraw(xr);
  
  if (xr->data) {
    free(xr->data);
  }
  free(xr);

  return 0;
}
