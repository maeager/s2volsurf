/* objrange.c
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
 * $Id: objrange.c 86 2013-01-10 00:20:44Z barnesd $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "s2plot.h"
#include "libobj.c"

int main(int argc, char **argv) {
  
  if (argc < 2) {
    fprintf(stderr, "usage: %s objfilename[s]\n", argv[0]);
    return -1;
  }

  XYZ minP, maxP;

  int i;
  for (i = 1; i < argc; i++) {
    OBJ_STRUCT *obj; // = loadObj(argv[argc -i], "null", 1, 1, 1, 1);
    char *ext = argv[argc-i] + strlen(argv[argc-i]) - 3;
    if (!strncmp(ext, "obj", 3)) {
      obj = loadObj(argv[argc-i], "null", 1,1,1,1);
    } else if (!strncmp(ext, "stl", 3)) {
      obj = loadObjFromSTL(argv[argc-i], "null", 1,1,1,1);
    } else {
      fprintf(stderr, "Cannot load file '%s'.\n", argv[argc-i]);
      exit(-1);
    }
    if (i == 1) {
      minP = obj->minP;
      maxP = obj->maxP;
    } else {
      if (obj->minP.x < minP.x) {
	minP.x = obj->minP.x;
      }
      if (obj->minP.y < minP.y) {
	minP.y = obj->minP.y;
      }
      if (obj->minP.z < minP.z) {
	minP.z = obj->minP.z;
      }

      if (obj->maxP.x > maxP.x) {
	maxP.x = obj->maxP.x;
      }
      if (obj->maxP.y > maxP.y) {
	maxP.y = obj->maxP.y;
      }
      if (obj->maxP.z > maxP.z) {
	maxP.z = obj->maxP.z;
      }
      
    }
  }

  // offset all model parts
  XYZ mid = VectorAdd(minP, maxP);
  mid = VectorMul(mid, 0.5);
  fprintf(stderr, "OBJ vertex range: (%f,%f,%f) -> (%f,%f,%f)\n", minP.x, minP.y, minP.z,
	  maxP.x, maxP.y, maxP.z);
  fprintf(stderr, "OBJ vertex midpoint: (%f,%f,%f)\n", mid.x, mid.y, mid.z);

  return 0;
}
