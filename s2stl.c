/* s2stl.c
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
 * $Id: s2stl.c 78 2012-11-26 22:41:06Z barnesd $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "s2plot.h"
#include "libobj.c"

typedef struct {
  char *filename;
  char *label;
  float r, g, b, a;
} STL_SURF;

// LIST YOUR SURFACE/S HERE: {filename, textlabel, r, g, b, alpha}
// knot.stl is from: http://www.eng.nus.edu.sg
// 	curl -o knot.zip http://www.eng.nus.edu.sg/LCEL/RP/u21/download/STLFiles/knot.zip
//	unzip knot.zip
#define NSURF 1
STL_SURF surf[] = {
  {"knot.stl", "Knot", 1.0,0.0,1.0, 1.0}
};

OBJ_STRUCT *obj[NSURF];

// callback function to set lighting and draw geometry
void cb(int *kc, double *t) {

  // arrange for "headlamp" - light that is just above the camera
  XYZ campos;
  XYZ up;
  ss2qc(&campos, &up, NULL, 1);
  COLOUR amb = {0.0,0.0,0.0};
  COLOUR light = {0.5, 0.5, 0.5};
  XYZ newpos = VectorAdd(campos, up);
  ss2sl(amb, 1, &newpos, &light, 1);

  // draw surfaces
  int i;
  for (i = 0; i < NSURF; i++) {
    drawObj(obj[i]);
  }
}

// direct access to some internal S2PLOT settings in lieu of functions
// to set these globals...
extern float specularcolour[4];
extern float shininess[1];
extern float emission[4];

int main(int argc, char **argv) {
  
  specularcolour[0] = specularcolour[1] = specularcolour[2] = 0.4;
  shininess[0] = 0.6;

  s2opend("/?", argc, argv);
  
  ss2sbc(1., 1., 1.);
  ss2sfc(0., 0., 0.);

  XYZ minP, maxP;

  int i;
  for (i = 0; i < NSURF; i++) {
    obj[i] = loadObjFromSTL(surf[i].filename, 
			    surf[i].label,
			    surf[i].r, surf[i].g, surf[i].b, surf[i].a);
    if (!obj[i]) {
      fprintf(stderr, "Failed to open STL file '%s'.\n", surf[i].filename);
      return(-1);
    }
    if (i == 0) {
      minP = obj[i]->minP;
      maxP = obj[i]->maxP;
    } else {
      if (obj[i]->minP.x < minP.x) {
	minP.x = obj[i]->minP.x;
      }
      if (obj[i]->minP.y < minP.y) {
	minP.y = obj[i]->minP.y;
      }
      if (obj[i]->minP.z < minP.z) {
	minP.z = obj[i]->minP.z;
      }

      if (obj[i]->maxP.x > maxP.x) {
	maxP.x = obj[i]->maxP.x;
      }
      if (obj[i]->maxP.y > maxP.y) {
	maxP.y = obj[i]->maxP.y;
      }
      if (obj[i]->maxP.z > maxP.z) {
	maxP.z = obj[i]->maxP.z;
      }
    }
  }

  // offset all model parts
  XYZ mid = VectorAdd(minP, maxP);
  mid = VectorMul(mid, 0.5);
  for (i = 0; i < NSURF; i++) {
    translateObj(obj[i], mid);
  }
  minP = VectorSub(mid, minP);
  maxP = VectorSub(mid, maxP);

  // now find global min (x,y,z) and max (x,y,z) so we have a
  // properly scaled world coordinate box to draw in
  float min = minP.x;
  if (minP.y < min) {
    min = minP.y;
  }
  if (minP.z < min) {
    min = minP.z;
  }

  float max = maxP.x;
  if (maxP.y > max) {
    max = maxP.y;
  }
  if (maxP.z > max) {
    max = maxP.z;
  }

  s2swin(min,max, min,max, min,max);
  
  cs2scb(cb);
  ss2spt(1);
  ss2srm(SHADE_SPECULAR);

  s2show(1);

  return 0;
}
