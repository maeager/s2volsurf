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


/** objcentre.c Michael Eager 2016
 * Centre all listed objects to the origin
 * {0,0,0}. All obj files are parsed for the global bounding box.
 * The original files are moved to <file>-orig.obj and the new file
 * is created with all vertices translated by the new mid point.
 * usage:  objcentre {list of obj files} 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>        /* errno */

#include "s2plot.h"
#include "libobj.c"

#define CMD_SIZE 256
char cmd[CMD_SIZE];
void SaveObj ( char *filename_old,  char *filename_new, XYZ newCentre );


int main(int argc, char **argv) {
  
  if (argc < 2) {
    fprintf(stderr, "usage: %s objfilename[s]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  int numObjFiles = argc-1;
  XYZ minP, maxP;

  int iarg;
  for (iarg = 1; iarg < argc; iarg++) {
    OBJ_STRUCT *obj; // = loadObj(argv[argc -iarg], "null", 1, 1, 1, 1);
    char *ext = argv[argc-iarg] + strlen(argv[argc-iarg]) - 3;
    if (!strncmp(ext, "obj", 3)) {
      obj = loadObj(argv[argc-iarg], "null", 1,1,1,1);
    } else if (!strncmp(ext, "stl", 3)) {
      obj = loadObjFromSTL(argv[argc-iarg], "null", 1,1,1,1);
    } else {
      fprintf(stderr, "Cannot load file '%s'.\n", argv[argc-iarg]);fflush(stderr);
      exit(EXIT_FAILURE);
    }

    if (iarg == 1) {
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
	  maxP.x, maxP.y, maxP.z);fflush(stderr);
  fprintf(stderr, "OBJ vertex midpoint: (%f,%f,%f)\n", mid.x, mid.y, mid.z);fflush(stderr);
  char *origfilename = NULL; //(char*)malloc(100*sizeof(char));
  char *centredfilename = NULL;

  fprintf(stderr, "Translating obj centre to 0,0,0....");  fflush(stderr);
  for (iarg = 1; iarg < argc; iarg++) {
    
    char *ext = argv[argc-iarg] + strlen(argv[argc-iarg]) - 3;
    origfilename = (char*)realloc(origfilename, (strlen(argv[argc-iarg]) + 5)*sizeof(char));
    centredfilename = argv[argc-iarg];

    strlcpy(origfilename,argv[argc-iarg], strlen(argv[argc-iarg]) - 3);

    if (!strncmp(ext, "obj", 3)) {
      strlcpy(&origfilename[strlen(argv[argc-iarg]) - 4],"-orig.obj", 10);
      fprintf(stderr,"Moving '%s' to '%s'.\n",argv[argc-iarg], origfilename);fflush(stderr);
      sprintf( cmd, "mv %s %s", argv[argc-iarg], origfilename );
      int ret = system(cmd);
      if (ret == -1) fprintf(stderr,"System command failed.");
      else fprintf(stderr,"System command returned: %d\n",ret);
      fprintf(stderr,"Translating centre to 0,0,0\n");fflush(stderr);
      int res = SaveTranslatedObj(origfilename, centredfilename, mid);
      if (res == -1) exit(EXIT_FAILURE);
    } else {
      fprintf(stderr, "Cannot load file '%s'.\n", argv[argc-iarg]);fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
  }
  fprintf(stderr, "Old midpoint: (%f,%f,%f)\n", mid.x, mid.y, mid.z);
  fprintf(stderr, "New vertex range: (%f,%f,%f) -> (%f,%f,%f)\n",
          minP.x-mid.x, minP.y-mid.y, minP.z-mid.z, maxP.x-mid.x,
          maxP.y-mid.y, maxP.z-mid.z);fflush(stderr);

  fprintf(stdout, "BoundingBox %.6f %.6f %.6f %.6f %.6f %.6f \n", minP.x-mid.x,
          minP.y-mid.y, minP.z-mid.z, maxP.x-mid.x, maxP.y-mid.y,
         maxP.z-mid.z);
  free(origfilename);
  exit(EXIT_SUCCESS);
}




int SaveTranslatedObj(char *filename_old,  char *filename_new, XYZ newCentre){
  
  int lcount = 0,nverts=0;
  char line[401];
  //char cmd;
  float fx, fy, fz;

  FILE *fin = fopen(filename_old, "r");
  if (fin == NULL){
    fprintf(stderr,"SaveTranslatedObj: Original file %s could not be opened.\n",filename_old);
    return -1;
  }
  FILE *fout = fopen(filename_new, "w");
  if (fout == NULL){
    fprintf(stderr,"SaveTranslatedObj: Cannot write to file %s.\n", filename_new);
    return -1;
  }
  
  while (fgets(line, 400, fin)) {
    lcount++;
    char *result = NULL;
    char space[] = " ";
    //fprintf(stderr,"%s\n",line);fflush(stderr);
    //fprintf(stderr, "token = %s\n", result);
    if (!strncmp(line, "v ", 2)) {
      result = strtok(line, space);
      nverts++;
      
      fy = atof(strtok(NULL, space));
      fz = atof(strtok(NULL, space));
      fx = atof(strtok(NULL, space));
      fprintf(fout,"v %.6f %.6f %.6f 1.0\n",
              fy - newCentre.y,
              fz - newCentre.z,
              fx - newCentre.x);
    } else{
      fprintf(fout,"%s",line); //line contains trailing \0 so no need for \n
    }
    fflush(fout);
  }
  fclose (fout);
  fclose (fin);

  fprintf(stderr, "Read %d lines from '%s' to '%s'.\n", lcount, filename_old, filename_new);
  fprintf(stderr, "Translated\t%d vertices to file \n", nverts);
  fflush(stderr);
  
  return 0;
}
