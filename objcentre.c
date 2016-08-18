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

#include <stdio.h>        /* fopen,fclose,fflush */
#include <stdlib.h>
#include <math.h>
#include <string.h>       /* strlcpy   */ 
#include <errno.h>        /* errno */

#include "s2plot.h"
#include "libobj.c"

#define CMD_SIZE 1024
char cmd[CMD_SIZE];

int QUIET=0;

int SaveTranslatedObj(char *filename_old,  char *filename_new, XYZ newCentre){
  
  int lcount = 0,nverts=0;
  char line[401];
  //char cmd;
  float fx, fy, fz;

  FILE *fin = fopen(filename_old, "r");
  if (fin == NULL){
    fprintf(stderr,"SaveTranslatedObj: Original file %s could not be opened.\n",filename_old);fflush(stderr);
    return -1;
  }
  FILE *fout = fopen(filename_new, "w");
  if (fout == NULL){
    fprintf(stderr,"SaveTranslatedObj: Cannot write to file %s.\n", filename_new);fflush(stderr);
    return -1;
  }
  char *result = NULL;
  while (fgets(line, 400, fin)) {
    lcount++;
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

  if(QUIET==0){
    fprintf(stderr, "Read %d lines from '%s' to '%s'.\n", lcount, filename_old, filename_new);
    fprintf(stderr, "Translated\t%d vertices to file \n", nverts);
}
  fflush(stderr);
  
  return 0;
}


int main(int argc, char **argv) {

  int numObjFiles = argc-1;
  if(strcmp(argv[1],"-q") == 0) {
    QUIET=1;
    numObjFiles-=1;
  }
  if (argc < 2) {
    fprintf(stderr, "usage: %s objfilename[s]\n", argv[0]);
    exit(EXIT_FAILURE);
  }


  XYZ minP, maxP;
  
  int iarg;
  for (iarg = 1; iarg < numObjFiles; iarg++) {
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
    if(QUIET==0)
      fprintf(stderr, "Found %d Obj files.\n",numObjFiles);
    fflush(stderr);

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

  
  if ( fabsf((float)mid.x) < 0.00001f &&
       fabsf((float)mid.y) < 0.00001f &&
       fabsf((float)mid.z) < 0.00001f) {
    fprintf(stderr, "OBJ vertex midpoint already 0,0,0.\n");fflush(stderr);
    fprintf(stdout, "No Changes.\nBoundingBox: %.6f %.6f %.6f %.6f %.6f %.6f \n",
            minP.x, minP.y, minP.z,
            maxP.x, maxP.y, maxP.z);
    exit(EXIT_SUCCESS);
  }
  

       char *origfilename= (char*)malloc(512*sizeof(char));
       char *centredfilename = NULL; //= (char*)malloc(512*sizeof(char));
  int res=0,ret=0;
  fprintf(stderr, "Translating obj centre to origin....\n");  fflush(stderr);
  for (iarg = 1; iarg < numObjFiles; iarg++) {
    centredfilename = argv[argc-iarg];    
    char *ext = centredfilename + strlen(centredfilename) - 3;
    origfilename = (char*)realloc(origfilename, (strlen(centredfilename) + 5)*sizeof(char));
    strlcpy(origfilename,centredfilename, strlen(centredfilename) - 3);

    //    centredfilename = (char*)realloc(centredfilename, (strlen(argv[argc-iarg]))*sizeof(char));
    //    strlcpy(centredfilename,argv[argc-iarg], strlen(argv[argc-iarg])+1);
    

    if (!strncmp(ext, "obj", 3)) {
      strlcpy(&origfilename[strlen(centredfilename) - 4],"-orig.obj", 10);

      if(QUIET==0) {
        fprintf(stderr,"Moving '%s' to '%s'.\n",centredfilename, origfilename);
      }
      //Move file and force sync
      sprintf( cmd, "mv -f \"%s\" \"%s\" ", centredfilename, origfilename );

      if(QUIET==0) {
        fprintf(stderr,"Cmd: %s\n",cmd);
        fflush(stderr);
      }
      ret = system(cmd);
      if (ret == -1) fprintf(stderr,"System command failed.");
      else if (ret == 127) fprintf(stderr,"Execution of the shell failed.");
      //else fprintf(stderr,"System command returned: %d\n",ret);

      if(QUIET==0) {
        fprintf(stderr,"\nTranslating centre to 0,0,0\n");fflush(stderr);
      }
      res = SaveTranslatedObj(origfilename, centredfilename, mid);
      if (res == -1) exit(EXIT_FAILURE);
    } else {
      fprintf(stderr, "Cannot load file '%s'.\n", centredfilename); fflush(stderr);
      exit(EXIT_FAILURE);
    }
  }
  fprintf(stderr, "Old midpoint: (%f,%f,%f)\n", mid.x, mid.y, mid.z);
  fprintf(stderr, "New vertex range: (%f,%f,%f) -> (%f,%f,%f)\n",
          minP.x-mid.x, minP.y-mid.y, minP.z-mid.z, maxP.x-mid.x,
          maxP.y-mid.y, maxP.z-mid.z);fflush(stderr);
  fprintf(stdout, "Format: XMIN YMIN ZMIN XMAX YMAX ZMAX \n"
  fprintf(stdout, "BoundingBox: %.6f %.6f %.6f %.6f %.6f %.6f \n", minP.x-mid.x,
          minP.y-mid.y, minP.z-mid.z, maxP.x-mid.x, maxP.y-mid.y,
         maxP.z-mid.z);
  if (origfilename) free(origfilename);
  // if (centredfilename) free(centredfilename);
  exit(EXIT_SUCCESS);
}

