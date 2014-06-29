/* xrw2stl.c
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
 * $Id: xrw2stl.c 78 2012-11-26 22:41:06Z barnesd $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "s2plot.h"
#include "libxrw.c"
#include "libobj.c"

extern _S2ISOSURFACE *_s2_isosurfs;

void usage(char *exename) {
  fprintf(stderr, "usage: %s [options] -f xrwfilename -o stlfilename\n", exename);
  fprintf(stderr, "-s s1 s2 s3\t\t set data averaging along each axis (1 1 1)\n");
  fprintf(stderr, "-l val\t\t set surface value (default 0.5)\n");
}

int main(int argc, char *argv[]) {

  #define FNAMELEN 400
  char ifname[FNAMELEN+1], ofname[FNAMELEN+1];
  int haveifname = 0, haveofname = 0;

  int invert = 0;
  int stride[3] = {1,1,1}; //{3,3,1};
  float surf_level = 0.5;

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
    } else if (!strcmp(argv[ic], "-l")) {
      surf_level = atof(argv[++ic]);
    } else if (!strcmp(argv[ic], "-s")) {
      stride[0] = atoi(argv[++ic]);
      stride[1] = atoi(argv[++ic]);
      stride[2] = atoi(argv[++ic]);
      if (stride[0] < 1 || stride[1] < 1 || stride[2] < 1) {
	fprintf(stderr, "stride must be at least 1 in each direction!\n");
	exit(-1);
      }    }
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
  
  VOL_STRUCT *vol = Xraw2Xvol(xr, stride);
  if (!vol) {
    fprintf(stderr, "Failed to parse data volume.\n");
    return -1;
  }

  showXvol(vol);

  int nx, ny, nz;
  nx = vol->nx;
  ny = vol->ny;
  nz = vol->nz;

  // while we don't actually draw any S2PLOT graphics, we do need S2PLOT
  // to be initialised so we can create an isosurface to save.
  s2opend("/S2MONO",argc,argv);			/* Open the display */
  
  if (invert) {
    ss2sbc(1., 1.,1.);
    ss2sfc(0., 0., 0.);
  }

  // viewport: simply generate appropriate aspect ratio
  float sx = 1.0;
  float sy = vol->wdy*(float)ny / (vol->wdx*(float)nx);
  float sz = vol->wdz*(float)nz / (vol->wdx*(float)nx);
  s2svp(-sx,sx, -sy,sy, -sz,sz);

  // world coordinate system: put centre of volume at (0,0,0)
  // world coordinates are physical units
  float ddx = vol->wdx*(float)nx * 0.5;
  float ddy = vol->wdy*(float)ny * 0.5;
  float ddz = vol->wdz*(float)nz * 0.5;
  s2swin(-ddx,ddx, -ddy,ddy, -ddz,ddz);
  fprintf(stderr, "s2swin(%f,%f,%f,%f,%f,%f)\n", -ddx,ddx,-ddy,ddy,-ddz,ddz);

  int i;
  float tr[12];
  for (i=0;i<12;i++) {				/* Set-up transfrom matrix */
    tr[i] = 0.0;
  }
  tr[0] = 0; // -vol->wdx*(float)nx*0.5;
  tr[1] = vol->wdx;
  tr[4] = 0; //-vol->wdy*(float)ny*0.5;
  tr[6] = vol->wdy;
  tr[8] = 0; //-vol->wdz*(float)nz*0.5;
  tr[11]= vol->wdz;
  
  char trans = 'o';

  int sidx = ns2cis(vol->data, nx, ny, nz, 0, nx-1, 0, ny-1, 0, nz-1, 
		tr, surf_level, 1, trans, 1.0, 1.0, 1., 1.);
  
  ns2dis(sidx, 1);

  FILE *fout = fopen(ofname, "wb");
  
  char title[80];
  sprintf(title, "%% STL file created by xrw2stl\n");
  fwrite(title, 1, 80, fout);
  
  _S2ISOSURFACE *surf = _s2_isosurfs + sidx;  

  int ntr = surf->ntri;
  fwrite(&ntr, sizeof(ntr), 1, fout);

  STL_FACET_STRUCT facet;
  unsigned short uu = 0;
  int vx;
  for (i = 0; i < surf->ntri; i++) {
    XYZ *iP = surf->trivert + 3 * i;
    for (vx = 0; vx < 3; vx++) {
      facet.vx[vx*3+0] = iP[vx].x;
      facet.vx[vx*3+1] = iP[vx].y;
      facet.vx[vx*3+2] = iP[vx].z;
    }
    XYZ N = CalcNormal(iP[0], iP[1], iP[2]);
    facet.n[0] = N.x;
    facet.n[1] = N.y;
    facet.n[2] = N.z;
    
    fwrite(&facet, sizeof(facet), 1, fout);

    fwrite(&uu, sizeof(unsigned short), 1, fout);

  }
  
  fclose(fout);

  return 0;
}
  
