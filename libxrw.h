/* libxrw.h
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
 * $Id: libxrw.c 228 2014-05-11 10:08:11Z barnesd $
 *
 */

/* * * * * N.B. read/write only ints, not longs (32/64-bit diff) * * * * */

#if !defined(_LIBXRW_H)
#define _LIBXRW_H

#include <stdio.h>
#include <string.h>

#include "s2types.h"

typedef struct {
  char filename[400];       // filename
  int nx, ny, nz;           // data dimensions
  float wdx, wdy, wdz;      // voxel dimensions
  unsigned char *data;      // voxels
  unsigned char red[256];   // red component of LUT
  unsigned char green[256]; // green
   unsigned char blue[256];  // blue
} XRAW_STRUCT;
XRAW_STRUCT *loadXraw(char *fname);
XRAW_STRUCT *createXraw(char *fname, int nx, int ny, int nz, int borderwidth);
XRAW_STRUCT *preloadXraw(char *fname); // just metadata
XRAW_STRUCT *trimXraw(XRAW_STRUCT *xr, int trim[3]);
void showXraw(XRAW_STRUCT *xr);
int saveXraw(XRAW_STRUCT *xr);
void deleteXraw(XRAW_STRUCT *xr);

typedef struct {
  char filename[400];  // filename
  int nx, ny, nz;      // data dimensions  
  float wdx, wdy, wdz; // voxel sizes
  float ***data;       // floating point data [0,1]
  float red[256];      // red components [0,1]
  float green[256];    // green [0,1]
  float blue[256];     // blue [0,1]
  float bzero, bscale; // raw = bzero + normalised_data * bscale
} VOL_STRUCT;
VOL_STRUCT *Xraw2Xvol(XRAW_STRUCT *xrs, int stride[3]);
VOL_STRUCT *makePow2Xvol(VOL_STRUCT *xrv); 
void showXvol(VOL_STRUCT *xv);
void deleteXvol(VOL_STRUCT *xv);

VOL_STRUCT *loadUshortRaw(char *ifname, int nx, int ny, int nz);
VOL_STRUCT *loadTGAstack(char *basename, int startframe, int endframe, int stride);
VOL_STRUCT *rebinXvol(VOL_STRUCT *ivol, int stride[3]);
void normaliseXvol(VOL_STRUCT *vol); // normalise in measured min, max
void rangeNormaliseXvol(VOL_STRUCT *vol, float min, float max); // normalise to provided min,max
void tightenXvol(VOL_STRUCT *vol, float losigma, float hisigma);
void derivXvol(VOL_STRUCT *vol);
VOL_STRUCT *cloneXvol(VOL_STRUCT *vol);
void accumulateXvol(VOL_STRUCT *dest, VOL_STRUCT *src, float scale);
XRAW_STRUCT *Xvol2Xraw(VOL_STRUCT *xv);
void Xvol2png(VOL_STRUCT *xv, char *filename, int use_rows);
VOL_STRUCT *loadRaw(char *filename, int nx, int ny, int nz, int type);

#endif

