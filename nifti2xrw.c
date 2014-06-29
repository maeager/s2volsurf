/* nifti2xrw.c
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
 * $Id: nifti2xrw.c 224 2014-04-23 22:54:03Z barnesd $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <nifti1_io.h> 
#include "libxrw.c"

void usage(char *exename) {
    fprintf(stderr, "usage: %s [options] -f niftifilename -o xrwfilename\n", exename);
    fprintf(stderr, "options:\n");
    fprintf(stderr, "\t -N\t\t do not normalise\n");
    fprintf(stderr, "\t -R min max\t\t normalise to range min-max\n");
}

int main(int argc, char *argv[]) {

  #define FNAMELEN 400
  char ifname[FNAMELEN+1], ofname[FNAMELEN+1];
  int haveifname = 0, haveofname = 0;
  int normalise = 1;
  int range_normalise = 0;
  float range_min = 0, range_max = 0;

  if (argc < 5) {
    usage(argv[0]);
    exit(-1);
  }

  int ic = 1;
  while (ic < argc) {
    if (!strcmp(argv[ic], "-f")) {
      strncpy(ifname, argv[++ic], FNAMELEN);
      haveifname = 1;
    } else if (!strcmp(argv[ic], "-N")) {
      normalise = 0;
    } else if (!strcmp(argv[ic], "-o")) {
      strncpy(ofname, argv[++ic], FNAMELEN);
      haveofname = 1;
    } else if (!strcmp(argv[ic], "-R")) {
      range_normalise = 1;
      sscanf(argv[++ic], "%f", &range_min);
      sscanf(argv[++ic], "%f", &range_max);
    }
    ic++;
  }
  
  if (!haveifname) {
    usage(argv[0]);
    return -1;
  }

  nifti_image *nim = nifti_image_read(ifname, 1);
  if (!nim) {
    fprintf(stderr, "Couldn't open nifti image (fa.nii)\n");
    exit(-1);
  }
  nifti_image_infodump(nim);

  if (nim->ndim < 3) {
    fprintf(stderr, "nifti image has less than 3 dimensions\n");
    exit(-1);
  }

  // oh what a lovely hack:
  if ((nim->nz == 1) && (nim->nt > 1)) {
    nim->nz = nim->nt;
    nim->nt = 1;
  }

  if (nim->nx < 16 || nim->ny < 16 || nim->nz < 16) {
    fprintf(stderr, "nifti image dimensions a bit smallish\n");
    exit(-1);
  }
  if (nim->datatype != 2 && nim->datatype != 8 && nim->datatype != 16 && nim->datatype != 4 && nim->datatype != 512) {
    fprintf(stderr, "unsupported data type\n");
    exit(-1);
  }
  if (!nim->data) {
    fprintf(stderr, "data is null!");
    exit(-1);
  }
  VOL_STRUCT niftiv;
  strcpy(niftiv.filename, "nnnn.nnn");

  // input crop
  int blc[3], trc[3];
  blc[0] = blc[1] = blc[2] = 0;
  trc[0] = nim->nx - 1;
  trc[1] = nim->ny - 1;
  trc[2] = nim->nz - 1;

  int stride[3] = {1,1,1};

  int w, h;
  w = nim->nx;
  h = nim->ny;

  int i;

  // fixed for input crop
  niftiv.nx = (trc[0] - blc[0] + 1) / stride[0];
  niftiv.ny = (trc[1] - blc[1] + 1) / stride[1];
  niftiv.nz = (trc[2] - blc[2] + 1) / stride[2];

  niftiv.wdx = nim->dx * stride[0];
  niftiv.wdy = nim->dy * stride[1];
  niftiv.wdz = nim->dz * stride[2];
  for (i = 0; i < 256; i++) {
    niftiv.red[i] = niftiv.green[i] = niftiv.blue[i] = powf((float)i / (float)255, 0.0);
  }
  // copy data with stride
  niftiv.data = (float ***)malloc(niftiv.nx * sizeof(float **));
  if (!niftiv.data) {
    fprintf(stderr,"Failed to allocate float volume.\n");
    //return NULL;
    exit(-1);
  }
  int j, k, ii, ij, ik, ax, ay, az;
  float sum;

  // fill y and z planes in reverse directions so that home view in 
  // s2plot (and pdf) matches orientation xrw was exported from in 
  // OsiriX - use jtarg, ktarg to control
  int jtarg, ktarg;

  for (i = 0; i < niftiv.nx; i++) {
    ii = blc[0] + i * stride[0];
    // copy this frame in to this slice
    niftiv.data[i] = (float **)malloc(niftiv.ny * sizeof(float *));
    if (!niftiv.data[i]) {
      fprintf(stderr, "Failed to allocate row in float volume.\n");
      exit(-1);
    }
    for (j = 0; j < niftiv.ny; j++) {
      ij = blc[1] + j * stride[1];

      //jtarg = niftiv.ny-1-j;
      jtarg = j;

      niftiv.data[i][jtarg] = (float *)malloc(niftiv.nz * sizeof(float));
      if (!niftiv.data[i][jtarg]) {
        fprintf(stderr, "Failed to allocate column in float volume.\n");
	exit(-1);
      }
      for (k = 0; k < niftiv.nz; k++) {
        ik = blc[2] + k * stride[2];

	//ktarg = niftiv.nz-1-k;
	ktarg = k;

        // averaging
        sum = 0.0;
	short *deref_s;
	int *deref_i;
	unsigned short *deref_us;
	unsigned char *deref_uc;
        float *deref;
	for (ax = 0; ax < stride[0]; ax++) {
          for (ay = 0; ay < stride[1]; ay++) {
            for (az = 0; az < stride[2]; az++) {
	      switch (nim->datatype) {
	      case 2: // unsigned char
		deref_uc = (unsigned char *)nim->data;
		sum += (float)(deref_uc[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]) / 255.0;
		break;
	      case 4: // signed short 
		deref_s = (short *)nim->data;
		sum += (float)(deref_s[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]);
		break;
	      case 8: // integer
		deref_i = (int *)nim->data;
		sum += (float)(deref_i[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]);
		break;
	      case 512: // unsigned short
		deref_us = (unsigned short *)nim->data;
		sum += (float)(deref_us[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]);
		break;
	      case 16: // float
		deref = nim->data;
		sum += (float)(deref[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]); /* / 255.0; */
		break;
		
	      }
            }
          }
        }

        niftiv.data[i][jtarg][ktarg] = sum / (float)(stride[0] * stride[1] * stride[2]);

      }
    }
  }

  // normaliseXraw into [0,1]
  if (range_normalise) {
    fprintf(stderr, "RANGE Normalising...\n");
    rangeNormaliseXvol(&niftiv, range_min, range_max);
  } else if (normalise) {
    fprintf(stderr, "Normalising ...\n");
    normaliseXvol(&niftiv);
  }

  XRAW_STRUCT *xraw = Xvol2Xraw(&niftiv);
  if (haveofname) {
    sprintf(xraw->filename, ofname);
  } else {
    sprintf(xraw->filename, "%s.xrw", ifname);
  }
  showXraw(xraw);
  saveXraw(xraw);

  char matfn[256];
  sprintf(matfn, "%s.qmat", xraw->filename);
  FILE *QMATO = fopen(matfn, "w");
  sprintf(matfn, "%s.smat", xraw->filename);
  FILE *SMATO = fopen(matfn, "w");
  int mi, mj;
  for (mi = 0; mi < 4; mi++) {
    for (mj = 0; mj < 4; mj++) {
      fprintf(QMATO, "%f ", nim->qto_xyz.m[mi][mj]);
      fprintf(SMATO, "%f ", nim->sto_xyz.m[mi][mj]);
    }
    fprintf(QMATO, "\n");
    fprintf(SMATO, "\n");
  }
  fclose(QMATO);
  fclose(SMATO);


  exit(0);

   
  return 1;
}

