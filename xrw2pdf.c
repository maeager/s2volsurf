/* xrw2pdf.c
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
 * $Id: xrw2pdf.c 224 2014-04-23 22:54:03Z barnesd $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "s2plot.h"
#include "libxrw.c"

#define MAXNVOL 4

#include "libobj.c"
void ReadViewFile(char *name);
void writePRC(void);

// surfaces / structures
OBJ_STRUCT **objs = NULL;
int nobjs = 0;

// tracks
#define TRACK_CHUNK 100
#define VERT_CHUNK 1000
int *count = NULL;    // how long is each track?
XYZ **start = NULL;    // where does each track start?
long ntracks = 0;
extern int _s2_conelines;


int NVR = 0;

void ds2dvrx(int vrid, int force_reload, int axis, int force_draw);
int vidx[4], vidy[4], vidz[4];
int doPRC;
  int suppress_vr = 0;
float isolev[4];
COLOUR isocol[4];
float isoalf[4];

int do_slice_reveal = 0;
// count of numeric keys pressed
int numcount[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
void numcb(int *num) {
  if ((*num >= 0) && (*num <= 9)) {
    numcount[*num]++;
  }

  if (*num == 9) {
    do_slice_reveal = 1;
  } else if (*num == 0) {
    do_slice_reveal = -1;
  }
}

// standards are dmin 0 dmax 1 amin 0 amax 0.2
float dmin[4], dmax[4], amin[4], amax[4];
float scapower[MAXNVOL];

int doDeriv[MAXNVOL];

float dscale[MAXNVOL][3];

int do_3d_render = 0;
int tex3d[MAXNVOL];
int _3d_speed_up = 1;

//float glow_period_frames = 1.0 * 25.0; // if >0 this is the period of glow oscillation for the VR
float glow_period_frames = -1.0;

int _nx[MAXNVOL], _ny[MAXNVOL], _nz[MAXNVOL];
float _tr[MAXNVOL][12];

#define FNAMELEN 400
int viewfile = 0;
char viewfile_fname[FNAMELEN];
int invert = 0;



void cb(double *t, int *kc) {

  //float _r,_g,_b;
  //ss2qbc(&_r, &_g, &_b);
  //fprintf(stderr, "bg is %f %f %f\n", _r, _g, _b);



  static int beenhere = 0;
  static int framecount;
  if (!beenhere) {
    if (viewfile) {
      fprintf(stderr, "reading view file %s\n", viewfile_fname);
      ReadViewFile(viewfile_fname);
      // reading view file resets the bg color!  So reset it here!
      if (invert) {
	ss2sbc(1., 1.,1.);
	ss2sfc(0., 0., 0.);
      } else {
	fprintf(stderr, "normal white-on-black\n");
	ss2sbc(0., 0., 0.);
	ss2sfc(1., 1., 1.);
      }
    }
    framecount = -1;
    beenhere = 1;
  }
  framecount++;


  int wh; // loop over the volume(s) to draw
  for (wh = 0; wh < NVR; wh++) {
    
    s2scir(1000+wh*256, 1000+wh*256+255);
    
    if (isolev[wh] >= 0.0) {
      char thisgrp[255];
      sprintf(thisgrp, "%s-%.3f", "isosurface", isolev[wh]);
      pushVRMLname(thisgrp);
      ns2dis(vidx[wh], 0);
      pushVRMLname("ANON");
    } else 
    
    if (do_3d_render) {
      
     if (wh == 0) {

      // to get this right we need to interleave the draws from each volume
      // (also would need to do this for orthogonal 2d texture rendering but
      // that's impossible because the planes for two different volumes are 
      // not nec. parallel!)

      // so this 3d render draw is only done once in the big "wh" loop.
      // Let's calculate some bounds etc., then for each slice of volume
      // 0 (wh == 0) we will draw interleaved slices from the other volumes.

      // these are the vertices of a unit cube
      XYZ unitverts[] = {{0, 0, 0},
			 {1, 0, 0},
			 {0, 1, 0},
			 {1, 1, 0},
			 {0, 0, 1},
			 {1, 0, 1},
			 {0, 1, 1},
			 {1, 1, 1}};
      // these are the edges of the cube, made up by joining these unitverts:
      int edges[12][2] = {{0,1}, {0,2}, {1,3}, {2,3},
			  {4,5}, {4,6}, {5,7}, {6,7},
			  {0,4}, {1,5}, {2,6}, {3,7}};
      
      //unsigned int texid = 1;
      char itrans = 's';
      float ialpha = 0.4;
    
      // the world vertices of the data being displayed (= dataverts * tr)
      XYZ worldverts[MAXNVOL][8];

      // 1. get camera position and view direction in world coords
      XYZ campos, upvec, viewdir, right;
      ss2qc(&campos, &upvec, &viewdir, 1);
      Normalise(&viewdir);
      right = CrossProduct(viewdir, upvec);
      
      int zh, i;
      // loop in reverse order so that final values in on-going vars are for wh==0
      for (zh = NVR-1; zh >= 0; zh--) {
	// these would be args to the vol renderer once it's built into 
	// S2PLOT 
	int adim = _nx[zh];
	int bdim = _ny[zh];
	int cdim = _nz[zh];
	int a1 = 0;
	int a2 = adim-1;
	int b1 = 0;
	int b2 = bdim-1;
	int c1 = 0;
	int c2 = cdim-1;
      
	for (i = 0; i < 8; i++) {
	  // the vertices of the data being displayed
	  XYZ dataverts[8];
	  // calculate this data vertex
	  dataverts[i].x = a1 + unitverts[i].x * (a2-a1);
	  dataverts[i].y = b1 + unitverts[i].y * (b2-b1);
	  dataverts[i].z = c1 + unitverts[i].z * (c2-c1);
	  
	  // and this world vertex position
	  worldverts[zh][i].x = _tr[zh][0] + _tr[zh][1] * dataverts[i].x 
	    + _tr[zh][2] * dataverts[i].y + _tr[zh][3] * dataverts[i].z;
	  worldverts[zh][i].y = _tr[zh][4] + _tr[zh][5] * dataverts[i].x 
	    + _tr[zh][6] * dataverts[i].y + _tr[zh][7] * dataverts[i].z;
	  worldverts[zh][i].z = _tr[zh][8] + _tr[zh][9] * dataverts[i].x 
	    + _tr[zh][10] * dataverts[i].y + _tr[zh][11] * dataverts[i].z;
	  //- ns2vthpoint(worldverts[i], GREEN, 5.0);
	}
      }

      // 2. find indices of first and last vertices: first is that vertex
      //    which a plane normal to the viewdir crosses, travelling towards
      //    the centre of the cube, from the camera position.
      int near_vtx, far_vtx;
      float near_dist, far_dist;
      XYZ tmp;
      float thisdist;  
      near_vtx = far_vtx = -1;
      near_dist = 9e30;
      far_dist = -9e30;
      for (i = 0; i < 8; i++) {
	// and now its distance from the camera position measured along the
	// view direction
	tmp = VectorSub(campos, worldverts[wh][i]);
	thisdist = DotProduct(tmp, viewdir);
	if (thisdist < near_dist) {
	  near_vtx = i;
	  near_dist = thisdist;
	}
	if (thisdist > far_dist) {
	  far_vtx = i;
	  far_dist = thisdist;
	}
      }
      //- ns2vthpoint(worldverts[near_vtx], BLUE, 10.0);
      //- ns2vthpoint(worldverts[ far_vtx], RED, 10.0);
      
      // 3. step from near distance to far distance, and calculate the 
      //    bounds of each polygon slice (intersection of cube and plane).
      XYZ p1, p2;
      int plidx; // plane index
      float fracdist; // 0 to 1 (near to far)
      XYZ pip; // point-in-plane
      XYZ pipvd; // point-in-plane, but along viewdir: should be centred!
      PLANE theplane;
      double mu;
      XYZ pt, pt2;
      
      // and here we place up to 6 vertices for a sliced polygon
      int npolyverts;
      int polyverts[6]; // which edge?
      float polyfracs[6]; // how far along edge?
      // and this is the position angle of the vertex in the viewplane
      float polyangs[6];
      
      int j,k;
      float ang;
      
      float xpts[7], ypts[7], zpts[7];
      XYZ iP[6], iTC[6];
      
      int NPL; // the number of planes we will draw
      
      // now scale NPL by dot product of (nearvtx - farvtx) . viewdir
      // because this says what is the "depth" of planes...
      pt2 = VectorSub(worldverts[wh][near_vtx], worldverts[wh][far_vtx]);
      ang = DotProduct(pt2, viewdir);
      pt2 = VectorSub(worldverts[wh][0], worldverts[wh][7]); // diagonal
      ang /= Modulus(pt2);
      NPL = (int)(2.0 * ang * sqrt((long)_nx[0]*_nx[0]+_ny[0]*_ny[0]+_nz[0]*_nz[0]) );
      
      NPL = NPL / _3d_speed_up;

      // if you halve it, you end up with about nx (ny, nz) slices, but
      // some aliasing, so leave it at "double density"
      //NPL *= 0.5;
      //fprintf(stderr, "NPL = %d\n", NPL);
      
      // loop in reverse order so farthest planes added to list (and then
      // drawn) first.
      
      // (in frames)
#define SLICE_PERIOD (20*15)
      int sub = NPL;
      static int framebase;
      if (do_slice_reveal != 0) {
	// start slice reveal
	if ((do_slice_reveal == 1) || (do_slice_reveal == -1)) {
	  framebase = framecount;
	  do_slice_reveal *= 2;
	}
	
	// end slice reveal
	//else if ((do_slice_reveal == 2) || (do_slice_reveal == -2)) {
	
	if (framecount - framebase > SLICE_PERIOD) {
	  do_slice_reveal = 0;
	} else {
	  
	  // hack up forward / backward slicing orthogonal to view dir
	  float domax = (framecount - framebase);
	  sub = (int)((domax / SLICE_PERIOD) * NPL);
	  sub = (sub < 0) ? sub : ((sub > NPL-1) ? NPL-1 : sub);
	}
	//}
	
	sub = (do_slice_reveal < 0) ? NPL-sub : sub;
      }
      
      for (plidx = NPL; plidx > NPL-sub; plidx--) {
	
	for (zh = 0; zh < NVR; zh++) {
	  
	  fracdist = ((float)plidx -(float)zh/(float)NVR) / (float)(NPL+1);
	  //                      ^^^^^^^^^^^^^^^^^^^^^^^^^  here is offset per volume 
	  thisdist = near_dist + fracdist * (far_dist-near_dist);
	  
	  // point-in-plane along near_vtx to far_vtx line
	  pip = VectorSub(worldverts[wh][near_vtx], worldverts[wh][far_vtx]);
	  pip = VectorMul(pip, fracdist);
	  pip = VectorAdd(worldverts[wh][near_vtx], pip);
	  //- ns2vthpoint(pip, YELLOW, 2.0);
	  
	  // plane equation: for n={a,b,c}, the plane is n.p=-d, giving
	  // ax+by+cz+d = 0
	  // So all we do is calculate what d is:
	  theplane.a = viewdir.x;
	  theplane.b = viewdir.y;
	  theplane.c = viewdir.z;
	  theplane.d = -1. * DotProduct(viewdir, pip);
	  
	  // point-in-plane along viewdir
	  p2 = VectorAdd(campos, viewdir);
	  if (!LinePlane(campos, p2, theplane, &mu, &pipvd)) {
	    fprintf(stderr, "Viewdir doesn't intersect plane: impossible!!!\n");
	    exit(-1);
	  }
	  //- ns2vthpoint(pipvd, YELLOW, 2.0);
	  
	  npolyverts = 0;
	  
	  for (i = 0; i < 12; i++) {
	    p1 = worldverts[zh][edges[i][0]];
	    p2 = worldverts[zh][edges[i][1]];
	    if (LinePlane(p1, p2, theplane, &mu, &pt)) {
	      if ((mu >= 0) && (mu <= 1.)) {
		// get position angle of vertex
		pt2 = VectorSub(pipvd, pt);
		Normalise(&pt2);
		ang = atan2(DotProduct(pt2, upvec), DotProduct(pt2, right));
		// and insert in list
		j = 0;
		while ((j < npolyverts) && (polyangs[j] < ang)) {
		  j++;
		}
		k = npolyverts - 1;
		while (k >= j) {
		  polyverts[k+1] = polyverts[k];
		  polyfracs[k+1] = polyfracs[k];
		  polyangs[k+1] = polyangs[k];
		  k--;
		}
		k++;
		polyverts[k] = i;
		polyfracs[k] = mu;
		polyangs[k] = ang;
		npolyverts++;
		//- ns2vthpoint(pt, YELLOW, 4.0);
	      }
	    }
	  }
	  
	  // ok, we have the edges, fraction along those edges, and the 
	  // position angle in the view plane of each vertex of this poylgon.
	  // Now we need to draw the polygon in eg. clockwise order...
	  for (i = 0; i < npolyverts; i++) {
	    p1 = worldverts[zh][edges[polyverts[i]][0]];
	    p2 = worldverts[zh][edges[polyverts[i]][1]];
	    mu = polyfracs[i];
	    pt = VectorAdd(p1, VectorMul(VectorSub(p1, p2), mu));
	    xpts[i] = pt.x;
	    ypts[i] = pt.y;
	    zpts[i] = pt.z;
	    
	    if (i == 0) {
	      xpts[npolyverts] = pt.x;
	      ypts[npolyverts] = pt.y;
	      zpts[npolyverts] = pt.z;
	    }
	    
	    // here are the XYZ arrays for 3d texturing via ns2texpoly3d...
	    iP[i] = pt;
	    p1 = unitverts[edges[polyverts[i]][0]];
	    p2 = unitverts[edges[polyverts[i]][1]];
	    iTC[i] = VectorAdd(p1, VectorMul(VectorSub(p1, p2), mu));
	    
	  }
	  
	  if (ss2qrm() == WIREFRAME) {
	    //s2sci(2 + (plidx % 12));
	    s2sci(2 + (zh % MAXNVOL));
	    //fprintf(stderr, "s2sci zh = %d\n", zh);
	    s2line(npolyverts+1, xpts, ypts, zpts);
	  } else {
	    if (do_slice_reveal && (plidx == NPL-sub+1)) {
	      //s2sci(S2_PG_ORANGE);
	      //s2line(npolyverts+1, xpts, ypts, zpts);
	      ns2texpoly3d(iP, iTC, npolyverts, tex3d[zh], itrans, 1.0);	  
	    }
	    ns2texpoly3d(iP, iTC, npolyverts, tex3d[zh], itrans, ialpha);
	  }
	  
	}
      }
     } 
    } else {
      
      static int ctr = 0;
      ctr++;
      /* need to repair for multiple volumes before re-enabling this code
	 if (glow_period_frames > 0.) {
	 float phase = (float)ctr / glow_period_frames;
	 float amp_mod = 0.5 + 0.5 * sin(phase * M_PI * 2);
	 amp_mod = amp_mod * amp_mod * amp_mod;
	 fprintf(stderr, "amp_mod = %f\t", amp_mod);
	 ns2svrl(vidx[0], dmin[0], dmax[0], amp_mod * amin[0] * dscale[0], amp_mod * amax[0] * dscale[0]);
	 ns2svrl(vidy[0], dmin[0], dmax[0], amp_mod * amin[0] * dscale[1], amp_mod * amax[0] * dscale[1]);
	 ns2svrl(vidz[0], dmin[0], dmax[0], amp_mod * amin[0] * dscale[2], amp_mod * amax[0] * dscale[2]);
	 }
      */
      // ds2dvrx draws the given volume rendering along the given axis, provided
      // the given axis is the appropriate one to draw, OR doPRC is true.
      if ((!suppress_vr) && !(numcount[1] % 2))  {
	pushVRMLname("VRSET1");
	ds2dvrx(vidx[wh], glow_period_frames > 0., 1, doPRC);
	pushVRMLname("VRSET2");
	ds2dvrx(vidy[wh], glow_period_frames > 0., 2, doPRC);
	pushVRMLname("VRSET3");
	ds2dvrx(vidz[wh], glow_period_frames > 0., 3, doPRC);
	pushVRMLname("ANON");
      }
    }
  }
    
  float ddx, ddy, ddz;
  float tx, ty, tz;
  s2qwin(&tx, &ddx, &ty, &ddy, &tz, &ddz);
  if ((nobjs > 0) && ((*kc+1) % 2)) {
    //fprintf(stderr, "WARNING: OBJS NO LONGER LIKELY TO HAVE CORRECT WCS!!!\n");
    // (fixed I think by below offset, reoffset)
    s2swin(0,(ddx-tx), 0,(ddy-ty), 0,(ddz-tz));
    int i;
    for (i = 0; i < nobjs; i++) {
      drawObj(objs[i]);
    }
    s2swin(tx,ddx, ty,ddy, tz,ddz);
  }

  if (start && numcount[2] % 2) {
    // draw tracks
    pushVRMLname("TRACKS");
    COLOUR col1 = {1., 1., 1.};
    //HSV hsv1 = {0., 0.5, 1.};

    float xpts[2], ypts[2], zpts[2];
    int i;
    for (i = 0; i < ntracks; i++) {
      int j;
      XYZ diff;
      for (j = 0; j < count[i]-1; j++) {
	
	// col1 standard rgb direction encoding (per segment)
	diff = VectorSub(start[i][j], start[i][j+1]);
	Normalise(&diff);
	col1.r = fabs(diff.x);
	col1.g = fabs(diff.y);
	col1.b = fabs(diff.z);

#if (0)
	// hsv1: try mapping total track length to saturation
	//       y-component effect track y-component to hue
	//       value to 1 for now
	diff = VectorSub(start[i][0], start[i][count[i]-1]);
	Normalise(&diff);
	hsv1.h = fabs(diff.y) * 360.;
	hsv1.s = fabs(diff.z);
	hsv1.v = fabs(diff.x);
	if (hsv1.s > 1.) {
	  hsv1.s = 1.;
	}
	
	if (hsv1.v > 1.) {
	  hsv1.v = 1.;
	}
	col1 = HSV2RGB(hsv1);
#endif
	
	if (1) {
	  ns2vthcline(start[i][j], start[i][j+1], col1, col1, 2.8);
	  
	} else {
	  
	  _s2_conelines = 1;
	  xpts[0] = start[i][j].x;
	  xpts[1] = start[i][j+1].x;
	  ypts[0] = start[i][j].y;
	  ypts[1] = start[i][j+1].y;
	  zpts[0] = start[i][j].z;
	  zpts[1] = start[i][j+1].z;
	  s2scr(11111, col1.r, col1.g, col1.b);
	  s2sci(11111);
	  s2slw(0.2);
	  s2line(2, xpts, ypts, zpts);
	}
      }
    }
    
    pushVRMLname("ANON");
  }
  
}

int kcb(unsigned char *key) {
  // override processing of 'P' to write PRC so we can make sure the user knows
  // to have prc enabled in the vol renering code
  if ((*key == 'P') && !doPRC) {
    fprintf(stderr, "TOGGLE VR MODE TO ALL SLICE SETS WITH 'B' THEN TRY AGAIN!!!\n");
    return 1;
  } else if (*key == 'B') {
    doPRC = !doPRC;
  }
  return 0; // not consumed
}


void usage(char *exename) {
  fprintf(stderr, "usage: %s [options] -f xrwfilename\n", exename);
  fprintf(stderr, " * * * per volume (<= 4) options\n");
  fprintf(stderr, "-s s1 s2 s3    set data averaging along each axis (1 1 1)\n");
  fprintf(stderr, "-2             force power-of-2 sidelength\n");
  fprintf(stderr, "-d d1 d2       set data min, max (default 0.0 1.0)\n");
  fprintf(stderr, "-a a1 a2       set alpha min, max (default 0.0 0.2)\n");
  fprintf(stderr, "-v v1 v2 v3    set data/alpha scaling per axis\n");
  fprintf(stderr, "-c name        set colourmap name (default use xrw cmap)\n");
  fprintf(stderr, "-p p1          set scaling power for data->colourmap (default 0.0)\n");
  fprintf(stderr, "-z             calculate and render volume derivative\n");
  fprintf(stderr, "-I lev r g b a isosurface level [0,1] (instead of volren)\n");
  fprintf(stderr, "-M filename    use matrix stored in filename\n");

  fprintf(stderr, " * * * surfaces and tracks\n");
  fprintf(stderr, "-o objfname label r g b a  draw surface in named .obj file\n");
  fprintf(stderr, "-Rx r1         radial blow-apart by factor\n");
  fprintf(stderr, "-Vx            calculate and use vertex normals\n");
  fprintf(stderr, "-Tf            (ascii) track name file pattern\n");
  fprintf(stderr, "-Ts            track file skip factor (default 1000)\n");

  fprintf(stderr, " * * * global options\n");
  fprintf(stderr, "-SV            suppress volume rendering\n");
  fprintf(stderr, "-DX            clip data outside data min, max\n");
  fprintf(stderr, "-BX n          zero n boundary planes each side\n");
  fprintf(stderr, "-i             invert (paper colour)\n");
  fprintf(stderr, "-x             don't draw box\n");
  fprintf(stderr, "-Vf filename   read view from filename\n");
  fprintf(stderr, "-NN            don't over-ride number keys for views\n");

  fprintf(stderr, " * * * following options defy PDF output\n");
  fprintf(stderr, "-3d            use 3-d rendering\n");
  fprintf(stderr, "-Sp            speed-up for 3-d rendering (default = 1)\n");
}

int main(int argc, char *argv[]) {

  int i; 					/* Loop variable */

  // per-volume options
  char ifname[MAXNVOL][FNAMELEN+1];
  int havefname[MAXNVOL];
  int stride[MAXNVOL][3];
  int pow2[MAXNVOL];
  int use_xrw_cmap[MAXNVOL];
  char cmapname[MAXNVOL][FNAMELEN+1];
  int matrix_supplied[MAXNVOL];
  char matrix_fname[MAXNVOL][FNAMELEN];

  float matrix[4][4];

  for (i = 0; i < MAXNVOL; i++) {
    havefname[i] = 0;
    stride[i][0] = stride[i][1] = stride[i][2] = 1;
    pow2[i] = 0;
    dmin[i] = 0.0;
    dmax[i] = 1.0;
    amin[i] = 0.0;
    amax[i] = 0.2;
    dscale[i][0] = dscale[i][1] = dscale[i][2] = 1.0;
    scapower[i] = 1.0; /// argh this used to be 0.0 (how silly!)
    doDeriv[i] = 0;
    use_xrw_cmap[i] = 1;
    matrix_supplied[i] = 0;
    tex3d[i] = -1;
    isolev[i] = -1.;
    isocol[i].r = isocol[i].g = isocol[i].b = 0.0;
    isoalf[i] = 0.0;
    vidx[i] = vidy[i] = vidz[i] = -1;
  }

  int drawbox = 1;
  float blow_apart = 1.0;
  int vxnormals = 0;

  int track_fpatt_supplied = 0;
  char track_fpatt[FNAMELEN];
  int track_skip = 1000;

  int clip_data = 0;
  int zero_boundaries = 0;

  int normal_numbers = 0;

  if (argc < 3) {
    //fprintf(stderr, "usage: %s xrwfilename\n", argv[0]);
    usage(argv[0]);
    return -1;
  }

  int ic = 1;
  while (ic < argc) {
    if (!strcmp(argv[ic], "-f")) {
      NVR++;
      if (NVR > 4) {
	fprintf(stderr, "Too many volumes!\n");
	return -1;
      }
      strncpy(ifname[NVR-1], argv[++ic], FNAMELEN);
      havefname[NVR-1] = 1;
    } else if (!strcmp(argv[ic], "-s")) {
      stride[NVR-1][0] = atoi(argv[++ic]);
      stride[NVR-1][1] = atoi(argv[++ic]);
      stride[NVR-1][2] = atoi(argv[++ic]);
      if (stride[NVR-1][0] < 1 || stride[NVR-1][1] < 1 || stride[NVR-1][2] < 1) {
	fprintf(stderr, "stride must be at least 1 in each direction!\n");
	exit(-1);
      }
    } else if (!strcmp(argv[ic], "-2")) {
      pow2[NVR-1] = 1;
    } else if (!strcmp(argv[ic], "-d")) {
      dmin[NVR-1] = atof(argv[++ic]);
      dmax[NVR-1] = atof(argv[++ic]);
    } else if (!strcmp(argv[ic], "-a")) {
      amin[NVR-1] = atof(argv[++ic]);
      amax[NVR-1] = atof(argv[++ic]);
    } else if (!strcmp(argv[ic], "-v")) {
      dscale[NVR-1][0] = atof(argv[++ic]);
      dscale[NVR-1][1] = atof(argv[++ic]);
      dscale[NVR-1][2] = atof(argv[++ic]);
    } else if (!strcmp(argv[ic], "-c")) {
      use_xrw_cmap[NVR-1] = 0;
      strncpy(cmapname[NVR-1], argv[++ic], FNAMELEN);
      fprintf(stderr, "cmapname = %s\n", cmapname[NVR-1]);
    } else if (!strcmp(argv[ic], "-p")) {
      scapower[NVR-1] = atof(argv[++ic]);
    } else if (!strcmp(argv[ic], "-z")) {
      doDeriv[NVR-1] = 1;
    } else if (!strcmp(argv[ic], "-M")) {
      matrix_supplied[NVR-1] = 1;
      strncpy(matrix_fname[NVR-1], argv[++ic], FNAMELEN);
    } else if (!strcmp(argv[ic], "-I")) {
      isolev[NVR-1] = atof(argv[++ic]);
      float r, g, b, a;
      r = atof(argv[++ic]);
      g = atof(argv[++ic]);
      b = atof(argv[++ic]);
      a = atof(argv[++ic]);
      isocol[NVR-1].r = r;
      isocol[NVR-1].g = g;
      isocol[NVR-1].b = b;
      isoalf[NVR-1] = a;

    } else if (!strcmp(argv[ic], "-o")) {
      nobjs++;
      char objname[FNAMELEN+1];
      char label[OBJ_STRUCT_LABEL_LEN+1];
      strncpy(objname, argv[++ic], FNAMELEN);
      strncpy(label, argv[++ic], OBJ_STRUCT_LABEL_LEN);
      float r, g, b, a;
      r = atof(argv[++ic]);
      g = atof(argv[++ic]);
      b = atof(argv[++ic]);
      a = atof(argv[++ic]);
      objs = (OBJ_STRUCT **)realloc(objs, nobjs * sizeof(OBJ_STRUCT *));
      if (!strncmp(objname+strlen(objname)-3, "obj", 3)) {
	objs[nobjs-1] = loadObj(objname, label, r, g, b, a);
      } else if (!strncmp(objname+strlen(objname)-3, "stl", 3)) {
	objs[nobjs-1] = loadObjFromSTL(objname, label, r, g, b, a);
      } else {
	fprintf(stderr, "unknown input surface format: %s\n", objname);
	exit(-1);
      }
    } else if (!strcmp(argv[ic], "-i")) {
      invert = 1;
    } else if (!strcmp(argv[ic], "-x")) {
      drawbox = 0;
    } else if (!strcmp(argv[ic], "-Rx")) {
      blow_apart = atof(argv[++ic]);
    } else if (!strcmp(argv[ic], "-Vx")) {
      vxnormals = 1;
    } else if (!strcmp(argv[ic], "-SV")) {
      suppress_vr = 1;
    } else if (!strcmp(argv[ic], "-Tf")) {
      track_fpatt_supplied = 1;
      strncpy(track_fpatt, argv[++ic], FNAMELEN);
    } else if (!strcmp(argv[ic], "-Ts")) {
      track_skip = atoi(argv[++ic]);
    } else if (!strcmp(argv[ic], "-3d")) {
      do_3d_render = 1;
    } else if (!strcmp(argv[ic], "-Sp")) {
      _3d_speed_up = atoi(argv[++ic]);
      fprintf(stderr, "speed-up factor %d\n", _3d_speed_up);
    } else if (!strcmp(argv[ic], "-DX")) {
      clip_data = 1;
    } else if (!strcmp(argv[ic], "-BX")) {
      zero_boundaries = atoi(argv[++ic]);
    } else if (!strcmp(argv[ic], "-Vf")) {
      viewfile = 1;
      strncpy(viewfile_fname, argv[++ic], FNAMELEN);
    } else if (!strcmp(argv[ic], "-NN")) {
      normal_numbers = 1;
    }
    ic++;
  }
  
  if (!havefname[0]) {
    usage(argv[0]);
    return -1;
  }
    

  if (track_fpatt_supplied) {
    XYZ *verts = NULL; // store all vertices here for all tracks
    long nverts = 0;
    long lenverts = 0;
    
    long lentracks = 0;
    
    int idx = 0;
    char infname[FNAMELEN];
    
    sprintf(infname, track_fpatt, idx);
    FILE *fin;
    fin = fopen(infname, "r");
    XYZ P;
    while (fin) {
      int nv = 0;
      while (fscanf(fin, "%lf %lf %lf", &P.x, &P.y, &P.z) == 3) {
	if (nverts >= lenverts) {
	  lenverts += VERT_CHUNK;
	  verts = (XYZ *)realloc(verts, lenverts * sizeof(XYZ));
	}
	verts[nverts] = P;
	nverts++;
	nv++;
      }
      
      if (ntracks >= lentracks) {
	lentracks += TRACK_CHUNK;
	count = (int *)realloc(count, lentracks * sizeof(int));
      }
      count[ntracks] = nv;
      
      fclose(fin);
      idx += track_skip;
      ntracks += 1;
      sprintf(infname, track_fpatt, idx);
      fin = fopen(infname, "r");
    }
    
    // calculate start offsets
    start = (XYZ **)malloc(ntracks * sizeof(XYZ *));
    int i;
    long accum = 0;
    for (i = 0; i < ntracks; i++) {
      start[i] = verts + accum;
      accum += count[i];
    }
    
    fprintf(stderr, "loaded %ld track files, lengths are:\n", ntracks);
    for (i = 0; i < ntracks; i++) {
      fprintf(stderr, "%7d ", count[i]);
    }    
  }


  float wx1,wx2,wy1,wy2,wz1,wz2;


  int wh;
  // major loop over volume/s
  for (wh = 0; wh < NVR; wh++) {
    
    XRAW_STRUCT *xr = loadXraw(ifname[wh]);
    if (!xr) {
      fprintf(stderr, "Failed to open or read '%s'.\n", ifname[wh]);
      return -1;
    }
    
    showXraw(xr);
    fprintf(stdout, "- - - - - - - - - - - - - - - - - - - - - - - - -\n");
    
    VOL_STRUCT *vol = Xraw2Xvol(xr, stride[wh]);
    if (!vol) {
      fprintf(stderr, "Failed to parse data volume.\n");
      return -1;
    }

    if (pow2[wh]) {
      VOL_STRUCT *vol2 = makePow2Xvol(vol);
      VOL_STRUCT *volbak = vol;
      vol = vol2;
      deleteXvol(volbak);
    }
    
    if (doDeriv[wh]) {
      // derivative of volume
      fprintf(stderr, " * * * DERIVATIVE OF VOLUME * * *\n");
      derivXvol(vol);
    }

    showXvol(vol);
    
    int nx, ny, nz;
    nx = vol->nx;
    ny = vol->ny;
    nz = vol->nz;
    
    if (wh == 0) {
      //s2opend("/S2MONO",argc,argv);			/* Open the display */
      s2opend("/?",argc,argv);			/* Open the display */
      if (invert) {
	ss2sbc(1., 1.,1.);
	ss2sfc(0., 0., 0.);
      } else {
	//fprintf(stderr, "normal white-on-black\n");
	ss2sbc(0., 0., 0.);
	ss2sfc(1., 1., 1.);
      }
      // viewport: simply generate appropriate aspect ratio
      float sx = 1.0;
      float sy = vol->wdy*(float)ny / (vol->wdx*(float)nx);
      float sz = vol->wdz*(float)nz / (vol->wdx*(float)nx);
      s2svp(-sx,sx, -sy,sy, -sz,sz);
    }

    // world coords: depends on whether a matrix was supplied
    float tr[12];
    for (i=0;i<12;i++) {				/* Set-up transfrom matrix */
      tr[i] = 0.0;
    }
        
    if (matrix_supplied[wh]) {
      FILE *MATI = fopen(matrix_fname[wh], "r");
      int mi, mj;
      for (mi = 0; mi < 4; mi++) {
	fscanf(MATI, "%f%f%f%f", &(matrix[mi][0]), &(matrix[mi][1]),
	       &(matrix[mi][2]), &(matrix[mi][3]));
      }
      fclose(MATI);
      
      // correct for stride
      for (mi = 0; mi < 3; mi++) {
	for (mj = 0; mj < 3; mj++) {
	  matrix[mi][mj] *= stride[wh][mj];
	}
      }
      
      for (mi = 0; mi < 4; mi++) {
	for (mj = 0; mj < 4; mj++) {
	  fprintf(stderr, "%f ", matrix[mi][mj]);
	}
	fprintf(stderr, "\n");
      }
      
      if (wh == 0) {

	XYZ ijk_trc, uvw_blc, uvw_trc;
	ijk_trc.x = nx-1;
	ijk_trc.y = ny-1;
	ijk_trc.z = nz-1;
	uvw_blc.x = matrix[0][3];
	uvw_blc.y = matrix[1][3];
	uvw_blc.z = matrix[2][3];
	uvw_trc.x = matrix[0][0] * ijk_trc.x /*+ matrix[0][1] * ijk_trc.y + matrix[0][2] * ijk_trc.z*/ + matrix[0][3];
	uvw_trc.y = /*matrix[1][0] * ijk_trc.x +*/ matrix[1][1] * ijk_trc.y /*+ matrix[1][2] * ijk_trc.z*/ + matrix[1][3];
	uvw_trc.z = /*matrix[2][0] * ijk_trc.x + matrix[2][1] * ijk_trc.y +*/ matrix[2][2] * ijk_trc.z + matrix[2][3];
	
	fprintf(stderr, "blc = %f,%f,%f; trc = %f,%f,%f\n", 
		uvw_blc.x, uvw_blc.y, uvw_blc.z,
		uvw_trc.x, uvw_trc.y, uvw_trc.z);
	
	float tswp;
	if (uvw_blc.x > uvw_trc.x) {
	  tswp = uvw_blc.x;
	  uvw_blc.x = uvw_trc.x;
	  uvw_trc.x = tswp;
	}
	if (uvw_blc.y > uvw_trc.y) {
	  tswp = uvw_blc.y;
	  uvw_blc.y = uvw_trc.y;
	  uvw_trc.y = tswp;
	}
	if (uvw_blc.z > uvw_trc.z) {
	  tswp = uvw_blc.z;
	  uvw_blc.z = uvw_trc.z;
	  uvw_trc.z = tswp;
	}
	
	s2swin(uvw_blc.x,uvw_trc.x, uvw_blc.y,uvw_trc.y, uvw_blc.z,uvw_trc.z);
      }
      
      
      // x
      tr[0] = matrix[0][3]; // offset
      tr[1] = matrix[0][0]; // dx
      tr[2] = matrix[0][1]; // dy
      tr[3] = matrix[0][2]; // dz
      
      // y
      tr[4] = matrix[1][3]; // offset
      tr[5] = matrix[1][0]; // dx
      tr[6] = matrix[1][1]; // dy
      tr[7] = matrix[1][2]; // dz
      
      // z
      tr[8] = matrix[2][3]; // offset
      tr[9] = matrix[2][0]; // dx
      tr[10] = matrix[2][1]; // dy
      tr[11] = matrix[2][2]; // dz
      
    } else {
      
      // world coordinate system: put centre of volume at (0,0,0)
      // world coordinates are physical units
      float ddx = vol->wdx*(float)nx * 0.5;
      float ddy = vol->wdy*(float)ny * 0.5;
      float ddz = vol->wdz*(float)nz * 0.5;

      if (wh == 0) {
	s2swin(-ddx,ddx, -ddy,ddy, -ddz,ddz);
      }
      
      tr[0] = -vol->wdx*(float)nx*0.5;
      tr[1] = vol->wdx;
      tr[4] = -vol->wdy*(float)ny*0.5;
      tr[6] = vol->wdy;
      tr[8] = -vol->wdz*(float)nz*0.5;
      tr[11]= vol->wdz;
    }
    
    s2qwin(&wx1,&wx2,&wy1,&wy2,&wz1,&wz2);
    fprintf(stderr, "s2swin(%f,%f,%f,%f,%f,%f)\n", wx1,wx2, wy1,wy2, wz1,wz2);
  
    if (wh == 0) {
      // set up lights
      COLOUR ambient = {0.7,0.7,0.7};
      XYZ lightpos[2];
      COLOUR lightcol[2];
      lightpos[0].x = wx1 + 0.5 * (wx1-wx2);
      lightpos[0].y = wy1 + 0.5 * (wy1-wy2);
      lightpos[0].z = wz1 + 0.5 * (wz1-wz2);
      lightpos[1].x = wx2 + 0.5 * (wx2-wx1);
      lightpos[1].y = wy2 + 0.5 * (wy2-wy1);
      lightpos[1].z = wz2 + 0.5 * (wz2-wz1);
      lightcol[0].r = lightcol[0].g = lightcol[0].b = 1.0;
      lightcol[1].r = lightcol[1].g = lightcol[1].b = 1.0;
      ss2sl(ambient, 2, lightpos, lightcol, 1);
      
      if (drawbox) {
	s2box("BCDE",0,0,"BCDE",0,0,"BCDE",0,0);		/* Draw coord box */  
      }
    }
  
    char trans = 's';			/* Transparency type */
    
    
    // install colour LUT
    s2scir(1000+wh*256, 1000+wh*256+255);
    if (use_xrw_cmap[wh]) {
      for (i = 0; i < 256; i++) {
	s2scr(1000+wh*256+i, vol->red[i], vol->green[i], vol->blue[i]);
      }
    } else {
      s2icm(cmapname[wh], 1000+wh*256, 1000+wh*256+255);
    }

    // brighten the colourmap because we have to use 's' (absorption trans) mode
    // to match look of PDF render
    float r, g, b;
    for (i = 0; i < 256; i++) {
      s2qcr(1000+wh*256+i, &r, &g, &b);
      r = powf(r, 0.2);
      g = powf(g, 0.2);
      b = powf(b, 0.2);
      s2scr(1000+wh*256+i, r, g, b);
    }
    
    // clip data?
    if (clip_data) {
      int j, k;
      for (i = 0; i < nx; i++) {
	for (j = 0; j < ny; j++) {
	  for (k = 0; k < nz; k++) {
	    // this should work because 2-d texture creation (internal to 
	    // _s2_priv_load_vr ... clips at min (sets alpha = 0.0) while
	    // 3-d texture creation (below) clips at min as well
	    if ((vol->data[i][j][k] < dmin[wh]) || (vol->data[i][j][k] > dmax[wh])) {
	      vol->data[i][j][k] = 0.0;
	    }
	  }
	}
      }
    }
    
    // zero boundary planes
    if (zero_boundaries) {
      int j, k;
      int bd;
      for (bd = 0; bd < zero_boundaries; bd++) {
	for (i = 0; i < nx; i++) {
	  for (j = 0; j < ny; j++) {
	    for (k = 0; k < nz; k++) {
	      if ((i == bd) || (i == nx-1-bd) ||
		  (j == bd) || (j == ny-1-bd) ||
		  (k == bd) || (k == nz-1-bd)) {
		vol->data[i][j][k] = 0.0;
	      }
	    }
	  }
	}
      }
    }
    
    if (isolev[wh] >= 0.0) {
      vidx[wh] = ns2cis(vol->data, nx, ny, nz, 0, nx-1, 0, ny-1, 0, nz-1, 
			tr, isolev[wh], 1, trans, isoalf[wh],
			isocol[wh].r, isocol[wh].g, isocol[wh].b);

    } else 

    if (do_3d_render) {
      // create and fill 3d texture
      tex3d[wh] = ss2c3dt(nx, ny, nz);
      
      // set globals 
      int j, k;
      _nx[wh] = nx;
      _ny[wh] = ny;
      _nz[wh] = nz;
      for (i = 0; i < 12; i++) {
	_tr[wh][i] = tr[i];
      }
      
      float ir, ig, ib, alf, frac;
      int idx1, idx2, cidx;
      long idx;
      int w,h,d;
      //float minaf =9e30, maxaf = -9e30;
      //s2qcir(&idx1, &idx2);
      idx1 = 1000 + wh*256;
      idx2 = 1000+wh*256+255;
      //      s2scir(1000+wh*256, 1000+wh*256+255);
      float denom_recip = 1.0 / (dmax[wh] - dmin[wh]);
      unsigned char *bits = (unsigned char *)ss2g3dt(tex3d[wh], &w, &h, &d);
      fprintf(stderr, "3d texture IS %d x %d x %d\n", w, h, d);
      //memset(bits, (unsigned char)0, (long)w*h*d*4);
      //fprintf(stderr, "memset done.\n");
      float ***volume = vol->data;
      for (i = 0; i < nx; i++) {
	for (j = 0; j < ny; j++) {
#pragma omp parallel for private (k,frac,cidx,alf,idx,ir,ig,ib)
	  for (k = 0; k < nz; k++) {
	    //if ((volume[i][j][k] > dmin[0]) && (volume[i][j][k] < dmax[0])) {
	    frac = (volume[i][j][k] - dmin[wh]) * denom_recip;
	    frac = (frac < 0.) ? 0. : ((frac > 1.) ? 1. : frac);
	    //frac = drand48()/RAND_MAX;
	    frac = powf(frac, scapower[wh]);
	    cidx = idx1 + (int)(frac * (float)(idx2-idx1));
	    cidx = (cidx < idx1) ? idx1 : ((cidx > idx2) ? idx2 : cidx);

	    //cidx = (idx1 + idx2) / 2;

	    //if (i == nx/2) {
	      //fprintf(stderr, "%8f    %8d   \n", frac, cidx);
	    //}
	    s2qcr(cidx, &ir, &ig, &ib);
	    alf = amin[wh] + (amax[wh] - amin[wh]) * (frac);
	    if (volume[i][j][k] < dmin[wh]) {
	      alf = 0.;
	    }
	    idx = (long)4 * ((long)i + (long)(j * nx) + (long)((long)k * nx * ny));
	    bits[idx + 0] = (int)(ir * 255.0);
	    bits[idx + 1] = (int)(ig * 255.0);
	    bits[idx + 2] = (int)(ib * 255.0);
	    bits[idx + 3] = (int)(alf* 255.0);
	    //if (bits[idx+3] < minaf) {
	    //  minaf = bits[idx+3];
	    //}
	    //if (bits[idx+3] > maxaf) {
	    //  maxaf = bits[idx+3];
	    //}
	    //}
	  }
	}
      }
      fprintf(stderr, "about to push texture...\n");
      ss2ptt(tex3d[wh]);
      fprintf(stderr, "pushed texture done.\n");
      //fprintf(stderr, "minaf = %f, maxaf = %f\n", minaf, maxaf);
      
    } else {
      vidx[wh] = ns2cvr(vol->data, nx, ny, nz, 0, nx-1, 0, ny-1, 0, nz-1, 
		       tr, trans, dmin[wh], dmax[wh],
		       amin[wh]*dscale[wh][0], amax[wh]*dscale[wh][0]);
      vidy[wh] = ns2cvr(vol->data, nx, ny, nz, 0, nx-1, 0, ny-1, 0, nz-1, 
		       tr, trans, dmin[wh], dmax[wh],
		       amin[wh]*dscale[wh][1], amax[wh]*dscale[wh][1]);
      vidz[wh] = ns2cvr(vol->data, nx, ny, nz, 0, nx-1, 0, ny-1, 0, nz-1, 
		       tr, trans, dmin[wh], dmax[wh],
		       amin[wh]*dscale[wh][2], amax[wh]*dscale[wh][2]);
    }
    
  }

  if (vxnormals) {
    fprintf(stderr, "calculating vertex normals...\n");
    for (i = 0; i < nobjs; i++) {
      calcVertexNormals(objs[i]);
      if (!objs[i]->norms) {
	fprintf(stderr, "huh?\n");
      }
    }
  }
  
  fprintf(stderr, "BLOW APART might not work due to multivolume support!\n");
  // blow apart the surfaces
  // origin of world space (= middle of displayed figure)
  XYZ origin;
  origin.x = 0.5 * (wx1+wx2); //ddx;
  origin.y = 0.5 * (wy1+wy2); //ddy;
  origin.z = 0.5 * (wz1+wz2); // ddz;  
  //s2qwin(&tx, &ddx, &ty, &ddy, &tz, &ddz);
  for (i = 0; i < nobjs; i++) {
    
    // position of surface (ie. mean vertex coord)
    XYZ geom = objs[i]->meanP;
    // offset
    XYZ off = VectorSub(geom, origin);
    
    // blowapart 1 = do nothing
    // blowapart >1 = expand away from origin
    // blowapart <1 = shrink towards origin
    XYZ sc_off = VectorMul(off, blow_apart);
    
    XYZ delta = VectorSub(off, sc_off);
    
    translateObj(objs[i], delta);
  }
  
  // install the callbacks
  cs2scb(cb);
  if (!normal_numbers) {
    cs2sncb(numcb);
  }
  cs2skcb(kcb);
  
  // set initial state: not in PRC mode
  doPRC = 0;

  //ss2txh(1);
  ss2srm(SHADE_FLAT);


  int stereo, fs, dome;
  ss2qsa(&stereo, &fs, &dome);
  if (stereo < 0) { // null device
    double t = 0.0;
    int kc = 0;
    doPRC = 1;
    cb(&t, &kc); // but what about billboards?
    writePRC();
  } else {
    s2show(1);					/* Show the S2PLOT window */
  }
  
  return 1;
}


