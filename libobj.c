/* libobj.c
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
 * $Id: libobj.c 86 2013-01-10 00:20:44Z barnesd $
 *
 */

#define OBJ_STRUCT_LABEL_LEN 80
typedef struct {
  char label[OBJ_STRUCT_LABEL_LEN+1];
  COLOUR col;
  float alpha;
  int nverts;
  XYZ *verts;
  int nfacets;
  int *facets; // 3 * int per facet = 3 * vertex indices
  int *facets_tcs; // 3 * int per facet = 3 * vertex texture coordinate indices
  XYZ minP, maxP; // min,max coordinates
  XYZ meanP; // mean (~central) coordinate
  int nnorms;
  XYZ *norms; // require as many as verts
  int nvtcs; // vertex texture coordinates
  XYZ *vtcs; // use only X,Y (as u,v)
} OBJ_STRUCT;

typedef struct {
  int nverts;
  int *vert_indices;
  float *vert_weights;
  float minW, maxW;
} OBJ_WGT_STRUCT;

OBJ_STRUCT *loadObj(char *filename, char *label, float r, float g, float b, float a);
OBJ_STRUCT *loadObjFromFS(char *filename, char *label, float r, float g, float b, float a);
OBJ_WGT_STRUCT *loadObjWgtFromFS(char *filename);
OBJ_STRUCT *loadObjFromSTL(char *filename, char *label, float r, float g, float b, float a);
void calcObjMinMaxMean(OBJ_STRUCT *obj);
void translateObj(OBJ_STRUCT *obj, XYZ vec);
void transformObj(OBJ_STRUCT *obj, float *m); // m[4c,3r] matrix
void drawObj(OBJ_STRUCT *obj); //, COLOUR c);
void drawObjWgts(OBJ_STRUCT *obj, OBJ_WGT_STRUCT *wgt, int wgted_vx_only,
		 int solid_col_idx);
void saveWgtedObj(OBJ_STRUCT *obj, OBJ_STRUCT *obj_norms, OBJ_WGT_STRUCT *wgt, char *fname);


OBJ_STRUCT *loadObj(char *filename, char *label, float r, float g, float b, float a) {
  
  OBJ_STRUCT *obj = (OBJ_STRUCT *)malloc(1 * sizeof(OBJ_STRUCT));

  strncpy(obj->label, label, OBJ_STRUCT_LABEL_LEN);
  obj->col.r = r;
  obj->col.g = g;
  obj->col.b = b;
  obj->alpha = a;

  obj->nverts = 0;
  obj->verts = NULL;
  obj->nfacets = 0;
  obj->facets = NULL;
  obj->facets_tcs = NULL;

  obj->nnorms = 0;
  obj->norms = NULL;

  obj->nvtcs = 0;
  obj->vtcs = NULL;

  XYZ one = {1,1,1};
  XYZ mone = {-1,-1,-1};
  XYZ zero = {0,0,0};
  obj->minP = mone;
  obj->maxP = one;
  obj->meanP = zero;

  //return obj;

  int lcount = 0;
  char line[401];
  //char cmd;
  //float f1, f2, f3;

  FILE *fin = fopen(filename, "r");
  while (fgets(line, 400, fin)) {
    lcount++;
#if (1)
    char *result = NULL;
    char space[] = " ";
    char slash[] = "/";
    char dslash[] = "//";
    result = strtok(line, space);
    //fprintf(stderr, "token = %s\n", result);
    if (!strcmp(result, "v")) {
      obj->nverts++;
      obj->verts = (XYZ *)realloc(obj->verts, obj->nverts * sizeof(XYZ));
      obj->verts[obj->nverts-1].y = atof(strtok(NULL, space));
      obj->verts[obj->nverts-1].z = atof(strtok(NULL, space));
      obj->verts[obj->nverts-1].x = atof(strtok(NULL, space));
    } else if (!strcmp(result, "vt")) {
      obj->nvtcs++;
      obj->vtcs = (XYZ *)realloc(obj->vtcs, obj->nvtcs * sizeof(XYZ));
      obj->vtcs[obj->nvtcs-1].x = atof(strtok(NULL, space));
      obj->vtcs[obj->nvtcs-1].y = atof(strtok(NULL, space));
    } else if (!strcmp(result, "vn")) {
      obj->nnorms++;
      obj->norms = (XYZ *)realloc(obj->norms, obj->nnorms * sizeof(XYZ));
      obj->norms[obj->nnorms-1].y = atof(strtok(NULL, space));
      obj->norms[obj->nnorms-1].z = atof(strtok(NULL, space));
      obj->norms[obj->nnorms-1].x = atof(strtok(NULL, space));
    } else if (!strcmp(result, "f")) {
      char *s1 = strtok(NULL, space);
      char *s2 = strtok(NULL, space);
      char *s3 = strtok(NULL, space);
      //fprintf(stderr, "%s : %s : %s\n", s1, s2, s3);
      if (strstr(s1, "//") && strstr(s2, "//") && strstr(s3, "//")) {
	static int ig_normal_msg = 0;
	if (!ig_normal_msg) {
	  fprintf(stderr, "* * * WARNING: vertex//normal format present: normals ignored!\n");
	  ig_normal_msg = 1;
	}
	obj->nfacets++;
	obj->facets = (int *)realloc(obj->facets, obj->nfacets * 3 * sizeof(int));
	//obj->facets_tcs = (int *)realloc(obj->facets_tcs, 
	//				 obj->nfacets * 3 * sizeof(int));
	obj->facets[(obj->nfacets-1)*3 + 0] = atoi(strtok(s1, dslash)) - 1;
	//obj->facets_tcs[(obj->nfacets-1)*3+0]= atoi(strtok(NULL, slash)) - 1;
	strtok(NULL, dslash);
	obj->facets[(obj->nfacets-1)*3 + 1] = atoi(strtok(s2, dslash)) - 1;
	//obj->facets_tcs[(obj->nfacets-1)*3+1]= atoi(strtok(NULL, slash)) - 1;
	strtok(NULL, dslash);
	obj->facets[(obj->nfacets-1)*3 + 2] = atoi(strtok(s3, dslash)) - 1;
	//obj->facets_tcs[(obj->nfacets-1)*3+2]= atoi(strtok(NULL, slash)) - 1;
	strtok(NULL, dslash);
	

      } else if (strstr(s1, "/") && strstr(s2, "/") && strstr(s3, "/")) {
	//fprintf(stderr, "parsing v/vt1\n");
	obj->nfacets++;
	obj->facets = (int *)realloc(obj->facets, obj->nfacets * 3 * sizeof(int));
	obj->facets_tcs = (int *)realloc(obj->facets_tcs, 
					 obj->nfacets * 3 * sizeof(int));
	obj->facets[(obj->nfacets-1)*3 + 0] = atoi(strtok(s1, slash)) - 1;
	obj->facets_tcs[(obj->nfacets-1)*3+0]= atoi(strtok(NULL, slash)) - 1;
	obj->facets[(obj->nfacets-1)*3 + 1] = atoi(strtok(s2, slash)) - 1;
	obj->facets_tcs[(obj->nfacets-1)*3+1]= atoi(strtok(NULL, slash)) - 1;
	obj->facets[(obj->nfacets-1)*3 + 2] = atoi(strtok(s3, slash)) - 1;
	obj->facets_tcs[(obj->nfacets-1)*3+2]= atoi(strtok(NULL, slash)) - 1;
      } else {
	obj->nfacets++;
	obj->facets = (int *)realloc(obj->facets, obj->nfacets * 3 * sizeof(int));
	obj->facets[(obj->nfacets-1)*3 + 0] = atoi(s1) - 1;
	obj->facets[(obj->nfacets-1)*3 + 1] = atoi(s2) - 1;
	obj->facets[(obj->nfacets-1)*3 + 2] = atoi(s3) - 1;
      }
    } else {
      //fprintf(stderr, "%s encountered (ignored)\n", result);
    }	

#else
    if (sscanf(line, "%c %f %f %f", &cmd, &f1, &f2, &f3) == 4) {
      switch(cmd) {
      case 'v':
	obj->nverts++;
	obj->verts = (XYZ *)realloc(obj->verts, obj->nverts * sizeof(XYZ));
	obj->verts[obj->nverts-1].y = f1;
	obj->verts[obj->nverts-1].z = f2;
	obj->verts[obj->nverts-1].x = f3;
	break;

      case 'f':
	obj->nfacets++;
	obj->facets = (int *)realloc(obj->facets, obj->nfacets * 3 * sizeof(int));
	obj->facets[(obj->nfacets-1)*3 + 0] = (int)(f1+0.5);
	obj->facets[(obj->nfacets-1)*3 + 1] = (int)(f2+0.5);
	obj->facets[(obj->nfacets-1)*3 + 2] = (int)(f3+0.5);
	break;

      default:
	fprintf(stderr, "Unknown obj file cmd '%c' at line %d, ignoring.\n", cmd, lcount);
	break;

      }
    } else {
      fprintf(stderr, "Could not parse obj file at line %d, ignoring.\n", lcount);
    }
#endif
  }
  
  fprintf(stderr, "Read %d lines from '%s'.\n", lcount, filename);
  fprintf(stderr, "\t%d vertices\n", obj->nverts);
  fprintf(stderr, "\t%d normals\n", obj->nnorms);
  fprintf(stderr, "\t%d texture coordinates\n", obj->nvtcs);
  fprintf(stderr, "\t%d facets\n", obj->nfacets);


  calcObjMinMaxMean(obj);

  return obj;
}

// load obj struct from *ASCII* freesurfer binary file, format descr here:
// http://wideman-one.com/gw/brain/fs/surfacefileformats.htm
// binary format was "awkward" to read because it had 3-byte ints: what is
// their endianness?!?!?!?!
OBJ_STRUCT *loadObjFromFS(char *filename, char *label, float r, float g, float b, float a) {
  OBJ_STRUCT *obj = (OBJ_STRUCT *)malloc(1 * sizeof(OBJ_STRUCT));

  strncpy(obj->label, label, OBJ_STRUCT_LABEL_LEN);
  obj->col.r = r;
  obj->col.g = g;
  obj->col.b = b;
  obj->alpha = a;

  obj->nverts = 0;
  obj->verts = NULL;
  obj->nfacets = 0;
  obj->facets = NULL;

  obj->norms = NULL;

  char line[401];
  float f1, f2, f3;
  int d1, d2, d3, d4;

  FILE *fin = fopen(filename, "r");

  fgets(line, 400, fin); // read comment
  fprintf(stderr, "read comment: %s\n", line);
  fgets(line, 400, fin); // read nvert, nface
  sscanf(line, "%d%d", &obj->nverts, &obj->nfacets);
  fprintf(stderr, "nverts = %d, nfacets = %d\n", obj->nverts, obj->nfacets);
  obj->verts = (XYZ *)malloc(obj->nverts * sizeof(XYZ));
  obj->facets = (int *)malloc(obj->nfacets * 3 * sizeof(int));
 
    
  // read verts
  int i;
  for (i = 0; i < obj->nverts; i++) {
    if (!fgets(line, 400, fin)) {
      fprintf(stderr, "failure reading file at vertex %d\n", i);
      exit(-1);
    }
    if (sscanf(line, "%f%f%f%d", &f1, &f2, &f3, &d4) != 4) {
      fprintf(stderr, "failure parsing file at vertex %d\n", i);
      exit(-1);
    }
    if (d4 != 0) {
      fprintf(stderr, "patch != 0\n");
      exit(-1);
    }
    //if (!(i % 100)) {
    //  fprintf(stderr, "x = %f\n", f1);
    //}
    obj->verts[i].x = f1;
    obj->verts[i].y = f2;
    obj->verts[i].z = f3;
  }

  // read facets
  for (i = 0; i < obj->nfacets; i++) {
    if (!fgets(line, 400, fin)) {
      fprintf(stderr, "failure reading file at facet %d\n", i);
      exit(-1);
    }
    if (sscanf(line, "%d%d%d%d", &d1, &d2, &d3, &d4) != 4) {
      fprintf(stderr, "failure parsing file at facet %d\n", i);
      exit(-1);
    }
    if (d4 != 0) {
      fprintf(stderr, "patch != 0\n");
      exit(-1);
    }
    obj->facets[i*3 + 0] = d1/*+1*/; //(int)(f1+0.5);
    obj->facets[i*3 + 1] = d2/*+1*/; //(int)(f2+0.5);
    obj->facets[i*3 + 2] = d3/*+1*/; //(int)(f3+0.5);
  }
  
  fprintf(stderr, "Read %d vertices (%d facets) from obj file.\n", obj->nverts, obj->nfacets);

  calcObjMinMaxMean(obj);

  return obj;
}

OBJ_WGT_STRUCT *loadObjWgtFromFS(char *filename) {
  OBJ_WGT_STRUCT *wgtobj = (OBJ_WGT_STRUCT *)malloc(1 * sizeof(OBJ_WGT_STRUCT));

  char line[401];
  wgtobj->minW = 9e99;
  wgtobj->maxW = -9e99;

  FILE *fin = fopen(filename, "r");

  fgets(line, 400, fin); // read (and ignore) "latency" line
  fgets(line, 400, fin); // read nvert
  sscanf(line, "%d", &wgtobj->nverts);
  fprintf(stderr, "nverts = %d\n", wgtobj->nverts);
  wgtobj->vert_indices = (int *)malloc(wgtobj->nverts * sizeof(int));
  wgtobj->vert_weights = (float *)malloc(wgtobj->nverts * sizeof(float));
  
  int i;
  for (i = 0; i < wgtobj->nverts; i++) {
    fgets(line, 400, fin);
    if (sscanf(line, "%d%f", &wgtobj->vert_indices[i], &wgtobj->vert_weights[i]) != 2) {
      fprintf(stderr, "Failed reading FS wgts file...\n");
      return NULL;
    }
    if (wgtobj->vert_weights[i] < wgtobj->minW) {
      wgtobj->minW = wgtobj->vert_weights[i];
    }
    if (wgtobj->vert_weights[i] > wgtobj->maxW) {
      wgtobj->maxW = wgtobj->vert_weights[i];
    }
  }

  return wgtobj;
}

typedef struct {
  float n[3];
  float vx[9];
  //unsigned short uu;
} /* __attribute__((packed)) */ STL_FACET_STRUCT;

#define VEPSI 0.0001

OBJ_STRUCT *loadObjFromSTL(char *filename, char *label, float r, float g, float b, float a) {

  OBJ_STRUCT *obj = (OBJ_STRUCT *)malloc(1 * sizeof(OBJ_STRUCT));

  strncpy(obj->label, label, OBJ_STRUCT_LABEL_LEN);
  obj->col.r = r;
  obj->col.g = g;
  obj->col.b = b;
  obj->alpha = a;

  obj->nverts = 0;
  obj->verts = NULL;
  obj->nfacets = 0;
  obj->facets = NULL;

  obj->norms = NULL;

  char title[81];
  
  FILE *fin = fopen(filename, "rb");
  if (!fin) {
    return NULL;
  }

  fread(title, sizeof(char), 80, fin);
  title[80] = '\0';

  unsigned int nf;
  fread(&nf, sizeof(nf), 1, fin);

  STL_FACET_STRUCT *stl_facet = (STL_FACET_STRUCT *)malloc(nf * sizeof(STL_FACET_STRUCT));
  int j;
  unsigned short uu;
  for (j = 0; j < nf; j++) {
    if (fread(&(stl_facet[j]), sizeof(STL_FACET_STRUCT), 1, fin) != 1) {
      fprintf(stderr, "failed reading facets\n");
      return NULL;
    }
    if (fread(&uu, sizeof(unsigned short), 1, fin) != 1) {
      fprintf(stderr, "failed reading unused short\n");
      return NULL;
    }
    if (uu) {
      //fprintf(stderr, "value of unused short is wrong!\n");
    }
  }
  
  unsigned char tmp;
  fread(&tmp, 1, 1, fin);
  if (!feof(fin)) {
    fprintf(stderr, "Not at end of file!\n");
  }
  
  int fidx, iv, sv;
  XYZ iP;
  obj->nfacets = nf;
  obj->facets = (int *)malloc(obj->nfacets * 3 * sizeof(int));

  for (fidx = 0; fidx < nf; fidx++) {
    // loop over input vertices
    for (iv = 0; iv < 3; iv++) {
      iP.x = stl_facet[fidx].vx[iv * 3 + 0];
      iP.y = stl_facet[fidx].vx[iv * 3 + 1];
      iP.z = stl_facet[fidx].vx[iv * 3 + 2];
#if (1)
      // find matching stored vertex
      sv = obj->nverts-1;
      while ((sv > -1) && (VectorLength(obj->verts[sv], iP) > VEPSI)) {
	sv--;
      }
      if (sv < 0) {
	sv = obj->nverts;
	obj->nverts++;
	obj->verts = (XYZ *)realloc(obj->verts, obj->nverts * sizeof(XYZ));
      }
      if (sv < 0 || sv >= obj->nverts) {
	fprintf(stderr, " * * * \n");
      }
      obj->verts[sv] = iP;
      
      obj->facets[fidx * 3 + iv] = sv/*+1*/;
#endif
    }
  }
  
  calcObjMinMaxMean(obj);

  return obj;
}

void calcObjMinMaxMean(OBJ_STRUCT *obj) {

  XYZ minP = {9e99,9e99,9e99};
  XYZ maxP = {-9e99,-9e99,-9e99};
  XYZ meanP = {0., 0., 0.};
  int i;

  for (i = 0; i < obj->nverts; i++) {

    if (obj->verts[i].x < minP.x) {
      minP.x = obj->verts[i].x;
    }
    if (obj->verts[i].x > maxP.x) {
      maxP.x = obj->verts[i].x;
    }
    
    if (obj->verts[i].y < minP.y) {
      minP.y = obj->verts[i].y;
    }
    if (obj->verts[i].y > maxP.y) {
      maxP.y = obj->verts[i].y;
    }
    
    if (obj->verts[i].z < minP.z) {
      minP.z = obj->verts[i].z;
    }
    if (obj->verts[i].z > maxP.z) {
      maxP.z = obj->verts[i].z;
    }

    meanP.x += obj->verts[i].x;
    meanP.y += obj->verts[i].y;
    meanP.z += obj->verts[i].z;

  }

  meanP.x /= (float)obj->nverts;
  meanP.y /= (float)obj->nverts;
  meanP.z /= (float)obj->nverts;

  obj->minP = minP;
  obj->maxP = maxP;
  obj->meanP = meanP;

  //fprintf(stderr, "Coordinate range is (%f,%f,%f) - (%f,%f,%f) [avg = (%f,%f,%f)]\n", 
  //	  obj->minP.x, obj->minP.y, obj->minP.z,
  //	  obj->maxP.x, obj->maxP.y, obj->maxP.z,
  //	  obj->meanP.x, obj->meanP.y, obj->meanP.z);
}

void calcVertexNormals(OBJ_STRUCT *obj) {

  if (obj->norms) {
    fprintf(stderr, "Already have normals. Not recalculating normals\n");
    return;
    free(obj->norms);
  }
  obj->norms = (XYZ *)malloc(obj->nverts * sizeof(XYZ));

  if (!obj->norms) {
    fprintf(stderr, "failed to allocate storage for normals!\n");
    exit(-1);
  }

  int i, j;
  int fcount;
  XYZ norm_acc, n;
  for (i = 0; i < obj->nverts; i++) {
    
    // find facets that use this vertex and accumulate the normals 
    fcount = 0;
    norm_acc.x = norm_acc.y = norm_acc.z = 0.0;

    for (j = 0; j < obj->nfacets; j++) {
      if ((obj->facets[j*3+0]/*-1*/ == i) ||
	  (obj->facets[j*3+1]/*-1*/ == i) ||
	  (obj->facets[j*3+2]/*-1*/ == i)) {
	fcount++;
	n = CalcNormal(obj->verts[obj->facets[j*3+0]/*-1*/],
		       obj->verts[obj->facets[j*3+1]/*-1*/],
		       obj->verts[obj->facets[j*3+2]/*-1*/]);
	// what weight should we add this with? ... vertex area?
	// use Heron's formula
	float a = VectorLength(obj->verts[obj->facets[j*3+0]/*-1*/],
			       obj->verts[obj->facets[j*3+1]/*-1*/]);
	float b = VectorLength(obj->verts[obj->facets[j*3+1]/*-1*/],
			       obj->verts[obj->facets[j*3+2]/*-1*/]);
	float c = VectorLength(obj->verts[obj->facets[j*3+2]/*-1*/],
			       obj->verts[obj->facets[j*3+0]/*-1*/]);
	float s = 0.5 * (a+b+c);
	float wgt = sqrt(s*(s-a)*(s-b)*(s-c));
	VectorMul(n, wgt);
	norm_acc = VectorAdd(norm_acc, n);
      }
    }
    //obj->norms[i] = VectorMul(norm_acc, -1);
    obj->norms[i] = norm_acc;
    if (fcount) {
      Normalise(&(obj->norms[i]));
    }

  }

}


void translateObj(OBJ_STRUCT *obj, XYZ vec) {
  int i;
  for (i = 0; i < obj->nverts; i++) {
    obj->verts[i] = VectorSub(vec, obj->verts[i]);
  }
  obj->minP = VectorSub(vec, obj->minP);
  obj->maxP = VectorSub(vec, obj->maxP);
}

void transformObj(OBJ_STRUCT *obj, float *m) {
  int i;
  XYZ in, out;
  for (i = 0; i < obj->nverts; i++) {
    in = obj->verts[i];
    out.x = m[0]*in.x + m[1]*in.y + m[2]*in.z + m[3];
    out.y = m[4]*in.x + m[5]*in.y + m[6]*in.z + m[7];
    out.z = m[8]*in.x + m[9]*in.y + m[10]*in.z + m[11];
    obj->verts[i] = out;
  }

  calcObjMinMaxMean(obj);  
}

void drawObjAsTexturedMesh(OBJ_STRUCT *obj, unsigned int texid) {
  
  ns2texmesh(obj->nverts, obj->verts,
	     obj->nnorms, obj->norms,
	     obj->nvtcs, obj->vtcs,
	     obj->nfacets, obj->facets, obj->facets_tcs,
	     texid, 'o', 1.0);

}

void drawObj(OBJ_STRUCT *obj) { // , COLOUR c) {

  int i, j;
  XYZ P[3];
  
  pushVRMLname(obj->label);

  for (i = 0; i < obj->nfacets; i++) {
    for (j = 0; j < 3; j++) {
      P[j] = obj->verts[obj->facets[i*3+j]];
    }
    if (obj->alpha < 0.999) {    
      if (obj->norms) {
	XYZ N[3];
	N[0] = obj->norms[obj->facets[i*3+0]];
	N[1] = obj->norms[obj->facets[i*3+1]];
	N[2] = obj->norms[obj->facets[i*3+2]];
	ns2vf3na(P, N, obj->col, 's', obj->alpha);
      } else {
	ns2vf3a(P, obj->col, 's', obj->alpha);
      }
    } else {
      if (obj->norms && 1) {
	//fprintf(stderr, "using normals...\n");
	XYZ N[3];
	N[0] = obj->norms[obj->facets[i*3+0]];
	N[1] = obj->norms[obj->facets[i*3+1]];
	N[2] = obj->norms[obj->facets[i*3+2]];
	ns2vf3n(P, N, obj->col);
      } else {
	ns2vf3(P, obj->col);
      }
    }
  }

  pushVRMLname("ANON");
}


void drawObjWgts(OBJ_STRUCT *obj, OBJ_WGT_STRUCT *wgt, int wgted_vx_only, 
		 int solid_col_idx) {

  int i, j, k;
  XYZ P[3];
  
  pushVRMLname(obj->label);

  int c1, c2, cidx;
  s2qcir(&c1, &c2);
  float r, g, b;

  COLOUR col;
  int draw_me;
  for (i = 0; i < obj->nfacets; i++) {
    draw_me = 0;
    col.r = col.g = col.b = 0.0;
    for (j = 0; j < 3; j++) {
      P[j] = obj->verts[obj->facets[i*3+j]/*-1*/];

      // do we have a wgt for this vertex?
      k = 0;
      while (k < wgt->nverts && wgt->vert_indices[k] != obj->facets[i*3+j]/*-1*/) {
	k++;
      }
      if (k < wgt->nverts) {
	//draw_me = 1;
	draw_me++;
	// map wgt [0,1] to colour in map
	float w = (wgt->vert_weights[k] - wgt->minW) / (wgt->maxW - wgt->minW);
	cidx = c1 + (int)(w * (c2-c1) + 0.5);
	if (cidx < c1) {
	  cidx = c1;
	}
	if (cidx > c2) {
	  cidx = c2;
	}
	s2qcr(cidx, &r, &g, &b);
	col.r += 2*r;
	col.g += 2*g;
	col.b += 2*b;
      } else {
	col.r += obj->col.r;
	col.g += obj->col.g;
	col.b += obj->col.b;
      }
    }
    col.r *= 0.3333;
    if (col.r > 1.) {
      col.r = 1.;
    }
    col.g *= 0.3333;
    if (col.g > 1.) {
      col.g = 1.;
    }
    col.b *= 0.3333;
    if (col.b > 1.) {
      col.b = 1.;
    }

    if (solid_col_idx > -1) {
      s2qcr(solid_col_idx, &r, &g, &b);
      col.r = r;
      col.g = g;
      col.b = b;
    }

    if (!wgted_vx_only || draw_me) {
      if (obj->alpha < 0.999) { 
	ns2vf3a(P, col, 's', obj->alpha);
      } else {	  
	ns2vf3(P, col);
      }
    }
  }

  pushVRMLname("ANON");
  
}


void saveWgtedObj(OBJ_STRUCT *obj, OBJ_STRUCT *obj_norms, OBJ_WGT_STRUCT *wgt, char *fname) {
  int i, j, k;

  FILE *fout = fopen(fname, "w");

  // 1. write vertices required
  for (i = 0; i < wgt->nverts; i++) {
    XYZ P = obj->verts[wgt->vert_indices[i]];
    fprintf(fout, "v %f %f %f\n", P.x, P.y, P.z);
  }

  // 1.5. write normals required
  for (i = 0; i < wgt->nverts; i++) {
    XYZ P = obj_norms->verts[wgt->vert_indices[i]];
    fprintf(fout, "vn %f %f %f\n", P.x, P.y, P.z);
  }

  // 2. go through all facets, and draw those that have wgts for all vertices
  for (i = 0; i < obj->nfacets; i++) {
    int draw_me = 0;
    int vidx[3];
    //XYZ P[3];
    for (j = 0; j < 3; j++) {
      //P[j] = obj->verts[obj->facets[i*3+j]/*-1*/];

      // do we have a wgt for this vertex?
      k = 0;
      while (k < wgt->nverts && wgt->vert_indices[k] != obj->facets[i*3+j]/*-1*/) {
	k++;
      }
      if (k < wgt->nverts) {
	//draw_me = 1;
	draw_me++;
	vidx[j] = k;
      }
    }
    if (draw_me == 3) {
      // write out this facet
      fprintf(fout, "f %d//%d %d//%d %d//%d\n", vidx[0]+1, vidx[0]+1,
	      vidx[1]+1, vidx[1]+1, vidx[2]+1, vidx[2]+1);
    }
  }

  fclose(fout);
  return;

#if (0)
 int j, k;
  XYZ P[3];
    pushVRMLname(obj->label);

  int c1, c2, cidx;
  s2qcir(&c1, &c2);
  float r, g, b;

  COLOUR col;
  int draw_me;
  for (i = 0; i < obj->nfacets; i++) {
    draw_me = 0;
    col.r = col.g = col.b = 0.0;
    for (j = 0; j < 3; j++) {
      P[j] = obj->verts[obj->facets[i*3+j]/*-1*/];

      // do we have a wgt for this vertex?
      k = 0;
      while (k < wgt->nverts && wgt->vert_indices[k] != obj->facets[i*3+j]/*-1*/) {
	k++;
      }
      if (k < wgt->nverts) {
	//draw_me = 1;
	draw_me++;
	// map wgt [0,1] to colour in map
	float w = (wgt->vert_weights[k] - wgt->minW) / (wgt->maxW - wgt->minW);
	cidx = c1 + (int)(w * (c2-c1) + 0.5);
	if (cidx < c1) {
	  cidx = c1;
	}
	if (cidx > c2) {
	  cidx = c2;
	}
	s2qcr(cidx, &r, &g, &b);
	col.r += 2*r;
	col.g += 2*g;
	col.b += 2*b;
      } else {
	col.r += obj->col.r;
	col.g += obj->col.g;
	col.b += obj->col.b;
      }
    }
    col.r *= 0.3333;
    if (col.r > 1.) {
      col.r = 1.;
    }
    col.g *= 0.3333;
    if (col.g > 1.) {
      col.g = 1.;
    }
    col.b *= 0.3333;
    if (col.b > 1.) {
      col.b = 1.;
    }

    if (solid_col_idx > -1) {
      s2qcr(solid_col_idx, &r, &g, &b);
      col.r = r;
      col.g = g;
      col.b = b;
    }

    if (!wgted_vx_only || draw_me) {
      if (obj->alpha < 0.999) { 
	ns2vf3a(P, col, 's', obj->alpha);
      } else {	  
	ns2vf3(P, col);
      }
    }
  }

  pushVRMLname("ANON");
#endif

}
