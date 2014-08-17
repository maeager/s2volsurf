/* s2texobj.c: demonstrate drawing textured OBJ
 *
 * David Barnes, August 2014
 */

#include "s2plot.h"
#include "libobj.c"

int main(int argc, char **argv) {

  s2opend("/?", argc, argv);
  
  OBJ_STRUCT *obj = loadObj("MU3248lowres_textured_pshop-resave.obj", "label", 1., 1., 1., 1.);

  calcObjMinMaxMean(obj);

  XYZ tran = VectorMul(VectorAdd(obj->minP, obj->maxP), 0.5);
  translateObj(obj, tran);

  float sx = 1.0;
  float sy = (obj->maxP.y - obj->minP.y) / (obj->maxP.x - obj->minP.x);
  float sz = (obj->maxP.z - obj->minP.z) / (obj->maxP.x - obj->minP.x);
  s2svp(-sx,sx, -sy,sy, -sz,sz);

  float m[16];
  int i;

  // scale - the only real point here is to shrink the model to be within the lighting environment!
  for (i = 0; i < 16; i++) {
    m[i] = 0.0;
  }
  m[0] = 0.02;
  m[5] = 0.02;
  m[10] = 0.02;
  transformObj(obj, m);

#define SC 1.4
  s2swin(SC*obj->minP.x, SC*obj->maxP.x, SC*obj->minP.y, SC*obj->maxP.y, SC*obj->minP.z, SC*obj->maxP.z);
  
  unsigned int texid = ss2lt("MU3248lowres_texturedImage10000.tga");

  //drawObj(obj);
  drawObjAsTexturedMesh(obj, texid);

  s2show(1);

  return 0;
}
