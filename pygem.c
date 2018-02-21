#include <Python.h>
#include <libgem.h>

static PyObject * dump_angle_line(PyObject *self, PyObject *args) {
  char *filename;
  float angle;
  int nb_points;
  struct line *output;
  struct point *point_min, *point_max;
  struct dataset *ds;
  double dispersion;
  int rnbp;
  
  if (!PyArg_ParseTuple(args,"fisf",&angle,&nb_points,&filename,&dispersion))
    return NULL;
  
  output=angle_line(angle,0,0);
  point_min=new_2D_point(-100,-100);
  point_max=new_2D_point(100,100);
  ds=boxed_line(output,point_min,point_max,nb_points);
  noise_dataset(ds,dispersion);
  dump_dataset(ds,filename);
  rnbp=ds->nb_points;
  free_point(point_min);
  free_point(point_max);
  free_line(output);
  free_dataset(ds);
  return Py_BuildValue("i",rnbp);
}

static PyMethodDef pygem_methods[] = {
  {"dump_angle_line",dump_angle_line,METH_VARARGS,"Dump angle line."},
  {NULL, NULL, 0, NULL} 
};

PyMODINIT_FUNC initpygem(void) {
  (void) Py_InitModule("pygem",pygem_methods);
}
