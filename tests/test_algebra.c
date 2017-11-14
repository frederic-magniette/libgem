/*
Copyright 2016-2018 Frédéric Magniette
This file is part of Libgem.

Libgem is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Libgem is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Libgem.  If not, see <http://www.gnu.org/licenses/>
*/

#include <libgem.h>

void test_gs() {
  double vals[3]={0.0,12.37,0};
  struct vector *v=new_valued_vector(3,vals);
  struct matrix *m=gram_schmidt(v);
  print_matrix(m);
  free_matrix(m);
  free_vector(v);
}


void test_proj() {
  struct dataset *ds=new_dataset_fromfile(2,"data/test_algebra.txt");
  struct graphics *gws=new_graphics(1,500,500,ds->point_min,ds->point_max,3);
  double vals[2]={0,1};
  struct vector *v=new_valued_vector(2,vals);
  struct matrix *m=gram_schmidt(v);
  struct point *ref=new_point(2);
  struct dataset *pds=ortho_projection_dataset(ds,m,v,ref,NULL);
  print_dataset(pds);
  plot_dataset(pds,gws);
  apply_graphics(gws);
  getchar();
}


int main() {

  test_proj();
  return 0;
}
