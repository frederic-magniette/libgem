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

#include "libgem.h"

//orthogonal projection of gram-schmidt based on v as first vector + random vectors + vectors normalization -> result is orthonormal base. result[1:] define plan wich is orthogonal to v
struct matrix *gram_schmidt(struct vector *v) {
  int i,j;
  struct vector **base=malloc(v->dim*sizeof(struct vector *));
  struct matrix *result=new_matrix(v->dim,v->dim);
  struct vector *tmp;
  struct vector *proj;
  double norme;
  base[0]=cp_vector(v);
  normalize_vector(base[0]);
  tmp=NULL;
  for(i=1;i<v->dim;i++) {
    do {
      if (tmp!=NULL)
        free_vector(tmp);
      tmp=new_random_vector(v->dim,1);
      for(j=0;j<i;j++) {
        proj=project_vector(tmp,base[j]);
        sub_vector(tmp,proj);
        free_vector(proj);
      }
      //print_vector(tmp);
      norme=norme_vector(tmp);
      //printf("norme=%f\n",norme);
      if (norme==0)
        continue;
      normalize_vector(tmp);
    } while (0);
    base[i]=tmp;
    tmp=NULL;
  }
  for(i=0;i<v->dim;i++) {
    set_col_matrix(result,base[i],i);
    free_vector(base[i]);
  }
  free(base);
  return result;
}

//project the points on a plan containing ref and orthogonal to v
struct point **orthonormal_projection(struct point **points,int nb_points,struct vector *v,struct point *ref) {
  struct matrix *base=gram_schmidt(v);
  struct point **result=malloc(nb_points*sizeof(struct point *));
  int i;
  for(i=0;i<nb_points;i++) {
    result[i]=ortho_project_point(points[i],base,v,ref); 
  }
  free_matrix(base);
  return result;
}
