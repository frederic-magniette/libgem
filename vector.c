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

//************************* VECTORS ****************************

struct vector *new_vector(int dim) {
  int i;
  struct vector *result=malloc(sizeof(struct vector));
  result->dim=dim;
  result->coords=malloc(sizeof(double)*dim);
  for(i=0;i<dim;i++)
    result->coords[i]=0;
  return result;
}

struct vector *new_valued_vector(int dim,double *values) {
  int i;
  struct vector *result=malloc(sizeof(struct vector));
  result->dim=dim;
  result->coords=malloc(sizeof(double)*dim);
  for(i=0;i<dim;i++)
    result->coords[i]=values[i];
  return result;
}

struct vector *new_random_vector(int dim,double max_value) {
  int i;
  struct vector *result=malloc(sizeof(struct vector));
  result->dim=dim;
  result->coords=malloc(sizeof(double)*dim);
  for(i=0;i<dim;i++)
    result->coords[i]=max_value*rand()/((double)RAND_MAX);
  return result;
}


void zero_vector(struct vector *v,double value) {
  int i;
  for(i=0;i<v->dim;i++)
    v->coords[i]=value;
}

int is_zero_vector(struct vector *v,double limit) {
  int i;
  for(i=0;i<v->dim;i++)
    if ((v->coords[i]>limit) || (v->coords[i]<0-limit))
      return 0;
  return 1;
}

struct vector *new_orthonormal_vector(int dim,int one_index) {
  int i;
  struct vector *result=malloc(sizeof(struct vector));
  result->dim=dim;
  result->coords=malloc(sizeof(double)*dim);
  for(i=0;i<dim;i++)
    if (i==one_index)
      result->coords[i]=1;
    else
      result->coords[i]=0;
  return result;
}

void free_vector(struct vector *v) {
  free(v->coords);
  free(v);
}

void print_vector(struct vector *v) {
  int i;
  printf("vector[%d]=",v->dim);
  for(i=0;i<v->dim;i++) {
    printf("v[%d]=%f ",i,v->coords[i]);
  }
  printf("\n");
}

void mult_lambda_vect(struct vector *v,double lambda) {
  int i;
  for(i=0;i<v->dim;i++) {
    v->coords[i]=v->coords[i]*lambda;
  }
}

//add src to dst
void add_vector(struct vector *dst,struct vector *src) {
  int i;
  for(i=0;i<dst->dim;i++) {
    dst->coords[i]=dst->coords[i]+src->coords[i];
  }
}

//substract src to dst
void sub_vector(struct vector *dst,struct vector *src) {
  int i;
  for(i=0;i<dst->dim;i++) {
    dst->coords[i]=dst->coords[i]-src->coords[i];
  }
}

//return the vector v1+v2
struct vector *sum_vector(struct vector *v1,struct vector *v2) {
  int i;
  struct vector *result=new_vector(v1->dim);
  for(i=0;i<v1->dim;i++) {
    result->coords[i]=v1->coords[i]+v2->coords[i];
  }
  return result;
}

//project vector source on vector support <u,v>/<u,u> u
struct vector *project_vector(struct vector *src,struct vector *support) {
  double factor;
  double denom;
  struct vector *result=cp_vector(support);
  denom=scalar_vector(support,support);
  if (denom==0) {
    return result;
  }
  factor=scalar_vector(support,src)/denom;
  mult_lambda_vect(result,factor);
  return result;
}

double coef_project_vector(struct vector *src,struct vector *support) {
  double denom;
  denom=scalar_vector(support,support);
  if (denom==0) {
    return 0;
  }
  return scalar_vector(support,src)/denom;
}

struct vector *cp_vector(struct vector *v) {
  int i;
  struct vector *result=new_vector(v->dim);
  for(i=0;i<v->dim;i++)
    result->coords[i]=v->coords[i];
  return result;
}

//return 1 if v1=v2 0 otherwise
int eq_vector(struct vector *v1,struct vector *v2) {
  int i;
  for(i=0;i<v1->dim;i++)
    if (v1->coords[i]!=v2->coords[i])
      return 0;
  return 1;
}
 
struct vector *vector_by2points(struct point *origin,struct point *dest) {
  int i;
  struct vector *result=new_vector(origin->dim);
  for(i=0;i<origin->dim;i++) {
    result->coords[i]=dest->coords[i]-origin->coords[i];
  }
  return result;
}

struct vector *vector_abs_by2points(struct point *origin,struct point *dest) {
  int i;
  struct vector *result=new_vector(origin->dim);
  for(i=0;i<origin->dim;i++) {
    result->coords[i]=fabs(dest->coords[i]-origin->coords[i]);
  }
  return result;
}

struct vector *vector_polarized_by2points(struct point *origin,struct point *dest) {
  int i;
  struct vector *result=vector_by2points(origin,dest);
  for(i=0;i<origin->dim;i++) {
    if (result->coords[i]>0)
      return result;
    if (result->coords[i]<0) {
      mult_lambda_vect(result,-1);
      return result;
    }
  }
  return result;
}

double scalar_vector(struct vector *v1,struct vector *v2) {
  int i;
  double result=0;
  for(i=0;i<v1->dim;i++) {
    result+=v1->coords[i]*v2->coords[i];
  }
  return result;
}

double norme_vector(struct vector *v) {
  return sqrt(scalar_vector(v,v));
}

void normalize_vector(struct vector *v) {
  double norme=norme_vector(v);
  if (norme==0)
    return;
  mult_lambda_vect(v,1/norme);
}


struct vector *weighted_mean_vector(struct dataset *ds,struct point *ref,struct weights *w) {
  int i,j;
  struct vector *result=new_vector(ds->dim);
  struct vector *vtemp;
  double sum_weights=0;

  for(i=0;i<ds->nb_points;i++) {
    sum_weights+=w->coefs[i];
  }
  if (sum_weights==0) {
    printf("empty weights\n");
    return result;
  }

  for(i=0;i<ds->nb_points;i++) {
    sum_weights+=w->coefs[i];
    vtemp=vector_polarized_by2points(ref,ds->points[i]);

    //sum the weighted value
    for(j=0;j<result->dim;j++) {
      result->coords[j]+=vtemp->coords[j]*w->coefs[i];
    }
    free_vector(vtemp);
  }
  //printf("sum_weights=%f\n",sum_weights);
  for(i=0;i<result->dim;i++) {
    result->coords[i]=result->coords[i]/sum_weights;
  }
  normalize_vector(result);
  return result;
}

void plot_vector(struct vector *v,struct point *origin,int id,struct graphics *gws) {
  struct point *dest=cp_point(origin);
  add_vector_point(dest,v);
  plot_point(dest,id,gws);
  free_point(dest);
}

int isnan_vector(struct vector *v) {
  int i;
  for(i=0;i<v->dim;i++) {
    if (isnan(v->coords[i]))
      return 1;
  }
  return 0;
}

struct vector *change_base_vector(struct vector *v,struct matrix *base) {
  return mult_vector_matrix(base,v);
}

double angle_between_2_vector(struct vector *v1,struct vector *v2) {
  double angle;
  double cosa;
  double sign;
  if ((v1->dim!=2)||(v2->dim!=2)) {
    printf("error : angle_between_2_vector works only in 2d but d(v1)=%d and d(v2)=%d\n",v1->dim,v2->dim);
    return -1;
  }
  //printf("norme(v1)=%f\n",norme_vector(v1));
  //printf("norme(v2)=%f\n",norme_vector(v2));
  //printf("<v1,v2>=%f\n",scalar_vector(v1,v2));

  sign=v1->coords[0]*v2->coords[1]-v1->coords[1]*v2->coords[0];
  cosa=scalar_vector(v1,v2)/(norme_vector(v1)*norme_vector(v2));
  //printf("cos a = %f\n",cosa);
  //printf("angle = %f rad\n",acos(cosa));
  angle=rad2deg(acos(cosa));
  //printf("non signed angle=%f deg \n",angle);
  //printf("sign=%f\n",sign);
  if (sign>0)
    return angle;
  else
    return -1*angle;
}

struct vector *cut_first_coord_vector(struct vector *v) {
  int i;
  struct vector *result=new_vector(v->dim-1);
  for(i=1;i<v->dim;i++)
    result->coords[i-1]=v->coords[i];
  return result;
}

