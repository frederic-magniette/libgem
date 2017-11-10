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

//************************** POINTS *******************************

struct point *new_point(int dim) {
  struct point *result=malloc(sizeof(struct point));
  int i;
  result->dim=dim;
  result->coords=malloc(dim*sizeof(double));
  for(i=0;i<dim;i++)
    result->coords[i]=0;
  return result;
}

struct point *new_valued_point(int dim,double *values) {
  struct point *result=new_point(dim);
  int i;
  for(i=0;i<dim;i++)
    result->coords[i]=values[i];
  return result;
}

struct point *new_random_point(int dim,struct point *pmin,struct point *pmax) {
  struct point *result=new_point(dim);
  int i;
  double range;
  for(i=0;i<dim;i++) {
    range=fabs(pmax->coords[i]-pmin->coords[i]);
    result->coords[i]=pmin->coords[i]+range/4+(rand()/((double)RAND_MAX)*range/2);
  }
  return result;
}

struct point *cp_point(struct point *p) {
  int i;
  struct point *result=new_point(p->dim);
  for(i=0;i<p->dim;i++)
    result->coords[i]=p->coords[i];
  return result;
}

void print_point(struct point *p) {
  int i;
  for(i=0;i<p->dim;i++) {
    printf("p[%d]=%f ",i,p->coords[i]);
  }
  printf("\n");
}

void dump_point(struct point *p,FILE *f) {
  int i;
  for(i=0;i<p->dim;i++) {
    fprintf(f,"%f ",p->coords[i]);
  }
  fprintf(f,"\n");
}


void free_point(struct point *p) {
  free(p->coords);
  free(p);
}

double dist_points(struct point *p1,struct point *p2) {
  struct vector *bp=vector_by2points(p1,p2);
  double result=norme_vector(bp);
  free_vector(bp);
  return result;
}

void plot_point(struct point *p,int id,struct graphics *gws) {
  int rref,gref,bref;
  FILE *f;
  char name[1024];
  int d;
  char title[1024];
  if (gws->gp) {
    get_plotfile(gws->gp,name);
    f=fopen(name,"w");
    for(d=0;d<p->dim;d++) {
      fprintf(f,"%lf ",p->coords[d]);
    }
    fprintf(f,"\n");
    fclose(f);
    sprintf(title,"point %d",id);
    plot_file(gws->gp,p->dim,name,0,1,title);
  }
  if (gws->sp) {
    attribute_color_splot(id,&rref,&gref,&bref);
    plot_real_point_splot(p,gws->sp,rref,gref,bref);
    /*printf("plotting point :");
    print_point(p);
    printf("with color %d,%d,%d\n",rref,gref,bref);*/
  }
}

void add_vector_point(struct point *p,struct vector *v) {
  int d;
  for(d=0;d<p->dim;d++) {
    p->coords[d]+=v->coords[d];
  }
} 

int is_valid_point(struct point *p,struct point *pmin,struct point *pmax) {
  int i;
  double val;
  double tolerance=0.000001;
  for(i=0;i<p->dim;i++) {
      val=p->coords[i];
      if (val<(pmin->coords[i]-tolerance)) {
        return 0;
      }
      
      if (val>(pmax->coords[i]+tolerance)) {
        return 0;
      }
  }
  return 1;
}

struct point *ortho_project_point(struct point *p,struct matrix *base,struct vector *dir,struct point *ref) {

  //project the point p in the plan orthogonal to dir and containing ref
  //return its coordinates in the projective base 

  struct vector *pr;
  struct vector *pr_proj;
  struct point *projection;
  struct point *result;
  struct vector *move;
  double longitudinal_dist;
  int i;

  //get the longitudinal distance between p and ref
  pr=vector_by2points(p,ref);
  printf("vector pr=");
  print_vector(pr);
  printf("norme=%f\n",norme_vector(pr));
  pr_proj=change_base_vector(pr,base);
  printf("proj vector=");
  print_vector(pr_proj);
  longitudinal_dist=pr_proj->coords[0];
  printf("long dist=%f\n",longitudinal_dist);
  free_vector(pr);
  free_vector(pr_proj);
  

  //build the projection point
  move=cp_vector(dir);
  normalize_vector(move);
  mult_lambda_vect(move,longitudinal_dist);
  projection=cp_point(p);
  add_vector_point(projection,move);
  free_vector(move);

  //convert coordinates to projective base
  for(i=0;i<ref->dim;i++) {
    projection->coords[i]-=ref->coords[i];
  }
  result=change_base_point(projection,base);
  result->coords[0]=longitudinal_dist;
  free_point(projection);

  //return the result
  return result;
}

struct point *cut_first_coord_point(struct point *p) {
  int i;
  struct point *result=new_point(p->dim-1);
  for(i=1;i<p->dim;i++)
    result->coords[i-1]=p->coords[i];
  return result;
}

struct point *change_base_point(struct point *p,struct matrix *base) {
  return mult_point_matrix(base,p);
}

struct point *barycenter_point(struct dataset *ds) {
  struct point *result=new_point(ds->dim);
  int i,j;
  for(i=0;i<ds->nb_points;i++) {
    for(j=0;j<ds->dim;j++) {
      result->coords[j]+=ds->points[i]->coords[j];
    }
  }
  for(j=0;j<ds->dim;j++) {
    result->coords[j]/=ds->nb_points;
  }
  return result;
}

struct point *weighted_barycenter_point(struct dataset *ds,struct weights *w) {
  struct point *result=new_point(ds->dim);
  double sum_w=0;
  int p,d;
  for(p=0;p<ds->nb_points;p++) {
    sum_w+=w->coefs[p];
    for(d=0;d<ds->dim;d++) {
      //printf("adding %lf to result[%d]\n",ds->points[p]->coords[d]*w->coefs[p],d);
      result->coords[d]+=ds->points[p]->coords[d]*w->coefs[p];
    }
  }
  if (sum_w==0) {
    free_point(result);
    printf("error : null weights for barycenter\n");
    return new_point(ds->dim);
  }
  for(d=0;d<ds->dim;d++) {
    result->coords[d]/=sum_w;
  }
  return result;
}
