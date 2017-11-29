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

//************************* LINES ******************************

struct line *new_line(struct point *ref,struct vector *dir_vect) {
  struct line *result=malloc(sizeof(struct line));
  result->dim=ref->dim;
  result->ref=ref;
  result->dir_vect=dir_vect;
  return result;
}

struct line *cp_line(struct line *l) {
  struct line *result=malloc(sizeof(struct line));
  result->dim=l->dim;
  result->ref=cp_point(l->ref);
  result->dir_vect=cp_vector(l->dir_vect);
  return result;
}

void free_line(struct line *l) {
  free_vector(l->dir_vect);
  free_point(l->ref);
  free(l);
}

void print_line(struct line *l) {
  printf("line(dim=%d): reference :",l->dim);
  print_point(l->ref);
  printf("director vector:");
  print_vector(l->dir_vect);
}

double dist_line_point(struct line *l,struct point *p) {
  struct vector *bp;
  struct vector *u;
  double normu;
  double lambda;
  double result;

  if (l->dim==1) {
    return fabs(l->ref->coords[0]-p->coords[0]);
  }
  if (is_zero_vector(l->dir_vect,0.000001)) {
    return dist_points(l->ref,p);
  }

  bp=vector_by2points(l->ref,p);
  u=cp_vector(l->dir_vect);
  normu=norme_vector(u);
  lambda=scalar_vector(bp,u)*(-1)/(normu*normu);
  
  
  mult_lambda_vect(u,lambda);
  add_vector(bp,u);
  result=norme_vector(bp);
  free_vector(u);
  free_vector(bp);

  if (isnan(result)) {
    printf("ERROR : distance is nan\n");
    print_line(l);
    print_point(p);
  }

  return result;
}

//return 1 if l1 > l2, -1 if l1 < l2, 0 otherwise
int compare_lines(struct line *l1,struct line *l2) {
  int i;
  //printf("l1dim=%d\n",l1->dim);
  for(i=0;i<l1->dim;i++) {
    if (l1->dir_vect->coords[i]>l2->dir_vect->coords[i]) {
      //printf("compare positive\n");
      return 1;
    }
    if (l1->dir_vect->coords[i]<l2->dir_vect->coords[i]) {
      //printf("compare negative\n");
      return -1;
    }
  }
  //printf("compare null\n");
  return 0;
}

//sort a lines array based on compare_lines 
void sort_line_array(struct line **la,int la_size) {
  struct line *tmp;
  int i;
  int finish=0;
  while(finish==0) {
    finish=1;
    for(i=1;i<la_size;i++) {
      //printf("comparing %d and %d\n",i-1,i);
      if (compare_lines(la[i-1],la[i])==-1) {
	//printf("inverting %d and %d\n",i-1,i);
	tmp=la[i-1];
	la[i-1]=la[i];
	la[i]=tmp;
	finish=0;
      }
    }
  }
}

void dump_vect_line(struct line *l,char *filename) {
  int i;
  FILE *f=fopen(filename,"w");
  
  //dump p
  for(i=0;i<l->dim;i++) {
    fprintf(f,"%lf ",l->ref->coords[i]);
  }
      fprintf(f,"\n");

  //dump p + v
  for(i=0;i<l->dim;i++) {
    fprintf(f,"%lf ",l->ref->coords[i]+l->dir_vect->coords[i]);
  }
  fprintf(f,"\n");
  fclose(f);
}

struct point *get_coord_line(struct line *l,double val,int axis) {
  double lambda;
  int i;
  struct point *result;
  if (l->dir_vect->coords[axis]==0)
    return cp_point(l->ref);
  lambda=(val-l->ref->coords[axis])/l->dir_vect->coords[axis];
  //printf("lambda=%f\n",lambda);
  result=new_point(l->dim);
  for(i=0;i<result->dim;i++) {
    if (i==axis) {
      result->coords[i]=val;
    } else {
      result->coords[i]=lambda*l->dir_vect->coords[i]+l->ref->coords[i];
    }
  }
  return result;
}

int best_axis_line(struct line *l) {
  int i;
  int result=0;
  double bestvalue=0.0;
  double cur;
  for(i=0;i<l->dim;i++) {
    cur=fabs(l->dir_vect->coords[i]);
    if (cur>bestvalue) {
      result=i;
      bestvalue=cur;
    }
  }
  return result;
}

struct point **boxed_line(struct line *l,struct point *point_min,struct point *point_max,int nb_steps) {
  int j,k;
  double ref;
  double refmin;
  double refmax;
  struct point *p;
  int best_axis;
  struct point **result;

  //allocate and initialize result
  result=malloc(nb_steps*sizeof(struct point *));
  for(j=0;j<nb_steps;j++) {
    result[j]=NULL;
  }

  //select output axis
  best_axis=best_axis_line(l);
  //printf("best axis is %d\n",best_axis);

  refmin=point_min->coords[best_axis];
  refmax=point_max->coords[best_axis];
  j=0;
  ref=refmin;
  for(k=0;k<nb_steps;k++) {
    p=get_coord_line(l,ref,best_axis);

    if (is_valid_point(p,point_min,point_max)) {
      result[j]=p;
      j++;
    } else {
      result[j]=NULL;
      free_point(p);
    }
    ref+=(refmax-refmin)/(nb_steps-1);
  }
  return result;
}

void dump_boxed_line(struct line *l,char *filename,struct point *point_min,struct point *point_max,int nb_steps) {
  int i,j;
  struct point **points=boxed_line(l,point_min,point_max,nb_steps);
  FILE *f=fopen(filename,"w");
  for(i=0;i<nb_steps;i++) {
    if (points[i]!=NULL) {
      for(j=0;j<l->dim;j++) {
        fprintf(f,"%lf ",points[i]->coords[j]);
      }
      fprintf(f,"\n");
      free_point(points[i]);
    }
  }
  free(points);
  fclose(f);
}

void dump_line(struct line *l,FILE *f) {
  int i;
  for(i=0;i<l->dim;i++) {
    fprintf(f,"%f ",l->ref->coords[i]);
  }
  for(i=0;i<l->dim;i++) {
    fprintf(f,"%f ",l->dir_vect->coords[i]);
  }  
  fprintf(f,"\n");
}

void plot_line(struct line *l,int id,struct graphics *gws) {
  char name[1024];
  char title[1024];
  int i;
  struct point **points;
  if (l->dim==2) {
    if (gws->sp) {
      points=boxed_line(l,gws->sp->point_min,gws->sp->point_max,gws->sp->w);
      for(i=0;i<gws->sp->w;i++) {
        if (points[i]!=NULL) {
          plot_point(points[i],id,gws);
          free_point(points[i]);
        }
      }
      free(points);
    }
  }
  if (l->dim==3) {
    if (gws->gp) {
      get_plotfile(gws->gp,name);
      dump_boxed_line(l,name,gws->point_min,gws->point_max,gws->w);
      sprintf(title,"line %d",id);
      plot_file(gws->gp,l->dim,name,1,1,title);
    }
  }
}

struct line *angle_line(int angle,double x,double y) {
  struct point *ref=new_point(2);
  struct vector *direct=new_vector(2);
  struct line *angline;
  ref->coords[0]=x;
  ref->coords[1]=y;
  direct->coords[0]=cos(deg2rad(angle));
  direct->coords[1]=sin(deg2rad(angle));
  angline=new_line(ref,direct);
  return angline;
}

struct line *random_line(int dim,struct point *pmin,struct point *pmax) {
  struct point *ref;
  struct vector *direct;
  double unif;
  int i;

  ref=new_point(dim);
  direct=new_vector(dim);

  //calculate the reference point
  for(i=0;i<dim;i++) {
    unif=rand()/((double)RAND_MAX);
    ref->coords[i]=pmin->coords[i]+unif*(pmax->coords[i]-pmin->coords[i]);
  }
  //printf("ref : \n");
  //print_point(ref);
    
  //calculate the directional vector
  for(i=0;i<dim;i++) {
    unif=rand()/((double)RAND_MAX);
    direct->coords[i]=unif*2000-1000;
  }
  normalize_vector(direct);
  //printf("directional vector : \n");
  //print_vector(direct);

  //output the line
  return new_line(ref,direct);
}



