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

//************************** CIRCLES *******************************


struct circle *new_circle(struct point *center,double radius) {
  struct circle *result=malloc(sizeof(struct circle));
  result->center=center;
  result->radius=radius;
  result->dim=center->dim;
  return result;
}

struct circle *cp_circle(struct circle *src) {
  struct circle *result=malloc(sizeof(struct circle));
  result->center=cp_point(src->center);
  result->radius=src->radius;
  result->dim=src->dim;
  return result;
}

void free_circle(struct circle *c) {
  free_point(c->center);
  free(c);
}

void print_circle(struct circle *c) {
  printf("circle(dim=%d): center :",c->dim);
  print_point(c->center);
  printf("radius:%f\n",c->radius);
}

double dist_circle_point(struct circle *c,struct point *p) {
  return fabs(dist_points(c->center,p)-c->radius);
}

void plot_circle(struct circle *c,int id,struct graphics *gws) {
  struct point *rp;
  struct point *p1=new_point(2);
  struct point *p2=new_point(2);
  double rx;
  double a=c->center->coords[0];
  double b=c->center->coords[1];
  double rad=c->radius;
  double sub;
  int x;
  if (gws->sp) {
    //plot the center
    p1->coords[0]=a;
    p1->coords[1]=b;
    //plot_point(p1,id,gws);

    for(x=0;x<gws->sp->w;x++) {
      //get the real coordinate if x
      rp=sdl2real_splot(x,0,gws->sp);
      rx=rp->coords[0];
      free_point(rp);
      sub=pow(rad,2)-pow(rx-a,2);
      if (sub<0)
        continue;
      p1->coords[0]=rx;
      p1->coords[1]=b+sqrt(sub);
      p2->coords[0]=rx;
      p2->coords[1]=b-sqrt(sub);
      //plot the point
      plot_point(p2,id,gws);
      plot_point(p1,id,gws);
    }
  }
  free_point(p1);
  free_point(p2);
}

void dump_boxed_circle(struct circle *c,char *filename,struct point *point_min,struct point *point_max,int nb_steps) {
  int i;
  double x;
  double xmin=point_min->coords[0];
  double xmax=point_max->coords[0];
  struct point *p1=new_point(c->dim);
  struct point *p2=new_point(c->dim);
  double a=c->center->coords[0];
  double b=c->center->coords[1];
  double rad=c->radius;
  double sub;
  FILE *f=fopen(filename,"w");
  
  if (xmin==xmax) {
    printf("error xmin=xmax=%f\n",xmin);
    fclose(f);
    return;
  }

 
  for(x=xmin;x<xmax;x+=(xmax-xmin)/nb_steps) {

    //calc the sqrt internal
    sub=pow(rad,2)-pow(x-a,2);
    if (sub<0)
      continue;
    
    //calc the high point
    p1->coords[0]=x;
    p1->coords[1]=b+sqrt(sub);
    for(i=0;i<c->dim;i++) {
      fprintf(f,"%lf ",p1->coords[i]);
    }
    fprintf(f,"\n");

    //calc the low point
    p2->coords[0]=x;
    p2->coords[1]=b-sqrt(sub);
    for(i=0;i<c->dim;i++) {
      fprintf(f,"%lf ",p2->coords[i]);
    }
    fprintf(f,"\n");
    
  }
  free_point(p1);
  free_point(p2);
  fclose(f);
}

void dump_circle(struct circle *c,FILE *f) {
  int i;
  for(i=0;i<c->dim;i++) {
    fprintf(f,"%f ",c->center->coords[i]);
  }
  fprintf(f,"%f ",c->radius);
  fprintf(f,"\n");
}



