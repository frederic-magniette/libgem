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

//*********************** GAUSSIAN LINES ***************

struct gline *new_gline(struct dataset *ds,struct weights *w,struct graphics *gws) {
  int p;
  struct gline *result=malloc(sizeof(struct gline));
  result->line=malloc(sizeof(struct line));

  //calculate the reference point
  result->line->ref=weighted_barycenter_point(ds,w);
  //print_point(result->line->ref);
  //printf("plot ref");
  //plot_point(result->line->ref,2,gws);
  //apply_graphics(gws);
  //getchar();

  result->line->dir_vect=weighted_mean_vector(ds,result->line->ref,w);
  result->line->dim=result->line->ref->dim;

  result->pgauss=new_distrib(CWGAUSS);
  for(p=0;p<ds->nb_points;p++) {
    //printf("adding distance %f\n",dist_line_point(result->line,ds->points[p]));
    add_data_distrib(result->pgauss,dist_line_point(result->line,ds->points[p]),w->coefs[p]);
  }

  return result;
}

struct gline *cp_gline(struct gline *gl) {
  struct gline *result=malloc(sizeof(struct gline));
  result->line=cp_line(gl->line);
  result->pgauss=cp_distrib(gl->pgauss);
  return result;
}

void print_gline(struct gline *gl) {
  print_line(gl->line);
  print_distrib(gl->pgauss);
}

void free_gline(struct gline *gl) {
  free_line(gl->line);
  free_distrib(gl->pgauss);
  free(gl);
}


double belonging_proba_gline(struct gline *gl,struct point *p) {
  return normalized_likelyhood_distrib(gl->pgauss,dist_line_point(gl->line,p));
}

double dist_gline(struct gline *gl,struct point *p) {
  return dist_line_point(gl->line,p);
}

void dump_boxed_gline(struct gline *gl,char *filename,struct point *point_min,struct point *point_max,int nb_steps) {
  dump_boxed_line(gl->line,filename,point_min,point_max,nb_steps);
}

void dump_gline(struct gline *gl,FILE *f) {
  dump_line(gl->line,f);
}

void plot_gline(struct gline *gl,int id,struct graphics *gws) {
  plot_line(gl->line,id,gws);
}

void plot_field_gline(struct gline *gl,int id,struct graphics *gws) {
  int rref,gref,bref;
  double distance;
  int x,y;
  int r,g,b;
  double proba;
  struct point *pt;
  double xmin,xmax;
  char title[1024];
  char name[1024];

  
  if (gl->line->dim==1) {
    if (gws->gp) {
      xmin=gws->point_min->coords[0];
      xmax=gws->point_max->coords[0];
      get_plotfile(gws->gp,name);
      dump_distrib(gl->pgauss,gl->line->ref->coords[0],xmin,xmax,100,name);
      sprintf(title,"distrib %d",id);
      plot_file(gws->gp,2,name,1,1,title);
    }
  }

    
  if (gl->line->dim==2) {
    if (gws->sp) {
      //plot the field
      for(x=0;x<gws->w;x++) {
        for(y=0;y<gws->h;y++) {
          
          //calc the proba
          pt=sdl2real_splot(x,y,gws->sp);
          distance=dist_gline(gl,pt);
          free_point(pt);
          proba=log(1+normalized_likelyhood_distrib(gl->pgauss,distance)*sqrt(2*PI));
          
          //calc the color : proportional and thresholded
          get_pixel_splot(gws->sp,x,y,&r,&g,&b);
          attribute_color_splot(id,&rref,&gref,&bref);
          r+=(int)(proba*rref);
          g+=(int)(proba*gref);
          b+=(int)(proba*bref);
          
          set_pixel_splot(gws->sp,x,y,r,g,b);
        }
      }
    }
  }
}


int same_gline(struct gline *gl1,struct gline *gl2,double distlim) {
  double dist_center;
  double dist_plusdir;
  struct point *p=cp_point(gl2->line->ref);
  struct vector *dir=cp_vector(gl2->line->dir_vect);
  mult_lambda_vect(dir,distlim*10);
  add_vector_point(p,dir);
  dist_center=dist_line_point(gl1->line,gl2->line->ref);
  dist_plusdir=dist_line_point(gl1->line,p);
  free_point(p);
  free_vector(dir);
  if (dist_center>distlim)
    return 0;
  if (dist_plusdir>distlim)
    return 0;
  return 1;
}
