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

struct spiral *new_spiral(struct line *support,double radius,double angular_offset,double angular_speed) {
  struct spiral *result=malloc(sizeof(struct spiral));
  result->dim=support->dim;
  result->support=support;
  result->proj_base=gram_schmidt(support->dir_vect);
  result->radius=radius;
  result->angular_offset=angular_offset;
  result->angular_speed=angular_speed;
  return result;
}

struct spiral *new_based_spiral(struct line *support,double radius,double angular_offset,double angular_speed,struct matrix *base) {
  struct spiral *result=malloc(sizeof(struct spiral));
  result->dim=support->dim;
  result->support=support;
  result->proj_base=base;
  result->radius=radius;
  result->angular_offset=angular_offset;
  result->angular_speed=angular_speed;
  return result;
}

void free_spiral(struct spiral *s) {
  if (!s)
    return;
  free_line(s->support);
  free(s);
}

void print_spiral(struct spiral *s) {
  printf("spiral[%d]=\nsupport=",s->dim);
  print_line(s->support);
  printf("radius=%f angular speed=%f angular offset=%f\n",s->radius,s->angular_speed,s->angular_offset);
  printf("projection base=\n");
  print_matrix(s->proj_base);
}

void dump_boxed_spiral(struct spiral *s,char *filename,struct point *point_min,struct point *point_max,int nb_steps) {
  int i;
  double angle;
  double long_distance;
  struct vector *v;
  struct vector *norm_sup;

  //build a reference line of points
  struct point **points=boxed_line(s->support,point_min,point_max,nb_steps);
  FILE *f=fopen(filename,"w");
  norm_sup=cp_vector(s->support->dir_vect);
  normalize_vector(norm_sup);
  //add the radius*angle to points
  for(i=0;i<nb_steps;i++) {
    if (points[i]!=NULL) {
      //printf("point %d\n",i);
      v=vector_by2points(points[i],s->support->ref);
      long_distance=coef_project_vector(v,norm_sup);
      free_vector(v);
      angle=fmod(s->angular_offset+s->angular_speed*long_distance,360);
      v=col_vector_matrix(s->proj_base,1);
      mult_lambda_vect(v,s->radius*cos(deg2rad(angle)));
      add_vector_point(points[i],v);
      free_vector(v);
      v=col_vector_matrix(s->proj_base,2);
      mult_lambda_vect(v,s->radius*sin(deg2rad(angle)));
      add_vector_point(points[i],v);

      free_vector(v);
      dump_point(points[i],f);
      //printf("for point %d dist to support=%f\n",i,dist_line_point(s->support,points[i]));
    }
  }

  free(points);
  fclose(f);
}

void plot_spiral(struct spiral *s,int id,struct graphics *gws) {
  char name[1024];
  char title[1024];
  get_plotfile(gws->gp,name);
  dump_boxed_spiral(s,name,gws->point_min,gws->point_max,50);
  sprintf(title,"spiral");
  plot_file(gws->gp,3,name,1,1,title);
  /*sprintf(title,"support");
  plot_line(s->support,0,gws);
  plot_point(s->support->ref,0,gws);
  get_plotfile(gws->gp,name);
  dump_base_matrix(s->proj_base,s->support->ref,name);
  sprintf(title,"ortho base");
  plot_file(gws->gp,3,name,0,1,title);*/
}
