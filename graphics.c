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

struct graphics *new_graphics(int dim,int w,int h,struct point *point_min,struct point*point_max) {
  struct graphics *result=malloc(sizeof(struct graphics));
  result->point_min=cp_point(point_min);
  result->point_max=cp_point(point_max);
  result->dim=dim;
  result->w=w;
  result->h=h;
  result->gp=NULL;
  result->sp=NULL;
  switch(dim) {
  case 1:
    printf("graphics : using gnuplot for 1d plots (histos)\n");
    result->gp=new_gplot(w,h);
    break; 
  case 2:
    printf("graphics : using SDL for 2d plots\n");
    result->sp=new_splot(w,h,result->point_min,result->point_max);
    clear_splot(result->sp);
    apply_splot(result->sp);
    break;
  case 3:
    printf("graphics : using gnuplot for 3d plots\n");
    result->gp=new_gplot(w,h);
    break;
  default:
    printf("graphics : unable to represent because dimension is %d\n",dim);
  }
  return result;
}

void free_graphics(struct graphics *gws) {
  if (gws==NULL)
    return;
  free_point(gws->point_min);
  free_point(gws->point_max);
  free_gplot(gws->gp);
  free_splot(gws->sp);
  free(gws);
}

void apply_graphics(struct graphics *gws) {
  if (gws==NULL)
    return;
  if (gws->sp)
    apply_splot(gws->sp);
}



  
