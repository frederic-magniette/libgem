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

int main() {
  double min[2]={0.0,0.0};
  double max[2]={10.0,10.0};
  double center[2]={5.0,5.0};
  double dir[2]={1.0,1.0};
  struct point *point_min=new_valued_point(2,min);
  struct point *point_max=new_valued_point(2,max);
  struct graphics *gws=new_graphics(2,600,600,point_min,point_max,1);
  struct vector *dir_vect=new_valued_vector(2,dir);
  struct point *ref=new_valued_point(2,center);
  struct line *l=new_line(ref,dir_vect);
  struct circle *c=new_circle(ref,3);
  clear_splot(gws->sp);
  plot_line(l,1,gws);
  plot_circle(c,2,gws);
  plot_point(ref,0,gws);
  apply_splot(gws->sp);
  printf("hit a key to finish\n");
  getchar();
  return 1;
}
