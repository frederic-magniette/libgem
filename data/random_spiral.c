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

int main(int argc,char **argv) {
  int nb_steps;
  char *output_name;
  double seed;
  int i;
  struct point *point_min;
  struct point *point_max;
  struct line *support;
  double radius;
  double angular_speed;
  double angular_offset;
  struct spiral *s;
  double length;

  //check first parameters
  if (argc!=8) {
    printf("usage %s length nb_steps radius angular_speed angular_offset output_filename seed\n",argv[0]);
    exit(1);
  }


  //read first parameters
  length=atof(argv[1]);
  nb_steps=atoi(argv[2]);
  radius=atof(argv[3]);
  angular_speed=atof(argv[4]);
  angular_offset=atof(argv[5]);
  output_name=argv[6];
  seed=atof(argv[7]);

  //init
  srand(seed);
  point_min=new_point(3);
  point_max=new_point(3);
  for(i=0;i<3;i++) {
    point_max->coords[i]+=(length+radius);
  }

  //create the spiral
  support=random_line(3,point_min,point_max);
  dump_boxed_line(support,"test_spiral_support.txt",point_min,point_max,nb_steps);
  s=new_spiral(support,radius,angular_offset,angular_speed);
  print_spiral(s);
  dump_boxed_spiral(s,output_name,point_min,point_max,nb_steps);
  struct graphics *gws=new_graphics(3,500,500,point_min,point_max,1);
  plot_spiral(s,0,gws);
  getchar();
  
  return 0;
}
