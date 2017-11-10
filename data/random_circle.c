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
  int dim;
  int nb_steps;
  char *output_name;
  double seed;
  int i;
  int fixargs=4;
  double unif;
  struct point *point_min;
  struct point *point_max;
  struct point *center;
  struct circle *output;
  double scale;
  double radius;

  //check first parameters
  if (argc<fixargs+1) {
    printf("usage %s dim nb_steps output_filename seed min[dim] max[dim]\n",argv[0]);
    exit(1);
  }


  //read first parameters
  dim=atoi(argv[1]);
  printf("dim=%d\n",dim);
  nb_steps=atoi(argv[2]);
  printf("nb_steps=%d\n",nb_steps);
  output_name=argv[3];
  seed=atof(argv[4]);
  printf("seed=%f\n",seed);
  if (argc<1+fixargs+2*dim) {
    printf("usage %s dim min[dim] max[dim]\n",argv[0]);
    exit(1);
  }

  //init
  srand(seed);
  point_min=new_point(dim);
  point_max=new_point(dim);
  center=new_point(dim);

  //read parameters
  for(i=0;i<dim;i++) {
    point_min->coords[i]=atoi(argv[1+fixargs+i]);
  }
  for(i=0;i<dim;i++) {
    point_max->coords[i]=atoi(argv[1+fixargs+dim+i]);
  }
  printf("point min :\n");
  print_point(point_min);
  printf("point max :\n");
  print_point(point_max);
  scale=dist_points(point_min,point_max);

  //get a random radius
  unif=rand()/((double)RAND_MAX);
  radius=unif*scale/2;
  printf("radius : %f\n",radius);

  //calculate the center point
  for(i=0;i<dim;i++) {
    unif=rand()/((double)RAND_MAX);
    center->coords[i]=point_min->coords[i]+radius+unif*(point_max->coords[i]-point_min->coords[i]-2*radius);
  }
  printf("center : \n");
  print_point(center);

  //output the line
  output=new_circle(center,radius);
  dump_boxed_circle(output,output_name,point_min,point_max,nb_steps);
  
  return 0;
}
