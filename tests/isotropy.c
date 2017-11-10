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

//check isotropy in 2D

int main(int argc,char **argv) {

  struct gem_ws *result;
  int angle_iter;
  int angle_lines;
  char fname[4096];
  struct dataset *ds;
  double convcrit;
  double scalecrit;
  struct graphics *gws;
  int a1,a2;
  double noise_level;
  static double pmin[2];
  static double pmax[2];
  struct point *point_min;
  struct point *point_max;
  int size;
  FILE *f;
  
  if (argc!=4) {
    printf("usage : %s angle convcrit scalecrit\n",argv[0]);
    return -1;
  }

  
  //parameters
  angle_lines=atoi(argv[1]);
  convcrit=atof(argv[2]);
  scalecrit=atof(argv[3]);
  noise_level=0.1;
  size=100;

  //size of generation space
  pmin[0]=-1*size;
  pmin[1]=-1*size;
  pmax[0]=size;
  pmax[1]=size;
  point_min=new_valued_point(2,pmin);
  point_max=new_valued_point(2,pmax);
  gws=new_graphics(2,1000,1000,point_min,point_max);
  
  //for(angle_lines=10;angle_lines<=90;angle_lines+=10) {
    sprintf(fname,"/tmp/isotropy_%d.txt",angle_lines);
    f=fopen(fname,"w");

    //iterate over base angles
    for(angle_iter=0;angle_iter<180;angle_iter++) {
      
      a1=angle_iter;
      a2=(angle_iter+angle_lines)%360;
      printf("\n\na1=%d a2=%d\n",a1,a2);
      
      //create dataset
      ds=new_empty_dataset(2);
      add_angle_line_dataset(ds,a1,10,0,0,size);
      add_angle_line_dataset(ds,a2,10,0,0,size);
      noise_dataset(ds,noise_level);
      
      //exec the multifit
      result=multifit(ds,convcrit,scalecrit,0,0,gws);
      remove_degenerated_objects_gem(result,3);
      
      //check for anisotropy
      if (result->nb_objects!=2) {
        printf("anisotropy for angle %d : %d\n",angle_iter,result->nb_objects);
        print_degenerated_objects_gem(result);
      }
      
      //plot the result
      plot_gem(result,gws);
      plot_dataset(ds,gws);

      if (result->nb_objects!=2) {
        dump_dataset(ds,"/tmp/problem.txt");
        getchar();
      }
      
      //printf("%f %d\n",result->global_dist,result->nb_objects);
      fprintf(f,"%d %f %d\n",angle_iter,result->global_dist,result->nb_objects);
      //getchar();
      
      //free data
      free_gem(result);
      free_dataset(ds);
    }
    fclose(f);
    //}
  return 1;
}

