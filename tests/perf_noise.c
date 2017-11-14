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

  struct gem_ws *result;
  struct dataset *ds;
  struct graphics *gws=NULL;
  double noise_level;
  double scalecrit,convcrit;
  static double pmin[2];
  static double pmax[2];
  struct point *point_min;
  struct point *point_max;
  int size;
  int nbp;
  int rep;
  FILE *f;
  int nberr,nbrep;
  int seed;
  
  if (argc!=5) {
    printf("usage : %s convcrit scalecrit output[detail|summary|none] seed\n",argv[0]);
    return -1;
  }

  
  //parameters
  convcrit=atof(argv[1]);
  scalecrit=atof(argv[2]);
  seed=atoi(argv[4]);
  noise_level=0.1;
  int nbppl=13;
  size=100;
  if (!strcmp(argv[3],"detail"))
    nbrep=1;
  else
    nbrep=1000;

  srand(seed);
  //size of generation space
  pmin[0]=-1*size;
  pmin[1]=-1*size;
  pmax[0]=size;
  pmax[1]=size;
  point_min=new_valued_point(2,pmin);
  point_max=new_valued_point(2,pmax);
  int o=param2verb(argv[3]);
  if (o!=0)
    gws=new_graphics(2,500,500,point_min,point_max,o);
  
  f=fopen("/tmp/perf_noise.txt","w");
  fprintf(f,"0.0 0.0\n");
  for(nbp=1;nbp<=nbppl;nbp++) {
    printf("iteration with %d points of noise\n",nbp);
    nberr=0;
    for(rep=0;rep<nbrep;rep++) {
      //create dataset
      ds=new_empty_dataset(2);
      add_angle_line_dataset(ds,10,nbppl,0,0,size);
      add_angle_line_dataset(ds,50,nbppl,0,0,size);
      noise_dataset(ds,noise_level);
      add_random_noise_dataset(ds,nbp);
      
      //exec the multifit
      result=multifit(ds,convcrit,scalecrit,0,0,gws);
      remove_degenerated_objects_gem(result,15*ds->nb_points/100);
      remove_dup_objects_gem(result,result->scale/20);
      
      //check for error
      if (result->nb_objects!=2) {
	printf("error for %d points of noise\n",nbp);
	nberr++;
      }

      //plot the result
      if (!strcmp(argv[3],"detail")||!strcmp(argv[3],"summary")) {
	plot_gem(result,gws);
	plot_dataset(ds,gws);
        if (!strcmp(argv[3],"detail"))
          getchar();  
      }
      
      //free data
      free_gem(result);
      free_dataset(ds);
    }
    fprintf(f,"%f %f\n",(double)nbp/(double)(2*nbppl+nbp)*100.0,((double)nberr/(double)nbrep)*100.0);
    fflush(f);
  }
  fclose(f);
  return 1;
}

