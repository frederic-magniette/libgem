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
#include <sys/time.h>

//check performances

int main(int argc,char **argv) {

  struct gem_ws *result;
  struct dataset *ds;
  double convcrit;
  int scalecrit;
  //struct graphics *gws;
  double noise_level;
  int size;
  FILE *f;
  int i,j;
  int nb_lines;
  struct timeval before;
  struct timeval after;
  double time;
  int nb_err;
  int nb_iter=5000;
  int seed;
  char fname[4096];
  struct distrib *convtime=new_distrib(WGAUSS);
  
  if (argc!=3) {
    printf("usage : %s seed convcrit\n",argv[0]);
    return -1;
  }

  
  //parameters
  seed=atoi(argv[1]);
  convcrit=atof(argv[2]);
  srand(seed);
  noise_level=0.1;
  size=100;

  for(scalecrit=20;scalecrit<60;scalecrit+=10) {
    sprintf(fname,"/tmp/perf_%d.txt",scalecrit);
    f=fopen(fname,"w");
    
    for(nb_lines=1;nb_lines<4;nb_lines++) {
      reinit_distrib(convtime);
      nb_err=0;
      
      for(i=0;i<nb_iter;i++) {
        
        printf("scalecrit %d iter %d for nb_lines=%d\n",scalecrit,i,nb_lines);
        
        //create dataset
        ds=new_empty_dataset(2);
        for(j=0;j<nb_lines;j++)
          add_random_line_dataset(ds,10,size);
        noise_dataset(ds,noise_level);
        
        //exec the multifit
        gettimeofday(&before,NULL);
        result=multifit(ds,convcrit,scalecrit,0,0,NULL);
        remove_degenerated_objects_gem(result,3);
        gettimeofday(&after,NULL);
        time=(after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec);
        add_data_distrib(convtime,time,1);
        //printf("nb_objects=%d\n",result->nb_objects);
        //printf("global distance = %f\n",result->global_dist);
        
        //check for error
        if (result->nb_objects!=nb_lines) {
          nb_err++;
        }
        
        //plot_gem(result,gws);
        //getchar();
        
        //free data
        free_gem(result);
        free_dataset(ds);
      }
      printf("for %d lines error rate=%f pc mean time=%f usec stdev=%f\n",nb_lines,(double)nb_err/(double)nb_iter*100.0,mean_distrib(convtime),stdev_distrib(convtime));
      fprintf(f,"%d %f %f %f\n",nb_lines,(double)nb_err/(double)nb_iter*100.0,mean_distrib(convtime),stdev_distrib(convtime));
    }
    fclose(f);
  }
  return 1;
}

