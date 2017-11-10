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

//check performances with limited error

int main(int argc,char **argv) {

  struct gem_ws *result;
  struct dataset *ds;
  double convcrit;
  int scalecrit;
  struct graphics *gws;
  double noise_level;
  static double pmin[2];
  static double pmax[2];
  struct point *point_min;
  struct point *point_max;
  int size;
  FILE *f;
  int i,j,k;
  int nb_lines;
  struct timeval before;
  struct timeval after;
  double time;
  int nb_err;
  int nb_iter;
  int seed;
  char fname[4096];
  struct distrib *convtime=new_distrib(WGAUSS);
  int maxlines=3;
  int angles[maxlines];
  int angres;
  int ok;
  struct line *rline;
  struct point *ref;
  
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
  nb_iter=1000;

  //size of generation space
  pmin[0]=-1*size;
  pmin[1]=-1*size;
  pmax[0]=size;
  pmax[1]=size;
  point_min=new_valued_point(2,pmin);
  point_max=new_valued_point(2,pmax);
  gws=new_graphics(2,500,500,point_min,point_max);


  for(scalecrit=20;scalecrit<60;scalecrit+=10) {
    sprintf(fname,"/tmp/cer_%d.txt",scalecrit);
    f=fopen(fname,"w");
    angres=20;
    
    for(nb_lines=2;nb_lines<=maxlines;nb_lines++) {
      reinit_distrib(convtime);
      nb_err=0;
      
      for(i=0;i<nb_iter;i++) {
        
        printf("scalecrit %d iter %d for nb_lines=%d\n",scalecrit,i,nb_lines);
        
        //create dataset
        ds=new_empty_dataset(2);
        for(j=0;j<nb_lines;j++) {
          do {
            ok=1;
            angles[j]=(int)(rand()/((double)RAND_MAX)*180);
            printf("new angle : %d\n",angles[j]);
            for(k=0;k<j;k++) {
              double diffa=abs(angles[j]-angles[k]);
              printf("diffa=%f\n",diffa);
              if (diffa<angres) {
                printf("incompatible with %d=%d\n",k,angles[k]);
                ok=0;
              }
              if ((diffa>180-angres)&&(diffa<180+angres)) {
                printf("incompatible with %d=%d\n",k,angles[k]);
                ok=0;
              }
            }
          } while(ok==0);
          ref=new_random_point(2,point_min,point_max);
          rline=angle_line(angles[j],ref->coords[0],ref->coords[1]);
          free_point(ref);
          add_line_dataset(ds,rline,10,100);
          free_line(rline);
        }
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
          plot_gem(result,gws);
          getchar();
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

