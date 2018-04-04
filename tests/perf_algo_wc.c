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
  double noise_level;
  int size;
  FILE *f;
  int i,j,k;
  int nb_lines;
  int nb_iter;
  int seed;
  double scalecrit;
  double pc;
  struct point *pmin;
  struct point *pmax;
  //struct point *gmin;
  //struct point *gmax;
  float uerr;
  float serr;
  //struct graphics *gws;
  int stop;
  float abl;
  float angres;

  //error histograms
  #define NB_BINS 8
  int histo[NB_BINS];
  struct line *targets[NB_BINS];
  
  if (argc!=3) {
    printf("usage : %s seed convcrit\n",argv[0]);
    return -1;
  }

  
  //parameters
  seed=atoi(argv[1]);
  convcrit=atof(argv[2]);
  srand(seed);
  noise_level=0.1;
  size=80;
  scalecrit=150;
  nb_iter=5000;
  pmin=new_2D_point(-50,-50);
  pmax=new_2D_point(50,50);
  //gmin=new_2D_point(-200,-200);
  //gmax=new_2D_point(200,200);
  //gws=new_graphics(2,1000,1000,gmin,gmax,param2verb("summary"));
  f=fopen("/tmp/perf_algo_wd.txt","w");
  uerr=0.0;
  angres=7.0;

  //initialize error histogram
  for(nb_lines=1;nb_lines<=NB_BINS;nb_lines++)
    histo[nb_lines-1]=0;

    
  for(nb_lines=1;nb_lines<=NB_BINS;nb_lines++) {
    printf("nb_lines=%d\n",nb_lines);
    for(i=0;i<nb_iter;i++) {
      printf("iter %d\n",i);
      
      //create dataset
      ds=new_empty_dataset(2);
      for(j=0;j<nb_lines;j++) {
        do {
          targets[j]=random_line(2,pmin,pmax);
          stop=1;
          for(k=0;k<j;k++) {
            abl=angle_between_line(targets[k],targets[j]);
            if ((abl<angres) || (abl>180-angres)) {
              printf("abl=%f\n",abl);
              stop=0;
            }
          }
          if (stop==0)
            free_line(targets[j]);
        } while(stop==0);
        add_line_dataset(ds,targets[j],20,size);
      }
      noise_dataset(ds,noise_level);
      
      //exec the multifit
      result=multifit(ds,convcrit,scalecrit,0,0,NULL);
      remove_degenerated_objects_gem(result,3);
      remove_dup_objects_gem(result,0.1);
      
      //check for error
      if (result->nb_objects!=nb_lines) {
        printf("waiting for %d, getting %d\n",nb_lines,result->nb_objects);
        for(j=0;j<nb_lines;j++) 
          for(k=0;k<j;k++)
            printf("angle between line %d and line %d : %f\n",j,k,angle_between_line(targets[k],targets[j]));
        histo[nb_lines-1]++;
        //plot_gem(result,gws);
        //getchar();
      }
      
      //free data
      free_gem(result);
      free_dataset(ds);
      for(j=0;j<nb_lines;j++) {
        free_line(targets[j]);
      }
    }

    pc=(float)histo[nb_lines-1]/(float)nb_iter*100.0;
    if (nb_lines==2)
      uerr=pc;
    serr=uerr*((float)nb_lines)*((float)nb_lines-1)/2.0;
    fprintf(f,"%d %d %f %f %f\n",nb_lines,histo[nb_lines-1],pc,serr,pc-serr);
    printf("%d %d %f %f %f\n",nb_lines,histo[nb_lines-1],pc,serr,pc-serr);

  }
 
  fclose(f);
  return 1;
}

